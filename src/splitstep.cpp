
#include <splitstep.hpp>
#include <omp.h>

#define EIGEN_DONT_VECTORIZE
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_NO_DEBUG

// #define SLICE_NUMBER 0

using namespace std;
using namespace Eigen;

SplitStep::SplitStep(Options &o) : opt(o) {}

void SplitStep::assignMatrixData(shared_ptr<MatrixData> d){
	w = d;
	setVariables();
}

// void SplitStep::setOptions(const Options &externaloptions){
// 	opt = externaloptions;
// }

void SplitStep::setVariables(){

  	x_axis.resize(w->meta.grid[0]);
  	y_axis.resize(w->meta.grid[1]);
  	for(int i=0;i<w->meta.grid[0];i++){
  		x_axis[i]=-w->meta.initCoord[0]+i*real(w->meta.initSpacing[0]);
  	}
  	for(int j=0;j<w->meta.grid[1];j++){
  		y_axis[j]=-w->meta.initCoord[1]+j*real(w->meta.initSpacing[1]);
  	}

  	VectorXcd X = VectorXcd(w->meta.grid[0]);
  	VectorXcd Y = VectorXcd(w->meta.grid[1]);
	for(int i = 0;i<w->meta.grid[0];i++){
		X(i) = complex<double>(x_axis[i],0.0);
	}
	for(int j = 0;j<w->meta.grid[1];j++){
		Y(j) = complex<double>(y_axis[j],0.0);
	}

	PotentialGrid = MatrixXcd::Zero(w->meta.grid[0],w->meta.grid[1]);
   	for(int i = 0; i < w->meta.grid[0]; i++){
   		for(int j = 0; j < w->meta.grid[1]; j++){
			PotentialGrid(i,j) = ( half * opt.omega_x * opt.omega_x * X(i) * X(i) +  half * opt.omega_y * opt.omega_y * Y(j) * Y(j) );
		}
	}

	kprop = MatrixXcd(w->meta.grid[0],w->meta.grid[1]);
	kpropITP = MatrixXcd(w->meta.grid[0],w->meta.grid[1]);
	kprop_x = MatrixXcd(w->meta.grid[0],w->meta.grid[1]);
	kprop_y = MatrixXcd(w->meta.grid[0],w->meta.grid[1]);
	kprop_x_strang = MatrixXcd(w->meta.grid[0],w->meta.grid[1]);
	kprop_y_strang = MatrixXcd(w->meta.grid[0],w->meta.grid[1]);
	kspace.resize(2);


	
	vector<vector<double>> kspace;
	vector<double> Kmax(2);
	Kmax[0] = M_PI / w->meta.spacing[0];
	Kmax[1] = M_PI / w->meta.spacing[1];
	vector<double> deltaK(2);
	deltaK[0] = Kmax[0] / (w->meta.grid[0] / 2.0);
	deltaK[1] = Kmax[1] / (w->meta.grid[1] / 2.0);

	double Tprinted = ( Kmax[0] * Kmax[0] + Kmax[1] * Kmax[1] );
	cout << "k squared = " << Tprinted << endl;
	Tprinted *= - 0.5 * opt.RTE_step;
	cout << "T = " << Tprinted << endl;

	double Vprinted = - opt.g * abs2(w->wavefunction[0](w->meta.grid[0]/2,w->meta.grid[1]/2)) * opt.RTE_step;
	cout << "V = " << Vprinted << endl;


	kspace.resize(2);
	for(int d = 0; d < 2; d++){
		// set k-space
		kspace[d].resize(w->meta.grid[d]);
		for(int i = 0; i <= w->meta.grid[d]/2; i++){
			// kspace[d][i] = (M_PI / rmax[d]) * (double)i;
			kspace[d][i] = deltaK[d] * (double)i;
		}
		for(int i = (w->meta.grid[d]/2)+1; i < w->meta.grid[d]; i++){
			// kspace[d][i] = -(M_PI / rmax[d]) * (double)(w->meta.grid[d] - i);
			kspace[d][i] = - deltaK[d] * (double)(w->meta.grid[d] - i);
		}
	}

	// plotVector("kspace", kspace[0], kspace[1], opt);

	#pragma omp parallel for
	for(int x = 0; x < w->meta.grid[0]; x++){
	    for(int y = 0; y < w->meta.grid[1]; y++){
	      	double T = - 0.5 * (kspace[0][x]*kspace[0][x] + kspace[1][y]*kspace[1][y]) * opt.RTE_step;
      		kprop(x,y) = complex<double>(cos(T),sin(T)) / complex<double>((double)(w->meta.grid[0]*w->meta.grid[1]),0.0);	    
	    }
	}

	cout << "Is Dimensionless " << w->meta.isDimensionless << endl;
	plotDataToPngEigen("kprop"+to_string(w->meta.steps), kprop, opt);

	#pragma omp parallel for
	for(int x = 0; x < w->meta.grid[0]; x++){
	    for(int y = 0; y < w->meta.grid[1]; y++){
	      	double T = - 0.5 * (kspace[0][x]*kspace[0][x] + kspace[1][y]*kspace[1][y]) * opt.RTE_step;
      		kpropITP(x,y) = complex<double>(exp(T),0.0) / complex<double>((double)(w->meta.grid[0]*w->meta.grid[1]),0.0);	    
	    }
	}

	vector<double> k_x_space(w->meta.grid[0]);
	vector<double> k_y_space(w->meta.grid[1]);

	for(int j = 0; j < w->meta.grid[0]; j++){
		int p = j - meta.grid[0]/2;
		k_x_space[j] = (double)p * M_PI / (w->meta.coord[0]);
	}

	for(int j = 0; j < w->meta.grid[1]; j++){
		int p = j - meta.grid[1]/2;
		k_y_space[j] = (double)p * M_PI / (w->meta.coord[1]);
	}

	#pragma omp parallel for
	for(int x = 0; x < w->meta.grid[0]; x++){
	    for(int y = 0; y < w->meta.grid[1]; y++){
	      	double T = - ( 0.5 * kspace[0][x]* kspace[0][x] + opt.omega_w.real() * y_axis[y] * kspace[0][x]) * opt.RTE_step;
      		kprop_x(x,y) = complex<double>(cos(T),sin(T))/* / complex<double>((double)(w->meta.grid[0]*w->meta.grid[1]),0.0)*/;	    
	    }
	}
	// plotDataToPngEigen("kprop_x", kprop_x, opt);

    #pragma omp parallel for
	for(int x = 0; x < w->meta.grid[0]; x++){
	    for(int y = 0; y < w->meta.grid[1]; y++){
	      	double T = - ( 0.5 * kspace[1][y] * kspace[1][y] - opt.omega_w.real() * x_axis[x] * kspace[1][y]) * opt.RTE_step;
      		kprop_y(x,y) = complex<double>(cos(T),sin(T)) /*/ complex<double>((double)(w->meta.grid[0]*w->meta.grid[1]),0.0)*/;	    
	    }
	}
	// plotDataToPngEigen("kprop_y", kprop_y, opt);

	#pragma omp parallel for
	for(int x = 0; x < w->meta.grid[0]; x++){
	    for(int y = 0; y < w->meta.grid[1]; y++){
	      	double T = - 0.5 * ( 0.5 * kspace[0][x]* kspace[0][x] + opt.omega_w.real() * y_axis[y] * kspace[0][x]) * opt.RTE_step;
      		kprop_y_strang(x,y) = complex<double>(cos(T),sin(T))/* / complex<double>((double)(w->meta.grid[0]*w->meta.grid[1]),0.0)*/;	    
	    }
	}
	// plotDataToPngEigen("kprop_x", kprop_x, opt);

    #pragma omp parallel for
	for(int x = 0; x < w->meta.grid[0]; x++){
	    for(int y = 0; y < w->meta.grid[1]; y++){
	      	double T = - 0.5 * ( 0.5 * kspace[1][y] * kspace[1][y] - opt.omega_w.real() * x_axis[x] * kspace[1][y]) * opt.RTE_step;
      		kprop_x_strang(x,y) = complex<double>(cos(T),sin(T)) /*/ complex<double>((double)(w->meta.grid[0]*w->meta.grid[1]),0.0)*/;	    
	    }
	}
	// plotDataToPngEigen("kprop_y", kprop_y, opt);

	// The following is very very unaesthetic, but necessary, because the fftw_plan in the eigenFFT implementation gets set
    // when the first fwd() or inv() gets called, but creating and destroying a plan is not threadsave, so to use #pragma omp parallel for later
    // one has to call beforehand to set the plan. only fixable, by implementing the FFT oneself.
    w->fft.Forward();
    w->fft.Backward();
    w->wavefunction[0].array() /= complex<double>((double)(w->meta.grid[0]*w->meta.grid[1]),0.0);
}

void SplitRot::timeStep(double delta_t){

	w->fft.Forward_X();
	w->wavefunction[0].array() *= kprop_y.array();
	w->fft.Backward_X();

	w->fft.Forward_Y();
	w->wavefunction[0].array() *= kprop_x.array();
	w->fft.Backward_Y();

	w->wavefunction[0].array() /= complex<double>((double)(w->meta.grid[0]*w->meta.grid[1]),0.0);

	#pragma omp parallel for
	for(int x = 0; x < w->meta.grid[0]; x++){
		for(int y = 0; y < w->meta.grid[1]; y++){
    		// complex<double> value = wavefctVec[i](0,x,y,z);
    		double V = - ( PotentialGrid(x,y).real() + opt.g * abs2(w->wavefunction[0](x,y)) ) * delta_t;
    		// potPlotGrid(0,x,y,0) = complex<double>(rotatingPotential(x,y,m) /*PotentialGrid(x,y).real()*/,0.0);
    		// potGrid(0,x,y,0) = complex<double>(cos(V),sin(V));
    		w->wavefunction[0](x,y) *= complex<double>(cos(V),sin(V));
		}
	}

	w->increment(delta_t, 1.0, 1.0);
}

void SplitTrap::timeStep(double delta_t){

	w->fft.Forward();
	w->wavefunction[0].array() *= kprop.array();
	w->fft.Backward(); 

	#pragma omp parallel for
	for(int x = 0; x < w->meta.grid[0]; x++){
		for(int y = 0; y < w->meta.grid[1]; y++){
    		// complex<double> value = wavefctVec[i](0,x,y,z);
    		double V = - ( PotentialGrid(x,y).real() + opt.g * abs2(w->wavefunction[0](x,y)) ) * delta_t;
    		// potPlotGrid(0,x,y,0) = complex<double>(rotatingPotential(x,y,m) /*PotentialGrid(x,y).real()*/,0.0);
    		// potGrid(0,x,y,0) = complex<double>(cos(V),sin(V));
    		w->wavefunction[0](x,y) *= complex<double>(cos(V),sin(V));
		}
	}

	w->increment(delta_t, 1.0, 1.0);
}

void SplitITP::timeStep(double delta_t){

	w->fft.Forward();
	w->wavefunction[0].array() *=  kpropITP.array();
	w->fft.Backward(); 

	#pragma omp parallel for
	for(int x = 0; x < w->meta.grid[0]; x++){
		for(int y = 0; y < w->meta.grid[1]; y++){
    		// complex<double> value = wavefctVec[i](0,x,y,z);
    		double V = - ( PotentialGrid(x,y).real() + opt.g * abs2(w->wavefunction[0](x,y)) ) * delta_t;
    		// potPlotGrid(0,x,y,0) = complex<double>(rotatingPotential(x,y,m) /*PotentialGrid(x,y).real()*/,0.0);
    		// potGrid(0,x,y,0) = complex<double>(cos(V),sin(V));
    		w->wavefunction[0](x,y) *= complex<double>(exp(V),0.0);
		}
	}

	w->increment(delta_t, 1.0, 1.0);
}

void SplitFree::timeStep(double delta_t){

	w->fft.Forward();
	w->wavefunction[0].array() *= kprop.array();
	w->fft.Backward(); 

	#pragma omp parallel for
	for(int x = 0; x < w->meta.grid[0]; x++){
		for(int y = 0; y < w->meta.grid[1]; y++){
    		// complex<double> value = wavefctVec[i](0,x,y,z);
    		double V = - opt.g * abs2(w->wavefunction[0](x,y)) * delta_t;
    		// potPlotGrid(0,x,y,0) = complex<double>(rotatingPotential(x,y,m) /*PotentialGrid(x,y).real()*/,0.0);
    		// potGrid(0,x,y,0) = complex<double>(cos(V),sin(V));
    		w->wavefunction[0](x,y) *= complex<double>(cos(V),sin(V));
		}
	}

	w->increment(delta_t, 1.0, 1.0);
}

void SplitRotStrang::timeStep(double delta_t){

	w->fft.Forward_X();
	w->wavefunction[0].array() *= kprop_x_strang.array();
	w->fft.Backward_X();

	w->fft.Forward_Y();
	w->wavefunction[0].array() *= kprop_y_strang.array();
	w->fft.Backward_Y();

	w->wavefunction[0].array() /= complex<double>((double)(w->meta.grid[0]*w->meta.grid[1]),0.0);

	#pragma omp parallel for
	for(int x = 0; x < w->meta.grid[0]; x++){
		for(int y = 0; y < w->meta.grid[1]; y++){
    		// complex<double> value = wavefctVec[i](0,x,y,z);
    		double V = - ( PotentialGrid(x,y).real() + opt.g * abs2(w->wavefunction[0](x,y)) ) * delta_t;
    		// potPlotGrid(0,x,y,0) = complex<double>(rotatingPotential(x,y,m) /*PotentialGrid(x,y).real()*/,0.0);
    		// potGrid(0,x,y,0) = complex<double>(cos(V),sin(V));
    		w->wavefunction[0](x,y) *= complex<double>(cos(V),sin(V));
		}
	}

	w->fft.Forward_Y();
	w->wavefunction[0].array() *= kprop_y_strang.array();
	w->fft.Backward_Y();

	w->fft.Forward_X();
	w->wavefunction[0].array() *= kprop_x_strang.array();
	w->fft.Backward_X();

	w->wavefunction[0].array() /= complex<double>((double)(w->meta.grid[0]*w->meta.grid[1]),0.0);

	w->increment(delta_t, 1.0, 1.0);
}
