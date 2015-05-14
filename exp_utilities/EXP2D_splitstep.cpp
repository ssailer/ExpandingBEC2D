
#include <EXP2D_splitstep.hpp>
#include <omp.h>

#define EIGEN_DONT_VECTORIZE
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_NO_DEBUG

// #define SLICE_NUMBER 0

using namespace std;
using namespace Eigen;

SplitStep::SplitStep(Options &o) : opt(o) {}

void SplitStep::assignMatrixData(MatrixData* &d){
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
  		x_axis[i]=-opt.min_x+i*real(w->meta.initSpacing[0]);
  	}
  	for(int j=0;j<w->meta.grid[1];j++){
  		y_axis[j]=-opt.min_y+j*real(w->meta.initSpacing[1]);
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
	kprop_x = MatrixXcd(w->meta.grid[0],w->meta.grid[1]);
	kprop_y = MatrixXcd(w->meta.grid[0],w->meta.grid[1]);
	kspace.resize(2);


	
	for(int d = 0; d < 2; d++){

		kspace[d].resize(w->meta.grid[d]);
		for(int i = 0; i <= w->meta.grid[d]/2; i++){
			kspace[d][i] = (M_PI / w->meta.coord[d]) * (double)i;
		}

		for(int i = (w->meta.grid[d]/2)+1; i < w->meta.grid[d]; i++){
			kspace[d][i] = -(M_PI / w->meta.coord[d]) * (double)(w->meta.grid[d] - i);
		}
	}

	plotVector("kspace", kspace[0], kspace[1], opt);

	#pragma omp parallel for
	for(int x = 0; x < w->meta.grid[0]; x++){
	    for(int y = 0; y < w->meta.grid[1]; y++){
	      	double T = - 0.5 * (kspace[0][x]*kspace[0][x] + kspace[1][y]*kspace[1][y]) * opt.RTE_step;
      		kprop(x,y) = complex<double>(cos(T),sin(T)) / complex<double>((double)(w->meta.grid[0]*w->meta.grid[1]),0.0);	    
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
	plotDataToPngEigen("kprop_x", kprop_x, opt);

    #pragma omp parallel for
	for(int x = 0; x < w->meta.grid[0]; x++){
	    for(int y = 0; y < w->meta.grid[1]; y++){
	      	double T = - ( 0.5 * kspace[1][y] * kspace[1][y] - opt.omega_w.real() * x_axis[x] * kspace[1][y]) * opt.RTE_step;
      		kprop_y(x,y) = complex<double>(cos(T),sin(T)) /*/ complex<double>((double)(w->meta.grid[0]*w->meta.grid[1]),0.0)*/;	    
	    }
	}
	plotDataToPngEigen("kprop_y", kprop_y, opt);
}

void SplitStep::timeStep(double delta_t){

	// w->fftForward();
	// w->wavefunction[0].array() *= kprop.array();

	w->fftForward_X();
	// plotDataToPngEigen("wfct_x_fft", w->wavefunction[0], opt);	
	w->wavefunction[0].array() *= kprop_y.array();
	w->fftBackward_X();

	w->fftForward_Y();
	w->wavefunction[0].array() *= kprop_x.array();
	w->fftBackward_Y();


	

	// MatrixXcd tmp = MatrixXcd(w->meta.grid[0],w->meta.grid[1]);
	// #pragma omp parallel for
	// for(int i = 0; i < w->meta.grid[0]; i++){
	// 	for(int j = 0; j < w->meta.grid[1]; j++){
	// 		tmp(j,i) = w->wavefunction[0](i,j);			
	// 	}
	// }
	// w->wavefunction[0] = tmp;

	// w->fftForward_Y();
	// // 	// plotDataToPngEigen("wfct_y_fft", w->wavefunction[0], opt);
	// w->wavefunction[0].array() *= kprop_y.array();
	// w->fftBackward_Y();

	// #pragma omp parallel for
	// for(int i = 0; i < w->meta.grid[0]; i++){
	// 	for(int j = 0; j < w->meta.grid[1]; j++){
	// 		tmp(j,i) = w->wavefunction[0](i,j);			
	// 	}
	// }
	// w->wavefunction[0] = tmp;



	// w->fftForward_X();
	// w->wavefunction[0].array() *= kprop_y.array();
	// w->fftBackward_X();

	// w->wavefunction[0].array() = w->wavefunction[0].transpose().array();

	// w->fftForward_X();
	// w->wavefunction[0].array() *= kprop_x.array();
	// w->fftBackward_X();

	// w->wavefunction[0].array() = w->wavefunction[0].transpose().array();




	w->wavefunction[0].array() /= complex<double>((double)(w->meta.grid[0]*w->meta.grid[1]),0.0);


	// w->fftBackward(); 

	// Vgrid.array() = exp( - complex<double>(opt.g,0.0) * w->wavefunction[0].conjugate().array() * w->wavefunction[0].array() * delta_t);
	// w->wavefunction[0].array() *= Vgrid.array();

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


	// // // plotDataToPngEigen("wfct_x_fft", w->wavefunction[0], opt);
	// w->fftForward_Y();
	// w->wavefunction[0].array() *= kprop_x.array();
	
	// w->fftBackward_Y();

	// w->fftForward_X();
	// w->wavefunction[0].array() *= kprop_y.array();
	// w->fftBackward_X();

	// w->wavefunction[0].array() = w->wavefunction[0].transpose().array();

	// w->fftForward_X();
	// w->wavefunction[0].array() *= kprop_x.array();
	// w->fftBackward_X();

	// w->wavefunction[0].array() = w->wavefunction[0].transpose().array();

	// w->fftForward_X();
	// w->wavefunction[0].array() *= kprop_y.array();
	// w->fftBackward_X();

	// w->fftForward_Y();
	// w->wavefunction[0].array() *= kprop_y.array();
	// w->fftBackward_Y();


	// w->fftForward_X();
	// w->wavefunction[0].array() *= kprop_x.array();
	// w->fftBackward_X();



	// w->fftBackward();

	// w->wavefunction[0].array() /= complex<double>((double)(w->meta.grid[0]*w->meta.grid[1]),0.0);

	w->increment(delta_t, 1.0, 1.0);
	// cout << "potential applied" << endl;

	
		
				// ComplexGrid::fft_unnormalized(wavefctVec[i], kgrid, true);
    
				// #pragma omp parallel for
				// for(int x = 0; x < opt.grid[1]; x++){
				// 	for(int y = 0; y < opt.grid[2]; y++){
				// 		for(int z = 0; z < opt.grid[3]; z++){				
				//    			kgrid(0,x,y,z) = kprop(0,x,y,z) * kgrid(0,x,y,z);
				//    		}
				// 	}
				// }

				// ComplexGrid::fft_unnormalized(kgrid, wavefctVec[i], false);

				// #pragma omp parallel for
				// for(int x = 0; x < opt.grid[1]; x++){
				// 	for(int y = 0; y < opt.grid[2]; y++){
				// 		for(int z = 0; z < opt.grid[3]; z++){	
				//     		// complex<double> value = wavefctVec[i](0,x,y,z);
				//     		double V = - ( /*PotentialGrid(x,y).real()*/ /*rotatingPotential(x,y,m)*//* +*/ opt.g * abs2(wavefctVec[i](0,x,y,z)) ) * timestepsize;
				//     		// potPlotGrid(0,x,y,0) = complex<double>(rotatingPotential(x,y,m) /*PotentialGrid(x,y).real()*/,0.0);
				//     		// potGrid(0,x,y,0) = complex<double>(cos(V),sin(V));
				//     		wavefctVec[i](0,x,y,z) = complex<double>(cos(V),sin(V)) * wavefctVec[i](0,x,y,z);
				//     	}
				// 	}
				// }
					
			

}


void SplitPotential::timeStep(double delta_t){

	w->fftForward();
	w->wavefunction[0].array() *= kprop.array();
	w->fftBackward(); 

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

void SplitFree::timeStep(double delta_t){

	w->fftForward();
	w->wavefunction[0].array() *= kprop.array();
	w->fftBackward(); 

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