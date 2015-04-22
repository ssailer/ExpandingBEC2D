
#include <EXP2D_splitstep.hpp>
#include <omp.h>

#define EIGEN_DONT_VECTORIZE
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_NO_DEBUG

// #define SLICE_NUMBER 0

using namespace std;
using namespace Eigen;

void SplitStep::cli(string name,int &slowestthread, vector<int> threadinfo, vector<int> stateOfLoops, int counter_max, double start)
{	
	for(int i = 0;i < stateOfLoops.size();i++){
		slowestthread = (stateOfLoops[slowestthread] <= stateOfLoops[i]) ? slowestthread : threadinfo[i];
	}
	
	if(fmod((float)stateOfLoops[slowestthread],(float)(counter_max/10))==0){
		int seconds, min, hour, total, expectedhour, expectedmin, expectedseconds;
		double totalstate = 0;
		double totalmaxpercent = (double)counter_max * (double)meta.samplesize / 100;
		for(int i = 0; i < meta.samplesize; i++){
			totalstate += stateOfLoops[i];
		}
		double totalPercent = totalstate/totalmaxpercent;

		int overallStepState = keeperOfTime.absoluteSteps + totalstate / meta.samplesize;

		total = omp_get_wtime() - start;

		overallStepState = (overallStepState == 0) ? 1 : overallStepState;

		int remainingSeconds = (total * opt.n_it_RTE / overallStepState) - total;


		
		hour = total / 3600;
		min = (total / 60) % 60;
		seconds = total % 60;
		expectedhour = (remainingSeconds / 3600);
		expectedmin = (remainingSeconds / 60) % 60;
		expectedseconds = remainingSeconds % 60;
		cout << "\r";
		cout << currentTime() <<  " " << name << " "
		 	 << std::setw(2) << std::setfill('0') << hour << ":"
			 << std::setw(2) << std::setfill('0') << min << ":"
			 << std::setw(2) << std::setfill('0') << seconds  << "    "
			 << std::setw(3) << std::setfill('0') << (int)totalPercent << "% "
			 << " threads: " << stateOfLoops.size()
			 << " remaining runtime: "
			 << std::setw(2) << std::setfill('0') << expectedhour << ":"
			 << std::setw(2) << std::setfill('0') << expectedmin << ":"
			 << std::setw(2) << std::setfill('0') << expectedseconds
			// cout << " | Slowest Thread: " << std::setw(3) << std::setfill('0') << (float)(stateOfLoops[slowestthread])/(float)(counter_max/10) << "% ";
			// for(int k = 0; k < stateOfLoops.size(); k++){
			// cout << k << "_" << threadinfo[k] << ": " << std::setw(3) << std::setfill('0') << (float)stateOfLoops[k]/((float)counter_max/100) << "% ";
		// }
		<< "    " << flush;
	}
}

void SplitStep::plot(const string name){

		plotDataToPng(name,"Control", wavefctVec[0],opt);

}

SplitStep::assignMatrixData(MatrixData* &d){
	w = d;
	setVariables();
}

// SplitStep::SplitStep(vector<ComplexGrid> &d,const MatrixData::MetaData &extMeta, const Options &externaloptions, int &extSLICE_NUMBER)
// {	
// 	setOptions(externaloptions);
// 	SLICE_NUMBER = extSLICE_NUMBER;
// 	meta = extMeta;
// 	// Both essential Variables
// 	wavefctVec.resize(opt.samplesize);
// 	for(int i = 0; i < opt.samplesize; i++){
// 		wavefctVec[i] = ComplexGrid(opt.grid[0], opt.grid[1], opt.grid[2], opt.grid[3]);
// 		wavefctVec[i] = d[i];
// 	}
	
  	

//   	// some constants used in computations to shorten stuff
// 	pi = M_PI;
//  	zero=complex<double>(0,0);
//  	half=complex<double>(0.5,0);
//  	one=complex<double>(1,0);
//  	two=complex<double>(2,0);
//  	four=complex<double>(4,0);
//  	six=complex<double>(6,0);
//  	i_unit=complex<double>(0,1);

//  	// setting up multithreading. Output to see what Eigen is doing.
// 	// omp_set_num_threads(omp_get_max_threads());
// 	cout << "Max Number of Threads in RTE: " << omp_get_max_threads() << endl;	
// 	// cout << "Eigenthreads: " << Eigen::nbThreads() << endl;

// 	// Using the setter function to initialize the stuff.
// 	RunSetup();

// }

void SplitStep::setOptions(const Options &externaloptions){
	opt = externaloptions;
}

void SplitStep::setVariables(){

		MatrixXcd kprop(w->meta.grid[0],w->meta.grid[1]);


	vector<vector<double>> kspace;	
	kspace.resize(3);
	for(int d = 0; d < 3; d++){
		kspace[d].resize(opt.grid[d+1]);
		for(int i = 0; i < opt.grid[d+1]/2; i++){
			kspace[d][i] = opt.klength[d]*2.0*sin( M_PI*((double)i)/((double)opt.grid[d+1]) );
			// kspace[d][i] = 2.0 * M_PI * i / (opt.grid[d+1]);
		}
		for(int i = opt.grid[d+1]/2; i < opt.grid[d+1]; i++){
			kspace[d][i] = opt.klength[d]*2.0*sin( M_PI*((double)(-opt.grid[d+1]+i))/((double)opt.grid[d+1]) );
			// kspace[d][i] = 2.0 * M_PI * (i - opt.grid[d+1]) / (opt.grid[d+1]);
		}
	}

	#pragma omp parallel for
	for(int x = 0; x < opt.grid[1]; x++){
	    for(int y = 0; y < opt.grid[2]; y++){

	      	double T = - 0.5 * (kspace[0][x]*kspace[0][x] + kspace[1][y]*kspace[1][y]) * opt.RTE_step; // / beta;	      
      		kprop(x,y) = complex<double>(cos(T),sin(T)) / complex<double>((double)(w->grid[0]*w->.grid[1]),0.0);	    
	    }
	}
}


void SplitStep::timeStep(double delta_t){


	w->fft_forward();
	w->data[0] *= kprop;
	w->fft_backward(); // YES?
	w->data[0] *= Vgrid;
		
				ComplexGrid::fft_unnormalized(wavefctVec[i], kgrid, true);
    
				#pragma omp parallel for
				for(int x = 0; x < opt.grid[1]; x++){
					for(int y = 0; y < opt.grid[2]; y++){
						for(int z = 0; z < opt.grid[3]; z++){				
				   			kgrid(0,x,y,z) = kprop(0,x,y,z) * kgrid(0,x,y,z);
				   		}
					}
				}

				ComplexGrid::fft_unnormalized(kgrid, wavefctVec[i], false);

				#pragma omp parallel for
				for(int x = 0; x < opt.grid[1]; x++){
					for(int y = 0; y < opt.grid[2]; y++){
						for(int z = 0; z < opt.grid[3]; z++){	
				    		// complex<double> value = wavefctVec[i](0,x,y,z);
				    		double V = - ( /*PotentialGrid(x,y).real()*/ /*rotatingPotential(x,y,m)*//* +*/ opt.g * abs2(wavefctVec[i](0,x,y,z)) ) * timestepsize;
				    		// potPlotGrid(0,x,y,0) = complex<double>(rotatingPotential(x,y,m) /*PotentialGrid(x,y).real()*/,0.0);
				    		// potGrid(0,x,y,0) = complex<double>(cos(V),sin(V));
				    		wavefctVec[i](0,x,y,z) = complex<double>(cos(V),sin(V)) * wavefctVec[i](0,x,y,z);
				    	}
					}
				}
					
			

}