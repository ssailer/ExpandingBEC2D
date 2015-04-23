
#include <EXP2D_splitstep.hpp>
#include <omp.h>

#define EIGEN_DONT_VECTORIZE
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_NO_DEBUG

// #define SLICE_NUMBER 0

using namespace std;
using namespace Eigen;

void SplitStep::assignMatrixData(MatrixData* &d){
	w = d;
	setVariables();
}

void SplitStep::setOptions(const Options &externaloptions){
	opt = externaloptions;
}

void SplitStep::setVariables(){

	kprop = MatrixXcd(w->meta.grid[0],w->meta.grid[1]);
	
	for(int d = 0; d < 2; d++){

		kspace[d].resize(w->meta.grid[d]);
		for(int i = 0; i <= w->meta.grid[d]/2; i++){
			kspace[d][i] = (M_PI / w->meta.coord[d]) * (double)i;
		}

		for(int i = (w->meta.grid[d]/2)+1; i < w->meta.grid[d]; i++){
			kspace[d][i] = -(M_PI / w->meta.coord[d]) * (double)(w->meta.grid[d] - i);
		}
	}

	#pragma omp parallel for
	for(int x = 0; x < w->meta.grid[0]; x++){
	    for(int y = 0; y < w->meta.grid[1]; y++){
	      	double T = - 0.5 * (kspace[0][x]*kspace[0][x] + kspace[1][y]*kspace[1][y]) * opt.RTE_step;
      		kprop(x,y) = complex<double>(cos(T),sin(T)) / complex<double>((double)(w->meta.grid[0]*w->meta.grid[1]),0.0);	    
	    }
	}
}


void SplitStep::timeStep(double delta_t){


	w->fftForward();
	w->wavefunction[0].array() *= kprop.array();
	w->fftBackward(); // YES?
	Vgrid.array() = exp( - complex<double>(opt.g,0.0) * w->wavefunction[0].conjugate().array() * w->wavefunction[0].array() * delta_t);
	w->wavefunction[0].array() *= Vgrid.array();
		
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