#include <EXP2D_evaluation.h>

using namespace std;
using namespace Eigen;


 Eval::Eval() {};
 Eval::~Eval() {};

void Eval::saveData(vector<MatrixXcd> &wavefctVec,Options &externalopt,int &external_snapshot_time){
		opt = externalopt;
		snapshot_time = external_snapshot_time;
		PsiVec.resize(wavefctVec.size());

		#pragma parallel for
		for(int k = 0; k < wavefctVec.size(); k++){
			PsiVec[k] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

			for(int i = 0; i < opt.grid[1]; i++){for(int j = 0; j < opt.grid[2]; j++){		
			PsiVec[k](0,i,j,0) = wavefctVec[k](i,j);}}
		}
}

void Eval::saveData(MatrixXcd &wavefct,Options &externalopt,int &external_snapshot_time){
		opt = externalopt;
		snapshot_time = external_snapshot_time;
		PsiVec.resize(1);

		PsiVec[0] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

		for(int i = 0; i < opt.grid[1]; i++){
			for(int j = 0; j < opt.grid[2]; j++){		
				PsiVec[0](0,i,j,0) = wavefct(i,j);
			}
		}		
}

void Eval::evaluateData(){
	vector<Observables> avResult(PsiVec.size());
	vortexLocationMap.resize(PsiVec.size());
	densityLocationMap.resize(PsiVec.size());
	vortexCoordinates.resize(PsiVec.size());
	densityCoordinates.resize(PsiVec.size());

	totalResult = Observables(3*opt.grid[1]);
		
	for(int k = 0; k < PsiVec.size(); k++){
		// avResult[k] = Observables(3*opt.grid[1]);
		// avResult[k] = calculator(PsiVec[k]);
		totalResult += calculator(PsiVec[k]);
		// avResult + evaluate(PsiVec[k]);
		findVortices(PsiVec[k],vortexLocationMap[k],vortexCoordinates[k]);
		findDensity(PsiVec[k],densityLocationMap[k],densityCoordinates[k]);

		for(int i = 0; i < opt.grid[1]; i++){
			for(int j = 0; j < opt.grid[2]; j++){
				vortexLocationMap[k](0,i,j,0) *= densityLocationMap[k](0,i,j,0);
			}
		}
	}




	
	// for(int k = 0; k < PsiVec.size(); k++){
	// 	 totalResult = avResult[k];
	// }

	totalResult /= PsiVec.size();
	

}

void Eval::plotData(){
	string filename = "Spectrum-" + to_string(snapshot_time); 
	plotspectrum(filename,totalResult);
	filename = "Vortices-" + to_string(snapshot_time);
	plotVortexLocationMap(filename,vortexLocationMap[0]);
	filename = "Control-Plot-" + to_string(snapshot_time);
	plotdatatopng(filename,PsiVec[0],opt);
	filename = "Density-" + to_string(snapshot_time);
	plotdatatopng(filename,densityLocationMap[0],opt);
	filename = "Density-Distribution-" + to_string(snapshot_time);
	plotVector(filename,x_dist,y_dist,opt);
	filename = "Angular-Dens" + to_string(snapshot_time);
	plotVector(filename,angularDensity,opt);
}

void Eval::findVortices(ComplexGrid &data, RealGrid &vortexLocationMap_local, vector<Coordinate<int32_t>> &vortexCoordinates){
	
	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2]; 
	
	RealGrid phase_grad_x(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	RealGrid phase_grad_y(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	vortexLocationMap_local = RealGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

	for(int x = 1; x < data.width() - 1; x++){
		for (int y = 1; y < data.height() - 1; y++){
			phase_grad_x(0,x,y,0) = ( arg(data(0,x+1,y,0))- arg(data(0,x-1,y,0)) ) / (2.0 * h_x );
			phase_grad_y(0,x,y,0) = ( arg(data(0,x,y+1,0))- arg(data(0,x,y+1,0)) ) / (2.0 * h_y );
		}
	}
	for(int x = 0; x < data.width(); x++){
		for (int y = 0; y < data.height(); y++){
			vortexLocationMap_local(0,x,y,0) = 0.0;
		}
	}
	for(int x = 2; x < data.width() - 2; x++){
		for (int y = 2; y < data.height() - 2; y++){
			vortexLocationMap_local(0,x,y,0) = ( phase_grad_y(0,x+1,y,0) - phase_grad_y(0,x-1,y,0) ) / (2.0 * h_x ) - ( phase_grad_x(0,x,y+1,0) - phase_grad_x(0,x,y-1,0) ) / (2.0 * h_y);
		}
	}

}

void Eval::findDensity(ComplexGrid &data, RealGrid &densityLocationMap_local, vector<Coordinate<int32_t>> &densityCoordinates){

	double lower_threshold = 10.; //abs2(data(0,opt.grid[1]/2,opt.grid[2]/2,0))*0.9;
	double upper_threshold = 20.;

	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2]; 	

	// RealGrid density_grad_x(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	// RealGrid density_grad_y(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);





	densityLocationMap_local = RealGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

	// for(int x = 1; x < data.width() - 1; x++){
	// 	for (int y = 1; y < data.height() - 1; y++){
	// 		density_grad_x(0,x,y,0) = ( arg(data(0,x+1,y,0))- arg(data(0,x-1,y,0)) ) / (2.0 * h_x );
	// 		density_grad_y(0,x,y,0) = ( arg(data(0,x,y+1,0))- arg(data(0,x,y+1,0)) ) / (2.0 * h_y );
	// 	}
	// }



	for(int i = 0; i < opt.grid[1]; i++){
	    for(int j = 0; j < opt.grid[2]; j++){
	    	if((abs2(data(0,i,j,0)) > lower_threshold)){
				densityLocationMap_local(0,i,j,0) = 1.;
				// densityCoordinates.push_back(data.make_coord(i,j,0));
			}else{
				densityLocationMap_local(0,i,j,0) = 0.;
			}
		}
	}

	vector<double> polarDensity; // first entry is 
	vector<int> phi;
	vector<double> radius;
	vector<Coordinate<int32_t>> cartesianCoordinates;
	// RealGrid conversionControl(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);


	for(int i = 0; i < opt.grid[1]; i++){
	    for(int j = 0; j < opt.grid[2]; j++){
				int x_shift = i - opt.grid[1]/2;
				int y_shift = j - opt.grid[2]/2;
				phi.push_back( atan2(x_shift,y_shift) * 360 / M_PI );
				radius.push_back(sqrt(x_shift*x_shift + y_shift*y_shift));
				polarDensity.push_back(abs2(data(0,i,j,0)));
				cartesianCoordinates.push_back(data.make_coord(i,j,0));
		}
	}
	
	int max_radius = (opt.grid[1] + opt.grid[2]) / 4;


	angularDensity.resize(36);
	vector<double>angularDensity_tmp(360);
	for(int i = 0; i < 360; i++){
		for(int j = 0; j < phi.size(); j++){
			if(phi[j] == i){
				if(radius[j] <= max_radius){
					angularDensity_tmp[i] += polarDensity[j];
				}
			}
		}
	}
	for(int i = 1; i < 11; i++){
		for(int j = 1; j < 37; j ++){
			angularDensity[j-1] +=angularDensity_tmp[i*j - 1];
		}
	}





	// x_dist.resize(opt.grid[1]);
	// y_dist.resize(opt.grid[2]);

	// for(int i = 0; i < opt.grid[1]; i++){
	// 	double sum = 0;
	// 	for(int j = 0; j < opt.grid[2]; j++){
	// 		sum += densityLocationMap_local(0,i,j,0);			
	// 	}
	// 	x_dist[i] = sum;
	// }
	// for(int j = 0; j < opt.grid[2]; j++){
	// 	double sum = 0;
	// 	for(int i = 0; i < opt.grid[1]; i++){
	// 		sum += densityLocationMap_local(0,i,j,0);
	// 	}
	// 	y_dist[j] = sum;
	// }











	double sum = 0;
	for(int k = 0; k < 10; k++){
		for(int x = 1; x < opt.grid[1]-1; x+=1){
			for(int y = 1; y < opt.grid[2]-1; y+=1){				
				for(int i = x-1; i <= x+1; i++){
					for(int j = y-1; j <= y+1; j++){
						sum += densityLocationMap_local(0,i,j,0);
					}
				}
				if((sum = 8) && (densityLocationMap_local(0,x,y,0) == 0)){ // now it is surround by stuff, and instelf zero, so we assume this is a vortex | this is a good place for a counter of vortices
						densityLocationMap_local(0,x,y,0) = 0; // 
					}
				else if(sum >= 5){ // Point is either half surrounded by stuff, or is itself stuff, so assume it is density, which we didn't catch before
					densityLocationMap_local(0,x,y,0) = 1.;
					// densityCoordinates.push_back(data.make_coord(i,j,0));
				}else{
					densityLocationMap_local(0,x,y,0) = 0.; //
				}
				if(sum > 9){cout << "ERROR: TO MUCH SUM" << endl;}
				sum = 0;
			}
		}
	}



}


Observables Eval::calculator(ComplexGrid data){
	
	Observables ares = Observables(3*opt.grid[1]);
	// R-Space
	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2]; 
	double raw_volume = h_x * opt.grid[1] * h_y * opt.grid[2];
	
	double threshold = abs2(data(0,opt.grid[1]/2,opt.grid[2]/2,0))*0.9;	
	
	for(int i = 0; i < opt.grid[1]-1; i++){
	    for(int j = 0; j < opt.grid[2]-1; j++){	    	    		
	      	ares.volume += h_x * h_y * (abs2(data(0,i,j,0))+abs2(data(0,i+1,j,0))+abs2(data(0,i,j+1,0))+abs2(data(0,i+1,j+1,0)))/4.0;
	    }
	}
	
	// K-Space
	ComplexGrid::fft(data, data);
	
	ArrayXd divisor(ares.number.size());
	divisor.setZero();
	
	vector<vector<double>> kspace;
	
	kspace.resize(2);
	for(int d = 0; d < 2; d++){
		// set k-space
		kspace[d].resize(opt.grid[d+1]);
		for (int i=0; i<opt.grid[d+1]/2; i++){
			// kspace[d][i] = opt.klength[d]*sin( M_PI*((double)i)/((double)opt.grid[d+1]) );
			kspace[d][i] = opt.klength[d]*((double)i)/((double)(opt.grid[d+1]/2));
		}
		for (int i=opt.grid[d+1]/2; i<opt.grid[d+1]; i++){
			// kspace[d][i] = opt.klength[d]*sin( M_PI*((double)(-opt.grid[d+1]+i))/((double)opt.grid[d+1]) );
			kspace[d][i] = opt.klength[d]*((double)(opt.grid[d+1]-i))/((double)opt.grid[d+1]/2);
		}
	}
	
	double kwidth2[3];
	
	for(int i = 0; i < 3; i++)
		kwidth2[i] = (opt.grid[i+1] == 1) ? 0 : opt.klength[i]*opt.klength[i];
	
	double index_factor = (ares.number.size() - 1) / sqrt(kwidth2[0] + kwidth2[1] + kwidth2[2]);
	
	for(int x = 0; x < data.width(); x++){
		for (int y = 0; y < data.height(); y++){
			for (int z = 0; z < data.depth(); z++){
				double k = sqrt(kspace[0][x]*kspace[0][x] + kspace[1][y]*kspace[1][y]);
				Coordinate<int32_t> c = data.make_coord(x,y,z);
				int index = index_factor * k;
				ares.k(index) += k;
				divisor(index)++;
				double number = abs2(data(0,c));
				ares.number(index) += number;
				ares.particle_count += number;
				ares.Ekin += number * k * k;
			}
		}
	}

	ares.particle_count /= raw_volume;

	
	ares.healing_length = 1 / sqrt(ares.particle_count * opt.g / ares.volume);
	
	#pragma omp parallel for schedule(guided,1)
	for(int l = 0; l < ares.number.size(); l++){
		if(divisor[l] == 0){
			divisor[l] = 1;
		}
	}

	ares.number /= divisor;
	ares.k /= divisor;
	
	
	return ares;
}
	
	