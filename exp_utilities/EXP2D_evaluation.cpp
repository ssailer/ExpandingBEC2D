#include <EXP2D_evaluation.h>

#define OBSERVABLES_DATA_POINTS_SIZE opt.grid[1]*opt.grid[2]

using namespace std;
using namespace Eigen;


 Eval::Eval() {};
 Eval::~Eval() {};

void Eval::saveData(vector<MatrixXcd> &wavefctVec,Options &externalopt,int &external_snapshot_time,string runname_external){
		runname = runname_external;
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

void Eval::saveData(MatrixXcd &wavefct,Options &externalopt,int &external_snapshot_time,string runname_external){
		runname = runname_external;
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

	vortexLocationMap.resize(PsiVec.size());
	densityLocationMap.resize(PsiVec.size());
	vortexCoordinates.resize(PsiVec.size());
	densityCoordinates.resize(PsiVec.size());

	totalResult = Observables(OBSERVABLES_DATA_POINTS_SIZE);
		
	for(int k = 0; k < PsiVec.size(); k++){
		findVortices(PsiVec[k],vortexLocationMap[k],vortexCoordinates[k]);
		findDensity(PsiVec[k],densityLocationMap[k],densityCoordinates[k]);
		totalResult += calculator(PsiVec[k],k);

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
	string filename = runname + "-Spectrum-" + to_string(snapshot_time); 
	plotspectrum(filename,totalResult);
	filename = runname + "-Vortices-" + to_string(snapshot_time);
	plotVortexLocationMap(filename,vortexLocationMap[0]);
	filename = runname + "-Control-Plot-" + to_string(snapshot_time);
	plotdatatopng(filename,PsiVec[0],opt);
	filename = runname + "-Density-" + to_string(snapshot_time);
	plotdatatopng(filename,densityLocationMap[0],opt);
	filename = runname + "-Density-Axial-Distribution-" + to_string(snapshot_time);
	plotVector(filename,x_dist,y_dist,opt);
	filename = runname + "-Angular-Dens" + to_string(snapshot_time);
	plotVector(filename,totalResult.angularDensity,opt);
}

void Eval::findVortices(ComplexGrid data, RealGrid &vortexLocationMap_local, vector<Coordinate<int32_t>> &vortexCoordinates){
	
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

void Eval::findDensity(ComplexGrid data, RealGrid &densityLocationMap_local, vector<Coordinate<int32_t>> &densityCoordinates){

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
				densityCoordinates.push_back(data.make_coord(i,j,0));
			}else{
				densityLocationMap_local(0,i,j,0) = 0.;
			}
		}
	}


	// angularDensity = polarDensity;

	// if(angularDensity.size() != phi.size())
	// 	cout << "ERROR: Angular Density index problems." << endl;

	x_dist.resize(opt.grid[1]);
	y_dist.resize(opt.grid[2]);

	for(int i = 0; i < opt.grid[1]; i++){
		double sum = 0;
		for(int j = 0; j < opt.grid[2]; j++){
			sum += densityLocationMap_local(0,i,j,0);			
		}
		x_dist[i] = sum;
	}
	for(int j = 0; j < opt.grid[2]; j++){
		double sum = 0;
		for(int i = 0; i < opt.grid[1]; i++){
			sum += densityLocationMap_local(0,i,j,0);
		}
		y_dist[j] = sum;
	}

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


Observables Eval::calculator(ComplexGrid data,int sampleindex){
	
	Observables ares = Observables(OBSERVABLES_DATA_POINTS_SIZE);
	// R-Space
	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2];
	double h[2];
	h[0] = h_x;
	h[1] = h_y; 
	// double raw_volume = h_x * opt.grid[1] * h_y * opt.grid[2];
	
	// double threshold = abs2(data(0,opt.grid[1]/2,opt.grid[2]/2,0))*0.9;

	ares.volume = h_x * h_y * densityCoordinates[sampleindex].size();
	


	// == Angular Density
	// vector<double> angularDensity;
	vector<double> phi;
	vector<double> polarDensity; // first entry is 
	vector<double> radius;
	vector<Coordinate<int32_t>> cartesianCoordinates;
	// RealGrid conversionControl(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);


	for(int i = 0; i < opt.grid[1]; i++){
	    for(int j = 0; j < opt.grid[2]; j++){
				int x_shift = i - opt.grid[1]/2;
				int y_shift = j - opt.grid[2]/2;
				phi.push_back( atan2(x_shift * h_x ,y_shift * h_y) * 180 / M_PI + 180);
				radius.push_back(sqrt(x_shift*x_shift * h_x*h_x + y_shift*y_shift *h_y*h_y));
				polarDensity.push_back(abs2(data(0,i,j,0)));
				cartesianCoordinates.push_back(data.make_coord(i,j,0));
		}
	}
	
	double max_radius = (opt.grid[1] + opt.grid[2]) / 4;

	vector<double>angularDensity_tmp(361);
	for(int i = 0; i <= 360; i++){
		for(int j = 0; j < phi.size(); j++){
			if(round(phi[j]) == i){
				if(radius[j] <= max_radius){
					angularDensity_tmp[i] += polarDensity[j];
				}
			}
		}
	}
	angularDensity_tmp[0] += angularDensity_tmp[360];
	angularDensity_tmp.pop_back();

	for(int i = 0; i < 360; i++){
		// vector<double>::iterator beginning = angularDensity_tmp.begin() + i;
		// vector<double>::iterator ending = beginning + 10;
		ares.angularDensity(i) = angularDensity_tmp[i]; //(accumulate(beginning,ending,0));
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
		// for (int32_t i = 0; i < kspace[d].size()/2; i++){
			// kspace[d][i] = opt.klength[d]*sin( M_PI*((double)i)/((double)opt.grid[d+1]) );
			// kspace[d][i] = opt.klength[d]*((double)i)/((double)(opt.grid[d+1]/2));
			kspace[d][i] = opt.klength[d] * 2 * M_PI  * ((double)i) / ((double)(opt.grid[d+1]*opt.grid[d+1]*h[d]));
		}
		for (int i=opt.grid[d+1]/2; i<opt.grid[d+1]; i++){
		// for (int32_t i = kspace[d].size()/2; i < kspace[d].size(); i++){
			// kspace[d][i] = opt.klength[d]*sin( M_PI*((double)(-opt.grid[d+1]+i))/((double)opt.grid[d+1]) );
			// kspace[d][i] = opt.klength[d]*((double)(opt.grid[d+1]-i))/((double)opt.grid[d+1]/2);
			kspace[d][i] = opt.klength[d] * 2 * M_PI  * ((double)(-opt.grid[d+1]+i)) / ((double)(opt.grid[d+1]*opt.grid[d+1]*h[d]));
		}
	}


	// DEFINITION OF KLENGTH FROM THE INTERNET! |||| ==>>     2*pi*i/(Nx*dx)
	// double kmax[2];

	// for(int i = 0; i < 2; i++){
	// 	kmax[i] = *max_element(kspace[i].begin(), kspace[i].end());
	// }
	
	double kwidth2[2];
	
	for(int i = 0; i < 2; i++)
		kwidth2[i] = (opt.grid[i+1] == 1) ? 0 : kspace[i][opt.grid[i+1]/2] * kspace[i][opt.grid[i+1]/2];
	
	double index_factor = (ares.number.size() - 1) / sqrt(kwidth2[0] + kwidth2[1]);

	for(int x = 0; x < data.width(); x++){
		for (int y = 0; y < data.height(); y++){
			for (int z = 0; z < data.depth(); z++){
				double k = sqrt(kspace[0][x]*kspace[0][x] + kspace[1][y]*kspace[1][y]);
				// Coordinate<int32_t> c = data.make_coord(x,y,z);
				int index = index_factor * k;
				// cout << k << "*" << index_factor << "=" << index << "/" << OBSERVABLES_DATA_POINTS_SIZE << endl;
				ares.k(index) += k;
				divisor(index)++;
				double number = abs2(data(0,x,y,z));
				ares.number(index) += number;
				ares.particle_count += number;
				ares.Ekin += number * k * k;
			}
		}
	}

	ares.density = ares.particle_count / ares.volume;

	// double density_integrated = 0;
	// for(int i = 0; i < opt.grid[1]-1; i++){
	//     for(int j = 0; j < opt.grid[2]-1; j++){	    	    		
	//       	density_integrated += h_x * h_y * (abs2(data(0,i,j,0))+abs2(data(0,i+1,j,0))+abs2(data(0,i,j+1,0))+abs2(data(0,i+1,j+1,0)))/4.0;
	//     }
	// }

	// cout << "Testmethods: " << density_integrated << "  " << ares.density << endl;

	
	// ares.healing_length = 1 / sqrt(ares.particle_count * opt.g / ares.volume);
	
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
	
	