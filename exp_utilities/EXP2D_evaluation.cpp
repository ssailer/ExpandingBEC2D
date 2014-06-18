#include <EXP2D_evaluation.h>

#define OBSERVABLES_DATA_POINTS_SIZE opt.grid[1]*opt.grid[2]
#define ANGULAR_AVERAGING_LENGTH 12

using namespace std;
using namespace Eigen;


Eval::Eval() {



};

Eval::~Eval() {};

void Eval::saveData(vector<MatrixXcd> &wavefctVec,Options &externalopt,int &external_snapshot_time,string runname_external){
		runname = runname_external;
		opt = externalopt;
		snapshot_time = external_snapshot_time;
		PsiVec.resize(wavefctVec.size());
		v_down = Vector<int32_t>(0,-1,0,opt.grid[1],opt.grid[2],opt.grid[3]);
		v_right = Vector<int32_t>(1,0,0,opt.grid[1],opt.grid[2],opt.grid[3]);
		v_up = Vector<int32_t>(0,1,0,opt.grid[1],opt.grid[2],opt.grid[3]);
		v_left = Vector<int32_t>(-1,0,0,opt.grid[1],opt.grid[2],opt.grid[3]);


		#pragma omp parallel for
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

	pres.vlist.clear();
	densityCoordinates.clear();

	// vortexLocationMap.resize(PsiVec.size());
	densityLocationMap.resize(PsiVec.size());
	// vortexCoordinates.resize(PsiVec.size());
	densityCoordinates.resize(PsiVec.size());
	phase = new RealGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	zeros = new RealGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);


	totalResult = Observables(OBSERVABLES_DATA_POINTS_SIZE);
		
	for(int k = 0; k < PsiVec.size(); k++){
		getDensity(PsiVec[k],densityLocationMap[k],densityCoordinates[k]);		
		totalResult += calculator(PsiVec[k],k);
	}
	getVortices(PsiVec[0],densityCoordinates[0]);


	
	// for(int k = 0; k < PsiVec.size(); k++){
	// 	 totalResult = avResult[k];
	// }

	totalResult /= PsiVec.size();
	

}

void Eval::plotData(){
	string filename = runname + "-Spectrum-" + to_string(snapshot_time); 
	plotspectrum(filename,totalResult);
	filename = runname + "-Vortices-" + to_string(snapshot_time);
	plotVortexList(filename,phase,pres,opt);	
	filename = runname + "-Control-Plot-" + to_string(snapshot_time);
	plotdatatopng(filename,PsiVec[0],opt);
	filename = runname + "-Density-" + to_string(snapshot_time);
	plotdatatopng(filename,densityLocationMap[0],opt);
	filename = runname + "-Density-Axial-Distribution-Gradient" + to_string(snapshot_time);
	plotVector(filename,x_dist_grad,y_dist_grad,opt);
	filename = runname + "-Angular-Dens" + to_string(snapshot_time);
	plotVector(filename,totalResult.angularDensity,opt);

	cout << "HERE COMES THE KILLED" << endl;
	// EXPERIMENTAL
	std::unordered_set<Coordinate<int32_t>> contour = trackContour(PsiVec[0],opt);
	filename = runname + "-Contour" + to_string(snapshot_time);
	plotContour(filename,phase,contour,opt);
	// EXPERIMENTAL END

	
	filename = runname + "-Observables" + ".dat";
	struct stat buffer;   
  	if(stat (filename.c_str(), &buffer) != 0){
  		ofstream datafile;
  		datafile.open(filename.c_str(), ios::out | ios::app);
  		datafile << std::left << std::setw(10) << "Timestep"
  						 << std::setw(10) << "X_max"
  						 << std::setw(10) << "Y_max"
  						 << std::setw(10) << "N"
  						 << std::setw(10) << "V"
  						 << std::setw(10) << "N/V"
  						 << std::setw(10) << "E_kin"
  				 << endl;
  		datafile.close();
  	} 

  	ofstream datafile(filename.c_str(), std::ios_base::out | std::ios_base::app);
	// datafile.open;
	datafile << std::left << std::setw(10) << snapshot_time
					 << std::setw(10) << opt.min_x * opt.stateInformation[0]
					 << std::setw(10) << opt.min_y * opt.stateInformation[1]
					 << std::setw(10) << totalResult.particle_count
					 << std::setw(10) << totalResult.volume
					 << std::setw(10) << totalResult.density
					 << std::setw(10) << totalResult.Ekin 
			 << endl;
	datafile.close();
}

void Eval::getVortices(const ComplexGrid &data, vector<Coordinate<int32_t>> &densityCoordinates){
	
	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2]; 
	
	// vortexLocationMap_local = RealGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);


	calc_fields(data,opt);
	pres.vlist.clear();
	find_vortices(phase,zeros,densityCoordinates,pres.vlist);

	cout << endl << "Vortices: " << endl;
	double number = 0;
	for(list<VortexData>::const_iterator it = pres.vlist.begin(); it != pres.vlist.end(); ++it){
		int x = it->x.x();
		int y = it->x.y();
		number += it->num_points;
		cout << " " << x << " " << y << "  " << abs2(PsiVec[0](0,x,y,0)) << " " << arg(PsiVec[0](0,x,y,0)) << endl;
	}
	cout << "Number of Vortices counted: " << number << "  " << endl;

	// calc_vortex_veloctities();
	// calc_vortex_discances();
	// calc_g2();

}

int Eval::get_phase_jump(const Coordinate<int32_t> &c, const Vector<int32_t> &v, const RealGrid *phase) 
{
	if(phase->at(0,c + v) + M_PI < phase->at(0,c))	// Phase ueberschreitet 0/2pi von unten
		return 1;
	else if(phase->at(0,c) + M_PI < phase->at(0,c + v))	// Phase ueberschreitet 0/2pi von oben
		return -1;
	else
		return 0;
}

void Eval::find_vortices(const RealGrid *phase, const RealGrid *zeros, vector<Coordinate<int32_t>> &densityCoordinates, list<VortexData> &vlist) 
{
	// Nullstellen zaehlen
	// vector< vector< vector<bool > > > checked(phase->width(), vector< vector<bool> >(phase->height(), vector<bool>(phase->depth(),false)));	// Welche felder schon ueberprueft wurden
	VortexData vortex;								// Charakteristika eines gefundenen Vortex

	for(vector<Coordinate<int32_t>>::const_iterator it = densityCoordinates.begin(); it != densityCoordinates.end(); ++it){

	// for (int z = 0; z < phase->depth(); z++)
	// {
	// 	for (int x = 0; x < phase->width(); x++)
	// 	{
	// 		for (int y = 0; y < phase->height(); y++)
	// 		{

				Coordinate<int32_t> c = *it; // phase->make_coord(x,y,z);
				Vector<int32_t> down = phase->make_vector(0,-1,0);
				Vector<int32_t> right = phase->make_vector(1,0,0);
				Vector<int32_t> up = phase->make_vector(0,1,0);
				Vector<int32_t> left = phase->make_vector(-1,0,0);
					
				int phase_winding = get_phase_jump(c, down, phase) + get_phase_jump(c+down, right, phase) + get_phase_jump(c+down+right, up, phase) + get_phase_jump(c+right, left, phase);
				// int mass_zeros = zeros->at(0,c) + zeros->at(0,c+down) + zeros->at(0,c+down+right) + zeros->at(0,c+right);
				
				// if(mass_zeros < 4 && mass_zeros >= 0)
				if(phase_winding != 0)
				{
					vortex.n = phase_winding;
					vortex.x = c + phase->make_vector(0.5, -0.5, 0);
					// cout << vortex.x << endl;
					vortex.points.clear();
					vortex.points.push_back(c);
					vortex.num_points = 1;
					vlist.push_back(vortex);					
				}
				/*if(get_vortex(phase->make_coord(x,y,z), phase, zeros, mass_zeros, checked, vortex)) // prueft auf Vortices
				{
					vlist.push_back(vortex);					// und er wird in der Liste zurueckgegeben
				}*/
	// 		}
	// 	}
	}
}




void Eval::calc_fields(const ComplexGrid &data, Options &opt)
{
	double zero_threshold = 0;//opt.N * 0.05 / (4. * opt.min_x * opt.stateInformation[0] * opt.min_y * opt.stateInformation[1]); ; //opt.N * 0.05 / data.width() / data.height() / data.depth();
	for(int x = 0; x < data.width(); x++)
	{
		for(int y = 0; y < data.height(); y++)
		{
			for(int z = 0; z < data.depth(); z++)
			{
				phase->at(0,x,y,z) = arg(data(0,x,y,z));
				if(abs2(data(0,x,y,z)) < zero_threshold)
					zeros->at(0,x,y,z) = 0.0;
				else
					zeros->at(0,x,y,z) = 1.0;
			}
		}
	}


}

inline void Eval::setDirection(int32_t &direction){
		switch(direction){
			case 0:
				direction = 6;
				break;
			case 1:
				direction = 0;
				break;
			case 2:
				direction = 0;
				break;
			case 3:
				direction = 2;
				break;
			case 4:
				direction = 2;
				break;
			case 5:
				direction = 4;
				break;
			case 6:
				direction = 4;
				break;
			case 7:
				direction = 6;
				break;
			default:
				cout << "Function setDirection() ERROR!" << endl;
				break;
		}
};

inline Coordinate<int32_t> Eval::nextClockwise(Coordinate<int32_t> &s, int32_t &direction){

		Coordinate<int32_t> c;

		switch(direction){
			case 0 :
				c = s + v_up;
				break;
			case 1 :
				c = s + v_right;
				break;
			case 2 :
				c = s + v_right;
				break;
			case 3 :
				c = s + v_down;
				break;
			case 4 :
				c = s + v_down;
				break;
			case 5 :
				c = s + v_left;
				break;
			case 6 :
				c = s + v_left;
				break;
			case 7 :
				c = s + v_up;
				break;
			default :
				cout << "Function nextClockwise() ERROR!" << endl;
		}
		return c; 

	};

void Eval::findInitialP(Coordinate<int32_t> &p,Coordinate<int32_t> &s, Coordinate<int32_t> *initial,const Options &opt){
	for(int x = p.x(); x < opt.grid[1]; x++){
		if(zeros->at(0,x,opt.grid[2]/2,0) > 0.0){
			p = zeros->make_coord(x,opt.grid[2]/2,0);
			// contour.push_back(p);
			s = p + v_left;
			initial[0] = p;
			initial[1] = s;
		}
	}

}

std::unordered_set<Coordinate<int32_t>> Eval::trackContour(const ComplexGrid &data, const Options &opt){

	std::unordered_set<Coordinate<int32_t>> contour;
	std::unordered_set<Coordinate<int32_t>> it;
	std::pair<std::unordered_set<Coordinate<int32_t>>::iterator,bool> ret;

	// IMPLEMENT OWN HASH FUNCTION FOR unordered_set<Coordinate<int32_t>>
	
	Coordinate<int32_t> s;
	Coordinate<int32_t> p = zeros->make_coord(0,opt.grid[2]/2,0);
	Coordinate<int32_t> initial[2];

	findInitialP(p,s,initial,opt);
	
	int direction = 0;
	int insert_counter = 0;
	bool singlepoint = false;
	bool stop = false;
	cout << endl;
	do{
		int32_t counter = 0;
		while(counter < 8){
			Coordinate<int32_t> c = nextClockwise(s,direction);
			if(zeros->at(0,c) > 0){
				p = c;
				setDirection(direction);
				ret = contour.insert(p);
				insert_counter++;
				cout << "Added Coordinate: " << p << endl;
				break;
			}
			s = c;
			counter++;
			if(counter == 8){
				singlepoint = true;
				break;
			}
		}

		if(singlepoint == true){
			cout << "Found single point, continuing the search." << endl;
			contour.clear();
			Coordinate<int32_t> p_old = p;
			p = p + v_right;
			findInitialP(p,p_old,initial,opt);
			break;
		}

		if(insert_counter >= 2 * contour.size()){
			cout << "Surrounded the contour two times. Stopping now." << endl;
			stop = true;
			break;
		}
		
		if((initial[0] != p) && (initial[1] != s)){
			cout << "Found initial conditions, seems like a good run." << endl;
			stop = true;(initial[0] != p) && (initial[1] != s)
			break;
		}

	}while(stop == false);
	cout << endl;

}


void Eval::getDensity(const ComplexGrid &data, RealGrid &densityLocationMap_local, vector<Coordinate<int32_t>> &densityCoordinates_local){

	double lower_threshold =  opt.N * 0.05 / (4. * opt.min_x * opt.stateInformation[0] * opt.min_y * opt.stateInformation[1]);  //abs2(data(0,opt.grid[1]/2,opt.grid[2]/2,0))*0.9;
	// double upper_threshold = 20.;
	// cout << lower_threshold << "+++" << endl;

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
	// }gi
	densityCounter = 0;
	densityCoordinates_local.clear();
	for(int i = 0; i < opt.grid[1]; i++){
	    for(int j = 0; j < opt.grid[2]; j++){
	    	if((abs2(data(0,i,j,0)) > lower_threshold)){
				densityLocationMap_local(0,i,j,0) = 1.;
				densityCoordinates_local.push_back(data.make_coord(i,j,0));
				densityCounter++;
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
		for(int j = 0; j < opt.grid[2]; j++){
			x_dist[i] += densityLocationMap_local(0,i,j,0);			
		}
	}
	for(int j = 0; j < opt.grid[2]; j++){
		for(int i = 0; i < opt.grid[1]; i++){
			y_dist[j] += densityLocationMap_local(0,i,j,0);
		}
	}

	x_dist_grad.resize(opt.grid[1]);
	y_dist_grad.resize(opt.grid[2]);
	for(int x = 1; x < x_dist.size() - 1; x++){
		x_dist_grad[x] = (x_dist[x+1] - x_dist[x-1] ) / (2.0 * h_x );
	}
	for(int y = 1; y < y_dist.size() - 1; y++){
		y_dist_grad[y] = (y_dist[y+1] - y_dist[y-1] ) / (2.0 * h_y );
	}
	x_dist_grad[0] = x_dist_grad[opt.grid[1]-1] = y_dist_grad[0] = y_dist_grad[opt.grid[2]] = 0.0;

	// double sum = 0;
	// for(int k = 0; k < 10; k++){
	// 	for(int x = 1; x < opt.grid[1]-1; x+=1){
	// 		for(int y = 1; y < opt.grid[2]-1; y+=1){				
	// 			for(int i = x-1; i <= x+1; i++){
	// 				for(int j = y-1; j <= y+1; j++){
	// 					sum += densityLocationMap_local(0,i,j,0);
	// 				}
	// 			}
	// 			if((sum = 8) && (densityLocationMap_local(0,x,y,0) == 0)){ // now it is surround by stuff, and instelf zero, so we assume this is a vortex | this is a good place for a counter of vortices
	// 					densityLocationMap_local(0,x,y,0) = 0; //
	// 					// IMPORTANT PART HERE: I assume, I found a vortex, so I'm adding it to the densityCoordinates_local, because I want the phase of this point checked in getVortices and find_Vortices!
	// 					// densityCoordinates_local.push_back(data.make_coord(x,y,0));
	// 				}
	// 			else if(sum >= 5){ // Point is either half surrounded by stuff, or is itself stuff, so assume it is density, which we didn't catch before
	// 				densityLocationMap_local(0,x,y,0) = 1.;
	// 				// densityCoordinates_local.push_back(data.make_coord(i,j,0));
	// 			}else{
	// 				densityLocationMap_local(0,x,y,0) = 0.; //
	// 			}
	// 			if(sum > 9){cout << "ERROR: TO MUCH SUM" << endl;}
	// 			sum = 0;
	// 		}
	// 	}
	// }



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

	ares.volume = h_x * h_y * densityCounter;
	for(int i = 0; i < opt.grid[1]-1; i++){
	    for(int j = 0; j < opt.grid[2]-1; j++){	    	    		
	      	ares.particle_count += h_x * h_y * abs2(data(0,i,j,0));
	    }
	}

	ares.density = ares.particle_count / ares.volume;

	// cout << "Testmethods: " << ares.volume << "  " << ares.particle_count << "  " << ares.density << endl;
	


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
	
	// double max_radius = (opt.grid[1] + opt.grid[2]) / 4;

	// vector<double>angularDensity_tmp(361);
	vector<vector<int>> index(361);
	int phiSize = phi.size();
	for(int j = 0; j < phiSize; j++){
		int i = round(phi[j]);
		index[i].push_back(j);
	}
	for(int j = 0; j < index[360].size(); j++){
		index[0].push_back(index[360][j]);
	}
	index.pop_back();
	for(int i = 0; i < 360; i++){
		double sum = 0;
		for(int k = -ANGULAR_AVERAGING_LENGTH; k <= ANGULAR_AVERAGING_LENGTH; k++){
			int l = i + k;
			l = (l < 0)        ? l + 360 : (l >= 360) ? l - 360 : l;			
			for(int m = 0; m < index[l].size(); m++){
				ares.angularDensity(i) += polarDensity[index[l][m]];
			}
		}					
	}
	ares.angularDensity /= (ANGULAR_AVERAGING_LENGTH*2 + 1);

	// angularDensity_tmp[0] += angularDensity_tmp[360];	

	// for(int i = 0; i < 360; i++){
	// 	// vector<double>::iterator beginning = angularDensity_tmp.begin() + i;
	// 	// vector<double>::iterator ending = beginning + 10;
	// 	ares.angularDensity(i) = angularDensity_tmp[i]; //(accumulate(beginning,ending,0));
	// }


	
	// K-Space
	ComplexGrid::fft(data, data);
	
	ArrayXd divisor(ares.number.size());
	divisor.setZero();
	
	vector<vector<double>> kspace;
	
	kspace.resize(2);
	for(int d = 0; d < 2; d++){
		// set k-space
		kspace[d].resize(opt.grid[d+1]);
		// for (int i=0; i<opt.grid[d+1]/2; i++){
		for(int i = 0; i < opt.grid[d+1]/2; i++){
		// for (int32_t i = 0; i < kspace[d].size()/2; i++){
			// kspace[d][i] = opt.klength[d]*sin( M_PI*((double)i)/((double)opt.grid[d+1]) );
			// kspace[d][i] = opt.klength[d]*((double)i)/((double)(opt.grid[d+1]/2));
			// kspace[d][i] = opt.klength[d] * 2 * M_PI  * ((double)i) / ((double)(opt.grid[d+1]*opt.grid[d+1]*h[d]));
			kspace[d][i] = opt.klength[d] * M_PI / ( (double)(opt.grid[d+1]/2 - i) * h[d] );
		}
		// for (int i=opt.grid[d+1]/2; i<opt.grid[d+1]; i++){
		for(int i = opt.grid[d+1]/2; i < opt.grid[d+1]; i++){
		// for (int32_t i = kspace[d].size()/2; i < kspace[d].size(); i++){
			// kspace[d][i] = opt.klength[d]*sin( M_PI*((double)(-opt.grid[d+1]+i))/((double)opt.grid[d+1]) );
			// kspace[d][i] = opt.klength[d]*((double)(opt.grid[d+1]-i))/((double)opt.grid[d+1]/2);
			// kspace[d][i] = opt.klength[d] * 2 * M_PI  * ((double)(-opt.grid[d+1]+i)) / ((double)(opt.grid[d+1]*opt.grid[d+1]*h[d]));
			kspace[d][i] = - opt.klength[d] * M_PI / ( (double)(i - opt.grid[d+1]/2 + 1) * h[d]);
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
				// ares.particle_count += number;
				ares.Ekin += number * k * k;
			}
		}
	}



	
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
	
	