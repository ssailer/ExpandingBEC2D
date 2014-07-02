#include <EXP2D_evaluation.h>

#define OBSERVABLES_DATA_POINTS_SIZE opt.grid[1]*opt.grid[2]
#define ANGULAR_AVERAGING_LENGTH 12

using namespace std;
using namespace Eigen;


Eval::Eval() {

		};

Eval::~Eval() {};

void Eval::saveData(vector<MatrixXcd> &wavefctVec,Options &external_opt,int &external_snapshot_time,string external_runname){
		runname = external_runname;
		opt = external_opt;
		snapshot_time = external_snapshot_time;
		PsiVec.resize(wavefctVec.size());

		#pragma omp parallel for
		for(int k = 0; k < wavefctVec.size(); k++){
			PsiVec[k] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
			for(int i = 0; i < opt.grid[1]; i++){
				for(int j = 0; j < opt.grid[2]; j++){		
					PsiVec[k](0,i,j,0) = wavefctVec[k](i,j);
				}
			}
		}

}

void Eval::saveData(MatrixXcd &wavefct,Options &external_opt,int &external_snapshot_time,string external_runname){
		runname = external_runname;
		opt = external_opt;
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
	contour.resize(PsiVec.size());

	densityLocationMap.resize(PsiVec.size());
	densityCoordinates.resize(PsiVec.size());
	phase = new RealGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	zeros = new RealGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	Contour tracker(opt);


	totalResult = Observables(OBSERVABLES_DATA_POINTS_SIZE);
		
	for(int k = 0; k < PsiVec.size(); k++){
		cout << endl << "Eval #" << k << endl;
		getDensity(PsiVec[k],densityLocationMap[k],densityCoordinates[k]);
		contour[k] = tracker.trackContour(densityLocationMap[k]);
		totalResult += calculator(PsiVec[k],k);		
	}	
	getVortices(PsiVec[0],densityCoordinates[0]);
	totalResult /= PsiVec.size();
}

void Eval::plotData(){
	string filename = runname + "-Spectrum-" + to_string(snapshot_time); 
	plotSpectrum(filename,totalResult);
	filename = runname + "-Vortices-" + to_string(snapshot_time);
	plotVortexList(filename,phase,pres,opt);	
	filename = runname + "-Control-Plot-" + to_string(snapshot_time);
	plotDataToPng(filename,PsiVec[0],opt);
	filename = runname + "-Density-" + to_string(snapshot_time);
	plotDataToPng(filename,densityLocationMap[0],opt);
	filename = runname + "-Density-Axial-Distribution-Gradient" + to_string(snapshot_time);
	plotVector(filename,x_dist_grad,y_dist_grad,opt);
	filename = runname + "-Angular-Dens" + to_string(snapshot_time);
	plotVector(filename,totalResult.angularDensity,opt);	
	filename = runname + "-Contour" + to_string(snapshot_time);
	plotContour(filename,PsiVec[0],contour[0],opt);

	
	filename = runname + "-Observables" + ".dat";
	struct stat buffer;   
  	if(stat (filename.c_str(), &buffer) != 0){
  		ofstream datafile;
  		datafile.open(filename.c_str(), ios::out | ios::app);
  		datafile << std::left << std::setw(10) << "Timestep"
  						 << std::setw(10) << "X_max"
  						 << std::setw(10) << "Y_max"
  						 << std::setw(10) << "A/R"
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
 					 << std::setw(10) << totalResult.aspectRatio 
					 << std::setw(10) << totalResult.particle_count
					 << std::setw(10) << totalResult.volume
					 << std::setw(10) << totalResult.density
					 << std::setw(10) << totalResult.Ekin
			 << endl;
	datafile.close();
}

void Eval::getVortices(ComplexGrid &data, vector<Coordinate<int32_t>> &densityCoordinates){
	
	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2]; 
	
	calc_fields(data,opt);
	pres.vlist.clear();
	find_vortices(densityCoordinates,pres.vlist);



	cout << "Vortices: " << endl;
	double number = 0;
	for(list<VortexData>::const_iterator it = pres.vlist.begin(); it != pres.vlist.end(); ++it){
		int x = it->x.x();
		int y = it->x.y();
		number += it->num_points;
		double sD = it->surroundDens;
		cout << " " << x << " " << y << "  " << abs2(PsiVec[0](0,x,y,0)) << " " << arg(PsiVec[0](0,x,y,0)) << " " << sD << endl;
	}
	cout << "Number of Vortices counted: " << number << "  " << endl;

	// calc_vortex_veloctities();
	// calc_vortex_discances();
	// calc_g2();
}

int Eval::get_phase_jump(const Coordinate<int32_t> &c, const Vector<int32_t> &v, const RealGrid *phase){
	if(phase->at(0,c + v) + M_PI < phase->at(0,c))	// Phase ueberschreitet 0/2pi von unten
		return 1;
	else if(phase->at(0,c) + M_PI < phase->at(0,c + v))	// Phase ueberschreitet 0/2pi von oben
		return -1;
	else
		return 0;
}

// bool Eval::checkRanges(Coordinate<int32_t> c,Coordinate<int32_t> d){
// 	if(c.x() > d.x()){
// 		if(c.y() > d.y()
// 	}
// }

void Eval::find_vortices(vector<Coordinate<int32_t>> &densityCoordinates, list<VortexData> &vlist) {

	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2];
	// Nullstellen zaehlen
	// vector< vector< vector<bool > > > checked(phase->width(), vector< vector<bool> >(phase->height(), vector<bool>(phase->depth(),false)));	// Welche felder schon ueberprueft wurden
	VortexData vortex;
					// Charakteristika eines gefundenen Vortex

	for(vector<Coordinate<int32_t>>::const_iterator it = densityCoordinates.begin(); it != densityCoordinates.end(); ++it){
		// if()

	// for (int z = 0; z < phase->depth(); z++)
	// {
		// for (int x = 0; x < phase->width(); x++)
	// 	{
			// for (int y = 0; y < phase->height(); y++)
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
		// }
		



	}

	// CHECK DENSITY AROUND VORTEX
	int max_radius = 10;
	for(list<VortexData>::iterator it = vlist.begin(); it != vlist.end(); ++it){
		vector<double> polarDensity;
		vector<double> radius;	
		for(int x_shift = - max_radius; x_shift < max_radius; x_shift++){
		    for(int y_shift = - max_radius; y_shift < max_radius; y_shift++){
				radius.push_back(sqrt(x_shift*x_shift * h_x*h_x + y_shift*y_shift *h_y*h_y));
				int x = it->points.front().x() + x_shift;
				int y = it->points.front().y() + y_shift;
				polarDensity.push_back(abs2(PsiVec[0](0,x,y,0)));
			}
		}
		double sum = 0;
		for(int i = 0; i < radius.size(); i++){
			if(radius[i] < max_radius){
				sum += polarDensity[i];
			}

		}
		it->surroundDens = sum;
	}

	// This Number is set at the start, maybe set this in run.cfg -> Options struct, or check how many got set inside the contour, if equal spacing vortices are used.
	const int NUMBER_OF_VORTICES = 4;
	vlist.sort([](VortexData &lhs, VortexData &rhs) {return lhs.surroundDens > rhs.surroundDens;});
	if(vlist.size() > NUMBER_OF_VORTICES){
		list<VortexData>::iterator it1 = vlist.begin();
		advance(it1,NUMBER_OF_VORTICES);
		vlist.erase(it1,vlist.end());
	}
}

void Eval::calc_fields(ComplexGrid &data, Options &opt){
	const double ZERO_THRESHOLD = 0;//opt.N * 0.05 / (4. * opt.min_x * opt.stateInformation[0] * opt.min_y * opt.stateInformation[1]); ; //opt.N * 0.05 / data.width() / data.height() / data.depth();
	for(int x = 0; x < data.width(); x++)
	{
		for(int y = 0; y < data.height(); y++)
		{
			for(int z = 0; z < data.depth(); z++)
			{
				if(abs2(data(0,x,y,z)) < ZERO_THRESHOLD)
					zeros->at(0,x,y,z) = 0.0;
				else
					zeros->at(0,x,y,z) = 1.0;
			}
		}
	}
}



// findLongestLines(){
// 	Coordinate<int32_t> c = PsiVec[0].make_coord(0,0,0);
// 	vector<lineData> xlines;
// 	vector<lineData> ylines;

// 	for(int32_t y = 0; y < opt.grid[2]; y++){
// 		c.y() = y;
// 		do{
// 			if(abs2(data(0,c)) > 0){
// 				lineData line;
// 				line.length = 1;
// 				line.start = c;
// 				c = c + v_right;

// 				do{	
// 					if(abs2(data(0,c)) > 0){
// 						c = c + v_right
// 					}
// 					if((abs2(data(0,c)) == 0) && abs2(data(0,c+v_right)) > 0){
// 						c = c + v_right + v_right;
// 					}
// 					line.stop = c;
// 					line.length++;
// 				}while(abs2(data(0,c)) > 0);
				

				
// 				xlines.push_back(line);
// 			}
// 		}
// 		c = c + v_up;
// 	}
// }

void Eval::getDensity(ComplexGrid &data, RealGrid &densityLocationMap_local, vector<Coordinate<int32_t>> &densityCoordinates_local){
	const double LOWER_THRESHOLD = opt.N * 0.05 / (4. * opt.min_x * opt.stateInformation[0] * opt.min_y * opt.stateInformation[1]);  //abs2(data(0,opt.grid[1]/2,opt.grid[2]/2,0))*0.9;
	// double upper_threshold = 20.;
	// cout << LOWER_THRESHOLD << "+++" << endl;

	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2]; 

	densityLocationMap_local = RealGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);


	densityCounter = 0;
	densityCoordinates_local.clear();
	for(int i = 0; i < opt.grid[1]; i++){
	    for(int j = 0; j < opt.grid[2]; j++){
	    	if((abs2(data(0,i,j,0)) > LOWER_THRESHOLD)){
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
	
	Observables obs = Observables(OBSERVABLES_DATA_POINTS_SIZE);
	// R-Space
	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2];
	double h[2];
	h[0] = h_x;
	h[1] = h_y; 
	// double raw_volume = h_x * opt.grid[1] * h_y * opt.grid[2];
	
	// double threshold = abs2(data(0,opt.grid[1]/2,opt.grid[2]/2,0))*0.9;

	obs.volume = h_x * h_y * densityCounter;
	for(int i = 0; i < opt.grid[1]-1; i++){
	    for(int j = 0; j < opt.grid[2]-1; j++){	    	    		
	      	obs.particle_count += h_x * h_y * abs2(data(0,i,j,0));
	    }
	}
	obs.density = obs.particle_count / obs.volume;

	// Aspect-Ratio
	vector<contourData> cData;

	for(c_set::iterator it = contour[sampleindex].begin(); it != contour[sampleindex].end(); ++it){
		contourData tmp;
		tmp.c = *it;
		int x_shift = tmp.c.x() - opt.grid[1]/2;
		int y_shift = tmp.c.y() - opt.grid[2]/2;
		tmp.phi = atan2(x_shift * h_x,y_shift * h_y) * 180 /M_PI + 180;
		tmp.r = sqrt(x_shift*x_shift * h_x*h_x + y_shift*y_shift *h_y*h_y);
		cData.push_back(tmp);
	}

	std::sort(cData.begin(),cData.end(),[](const contourData &lhs, const contourData &rhs) -> bool {return (lhs.phi < rhs.phi);});

	vector<double> cRadius(361);
	vector<int> divisor_counter(361);
	for(vector<contourData>::iterator it = cData.begin(); it != cData.end(); ++it){
		int index = round(it->phi);
		cRadius[index] += it->r;
		divisor_counter[index]++;
	}
	// cRadius.erase(cRadius.begin());
	// divisor_counter.erase(divisor_counter.begin());
	cRadius[0] += cRadius[360]; cRadius.pop_back();
	divisor_counter[0] += cRadius[360]; divisor_counter.pop_back();

	for(int i = 0; i < 360; i++){
		if(divisor_counter[i] == 0){
			divisor_counter[i] = 1;
		}

		cRadius[i] /= divisor_counter[i];
	}

	// string name = "cRadius_" + to_string(sampleindex) + "_" + to_string(sampleindex);
	// plotVector(name,cRadius,opt);

	vector<double> cDistance(180);
	for(int i = 0; i < 180; i++){
		cDistance[i] = fabs(cRadius[i] + cRadius[i+180]);
	}
	double tmp_ratio = 0;
	for(int i = 0; i < 89; i++){
		double tmp1 = cDistance[i] / cDistance[i+90];
		double tmp2 = cDistance[i+1] / cDistance[i+91];
		tmp_ratio = (tmp1 > tmp2) ? tmp1 : tmp2;
	}
	obs.aspectRatio = tmp_ratio; // zwischen 0 und 90 grad, also effektiv x und y richtung
	// FIXME replace this with a check for the max and min values, save the corresponding angles and check if they change (= overall rotation in the gas!)

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
				obs.angularDensity(i) += polarDensity[index[l][m]];
			}
		}					
	}
	obs.angularDensity /= (ANGULAR_AVERAGING_LENGTH*2 + 1);

	// angularDensity_tmp[0] += angularDensity_tmp[360];	

	// for(int i = 0; i < 360; i++){
	// 	// vector<double>::iterator beginning = angularDensity_tmp.begin() + i;
	// 	// vector<double>::iterator ending = beginning + 10;
	// 	obs.angularDensity(i) = angularDensity_tmp[i]; //(accumulate(beginning,ending,0));
	// }


	
	// K-Space
	ComplexGrid::fft(data, data);
	
	ArrayXd divisor(obs.number.size());
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
	
	double index_factor = (obs.number.size() - 1) / sqrt(kwidth2[0] + kwidth2[1]);

	for(int x = 0; x < data.width(); x++){
		for (int y = 0; y < data.height(); y++){
			for (int z = 0; z < data.depth(); z++){
				double k = sqrt(kspace[0][x]*kspace[0][x] + kspace[1][y]*kspace[1][y]);
				// Coordinate<int32_t> c = data.make_coord(x,y,z);
				int index = index_factor * k;
				// cout << k << "*" << index_factor << "=" << index << "/" << OBSERVABLES_DATA_POINTS_SIZE << endl;
				obs.k(index) += k;
				divisor(index)++;
				double number = abs2(data(0,x,y,z));
				obs.number(index) += number;
				// obs.particle_count += number;
				obs.Ekin += number * k * k;
			}
		}
	}



	
	// obs.healing_length = 1 / sqrt(obs.particle_count * opt.g / obs.volume);
	
	#pragma omp parallel for schedule(guided,1)
	for(int l = 0; l < obs.number.size(); l++){
		if(divisor[l] == 0){
			divisor[l] = 1;
		}
	}

	obs.number /= divisor;
	obs.k /= divisor;	
	
	return obs;
}
	
	