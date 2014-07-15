#include <EXP2D_evaluation.h>

#define OBSERVABLES_DATA_POINTS_SIZE opt.grid[1]*opt.grid[2]
#define ANGULAR_AVERAGING_LENGTH 12
#define NUMBER_OF_VORTICES 100
#define VORTEX_SURROUND_DENSITY_RADIUS 10
#define EDGE_RANGE_CHECK 10

using namespace std;
using namespace Eigen;


Eval::Eval() {};

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

void Eval::checkEdges(){
	int numberOfEdgePoints = 2 * EDGE_RANGE_CHECK * opt.grid[1] + 2 * EDGE_RANGE_CHECK * opt.grid[2] - 4 * EDGE_RANGE_CHECK * EDGE_RANGE_CHECK;
	double threshold = numberOfEdgePoints * opt.N * 0.05 / (4. * opt.min_x * opt.stateInformation[0] * opt.min_y * opt.stateInformation[1]);


	// cout << "Size of checkEdges::sum " << PsiVec.size() << endl;
	vector<double> sum(PsiVec.size());
	for(int k = 0; k < PsiVec.size(); k++){
		sum[k] = 0;
		for(int x = 0; x < EDGE_RANGE_CHECK; x++){
			for(int y = 0; y < opt.grid[2]; y++){
				sum[k] += abs2(PsiVec[k](0,x,y,0));
			}
		}
		for(int x = opt.grid[1] - EDGE_RANGE_CHECK; x < opt.grid[1]; x++){
			for(int y = 0; y < opt.grid[2]; y++){
				sum[k] += abs2(PsiVec[k](0,x,y,0));
			}
		}
		for(int y = 0; y < EDGE_RANGE_CHECK; y++){
			for(int x = 0; x < opt.grid[1]; x++){
				sum[k] += abs2(PsiVec[k](0,x,y,0));
			}
		}
		for(int y = opt.grid[2] - EDGE_RANGE_CHECK; y < opt.grid[2]; y++){
			for(int x = 0; x < opt.grid[1]; x++){
				sum[k] += abs2(PsiVec[k](0,x,y,0));
			}
		}
	}
	for(int k = 0; k < PsiVec.size(); k++){
		std::string error = "Gas reached the edges of the grid, with a density of " + to_string(sum[k]) + "/" + to_string(threshold) + ".  The simulation failed in step: ";
		// cout << "checkEdges[" << k << "] result: " << sum[k] << "/" << threshold << " with " << numberOfEdgePoints << " of points on the edge." << endl;
		if(sum[k] > threshold){			
			expException e(error);
			throw;
		}
	}
}

void Eval::evaluateData(){
	// cout << "evaluateData" << endl;	
	// cout << "checkEdges call: " << endl;
	// checkEdges();

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
		cout << "-getDensity" << endl;
		contour[k] = tracker.trackContour(densityLocationMap[k]);
		cout << "-trackContour" << endl;
		totalResult += calculator(PsiVec[k],k);
		cout << "-calculator" << endl;		
	}
	getVortices(PsiVec[0],densityCoordinates[0]);
	cout << endl << "-getVortices" << endl;
	totalResult /= PsiVec.size();

	
}

void Eval::evaluateDataITP(){

	pres.vlist.clear();
	densityCoordinates.clear();
	contour.resize(PsiVec.size());

	densityLocationMap.resize(PsiVec.size());
	densityCoordinates.resize(PsiVec.size());
	phase = new RealGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	zeros = new RealGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

	totalResult = Observables(OBSERVABLES_DATA_POINTS_SIZE);
		
	for(int k = 0; k < PsiVec.size(); k++){
		totalResult += calculator(PsiVec[k],k);
	}
	totalResult /= PsiVec.size();
}

void Eval::plotData(){
	std::string snapShotString = to_string(snapshot_time);
	std::stringstream ss;
	ss << std::setfill('0') << std::setw(4) << snapShotString;
	snapShotString = ss.str();

	string filename = runname + "-Control-Plot-" + snapShotString;
	plotDataToPngExpanding(filename,PsiVec[0],opt);

	filename = runname + "-Spectrum-" + snapShotString; 
	plotSpectrum(filename,totalResult);

	filename = runname + "-Vortices-" + snapShotString;
	plotVortexList(filename,phase,pres,opt);	

	filename = runname + "-Density-" + snapShotString;
	plotDataToPng(filename,densityLocationMap[0],opt);

	// filename = runname + "-Density-Axial-Distribution-Gradient-" + snapShotString;
	// plotVector(filename,x_dist_grad,y_dist_grad,opt);

	filename = runname + "-Angular-Dens-" + snapShotString;
	plotVector(filename,totalResult.angularDensity,opt);	

	filename = runname + "-Contour-" + snapShotString;
	plotContour(filename,PsiVec[0],contour[0],opt);

	
	filename = runname + "-Observables" + ".dat";
	struct stat buffer;   
  	if(stat (filename.c_str(), &buffer) != 0){
  		ofstream datafile;
  		datafile.open(filename.c_str(), ios::out | ios::app);
  		datafile << std::left << std::setw(10) << "Timestep"
  						 << std::setw(10) << "X_max"
  						 << std::setw(10) << "Y_max"
  						 << std::setw(10) << "R_max"
  						 << std::setw(10) << "R_min"
  						 << std::setw(10) << "R_max/R_min"
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
 					 << std::setw(10) << totalResult.r_max
 					 << std::setw(10) << totalResult.r_min
 					 << std::setw(10) << totalResult.r_max / totalResult.r_min  
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
	// cout << "calc_fields" << endl;
	pres.vlist.clear();
	findVortices(densityCoordinates,pres.vlist);

	// cout << "Vortices: " << endl;
	// double number = 0;
	// for(list<VortexData>::const_iterator it = pres.vlist.begin(); it != pres.vlist.end(); ++it){
	// 	int x = it->x.x();
	// 	int y = it->x.y();
	// 	number += it->num_points;
	// 	double sD = it->surroundDens;
	// 	cout << " " << x << " " << y << "  " << abs2(PsiVec[0](0,x,y,0)) << " " << arg(PsiVec[0](0,x,y,0)) << " " << sD << endl;
	// }
	// cout << "Number of Vortices counted: " << number << "  " << endl;

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

void Eval::findVortices(vector<Coordinate<int32_t>> &densityCoordinates, list<VortexData> &vlist) {

	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2];
	// Nullstellen zaehlen
	// vector< vector< vector<bool > > > checked(phase->width(), vector< vector<bool> >(phase->height(), vector<bool>(phase->depth(),false)));	// Welche felder schon ueberprueft wurden
	VortexData vortex;
					// Charakteristika eines gefundenen Vortex

	// cout << "findVortices_before iterator" << endl;



	for(vector<Coordinate<int32_t>>::const_iterator it = densityCoordinates.begin(); it != densityCoordinates.end(); ++it){

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


	// cout << "findVortices before denscheck" << endl;

	for(list<VortexData>::iterator it = vlist.begin(); it != vlist.end(); ++it){
		vector<double> polarDensity;
		vector<double> radius;	
		for(int x_shift = - VORTEX_SURROUND_DENSITY_RADIUS; x_shift < VORTEX_SURROUND_DENSITY_RADIUS; x_shift++){
		    for(int y_shift = - VORTEX_SURROUND_DENSITY_RADIUS; y_shift < VORTEX_SURROUND_DENSITY_RADIUS; y_shift++){
		    	int x = it->points.front().x() + x_shift;
				int y = it->points.front().y() + y_shift;
				radius.push_back(sqrt(x_shift*x_shift * h_x*h_x + y_shift*y_shift *h_y*h_y));
				polarDensity.push_back(abs2(PsiVec[0](0,x,y,0)));
			}
		}
		double sum = 0;
		double zeroDensity = 0;
		for(int i = 0; i < radius.size(); i++){
			if(radius[i] < VORTEX_SURROUND_DENSITY_RADIUS){
				sum += polarDensity[i];
			}
			if(radius[i] < 2){
				zeroDensity += polarDensity[i];
			}
		}
		it->zeroDensity = zeroDensity;
		it->surroundDens = sum;
	}

	// This Number is set at the start, maybe set this in run.cfg -> Options struct, or check how many got set inside the contour, if equal spacing vortices are used.
	 // cout << "findVortices before sorting" << endl;
	vlist.sort([](VortexData &lhs, VortexData &rhs) {return lhs.zeroDensity < rhs.zeroDensity;});
	if(vlist.size() > opt.vortexnumber){
		list<VortexData>::iterator it1 = vlist.begin();
		advance(it1,opt.vortexnumber);
		vlist.erase(it1,vlist.end());
	}
}

void Eval::calc_fields(ComplexGrid &data, Options &opt){
	double LOWER_THRESHOLD = opt.N * 0.05 / (4. * opt.min_x * opt.stateInformation[0] * opt.min_y * opt.stateInformation[1]); ; //opt.N * 0.05 / data.width() / data.height() / data.depth();
	for(int x = 0; x < data.width(); x++)
	{
		for(int y = 0; y < data.height(); y++)
		{
			for(int z = 0; z < data.depth(); z++)
			{	
				phase->at(0,x,y,z) = arg(data(0,x,y,z));
				if(abs2(data(0,x,y,z)) <= LOWER_THRESHOLD)
					zeros->at(0,x,y,z) = 0.0;
				else
					zeros->at(0,x,y,z) = 1.0;
			}
		}
	}
	// string name = "Test Zeros";
	// plotDataToPng(name, *zeros, opt);
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

void Eval::getDensity(ComplexGrid &data, RealGrid &densityLocationMap, vector<Coordinate<int32_t>> &densityCoordinates_local){
	double threshold = 10;//opt.N * 0.10 / (4. * opt.min_x * opt.stateInformation[0] * opt.min_y * opt.stateInformation[1]);  //abs2(data(0,opt.grid[1]/2,opt.grid[2]/2,0))*0.9;
	// double upper_threshold = 20.;
	// cout << "Threshold " << threshold << endl;

	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2]; 

	RealGrid densityLocationMap_local = RealGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);



	densityCounter = 0;
	densityCoordinates_local.clear();
	for(int i = VORTEX_SURROUND_DENSITY_RADIUS; i < opt.grid[1] - VORTEX_SURROUND_DENSITY_RADIUS; i++){
	    for(int j = VORTEX_SURROUND_DENSITY_RADIUS; j < opt.grid[2] - VORTEX_SURROUND_DENSITY_RADIUS; j++){
	    	if((abs2(data(0,i,j,0)) > threshold)){
   				densityLocationMap_local(0,i,j,0) = 1.;
				densityCoordinates_local.push_back(data.make_coord(i,j,0));
				densityCounter++;
			}else{
				densityLocationMap_local(0,i,j,0) = 0.;
			}
		}
	}

	// string testname = "ERROR_0-getDensity"+to_string(omp_get_wtime());
	// plotDataToPng(testname,densityLocationMap_local,opt);
	// testname = "ERROR_0-data"+to_string(omp_get_wtime());
	// plotDataToPng(testname,data,opt);




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

	densityLocationMap = densityLocationMap_local;
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
	obs.r_max = 0;
	obs.r_min = (opt.grid[1] * h_x >= opt.grid[2] * h_y) ? opt.grid[1] * h_x : opt.grid[2] * h_y;
	vector<contourData> cData;
	for(c_set::iterator it = contour[sampleindex].begin(); it != contour[sampleindex].end(); ++it){
		contourData tmp;
		tmp.c = *it;
		int x_shift = tmp.c.x() - opt.grid[1]/2;
		int y_shift = tmp.c.y() - opt.grid[2]/2;
		tmp.phi = atan2(x_shift * h_x,y_shift * h_y) * 180 /M_PI + 180;
		tmp.r = sqrt(x_shift*x_shift * h_x*h_x + y_shift*y_shift *h_y*h_y);
		cData.push_back(tmp);
		if(obs.r_max <= tmp.r){
			obs.r_max = tmp.r;
			obs.r_max_phi = tmp.phi;
		}
		if(obs.r_min >= tmp.r){
			obs.r_min = tmp.r;
			obs.r_min_phi = tmp.phi;
		}

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
			kspace[d][i] = opt.klength[d]/**opt.stateInformation[0]*/*2.0*sin( M_PI*((double)i)/((double)opt.grid[d+1]) );
			// kspace[d][i] = opt.klength[d]*((double)i)/((double)(opt.grid[d+1]/2));
			// kspace[d][i] = opt.klength[d] * 2 * M_PI  * ((double)i) / ((double)(opt.grid[d+1]*opt.grid[d+1]*h[d]));
			// kspace[d][i] = opt.klength[d] * M_PI / ( (double)(opt.grid[d+1]/2 - i) * h[d] );
		}
		// for (int i=opt.grid[d+1]/2; i<opt.grid[d+1]; i++){
		for(int i = opt.grid[d+1]/2; i < opt.grid[d+1]; i++){
		// for (int32_t i = kspace[d].size()/2; i < kspace[d].size(); i++){
			kspace[d][i] = opt.klength[d]/**opt.stateInformation[1]*/*2.0*sin( M_PI*((double)(-opt.grid[d+1]+i))/((double)opt.grid[d+1]) );
			// kspace[d][i] = opt.klength[d]*((double)(opt.grid[d+1]-i))/((double)opt.grid[d+1]/2);
			// kspace[d][i] = opt.klength[d] * 2 * M_PI  * ((double)(-opt.grid[d+1]+i)) / ((double)(opt.grid[d+1]*opt.grid[d+1]*h[d]));
			// kspace[d][i] = - opt.klength[d] * M_PI / ( (double)(i - opt.grid[d+1]/2 + 1) * h[d]);
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
	
	