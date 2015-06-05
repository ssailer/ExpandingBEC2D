#include <EXP2D_evaluation.h>

#define OBSERVABLES_DATA_POINTS_SIZE data.meta.grid[0]*data.meta.grid[1]
#define ANGULAR_AVERAGING_LENGTH 12
#define NUMBER_OF_VORTICES 100
#define DENSITY_CHECK_DISTANCE 5
#define EDGE_RANGE_CHECK 0.9

using namespace std;
using namespace Eigen;


Eval::Eval(MatrixData d,Options o, string runName) : data(d),  opt(o) , runname(runName) {
	// data = d;
	// opt = o;
	toPhysicalUnits(opt);
	data.convertToPhysicalUnits();
};

void Eval::process(){

	// pres.resize(data.wavefunction.size());

	vlist.resize(data.wavefunction.size());

	densityCoordinates.clear();

	contour.resize(data.wavefunction.size());

	densityCounter.resize(data.wavefunction.size());

	densityLocationMap.resize(data.wavefunction.size());
	densityCoordinates.resize(data.wavefunction.size());

	phase = MatrixXd::Zero(data.meta.grid[0],data.meta.grid[1]);

	Contour tracker(data.meta);

	totalResult = Observables(OBSERVABLES_DATA_POINTS_SIZE);

	cout << currentTime() <<  " Step: " << data.meta.steps << " Time : " << data.meta.time << " s ";		 

	getDensity();
	cout << "dens " ;
	// cout  << "Evaluating sample #: ";
	for(int k = 0; k < data.wavefunction.size(); k++){

		contour[k] = tracker.trackContour(densityLocationMap[k]);
		cout << "con " ;
		totalResult += calculator(data.wavefunction[k],k);
		cout << "calc " ;
		getVortices(data.wavefunction[k],densityCoordinates[k],vlist[k]);
		cout << "vort " ;
		// getVortexDistance(pres[k]);
		// cout << "-getVortexDistance" ;
	}	
	totalResult /= data.wavefunction.size();

}

void Eval::save(){

	string dirname = "runObservables";
    struct stat st;
    	if(stat(dirname.c_str(),&st) != 0){
        mkdir(dirname.c_str(),0755);
    }

	string filename = dirname + "/" + opt.runmode + "_Observables.dat";	
	
	struct stat buffer;   
  	if(stat (filename.c_str(), &buffer) != 0){
  		ofstream datafile;
  		datafile.open(filename.c_str(), ios::out | ios::app);
  		datafile << std::left << std::setw(15) << "Timestep" << ","
  						 << std::setw(15) << "Time" << ","
  						 << std::setw(15) << "X_max" << ","
  						 << std::setw(15) << "Y_max" << ","
					 	 << std::setw(15) << "Vortexnumber" << ","
  						 << std::setw(15) << "D_max" << ","
  						 << std::setw(15) << "D_min" << ","
  						 << std::setw(15) << "Rx" << ","
						 << std::setw(15) << "Ry" << ","
						 << std::setw(15) << "Rx/Ry" << ","
  						 << std::setw(15) << "D_max/D_min" << ","
  						 << std::setw(15) << "D_max Angle" << ","
  						 << std::setw(15) << "D_min Angle" << ","
  						 << std::setw(15) << "Ratio" << ","
  						 << std::setw(15) << "RatioAngle" << ","
  						 << std::setw(15) << "N" << ","
  						 << std::setw(15) << "V" << ","
  						 << std::setw(15) << "N/V" << ","
  						 << std::setw(15) << "E_kin" << ","
  						 << std::setw(15) << "n0"
  				 << endl;
  		datafile.close();
  	}

  	double n0 = 2 * (opt.N / M_PI) * (1 / (totalResult.Rx * totalResult.Ry)); 

  	ofstream datafile(filename.c_str(), std::ios_base::out | std::ios_base::app);
	// datafile.open;
	datafile << std::left << std::setw(15) << data.meta.steps << ","
					 << std::setw(15)  << data.meta.time << ","
					 << std::setw(15)  << opt.min_x * opt.stateInformation[0] << ","
					 << std::setw(15)  << opt.min_y * opt.stateInformation[1] << ","
					 << std::setw(15)  << opt.vortexnumber << ","
					 << std::setw(15)  << totalResult.r_max << ","
 					 << std::setw(15)  << totalResult.r_min << ","
 					 << std::setw(15)  << totalResult.Rx << ","
 					 << std::setw(15)  << totalResult.Ry << ","
 					 << std::setw(15)  << totalResult.Rx / totalResult.Ry << ","
 					 << std::setw(15)  << totalResult.r_max / totalResult.r_min   << ","
 					 << std::setw(15)  << totalResult.r_max_phi << ","
 					 << std::setw(15)  << totalResult.r_min_phi << ","
 					 << std::setw(15)  << totalResult.aspectRatio  << ","
 					 << std::setw(15)  << totalResult.aspectRatioAngle  << ","
					 << std::setw(15)  << totalResult.particle_count << ","
					 << std::setw(15)  << totalResult.volume << ","
					 << std::setw(15)  << totalResult.density << ","
					 << std::setw(15)  << totalResult.Ekin << ","
					 << std::setw(15)  << n0
			 << endl;
	datafile.close();
	cout << "save" << endl;
}



void Eval::getVortices(MatrixXcd &DATA, vector<Coordinate<int32_t>> &densityCoordinates,list<VortexData> &vlist){
	
	calc_fields(DATA,opt);

	vlist.clear();
	findVortices(densityCoordinates,vlist);

	// CHECK DENSITY AROUND VORTEX

	for(list<VortexData>::iterator it = vlist.begin(); it != vlist.end(); ++it){
		list<VortexData>::iterator it1 = it;
		
		for(++it1; it1 != vlist.end();){
			Vector<int32_t> distance = it->c - it1->c;
			if( distance.norm() <= 2.0){
				it1 = vlist.erase(it1);
			}else{
				++it1;
			}
		}

	}

	for(list<VortexData>::iterator it = vlist.begin(); it != vlist.end(); ++it){
		vector<double> polarDensity;
		vector<double> radius;	
		for(int x_shift = - DENSITY_CHECK_DISTANCE; x_shift < DENSITY_CHECK_DISTANCE; x_shift++){
		    for(int y_shift = - DENSITY_CHECK_DISTANCE; y_shift < DENSITY_CHECK_DISTANCE; y_shift++){
		    	int x = it->c.x() + x_shift;
				int y = it->c.y() + y_shift;
				radius.push_back(sqrt(x_shift*x_shift * data.meta.spacing[0]*data.meta.spacing[0] + y_shift*y_shift *data.meta.spacing[1]*data.meta.spacing[1]));
				polarDensity.push_back(abs2(data.wavefunction[0](x,y)));
			}
		}
		double sum = 0;
		double zeroDensity = 0;
		for(int i = 0; i < radius.size(); i++){
			if(radius[i] < DENSITY_CHECK_DISTANCE){
				sum += polarDensity[i];
			}
			if(radius[i] < 2){
				zeroDensity += polarDensity[i];
			}
		}
		it->zeroDensity = zeroDensity;
		it->surroundDens = sum;
	}

	// COMMENTED OUT TO GET INCREASING NUMBERS OF VORTICES

	// vlist.sort([](VortexData &lhs, VortexData &rhs) {return lhs.surroundDens > rhs.surroundDens;});
	
	// if(opt.initialRun == true){
		opt.vortexnumber = vlist.size();
	// 	cout << "Evaluation found " << opt.vortexnumber << " Vortices." << endl;
	// } else {
	// 	if(vlist.size() > opt.vortexnumber){
	// 		list<VortexData>::iterator it1 = vlist.begin();
	// 		advance(it1,opt.vortexnumber);
	// 		vlist.erase(it1,vlist.end());
	// 	}
	// }
}

int Eval::get_phase_jump(const Coordinate<int32_t> &c, const Vector<int32_t> &v){
	if(phase(c.x() + v.x(),c.y()+v.y()) + M_PI < phase(c.x(),c.y()))	// Phase ueberschreitet 0/2pi von unten
		return 1;
	else if(phase(c.x(),c.y()) + M_PI < phase(c.x() + v.x(),c.y()+v.y()))	// Phase ueberschreitet 0/2pi von oben
		return -1;
	else
		return 0;
}

void Eval::findVortices(vector<Coordinate<int32_t>> &densityCoordinates, list<VortexData> &vlist) {

	VortexData vortex;

	Vector<int32_t> down = Vector<int32_t>(0,-1,0,data.meta.grid[0],data.meta.grid[1],1);
	Vector<int32_t> right = Vector<int32_t>(1,0,0,data.meta.grid[0],data.meta.grid[1],1);
	Vector<int32_t> up = Vector<int32_t>(0,1,0,data.meta.grid[0],data.meta.grid[1],1);
	Vector<int32_t> left = Vector<int32_t>(-1,0,0,data.meta.grid[0],data.meta.grid[1],1);
	Vector<int32_t> rightdown = Vector<int32_t>(0.5, -0.5, 0,data.meta.grid[0],data.meta.grid[1],1);

	// #pragma omp parallel for
	for(int i = 0; i < densityCoordinates.size(); i++){
	// for(vector<Coordinate<int32_t>>::const_iterator it = densityCoordinates.begin(); it != densityCoordinates.end(); ++it){

		Coordinate<int32_t> c = densityCoordinates[i];
			
		int phase_winding = get_phase_jump(c, down) + get_phase_jump(c+down, right) + get_phase_jump(c+down+right, up) + get_phase_jump(c+right, left);

		if(phase_winding != 0){
			vortex.n = phase_winding;
			vortex.c = c +right; // /*+ rightdown*/;
			vlist.push_back(vortex);
		}
	}

	// cout << "Vortexnumber before surround dens check: " << vlist.size() << endl;


}

inline double Eval::norm(Coordinate<double> &a, Coordinate<double> &b, double &h_x, double &h_y){
	return sqrt( (a.x() - b.x()) * (a.x() - b.x()) * h_x * h_x + (a.y() - b.y()) * (a.y() - b.y()) * h_y * h_y);
}

// void Eval::getVortexDistance(list<VortexData &vlist){

// 	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
// 	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2];

// 	pres.histogram.resize(OBSERVABLES_DATA_POINTS_SIZE);
// 	pres.distance.resize(OBSERVABLES_DATA_POINTS_SIZE);
// 	for (list<VortexData>::iterator it = pres.vlist.begin(); it != pres.vlist.end(); ++it)
// 	{
// 		if(it->n != 0)
// 		{
// 			double shortest_distance = 65536.0;
// 			for(list<VortexData>::iterator oit = pres.vlist.begin(); oit != pres.vlist.end(); ++oit)
// 			{
// 				if((it != oit) && (oit->n!=0))
// 				{
// 					double distance = (it->x - oit->x).norm();
// 					double coordDistance = Eval::norm(it->x, oit->x,h_x,h_y);
// 					// cout << "Coordinate Distance: " << coordDistance << endl;
// 					// cout << "Grid Distance: " << distance << endl;
// 					if(distance < shortest_distance)
// 						shortest_distance = distance;
// 					pairDistanceHistogram(pres, distance, coordDistance);
// 				}
// 			}
// 			// if(shortest_distance != 65536.0)
// 			// {
// 			// 	it->pair_distance = shortest_distance;
// 			// 	//inc_pd_histogram(ares.pd_histogram_closest, shortest_distance);
// 			// 	av_pair_distance.average(it->pair_distance);
// 			// 	if(it->pair_distance > ares.max_pair_distance)
// 			// 		ares.max_pair_distance = it->pair_distance;
// 			// }
// 			// else
// 			// 	it->pair_distance = 0.0;
// 		}
// 		// else
// 		// 	it->pair_distance = 0.0;
// 	}
	
// 	// // mittlerer Abstand der Vortices
// 	// ares.pair_distance_all = av_pair_distance.av();
// 	// if(av_pair_distance.av() != 0)
// 	// 	ares.pair_distance_nonzero = av_pair_distance.av();
// }

// inline void Eval::pairDistanceHistogram(PathResults &pres, double &distance, double &coordDistance)
// {
// 	double max_distance = sqrt(opt.grid[1]*opt.grid[1] + opt.grid[2]*opt.grid[2] + opt.grid[3]*opt.grid[3]) / 2.0;
// 	int index = (pres.histogram.size() - 1) * distance / max_distance;
// 	pres.histogram[index] += 0.5;
// 	pres.distance[index] = coordDistance;

// }

int Eval::getVortexNumber(){
	return opt.vortexnumber;
}

void Eval::calc_fields(MatrixXcd &DATA, Options &opt){
	#pragma omp parallel for
	for(int x = 0; x < data.meta.grid[0]; x++)
	{
		for(int y = 0; y < data.meta.grid[1]; y++)
		{
			phase(x,y) = arg(DATA(x,y));
		}
	}
}

void Eval::erosion(MatrixXi &d){
	MatrixXi tmp_map = MatrixXi::Zero(data.meta.grid[0],data.meta.grid[1]);

	for(int i = DENSITY_CHECK_DISTANCE; i < data.meta.grid[0] - DENSITY_CHECK_DISTANCE; i++){
	 	for(int j = DENSITY_CHECK_DISTANCE; j < data.meta.grid[1] - DENSITY_CHECK_DISTANCE; j++){
			int sum = d.block<3,3>(i-1,j-1).sum();
			if(sum == 9){
				tmp_map(i,j) = 1;
			} 
		}
	}
	d = tmp_map;

}

void Eval::dilation(MatrixXi &d){
	MatrixXi tmp_map = MatrixXi::Zero(data.meta.grid[0],data.meta.grid[1]);

	for(int i = DENSITY_CHECK_DISTANCE; i < data.meta.grid[0] - DENSITY_CHECK_DISTANCE; i++){
	 	for(int j = DENSITY_CHECK_DISTANCE; j < data.meta.grid[1] - DENSITY_CHECK_DISTANCE; j++){
			int sum = d.block<3,3>(i-1,j-1).sum();
			if(sum != 0){
				tmp_map(i,j) = 1;
			} 
		}
	}
	d = tmp_map;
}

void Eval::floodFillUtil(MatrixXi &checkedCounter ,MatrixXi &mask ,MatrixXi &dens, int x, int y, int prevC, int newC)
{
    // Base cases
    if (x < 0 || x >= data.meta.grid[0] || y < 0 || y >= data.meta.grid[1])
        return;
    if (dens(x,y) != prevC)
        return;
 
    // Replace the color at (x, y)
    if(checkedCounter(x,y) == 1)
    	return;
    mask(x,y) = newC;
    checkedCounter(x,y) = 1;
    // Recur for north, east, south and west
    // if(checkedCounter(x+1,y) == 0)
    	floodFillUtil(checkedCounter,mask,dens, x+1, y, prevC, newC);
	// if(checkedCounter(x-1,y) == 0)
    	floodFillUtil(checkedCounter,mask,dens, x-1, y, prevC, newC);
    // if(checkedCounter(x,y+1) == 0)
    	floodFillUtil(checkedCounter,mask,dens, x, y+1, prevC, newC);
    // if(checkedCounter(x,y-1) == 0)
    	floodFillUtil(checkedCounter,mask,dens, x, y-1, prevC, newC);
}

MatrixXi Eval::floodFill(MatrixXi &dens,int prevC, int newC){
	MatrixXi tmpD = MatrixXi::Zero(data.meta.grid[0],data.meta.grid[1]);
	MatrixXi checkedCounter = MatrixXi::Zero(data.meta.grid[0],data.meta.grid[1]);
	floodFillUtil(checkedCounter,tmpD,dens,1,1,prevC,newC); // Start at (1,1) and replace prevC of dens with newC's in tmpD, create a mask for the zero regions
	return tmpD;
}

void Eval::fillHoles(MatrixXi &dens, MatrixXi &mask){
	for(int i = 0; i < data.meta.grid[0]; i++){
		for(int j = 0; j < data.meta.grid[1]; j++){
			if(mask(i,j) != 1){
				dens(i,j) = 1;
			}
		}
	}
}


void Eval::getDensity(){

	double maximum = 0;
	for(int k = 0; k < data.wavefunction.size(); k++){
		for(int i = 0; i < data.meta.grid[0]; i++){
			for(int j = 0; j < data.meta.grid[1]; j++){
				double value = abs2(data.wavefunction[k](i,j));
				maximum = ( value > maximum) ? value : maximum;
			}
		}
	}
	double threshold = maximum * 0.05;

	for(int k = 0; k < data.wavefunction.size(); k++){		
		densityLocationMap[k] = MatrixXi::Zero(data.meta.grid[0],data.meta.grid[1]);	
		densityCounter[k] = 0;
		densityCoordinates[k].clear();
		for(int i = DENSITY_CHECK_DISTANCE; i < data.meta.grid[0] - DENSITY_CHECK_DISTANCE; i++){
			for(int j = DENSITY_CHECK_DISTANCE; j < data.meta.grid[1] - DENSITY_CHECK_DISTANCE; j++){
			    if((abs2(data.wavefunction[k](i,j)) > threshold)){
	   				densityLocationMap[k](i,j) = 1;
	   				// Coordinate<int32_t> tmpCoord = Coordinate<int32_t>(i,j,0,data.meta.grid[0],data.meta.grid[1],1);
					// densityCoordinates[k].push_back(tmpCoord);
					// densityCounter[k]++;
				}
			}
		}

		
		// MatrixXf densMapGrad = MatrixXf::Zero(data.meta.grid[0],data.meta.grid[1]);
		// for(int i = DENSITY_CHECK_DISTANCE; i < data.meta.grid[0] - DENSITY_CHECK_DISTANCE; i++){
		// 	for(int j = DENSITY_CHECK_DISTANCE; j < data.meta.grid[1] - DENSITY_CHECK_DISTANCE; j++){
		// 		// smoothing(densityLocationMap[k],i,j);
		// 		for(int m = i - DENSITY_CHECK_DISTANCE; m < i + DENSITY_CHECK_DISTANCE; m++){
		// 			for(int n = j - DENSITY_CHECK_DISTANCE; n < j + DENSITY_CHECK_DISTANCE; n++){
		// 				densMapGrad(i,j) += densityLocationMap[k](m,n);
		// 			}		
		// 		}
		// 		densMapGrad(i,j) /= (DENSITY_CHECK_DISTANCE * DENSITY_CHECK_DISTANCE * 4);
		// 	}
		// }

				// if(densMapGrad(i,j) < 0.9 && densMapGrad(i,j) > 0.1){
		erosion(densityLocationMap[k]);
		dilation(densityLocationMap[k]);
		// dilation(densityLocationMap[k]);	

		MatrixXi zeroMask = floodFill(densityLocationMap[k],0,1);

		fillHoles(densityLocationMap[k],zeroMask);

		for(int i = DENSITY_CHECK_DISTANCE; i < data.meta.grid[0] - DENSITY_CHECK_DISTANCE; i++){
			for(int j = DENSITY_CHECK_DISTANCE; j < data.meta.grid[1] - DENSITY_CHECK_DISTANCE; j++){
			    if(densityLocationMap[k](i,j) == 1){
	   				Coordinate<int32_t> tmpCoord = Coordinate<int32_t>(i,j,0,data.meta.grid[0],data.meta.grid[1],1);
					densityCoordinates[k].push_back(tmpCoord);
					densityCounter[k]++;
				}
			}
		}
	}

	// string testname = "ERROR_0-getDensity"+to_string(omp_get_wtime());
	// plotDataToPng(testname,densityLocationMap[0],opt);
	// testname = "ERROR_0-data"+to_string(omp_get_wtime());
	// plotDataToPngEigen(testname,data.wavefunction[0],opt);




	// angularDensity = polarDensity;

	// if(angularDensity.size() != phi.size())
	// 	cout << "ERROR: Angular Density index problems." << endl;



	// x_dist.resize(data.meta.grid[0]);
	// y_dist.resize(data.meta.grid[1]);

	// for(int i = 0; i < data.meta.grid[0]; i++){
	// 	for(int j = 0; j < data.meta.grid[1]; j++){
	// 		x_dist[i] += densMap(i,j);
	// 		y_dist[j] += densMap(i,j);			
	// 	}
	// }


	// x_dist_grad.resize(data.meta.grid[0]);
	// y_dist_grad.resize(data.meta.grid[1]);
	// for(int x = 1; x < x_dist.size() - 1; x++){
	// 	x_dist_grad[x] = (x_dist[x+1] - x_dist[x-1] ) / (2.0 * data.meta.spacing[0] );
	// }
	// for(int y = 1; y < y_dist.size() - 1; y++){
	// 	y_dist_grad[y] = (y_dist[y+1] - y_dist[y-1] ) / (2.0 * data.meta.spacing[1] );
	// }
	// x_dist_grad[0] = x_dist_grad[data.meta.grid[0]-1] = y_dist_grad[0] = y_dist_grad[data.meta.grid[1]-1] = 0.0;



	// double sum = 0;
	// for(int k = 0; k < 10; k++){
	// 	for(int x = 1; x < opt.grid[0]-1; x+=1){
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

	// densityLocationMap = densityLocationMap_local;
}

void Eval::aspectRatio(Observables &obs, int &sampleindex){
	double h_x = data.meta.spacing[0];
	double h_y = data.meta.spacing[1];

	// Aspect-Ratio
	obs.r_max = 0;
	obs.r_min = (data.meta.grid[0] * h_x >= data.meta.grid[1] * h_y) ? data.meta.grid[0] * h_x : data.meta.grid[1] * h_y;
	vector<contourData> cData;
	for(c_set::iterator it = contour[sampleindex].begin(); it != contour[sampleindex].end(); ++it){
		contourData tmp;
		tmp.c = *it;
		int x = (tmp.c.x() - data.meta.grid[0]/2) * h_x;
		int y = (tmp.c.y() - data.meta.grid[1]/2) * h_y;
		tmp.phi = atan2(y,x) * 180 /M_PI + 180;
		tmp.r = sqrt(x*x  + y*y);
		cData.push_back(tmp);
		// if(obs.r_max <= tmp.r){
		// 	obs.r_max = tmp.r;
		// 	obs.r_max_phi = tmp.phi;
		// }
		// if(obs.r_min >= tmp.r){
		// 	obs.r_min = tmp.r;
		// 	obs.r_min_phi = tmp.phi;
		// }

	}

	std::sort(cData.begin(),cData.end(),[](const contourData &lhs, const contourData &rhs) -> bool {return (lhs.phi < rhs.phi);});

	vector<double> cRadius(361);
	vector<int> divisor_counter(361);
	for(vector<contourData>::const_iterator it = cData.begin(); it != cData.end(); ++it){
		int index = round(it->phi);
		cRadius[index] += it->r;
		divisor_counter[index]++;
	}
	// cRadius.erase(cRadius.begin());
	// divisor_counter.erase(divisor_counter.begin());
	cRadius[0] += cRadius[360]; cRadius.pop_back();
	divisor_counter[0] += divisor_counter[360]; divisor_counter.pop_back();


	for(int i = 0; i < 360; i++){
		if(divisor_counter[i] != 0){
			cRadius[i] /= divisor_counter[i];
		}
	}
	for(int i = 0; i < 360; i++){
		checkNextAngles(cRadius,i);
	}

	obs.Rx = cRadius[0];
	obs.Ry = cRadius[90];
	// string name = "cRadius_" + to_string(sampleindex) + "_" + to_string(sampleindex);
	// plotVector(name,cRadius,opt);

	vector<double> cDistance(180);
	for(int i = 0; i < 180; i++){
		cDistance[i] = fabs(cRadius[i] + cRadius[i+180]);
	}
	vector<double>::iterator maxDistance = std::max_element(cDistance.begin(), cDistance.end());
	obs.r_max = *maxDistance;
	obs.r_max_phi = std::distance(cDistance.begin(), maxDistance);

	vector<double>::iterator minDistance = std::min_element(cDistance.begin(), cDistance.end());
	obs.r_min = *minDistance;
	obs.r_min_phi = std::distance(cDistance.begin(), minDistance);

	vector<double> tmp_ratio(90);
	// double tmp_ratio = 0;
	for(int i = 0; i < 90; i++){
		if(cDistance[i+90] >= 0.0){
			tmp_ratio[i] = cDistance[i] / cDistance[i+90];
			obs.fixedAspectRatio(i) = tmp_ratio[i];
		} else {
			cout << "WARNING: Aspect-Ratio: Calculated Distance smaller than zero!" << endl;
		}
		// double tmp1 = cDistance[i] / cDistance[i+90];
		// double tmp2 = cDistance[i+1] / cDistance[i+91];
		// tmp_ratio = (tmp1 > tmp2) ? tmp1 : tmp2;
	}
	vector<double>::iterator maxElement = std::max_element(tmp_ratio.begin(),tmp_ratio.end());
	int maxAspectRatioIndex = std::distance(tmp_ratio.begin(), maxElement);
	obs.aspectRatioAngle = maxAspectRatioIndex;
	obs.aspectRatio = *maxElement; // zwischen 0 und 90 grad, also effektiv x und y richtung
	// FIXME replace this with a check for the max and min values, save the corresponding angles and check if they change (= overall rotation in the gas!)
}

Observables Eval::calculator(MatrixXcd DATA,int sampleindex){
	
	Observables obs = Observables(OBSERVABLES_DATA_POINTS_SIZE);
	// R-Space
	double h_x = data.meta.spacing[0];
	double h_y = data.meta.spacing[1];
	double h[2];
	double x_max = data.meta.coord[0];
	double y_max = data.meta.coord[1];
	vector<double> rmax(2);
	rmax[0] = data.meta.coord[0];
	rmax[1] = data.meta.coord[1];
	h[0] = data.meta.spacing[0];
	h[1] = data.meta.spacing[1]; 
	// double raw_volume = h_x * opt.grid[1] * h_y * opt.grid[2];
	
	// double threshold = abs2(DATA(0,opt.grid[1]/2,opt.grid[2]/2,0))*0.9;

	// cout << "DensityCounter " << sampleindex << " : " << densityCounter[sampleindex] << endl;

	obs.volume = data.meta.spacing[0] * data.meta.spacing[1] * densityCounter[sampleindex];
	double sum = 0;
	#pragma omp parallel for reduction(+:sum)
	for(int i = 0; i < data.meta.grid[0]; i++){
	    for(int j = 0; j < data.meta.grid[1]; j++){	    	    		
	      	sum += abs2(DATA(i,j));

	    }
	}
	obs.particle_count = sum;
	obs.particle_count *= data.meta.spacing[0] * data.meta.spacing[1];
	obs.density = obs.particle_count / obs.volume;

	aspectRatio(obs,sampleindex);

	// == Angular Density
	// vector<double> angularDensity;
	vector<double> phi;
	vector<double> polarDensity; // first entry is 
	vector<double> radius;
	vector<Coordinate<int32_t>> cartesianCoordinates;

	ArrayXd divisor2(obs.number.size());

	double r_index_factor = (obs.radialDensity.size() -1) / sqrt(data.meta.coord[0] * data.meta.coord[0] + data.meta.coord[1] * data.meta.coord[1]);
	// cout << "INDEX_FACTOR " << r_index_factor << endl;
	for(int i = 0; i < data.meta.grid[0]; i++){
	    for(int j = 0; j < data.meta.grid[1]; j++){
				int x_shift = i - data.meta.grid[0]/2;
				int y_shift = j - data.meta.grid[1]/2;
				phi.push_back( atan2(x_shift * data.meta.spacing[0] ,y_shift * data.meta.spacing[1]) * 180 / M_PI + 180);
				double r_tmp = sqrt(x_shift*x_shift * data.meta.spacing[0]*data.meta.spacing[0] + y_shift*y_shift *data.meta.spacing[1]*data.meta.spacing[1]);

				

				radius.push_back(r_tmp);
				double dens_tmp = abs2(DATA(i,j));
				polarDensity.push_back(dens_tmp);
				Coordinate<int32_t> tmpCoord =  Coordinate<int32_t>(i,j,0,data.meta.grid[0],data.meta.grid[1],1);
				cartesianCoordinates.push_back(tmpCoord);

				int index = r_index_factor * r_tmp;
				// cout << " radius " << r_tmp << " / " << data.meta.coord[0] << " index " << index << " / " << obs.radialDensity.size() << endl;
				// if(index >= obs.number.size()){
				// 	cout << "index " << index << " r_tmp " << r_tmp << " " << data.meta.spacing[0] << " " << endl;
				// }
				divisor2(index) += 1.0;
				obs.r(index) += r_tmp;
				obs.radialDensity(index) += dens_tmp;
		}
	}
	obs.r /= divisor2;
	obs.radialDensity /= divisor2;

	

	vector<vector<int>> index(361);
	for(int j = 0; j < phi.size(); j++){
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



	// vector<int> edges;
	// checkResizeCondition(edges);

	// MatrixXcd smallData = DATA.block(edges[0],edges[1],edges[2],edges[3]); // MatrixXcd((int)*x_minmax.second - (int)*x_minmax.first,(int)*y_minmax.second - (int)*y_minmax.first);

	// plotDataToPngEigen("smallData",smallData,opt);

	// K-Space
	// ComplexGrid::fft(DATA, DATA);
	// data.fft.Forward(smallData);
	data.fft.Forward(DATA);

	DATA /= sqrt(DATA.cols() * DATA.rows());
	// smallData /= sqrt(smallData.cols() * smallData.rows());

	// plotDataToPngEigen("smallData_kspace",smallData,opt);

	// DATA = MatrixXcd::Zero(data.meta.grid[0],data.meta.grid[1]);
	// DATA.block((int)*x_minmax.first,(int)*y_minmax.first,diff_x,diff_y) = smallData;

	// double sum2 = 0;
	// #pragma omp parallel for reduction(+:sum2)
	// for(int i = 0; i < data.meta.grid[0]; i++){
	//     for(int j = 0; j < data.meta.grid[1]; j++){	    	    		
	//       	sum2 += abs2(DATA(i,j));

	//     }
	// }
	// cout << "kSUM : " << sum2 << endl;
	
	ArrayXd divisor(obs.number.size());
	divisor.setZero();
	
	vector<vector<double>> kspace;

	kspace.resize(2);
	for(int d = 0; d < 2; d++){
		// set k-space
		kspace[d].resize(data.meta.grid[d]);
		for(int i = 0; i <= data.meta.grid[d]/2; i++){
			kspace[d][i] = (M_PI / rmax[d]) * (double)i;
		}
		for(int i = (data.meta.grid[d]/2)+1; i < data.meta.grid[d]; i++){
			kspace[d][i] = -(M_PI / rmax[d]) * (double)(data.meta.grid[d] - i);
		}
	}


	// DEFINITION OF KLENGTH FROM THE INTERNET! |||| ==>>     2*pi*i/(Nx*dx)
	// double kmax[2];

	// for(int i = 0; i < 2; i++){
	// 	kmax[i] = *max_element(kspace[i].begin(), kspace[i].end());
	// }
	
	double kwidth2[2];

	for(int i = 0; i < 2; i++)
		kwidth2[i] = (data.meta.grid[i] == 1) ? 0 : kspace[i][data.meta.grid[i]/2] * kspace[i][data.meta.grid[i]/2];
	
	double index_factor = (obs.number.size() - 1) / sqrt(kwidth2[0] + kwidth2[1]);

	for(int x = 0; x < data.meta.grid[0]; x++){
		for (int y = 0; y < data.meta.grid[1]; y++){
				double k = sqrt(kspace[0][x]*kspace[0][x] + kspace[1][y]*kspace[1][y]);
				int index = index_factor * k;
				obs.k(index) += k;
				divisor(index)++;
				double number = abs2(DATA(x,y));
				obs.number(index) += number;
				obs.Ekin += number * k * k;
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

bool Eval::checkResizeCondition(vector<int> &edges){

	edges.resize(4);

	vector<double> x_tmp(densityCoordinates[0].size());
	vector<double> y_tmp(densityCoordinates[0].size());

	for(int i = 0; i < densityCoordinates[0].size(); i++){
		x_tmp[i] = densityCoordinates[0][i].x();
		y_tmp[i] = densityCoordinates[0][i].y();
	}
	auto x_minmax = std::minmax_element (x_tmp.begin(),x_tmp.end());
	auto y_minmax = std::minmax_element (y_tmp.begin(),y_tmp.end());

	edges[0] = (int)*x_minmax.first;
	edges[1] = (int)*y_minmax.first;

	edges[2] = (int)*x_minmax.second - (int)*x_minmax.first;
	edges[3] = (int)*y_minmax.second - (int)*y_minmax.first;

	return (edges[2] >= (data.meta.grid[0] * EDGE_RANGE_CHECK) || edges[3] >= (data.meta.grid[1] * EDGE_RANGE_CHECK));
}		

void Eval::checkNextAngles(vector<double> &r, int &i){
		int l_index = i-1;
		int r_index = i;
		while(cyclicReadout(r,r_index) == 0.0) {
			r_index++;
		}
		if(r_index != i){
			double average = (cyclicReadout(r,l_index) + cyclicReadout(r,r_index))/2.0;
			for(int k = l_index+1;k < r_index;k++){
				cyclicAssignment(r,k,average);
			}
		}
	}

void Eval::cyclicAssignment(vector<double> &r, int i, double rvalue){
		if(i < 0){  cyclicAssignment(r,i+360,rvalue); }
		else if(i > 359){ cyclicAssignment(r,i-360,rvalue); }
		else{  r[i] = rvalue; }
	}

double Eval::cyclicReadout(vector<double> &r, int i){
		if(i < 0){ return cyclicReadout(r,i+360); }
		else if(i > 359){ return cyclicReadout(r,i-360); }
		else{ return r[i]; }
	}

// Observables Eval::calculatorITP(ComplexGrid data,int sampleindex){
	
// 	Observables obs = Observables(OBSERVABLES_DATA_POINTS_SIZE);
// 	// R-Space
// 	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
// 	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2];
// 	double h[2];
// 	h[0] = h_x;
// 	h[1] = h_y; 
// 	// double raw_volume = h_x * opt.grid[1] * h_y * opt.grid[2];
	
// 	// double threshold = abs2(data(0,opt.grid[1]/2,opt.grid[2]/2,0))*0.9;

// 	obs.volume = h_x * h_y * densityCounter[sampleindex];
// 	for(int i = 0; i < opt.grid[1]; i++){
// 	    for(int j = 0; j < opt.grid[2]; j++){	    	    		
// 	      	obs.particle_count += abs2(data(0,i,j,0));
	      	
// 	    }
// 	}
// 	obs.particle_count *= h_x * h_y;
// 	obs.density = obs.particle_count / obs.volume;

// 	// K-Space
// 	ComplexGrid::fft(data, data);
	
// 	ArrayXd divisor(obs.number.size());
// 	divisor.setZero();
	
// 	vector<vector<double>> kspace;
	
// 	kspace.resize(2);
// 	for(int d = 0; d < 2; d++){
// 		// set k-space
// 		kspace[d].resize(opt.grid[d+1]);
// 		// for (int i=0; i<opt.grid[d+1]/2; i++){
// 		for(int i = 0; i < opt.grid[d+1]/2; i++){
// 		// for (int32_t i = 0; i < kspace[d].size()/2; i++){
// 			kspace[d][i] = opt.klength[d]/**opt.stateInformation[0]*/*2.0*sin( M_PI*((double)i)/((double)opt.grid[d+1]) );
// 			// kspace[d][i] = opt.klength[d]*((double)i)/((double)(opt.grid[d+1]/2));
// 			// kspace[d][i] = opt.klength[d] * 2 * M_PI  * ((double)i) / ((double)(opt.grid[d+1]*opt.grid[d+1]*h[d]));
// 			// kspace[d][i] = opt.klength[d] * M_PI / ( (double)(opt.grid[d+1]/2 - i) * h[d] );
// 		}
// 		// for (int i=opt.grid[d+1]/2; i<opt.grid[d+1]; i++){
// 		for(int i = opt.grid[d+1]/2; i < opt.grid[d+1]; i++){
// 		// for (int32_t i = kspace[d].size()/2; i < kspace[d].size(); i++){
// 			kspace[d][i] = opt.klength[d]/**opt.stateInformation[1]*/*2.0*sin( M_PI*((double)(-opt.grid[d+1]+i))/((double)opt.grid[d+1]) );
// 			// kspace[d][i] = opt.klength[d]*((double)(opt.grid[d+1]-i))/((double)opt.grid[d+1]/2);
// 			// kspace[d][i] = opt.klength[d] * 2 * M_PI  * ((double)(-opt.grid[d+1]+i)) / ((double)(opt.grid[d+1]*opt.grid[d+1]*h[d]));
// 			// kspace[d][i] = - opt.klength[d] * M_PI / ( (double)(i - opt.grid[d+1]/2 + 1) * h[d]);
// 		}
// 	}


// 	// DEFINITION OF KLENGTH FROM THE INTERNET! |||| ==>>     2*pi*i/(Nx*dx)
// 	// double kmax[2];

// 	// for(int i = 0; i < 2; i++){
// 	// 	kmax[i] = *max_element(kspace[i].begin(), kspace[i].end());
// 	// }
	
// 	double kwidth2[2];

// 	for(int i = 0; i < 2; i++)
// 		kwidth2[i] = (opt.grid[i+1] == 1) ? 0 : kspace[i][opt.grid[i+1]/2] * kspace[i][opt.grid[i+1]/2];
	
// 	double index_factor = (obs.number.size() - 1) / sqrt(kwidth2[0] + kwidth2[1]);

// 	for(int x = 0; x < data.width(); x++){
// 		for (int y = 0; y < data.height(); y++){
// 			for (int z = 0; z < data.depth(); z++){
// 				double k = sqrt(kspace[0][x]*kspace[0][x] + kspace[1][y]*kspace[1][y]);
// 				// Coordinate<int32_t> c = data.make_coord(x,y,z);
// 				int index = index_factor * k;
// 				// cout << k << "*" << index_factor << "=" << index << "/" << OBSERVABLES_DATA_POINTS_SIZE << endl;
// 				obs.k(index) += k;
// 				divisor(index)++;
// 				double number = abs2(data(0,x,y,z));
// 				obs.number(index) += number;
// 				// obs.particle_count += number;
// 				obs.Ekin += number * k * k;
// 			}
// 		}
// 	}
	
// 	#pragma omp parallel for schedule(guided,1)
// 	for(int l = 0; l < obs.number.size(); l++){
// 		if(divisor[l] == 0){
// 			divisor[l] = 1;
// 		}
// 	}

// 	obs.number /= divisor;
// 	obs.k /= divisor;	
	
// 	return obs;
// }

// void Eval::evaluateDataITP(){

// 	// pres.vlist.clear();
// 	densityCoordinates.clear();
// 	contour.resize(PsiVec.size());
// 	densityCounter.resize(PsiVec.size());

// 	densityLocationMap.resize(PsiVec.size());
// 	densityCoordinates.resize(PsiVec.size());
// 	phase.resize(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
// 	zeros.resize(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

// 	totalResult = Observables(OBSERVABLES_DATA_POINTS_SIZE);
		
// 	// for(int k = 0; k < PsiVec.size(); k++){
// 		getDensity(PsiVec[0],densityLocationMap[0],densityCoordinates[0],densityCounter[0]);
// 		totalResult = calculatorITP(PsiVec[0],0);
// 	// }
// 	// totalResult /= PsiVec.size();
// }

// void Eval::CombinedSpectrum(){
// 	string dirname = "CombinedRunPlots";
//     struct stat st;
//     	if(stat(dirname.c_str(),&st) != 0){
//         mkdir(dirname.c_str(),0755);
//     }
    
// 	std::string snapShotString = to_string(snapshot_time);
// 	std::stringstream ss;
// 	ss << std::setfill('0') << std::setw(5) << snapShotString;
// 	snapShotString = ss.str();
	
// 	runname = dirname + "/" + runname;

// 	string plotname = runname + "-Spectrum-" + snapShotString;
// 	string title = "Spectrum " + snapShotString; 
// 	plotSpectrum(plotname,title, totalResult);

// }

// void Eval::plot(){
// 	string dirname = "runPlots";
//     struct stat st;
//     	if(stat(dirname.c_str(),&st) != 0){
//         mkdir(dirname.c_str(),0755);
//     }
    
// 	std::string snapShotString = to_string(data.meta.steps);
// 	std::stringstream ss;
// 	ss << std::setfill('0') << std::setw(5) << snapShotString;
// 	snapShotString = ss.str();
	
// 	runname = dirname + "/" + runname;

// 	vector<double> Xexpanding(opt.grid[1]);
// 	vector<double> Yexpanding(opt.grid[2]);
// 	double b_x = opt.min_x * opt.stateInformation[0];
// 	double b_y = opt.min_y * opt.stateInformation[1];
// 	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
// 	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2];
// 	for(int i = 0; i < opt.grid[1]; i++){
// 		Xexpanding[i] = -b_x + h_x * i;
// 	}
// 	for(int i = 0; i < opt.grid[2]; i++){
// 		Yexpanding[i] = -b_y + h_y * i;
// 	}

// 	vector<double> ranges(2);
// 	complex<double> tmp3 = complex<double>(opt.RTE_step * opt.n_it_RTE,0.0);
// 	ranges[0] = opt.min_x * real(sqrt(complex<double>(1.0,0.0)+opt.exp_factor*opt.dispersion_x*opt.dispersion_x*tmp3*tmp3));
// 	ranges[1] = opt.min_y * real(sqrt(complex<double>(1.0,0.0)+opt.exp_factor*opt.dispersion_y*opt.dispersion_y*tmp3*tmp3));
	

	// string plotname = "Control-Plot-" + snapShotString;
	// string title = "Density " + snapShotString + " " + to_string(data.meta.time);
	// plotDataToPngEigen(plotname,data.wavefunction[0],opt);

	// if(opt.runmode.compare(1,1,"1") == 0){
	// 	title = "Density " + snapShotString;
	// 	plotname = "ExpandingFrame-" + snapShotString;
	// 	plotWithExpandingFrame(plotname,title,PsiVec[0],ranges,Xexpanding,Yexpanding,opt);
	// }

	// plotname = "Spectrum-" + snapShotString;
	// title = "Spectrum " + snapShotString; 
	// plotSpectrum(plotname,title,totalResult);

	// plotname = "Radial-Density-" + snapShotString;
	// title = "Radial-Density " + snapShotString; 
	// plotRadialDensity(plotname,title,totalResult);

	// if(pres[0].vlist.size() >= 0){
	// 	plotname = "PairDistance" + snapShotString;
	// 	title = "PairDistance " + snapShotString;
	// 	plotPairDistance(plotname,title,pres[0]);
	// }

	// plotname = "Vortices-" + snapShotString;
	// title = "Vortices " + snapShotString;
	// plotVortexList(plotname,title,phase,pres[0],opt);	

	// plotname = "Density-" + snapShotString;
	// title = "Density " + snapShotString;
	// plotDataToPng(plotname,title,densityLocationMap[0],opt);

	// plotname = "Density-Axial-Distribution-Gradient-" + snapShotString;
	// title = "Density " + snapShotString;
	// plotVector(plotname,title,x_dist_grad,y_dist_grad,opt);

	// plotname = "Angular-Dens-" + snapShotString;
	// title = "Angular Density " + snapShotString;
	// plotVector(plotname,title,totalResult.angularDensity,opt);	

	// plotname = "Contour-" + snapShotString;
	// title = "Contour " + snapShotString + " " + to_string(opt.t_abs.real());
	// plotContour(plotname,title,PsiVec[0],contour[0],opt);


// }

// void Eval::saveData2DSlice(vector<ComplexGrid> &wavefctVec, Options & external_opt, int external_snapshot_time, string external_runname, int sliceNumber){
// 	runname = external_runname;
// 	opt = external_opt;
// 	snapshot_time = external_snapshot_time;
// 	PsiVec.resize(wavefctVec.size());
// 	#pragma omp parallel for
// 	for(int k = 0; k < wavefctVec.size(); k++){
// 		PsiVec[k] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
// 		for(int i = 0; i < opt.grid[1]; i++){
// 			for(int j = 0; j < opt.grid[2]; j++){
// 				PsiVec[k](0,i,j,0) = wavefctVec[k](0,i,j,sliceNumber) / complex<double>(opt.Ag,0.0);
// 			}
// 		}
// 	}
// 	convertFromDimensionless();
// }

// void Eval::saveData(MatrixXcd &wavefct,Options &external_opt,int external_snapshot_time,string external_runname){
// 	runname = external_runname;
// 	opt = external_opt;
// 	snapshot_time = external_snapshot_time;
// 	PsiVec.resize(1);

// 	PsiVec[0] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	
// 	for(int i = 0; i < opt.grid[1]; i++){
// 		for(int j = 0; j < opt.grid[2]; j++){		
// 			PsiVec[0](0,i,j,0) = wavefct(i,j) / complex<double>(opt.Ag,0.0);
// 		}
// 	}
// 	convertFromDimensionless();		
// }

// void Eval::saveDataFromEval(Options &external_opt,int &external_snapshot_time,string &external_runname,vector<Eval> &extEval){
// 	runname = external_runname;
// 	opt = external_opt;
// 	snapshot_time = external_snapshot_time;

// 	int numberOfSamples = extEval.size();


// 	totalResult = Observables(extEval[0].totalResult.number.size());
// 	contour.resize(numberOfSamples*opt.samplesize);
// 	pres.resize(numberOfSamples*opt.samplesize);
	
// 	for(int k = 0; k < numberOfSamples; k++){
// 		if(k == 0){
// 			totalResult = extEval[k].totalResult;
// 		}else{
// 			totalResult += extEval[k].totalResult;
// 		}
// 		for(int i = 0; i < opt.samplesize;i++){
// 			contour[i+k*opt.samplesize] = extEval[k].contour[i];
// 			pres[i+k*opt.samplesize].vlist = extEval[k].pres[i].vlist;
// 		}
// 	}
// 	totalResult /= numberOfSamples;
// 	CombinedEval();
// 	CombinedSpectrum();
// }

// void Eval::CombinedEval(){

// 	string dirname = "CombinedRunObservables";
//     struct stat st;
//     	if(stat(dirname.c_str(),&st) != 0){
//         mkdir(dirname.c_str(),0755);
//     }

// 	string filename = dirname + "/" + runname + "_Observables.dat";	
	
// 	struct stat buffer;
//   	if(stat (filename.c_str(), &buffer) != 0){
//   		ofstream datafile;
//   		datafile.open(filename.c_str(), ios::out | ios::app);
//   		datafile << std::left << std::setw(15) << "Timestep"
//   						 << std::setw(15) << "X_max"
//   						 << std::setw(15) << "Y_max"
//   						 << std::setw(15) << "D_max"
//   						 << std::setw(15) << "D_min"
//   						 << std::setw(15) << "Rx"
// 						 << std::setw(15) << "Ry"
//   						 << std::setw(15) << "D_max/D_min"
//   						 << std::setw(15) << "D_max Angle"
//   						 << std::setw(15) << "D_min Angle"
//   						 << std::setw(15) << "Ratio"
//   						 << std::setw(15) << "RatioAngle"
//   						 << std::setw(15) << "N"
//   						 << std::setw(15) << "V"
//   						 << std::setw(15) << "N/V"
//   						 << std::setw(15) << "E_kin"
//   				 << endl;
//   		datafile.close();
//   	} 

//   	ofstream datafile(filename.c_str(), std::ios_base::out | std::ios_base::app);
// 	// datafile.open;
// 	datafile << std::left << std::setw(15) << snapshot_time
// 					 << std::setw(15) << opt.min_x * opt.stateInformation[0]
// 					 << std::setw(15) << opt.min_y * opt.stateInformation[1]
//  					 << std::setw(15) << totalResult.r_max
//  					 << std::setw(15) << totalResult.r_min
//  					 << std::setw(15) << totalResult.Rx
//  					 << std::setw(15) << totalResult.Ry
//  					 << std::setw(15) << totalResult.r_max / totalResult.r_min  
//  					 << std::setw(15) << totalResult.r_max_phi
//  					 << std::setw(15) << totalResult.r_min_phi
//  					 << std::setw(15) << totalResult.aspectRatio 
//  					 << std::setw(15) << totalResult.aspectRatioAngle 
// 					 << std::setw(15) << totalResult.particle_count
// 					 << std::setw(15) << totalResult.volume
// 					 << std::setw(15) << totalResult.density
// 					 << std::setw(15) << totalResult.Ekin
// 			 << endl;
// 	datafile.close();


// 	filename = dirname + "/" + runname + "_Observables.csv";	
	
//   	if(stat (filename.c_str(), &buffer) != 0){
//   		ofstream datafile1;
//   		datafile1.open(filename.c_str(), ios::out | ios::app);
//   		datafile1 << std::left << "," << "Timestep"
//   						 << "," << "X_max"
//   						 << "," << "Y_max"
//   						 << "," << "D_max"
//   						 << "," << "D_min"
//   						 << "," << "Rx"
// 						 << "," << "Ry"
//   						 << "," << "D_max/D_min"
//   						 << "," << "D_max Angle"
//   						 << "," << "D_min Angle"
//   						 << "," << "Ratio"
//   						 << "," << "RatioAngle"
//   						 << "," << "N"
//   						 << "," << "V"
//   						 << "," << "N/V"
//   						 << "," << "E_kin"
//   				 << endl;
//   		datafile1.close();
//   	} 

//   	ofstream datafile1(filename.c_str(), std::ios_base::out | std::ios_base::app);
// 	// datafile.open;
// 	datafile1 << std::left << "," << snapshot_time
// 					 << "," << opt.min_x * opt.stateInformation[0]
// 					 << "," << opt.min_y * opt.stateInformation[1]
//  					 << "," << totalResult.r_max
//  					 << "," << totalResult.r_min
//  					 << "," << totalResult.Rx
// 					 << "," << totalResult.Ry
//  					 << "," << totalResult.r_max / totalResult.r_min  
//  					 << "," << totalResult.r_max_phi
//  					 << "," << totalResult.r_min_phi
//  					 << "," << totalResult.aspectRatio 
//  					 << "," << totalResult.aspectRatioAngle 
// 					 << "," << totalResult.particle_count
// 					 << "," << totalResult.volume
// 					 << "," << totalResult.density
// 					 << "," << totalResult.Ekin
// 			 << endl;
// 	datafile1.close();

// 	filename = dirname + "/" + runname + "_Ratios.csv";	
	
//   	if(stat (filename.c_str(), &buffer) != 0){
//   		ofstream datafile2;
//   		datafile2.open(filename.c_str(), ios::out | ios::app);
//   		datafile2 << std::left << "," << "Timestep";
//   		for(int i = 0; i < totalResult.fixedAspectRatio.size(); i++){
//   			datafile2 << "," << std::left << i;
//   		}
//   		datafile2 << endl;
//   		datafile2.close();
//   	} 

//   	ofstream datafile2(filename.c_str(), std::ios_base::out | std::ios_base::app);
// 	// datafile2.open;
// 	datafile2 << std::left << "," << snapshot_time;
//   	for(int i = 0; i < totalResult.fixedAspectRatio.size(); i++){
//   		datafile2 << "," << totalResult.fixedAspectRatio(i);
//   	}
//   	datafile2 << endl;
// 	datafile2.close();
// }
	
	