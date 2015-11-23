#include <evaluation.h>

#define OBSERVABLES_DATA_POINTS_SIZE data->meta.grid[0]*data->meta.grid[1]
#define ANGULAR_AVERAGING_LENGTH 12
#define NUMBER_OF_VORTICES 100
#define DENSITY_CHECK_DISTANCE 2
#define EDGE_RANGE_CHECK 0.85

using namespace std;
using namespace Eigen;


Eval::Eval(shared_ptr<MatrixData> d,Options o) : data(d),  opt(o) {
	data = d;
	// opt = o;
	// cout << opt.Ag << "  " << opt.OmegaG << endl;
	// cout << data->meta.Ag << "  " << data->meta.OmegaG << endl;
	toPhysicalUnits(opt);
	data->convertToPhysicalUnits();
	// cout << "convertToPhysicalUnits()" << endl;

}

Eval::~Eval() {
	toDimensionlessUnits(opt);
	data->convertToDimensionlessUnits();
	// cout << "convertToDimensionlessUnits()" << endl;
}

void Eval::process(){

	// pres.resize(data->wavefunction.size());

	vlist.resize(data->wavefunction.size());

	contour.resize(data->wavefunction.size());

	densityCounter.resize(data->wavefunction.size());

	// densityLocationMap.resize(data->wavefunction.size());
	densityCoordinates.resize(data->wavefunction.size());

	phase = MatrixXd::Zero(data->meta.grid[0],data->meta.grid[1]);
	density = MatrixXd::Zero(data->meta.grid[0],data->meta.grid[1]);

	Contour tracker(data->meta);

	totalResult = Observables(OBSERVABLES_DATA_POINTS_SIZE);

	cout << currentTime() <<  " Step: " << data->meta.steps << " Time : " << data->meta.time << " s " << endl;		 


	// getDensity();
	for(int k = 0; k < data->wavefunction.size(); k++){
		calc_fields(data->wavefunction[k],opt);
		// contour[k] = tracker.trackContour(densityLocationMap[k]);
		totalResult += calculator(data->wavefunction[k],k);
		getVortices(data->wavefunction[k],densityCoordinates[k],vlist[k]);
		// cout << "-getVortexDistance" ;
	}	
	totalResult /= data->wavefunction.size();

}

void Eval::save(){

	string dirname = "runObservables";
    struct stat st;
    	if(stat(dirname.c_str(),&st) != 0){
        // mkdir(dirname.c_str(),0755); // LINUX 
        mkdir(dirname.c_str()); // WINDOWS
    }

	string filename = dirname + "/" + opt.runmode + "_Observables.dat";	

	int setwidth = 15;
	int setprec = 9;
	
	struct stat buffer;   
  	if(stat (filename.c_str(), &buffer) != 0){
  		ofstream datafile;
  		datafile.open(filename.c_str(), ios::out | ios::app);
  		datafile << std::left << std::setw(setwidth) << std::setprecision(setprec) << "Timestep" << ","
  						 	  << std::setw(setwidth) << std::setprecision(setprec) << "Time" << ","
  						 	  << std::setw(setwidth) << std::setprecision(setprec) << "X_max" << ","
  						 	  << std::setw(setwidth) << std::setprecision(setprec) << "Y_max" << ","
					 	 	  << std::setw(setwidth) << std::setprecision(setprec) << "Vortexnumber" << ","
					 	 	  << std::setw(setwidth) << std::setprecision(setprec) << "Alpha" << ","
  						 	  << std::setw(setwidth) << std::setprecision(setprec) << "Rx" << ","
						 	  << std::setw(setwidth) << std::setprecision(setprec) << "Ry" << ","
						 	  << std::setw(setwidth) << std::setprecision(setprec) << "Rx/Ry" << ","
  						 	  << std::setw(setwidth) << std::setprecision(setprec) << "E_Major" << ","
  						 	  << std::setw(setwidth) << std::setprecision(setprec) << "E_Minor" << ","
  						 	  << std::setw(setwidth) << std::setprecision(setprec) << "E_Major_Angle" << ","
  						 	  << std::setw(setwidth) << std::setprecision(setprec) << "E_Minor_Angle" << ","
  						 	  << std::setw(setwidth) << std::setprecision(setprec) << "E_Ratio" << ","
  						 	  << std::setw(setwidth) << std::setprecision(setprec) << "N" << ","
  						 	  << std::setw(setwidth) << std::setprecision(setprec) << "V" << ","
  						 	  << std::setw(setwidth) << std::setprecision(setprec) << "N/V" << ","
  						 	  << std::setw(setwidth) << std::setprecision(setprec) << "E_kin" << ","
  						 	  << std::setw(setwidth) << std::setprecision(setprec) << "n0"
  				 << endl;
  		datafile.close();
  	}

  	// Find n0 from fit, legacy code 

  	// double last_largest_element = 0;
   //  std::queue<double> largest_elements;
   //  for(int i = 0; i < data->meta.grid[0]; ++i){
   //  	for(int j = 0; j < data->meta.grid[1]; ++j){
   //  		const double val = density(i,j);
   //  		if(val > last_largest_element){
   //  			largest_elements.push(val);
   //  			last_largest_element = val;
   //  		}
   //  		if(largest_elements.size() > /*meta.grid[0]*meta.grid[1]*0.001*/ 10){
   //  			largest_elements.pop();
   //  		}
   //  	}
   //  }
   //  double n0 = 0;
   //  int le_size = largest_elements.size();
   //  while(!largest_elements.empty()){
   //  	n0 += largest_elements.front();
   //  	largest_elements.pop();
   //  }
    // cout << "n0 from totalresult = " << totalResult.n0 << endl;
    // cout << "n0 from calcu = " << n0 / le_size << endl;

    // END find n0

  	// double n0 = 2.0 * (totalResult.particle_count / M_PI) * (1 / (totalResult.r_max * totalResult.r_min)); 

  	ofstream datafile(filename.c_str(), std::ios_base::out | std::ios_base::app);
	// datafile.open;
	datafile << std::left << std::setw(setwidth) << std::setprecision(setprec) << data->meta.steps << ","
					  	  << std::setw(setwidth) << std::setprecision(setprec) << data->meta.time << ","
					  	  << std::setw(setwidth) << std::setprecision(setprec) << data->meta.coord[0] << ","
					  	  << std::setw(setwidth) << std::setprecision(setprec) << data->meta.coord[1] << ","
					  	  << std::setw(setwidth) << std::setprecision(setprec) << opt.vortexnumber << ","
					  	  << std::setw(setwidth) << std::setprecision(setprec) << totalResult.alpha << ","
 					  	  << std::setw(setwidth) << std::setprecision(setprec) << totalResult.Rx << ","
 					  	  << std::setw(setwidth) << std::setprecision(setprec) << totalResult.Ry << ","
 					  	  << std::setw(setwidth) << std::setprecision(setprec) << totalResult.Rx / totalResult.Ry << ","
 					  	  << std::setw(setwidth) << std::setprecision(setprec) << totalResult.r_max << ","
 					  	  << std::setw(setwidth) << std::setprecision(setprec) << totalResult.r_min << ","
 					  	  << std::setw(setwidth) << std::setprecision(setprec) << totalResult.r_max_phi << ","
 					  	  << std::setw(setwidth) << std::setprecision(setprec) << totalResult.r_min_phi << ","
 					  	  << std::setw(setwidth) << std::setprecision(setprec) << totalResult.aspectRatio  << ","
					  	  << std::setw(setwidth) << std::setprecision(setprec) << totalResult.particle_count << ","
					  	  << std::setw(setwidth) << std::setprecision(setprec) << totalResult.volume << ","
					  	  << std::setw(setwidth) << std::setprecision(setprec) << totalResult.density << ","
					  	  << std::setw(setwidth) << std::setprecision(setprec) << totalResult.Ekin << ","
					  	  << std::setw(setwidth) << std::setprecision(setprec) << totalResult.n0
			 << endl;
	datafile.close();
	cout << "\t Saved results to file." << endl;
}



void Eval::getVortices(MatrixXcd &DATA, vector<Coordinate<int32_t>> &densityCoordinates,list<VortexData> &vlist){
	
	

	vlist.clear();
	findVortices(densityCoordinates,vlist);

	// CHECK DENSITY AROUND VORTEX


	// erase vortices to close together
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
				radius.push_back(sqrt(x_shift*x_shift * data->meta.spacing[0]*data->meta.spacing[0] + y_shift*y_shift *data->meta.spacing[1]*data->meta.spacing[1]));
				polarDensity.push_back(abs2(data->wavefunction[0](x,y)));
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

	Vector<int32_t> down = Vector<int32_t>(0,-1,0,data->meta.grid[0],data->meta.grid[1],1);
	Vector<int32_t> right = Vector<int32_t>(1,0,0,data->meta.grid[0],data->meta.grid[1],1);
	Vector<int32_t> up = Vector<int32_t>(0,1,0,data->meta.grid[0],data->meta.grid[1],1);
	Vector<int32_t> left = Vector<int32_t>(-1,0,0,data->meta.grid[0],data->meta.grid[1],1);
	Vector<int32_t> rightdown = Vector<int32_t>(0.5, -0.5, 0,data->meta.grid[0],data->meta.grid[1],1);

	// #pragma omp parallel for
	for(int i = 0; i < densityCoordinates.size(); i++){
	// for(vector<Coordinate<int32_t>>::const_iterator it = densityCoordinates.begin(); it != densityCoordinates.end(); ++it){

		Coordinate<int32_t> c = densityCoordinates[i];

		if(c.x() > 1 && c.x() < data->meta.grid[0]-1){
			if(c.y() > 1 && c.y() < data->meta.grid[1]-1){
				int phase_winding = get_phase_jump(c, down) + get_phase_jump(c+down, right) + get_phase_jump(c+down+right, up) + get_phase_jump(c+right, left);

				if(phase_winding != 0){
					vortex.n = phase_winding;
					vortex.c = c;
					vlist.push_back(vortex);
				}
			}
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
	for(int x = 0; x < data->meta.grid[0]; x++)
	{
		for(int y = 0; y < data->meta.grid[1]; y++)
		{
			phase(x,y) = arg(DATA(x,y));
			density(x,y) = abs2(DATA(x,y));
		}
	}
}

int Eval::checkSum(MatrixXi &d,int &i, int &j){
	int radius = DENSITY_CHECK_DISTANCE;

	int sum = 0;
	int counter = 0;
	for(int x = 0; x <= 2 * radius; x++){
		for(int y = 0; y <= 2 * radius; y++){
			int k = x - radius;
			int l = y - radius;
			int distance = sqrt(k*k + l*l);
			if(distance < radius){
				int m = i + k;
				int n = j + l;
				sum += d(m,n);
				counter++;				
			}

		}
	}
	if(sum != 0){				
		if(counter == sum){
			return 1;
		}
		return 2;		
	}
	return 0;
}

void Eval::erosion(MatrixXi &d){
	MatrixXi tmp_map = MatrixXi::Zero(data->meta.grid[0],data->meta.grid[1]);

	for(int i = DENSITY_CHECK_DISTANCE; i < data->meta.grid[0] - DENSITY_CHECK_DISTANCE; i++){
	 	for(int j = DENSITY_CHECK_DISTANCE; j < data->meta.grid[1] - DENSITY_CHECK_DISTANCE; j++){
			if(checkSum(d,i,j) == 1){
				tmp_map(i,j) = 1;
			} 
		}
	}
	d = tmp_map;

}

void Eval::dilation(MatrixXi &d){
	MatrixXi tmp_map = MatrixXi::Zero(data->meta.grid[0],data->meta.grid[1]);

	for(int i = DENSITY_CHECK_DISTANCE; i < data->meta.grid[0] - DENSITY_CHECK_DISTANCE; i++){
	 	for(int j = DENSITY_CHECK_DISTANCE; j < data->meta.grid[1] - DENSITY_CHECK_DISTANCE; j++){
			if(checkSum(d,i,j) >= 1){
				tmp_map(i,j) = 1;
			} 
		}
	}
	d = tmp_map;
}



void Eval::floodFill(MatrixXi &dens){

	typedef struct {
			int32_t x;
			int32_t y;
	} Coord;

	std::stack<Coord> cStack;
	cStack.push(Coord{0,0});
	while(!cStack.empty()){
		Coord now = cStack.top();
		if(now.x < 0 || now.x >= data->meta.grid[0] || now.y < 0 || now.y >= data->meta.grid[1])
			cStack.pop();
		else if(dens(now.x,now.y) != 0)
			cStack.pop();
		else {
			dens(now.x,now.y) = 2;
			cStack.push(Coord{now.x+1,now.y});
			cStack.push(Coord{now.x-1,now.y});
			cStack.push(Coord{now.x,now.y+1});
			cStack.push(Coord{now.x,now.y-1});
		}
	}
}

void Eval::fillHoles(MatrixXi &dens){
	for(int i = 0; i < data->meta.grid[0]; i++){
		for(int j = 0; j < data->meta.grid[1]; j++){

			if(dens(i,j) != 2){
				dens(i,j) = 1;
			} else {
				dens(i,j) = 0;
			}

		}
	}
}

void Eval::smooth(MatrixXd &dens){
	MatrixXd tmp = dens;
	for(int i = DENSITY_CHECK_DISTANCE; i < data->meta.grid[0] - DENSITY_CHECK_DISTANCE; i++){
		for(int j = DENSITY_CHECK_DISTANCE; j < data->meta.grid[1] - DENSITY_CHECK_DISTANCE; j++){
			tmp(i,j) = dens(i-1,j) + dens(i+1,j) + dens(i,j-1) + dens(i,j+1) + dens(i,j);
			tmp(i,j) += dens(i-2,j) + dens(i+2,j) + dens(i,j-2) + dens(i,j+2);
			tmp(i,j) += dens(i-1,j-1) + dens(i+1,j+1) + dens(i+1,j-1) + dens(i-1,j+1);
		}
	}
	dens = tmp / 13.0;
}


vector<double> Eval::fitTF()
{
	lmfitter fit(density,data->meta);
	return fit.optimize();
}




void Eval::getDensity(){

	// smooth(density);

	double maximum = 0;
	for(int k = 0; k < data->wavefunction.size(); k++){
		for(int i = 0; i < data->meta.grid[0]; i++){
			for(int j = 0; j < data->meta.grid[1]; j++){
				double value = density(i,j);
				maximum = ( value > maximum) ? value : maximum;
			}
		}
	}
	double threshold = maximum * 0.01;

	for(int k = 0; k < data->wavefunction.size(); k++){

		
		densityLocationMap[k] = MatrixXi::Zero(data->meta.grid[0],data->meta.grid[1]);	
		densityCounter[k] = 0;
		densityCoordinates[k].clear();
		for(int i = DENSITY_CHECK_DISTANCE; i < data->meta.grid[0] - DENSITY_CHECK_DISTANCE; i++){
			for(int j = DENSITY_CHECK_DISTANCE; j < data->meta.grid[1] - DENSITY_CHECK_DISTANCE; j++){
			    if(density(i,j) > threshold){
	   				densityLocationMap[k](i,j) = 1;
				}
			}
		}
		
		floodFill(densityLocationMap[k]);
		fillHoles(densityLocationMap[k]);
		erosion(densityLocationMap[k]);
		dilation(densityLocationMap[k]);

		for(int i = DENSITY_CHECK_DISTANCE; i < data->meta.grid[0] - DENSITY_CHECK_DISTANCE; i++){
			for(int j = DENSITY_CHECK_DISTANCE; j < data->meta.grid[1] - DENSITY_CHECK_DISTANCE; j++){
			    if(densityLocationMap[k](i,j) == 1){
	   				Coordinate<int32_t> tmpCoord = Coordinate<int32_t>(i,j,0,data->meta.grid[0],data->meta.grid[1],1);
					densityCoordinates[k].push_back(tmpCoord);
					densityCounter[k]++;
				}
			}
		}


	}
}

vector<int> Eval::findMajorMinor(){
	vector<double> p = polarDensity();
	vector<double> halfaxis(360);
	int size = 40;
	int middle = (size-1) / 2;
	for(int i = 0; i < 360; ++i){
		for(int k = 0; k < size; ++k){
			int j = i - middle + k;
			if(j < 0) j += 360;
			if(j >= 360) j -= 360;
			halfaxis[i] += p[j];
		}
	}
	vector<double> halfaxis2(360);
	for(int i = 0; i < 360; ++i){
			int j = i + 180;
			if(j < 0) j += 360;
			if(j >= 360) j -= 360;
			halfaxis2[i] += halfaxis[j];
	}
	halfaxis = halfaxis2;
	vector<int> axis(2);
	vector<double>::iterator maxAngle = std::max_element(halfaxis.begin(), halfaxis.end());
	axis[0] = std::distance(halfaxis.begin(), maxAngle);

	vector<double>::iterator minangle = std::min_element(halfaxis.begin(), halfaxis.end());
	axis[1] = std::distance(halfaxis.begin(), minangle);
	return axis;
}

vector<double> Eval::polarDensity(){
	vector<double> pDensity(360);
	
	for(int x = 0; x < data->meta.grid[0]; ++x){
		for(int y = 0; y < data->meta.grid[1]; ++y){
			int x_shift = x - data->meta.grid[0]/2;
			int y_shift = y - data->meta.grid[1]/2;
			int i;
			if(x_shift != 0 && y_shift != 0)
				i = round(atan2(y_shift,x_shift) * 180 / M_PI);
			i = (i >= 0 ) ? i % 360 : ( 360 - abs ( i%360 ) ) % 360;
			pDensity[i]	+= densityLocationMap[0](x,y);
		}
	}
	// moving average
	vector<double> av(360);
	for(int i = 0; i < 360; ++i){
		int size = 5;
		int middle = (size-1)/2;
		double sum = 0;
		for(int j = 0; j < size; j++){
			int k = i + j - middle;
			if(k < 0) k += 360;
			if(k >= 360) k -= 360;
			sum += pDensity[k];
		}
		av[i] = sum / size;
	}
	return av;
}

Ellipse Eval::fitEllipse(c_set &c_data){

	// http://www.r-bloggers.com/fitting-an-ellipse-to-point-data/

	int numPoints = c_data.size();
	MatrixXd D1(numPoints,3);
	MatrixXd D2(numPoints,3);
	Matrix<double, 3, 3> S1;
	Matrix<double, 3, 3> S2;
	Matrix<double, 3, 3> S3;
	Matrix<double, 3, 3> T;
	Matrix<double, 3, 3> M;
	Matrix<double, 3, 3> C1;
	Matrix<double, 3, 1> a1;
	Matrix<double, 3, 1> a2;
	Matrix<double, 6, 1> f;
	MatrixXd temp;

	C1(0,2) = 0.5;
	C1(1,1) = -1.0;
	C1(2,0) = 0.5;

	// cerr << endl << "C1 " << C1 << endl;


	int i = 0;
	for(c_set::iterator it = c_data.begin(); it != c_data.end(); ++it){
	// for(int i = 0; i < numPoints; ++i){
		int x = it->x();
		int y = it->y();
		D1(i,0) = x * x;
		D1(i,1) = x * y;
		D1(i,2) = y * y;

		D2(i,0) = x;
		D2(i,1) = y;
		D2(i,2) = 1;
		i++;
	}
	S1 = D1.transpose() * D1;
	S2 = D1.transpose() * D2;
	S3 = D2.transpose() * D2;

	T = -1.0 * S3.inverse() * S2.transpose();

	M = S1 + S2 * T;
	M = C1 * M;

	EigenSolver<MatrixXd> es(M);
	for(int i = 0; i < 3; ++i){
		Vector3cd v = es.eigenvectors().col(i);
		complex<double> condition = v(1) * v(1) - 4.0 * v(0) * v(2);
		if(condition.imag() == 0 && condition.real() < 0){
			a1 = v.real();
		}
		// cerr << condition << endl;
		// cerr << endl << "A1 = " << a1 << endl;
	}

	a2 = T * a1;
	f(0) = a1(0);
	f(1) = a1(1);
	f(2) = a1(2);
	f(3) = a2(0);
	f(4) = a2(1);
	f(5) = a2(2);

	// if(f(0) < 0 && f(1) > 0 && f(2) < 0 && f(3) > 0 && f(4) > 0 && f(5) < 0 ){
	// 	f = -f;
	// }

	Matrix<double,2,2> A;
	A(0,0) = 2*f(0);
	A(0,1) = f(1);
	A(1,0) = f(1);
	A(1,1) = 2*f(2);
  	Matrix<double,2,1> b(-f(3), -f(4));
    Matrix<double,2,1> soln = A.inverse() * b;

  	double b2 = f(1) * f(1) / 4;
  
  	double center[2] = {soln(0), soln(1)}; 

  	double cond = f(0) * center[0] * center[0] + f(1) * center[0] * center[1] + f(2) * center[1] * center[1] + f(3) * center[0] + f(4) * center[1] + f(5);
  	if(cond >= 0.0){
  		f = -f;
		A(0,0) = 2*f(0);
		A(0,1) = f(1);
		A(1,0) = f(1);
		A(1,1) = 2*f(2);
	  	b = Matrix<double,2,1>(-f(3), -f(4));
	    soln = A.inverse() * b;	
  		b2 = f(1) * f(1) / 4; 	
  		center[0] = soln(0);
  		center[1] = soln(1);
  	}
  
    double num = 2 * (f(0) * f(4) * f(4) / 4 + f(2) * f(3) * f(3) / 4 + f(5) * b2 - f(1)*f(3)*f(4)/4 - f(0)*f(2)*f(5)); 
    double den1 = (b2 - f(0)*f(2));
    double den2 = sqrt((f(0) - f(2)) * (f(0) - f(2)) + 4*b2) ;
    double den3 = f(0) + f(2) ;
  
  	double axes[2] = {sqrt( num / (den1 * (den2 - den3))), sqrt( num / (den1 * (-den2 - den3)))};
  
  // calculate the angle of rotation 
    // double term = (f(0) - f(2)) / f(1);
    double y = /*fabs*/(f(1));
    double x = /*fabs*/(f(0) - f(2));

    // double angle = atan2( y , x ) / 2 + M_PI/2;
    double angle = atan(y/x)/2.0;
    // cerr << "                                 ";
    // cerr << " cond = " << cond;
    // cerr << " y = " << y << " x = " << x;
    // cerr << " a = " << angle * 180/M_PI;
    // if(angle < 0) {angle += (M_PI/2);}
    // if(x > 0 && y > 0) {angle += (M_PI/2); cerr << " 90 >0";}
    // cerr << " a = " << angle * 180/M_PI << endl;
    // else if(angle >= (M_PI / 2)) {angle -= M_PI;  cerr << " -180 ";}
    // cerr << " =  " << angle * 180/M_PI << endl;;
  // inline double vortex(int b, int y, int a, int x) //Vortex with phase [0,2*pi)          
  // {
  //         if(atan2(b-y,a-x)<0){ return 2*M_PI+atan2(b-y,a-x); } //atan2 is defined from [-pi,pi) so it needs to be changed to [0,2*pi)
  //     else{ return atan2(b-y,a-x); }        
  // }
  	
  	Ellipse eFit;
  // list(coef=f, center = center, major = max(semi.axes), minor = min(semi.axes), angle = unname(angle)) 
    eFit.coef = f;
    eFit.center.resize(2);
    eFit.center[0] = center[0];
    eFit.center[1] = center[1];
    eFit.major = (axes[0] > axes[1]) ? axes[0] : axes[1];
    eFit.minor = (axes[0] < axes[1]) ? axes[0] : axes[1];
    eFit.angle = angle;
    return eFit;

}

c_set Eval::generateContour(Ellipse &ellipse){
	c_set tmp;
	double t = 0;
	double delta = 0.0001;
	while(t < 2 * M_PI){
		int32_t x = ellipse.center[0] + ellipse.major * cos(t) * cos(ellipse.angle) - ellipse.minor * sin(t) * sin(ellipse.angle);
		int32_t y = ellipse.center[0] + ellipse.major * cos(t) * sin(ellipse.angle) + ellipse.minor * sin(t) * cos(ellipse.angle);	
		Coordinate<int32_t> c = Coordinate<int32_t>(x,y,0,data->meta.grid[0],data->meta.grid[1],1);
		pair<c_set::iterator,bool> test = tmp.insert(c);
		if(test.second == false){
			t += delta;
		} else {
			t += delta;
		}
	}
	return tmp;
}

c_set Eval::generateContour(vector<double>& params_tf){
	c_set tmp;
	double t = 0;
	double delta = 0.0001;
	while(t < 2 * M_PI){
		int32_t x = data->meta.grid[0]/2 + params_tf[1]/data->meta.spacing[0] * cos(t) * cos(params_tf[2]) - params_tf[3]/data->meta.spacing[0] * sin(t) * sin(params_tf[2]);
		int32_t y = data->meta.grid[1]/2 + params_tf[1]/data->meta.spacing[1] * cos(t) * sin(params_tf[2]) + params_tf[3]/data->meta.spacing[1] * sin(t) * cos(params_tf[2]);	
		Coordinate<int32_t> c = Coordinate<int32_t>(x,y,0,data->meta.grid[0],data->meta.grid[1],1);
		pair<c_set::iterator,bool> test = tmp.insert(c);
		if(test.second == false){
			t += delta;
		} else {
			t += delta;
		}
	}
	return tmp;
}

vector<Coordinate<int32_t>> Eval::generate_density_coordinates(vector<double>& params_tf)
{	
	std::vector<Coordinate<int32_t>> v;
	for(int i = 0; i < data->meta.grid[0]; i++){
		for(int j = 0; j < data->meta.grid[1]; j++){
			double i0 = -data->meta.coord[0] + data->meta.spacing[0] * i;
			double i1 = -data->meta.coord[1] + data->meta.spacing[1] * j;
			double value = params_tf[0] * (1 - (i0*i0)/(params_tf[1]*params_tf[1]) - (i1*i1)/(params_tf[3]*params_tf[3]) - params_tf[2] * i0 * i1);
			if(value >= 0.0){
				v.push_back(Coordinate<int32_t>(i,j,0,data->meta.grid[0],data->meta.grid[1],1));				
			}
		}
	}
	return v;
}

void Eval::aspectRatio(Observables &obs, int &sampleindex){

	double h_x = data->meta.spacing[0];
	double h_y = data->meta.spacing[1];



	ellipse = fitEllipse(contour[sampleindex]);
	// c_set cEllipse = generateContour(ellipse);

	// FROM ELLIPSE FIT 
	// obs.r_max = ellipse.major * h_x;
	// obs.r_min = ellipse.minor * h_y;

	// obs.r_max_phi = ellipse.angle * 180/M_PI;
	// obs.r_min_phi = ellipse.angle * 180/M_PI + 90;

	// obs.aspectRatio = obs.r_max / obs.r_min;

	// END ELLIPSE FIT

	// FROM THOMAS FERMI FIT

	vector<double> params_tf = fitTF();

	// c_set cEllipse = 
	contour[sampleindex] = generateContour(params_tf);
	densityCoordinates[sampleindex] = generate_density_coordinates(params_tf);
	densityCounter[sampleindex] = densityCoordinates[sampleindex].size();

	obs.n0 = params_tf[0];
	obs.r_max = params_tf[1];
	obs.r_min = params_tf[3];

	double tmp1 = - params_tf[2] * params_tf[1] * params_tf[1] * params_tf[3] * params_tf[3];
    double tmp2 = (params_tf[1] * params_tf[1] - params_tf[3] * params_tf[3]);
    // double at = atan2(tmp1,tmp2);
    double at = atan(tmp1/tmp2);
    // double at = atan(tmp1/tmp2);
    at *= ( 180 / M_PI ) / 2.0;
	obs.r_max_phi = at;
	obs.r_min_phi = at + 90;
	obs.aspectRatio = obs.r_max / obs.r_min;



	// END TF FIT

	// contour[sampleindex] = cEllipse;


	vector<contourData> cData;
	for(c_set::iterator it =  contour[sampleindex].begin(); it !=  contour[sampleindex].end(); ++it){
		contourData tmp;
		tmp.c = *it;
		double x = (tmp.c.x() - data->meta.grid[0]/2) * h_x;
		double y = (tmp.c.y() - data->meta.grid[1]/2) * h_y;
		tmp.phi = atan2(y,x) * 180 / M_PI + 180;
		tmp.r = sqrt(x*x  + y*y);
		cData.push_back(tmp);
	}

	std::sort(cData.begin(),cData.end(),[](const contourData &lhs, const contourData &rhs) -> bool {return (lhs.phi < rhs.phi);});

	vector<double> cRadius(361);
	vector<int> divisor_counter(361);
	for(vector<contourData>::const_iterator it = cData.begin(); it != cData.end(); ++it){
		int index = round(it->phi);
		cRadius[index] += it->r;
		divisor_counter[index]++;
	}
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

	// vector<int> axis(2);
	// axis = findMajorMinor();
	// obs.Rx = ellipse.major * h_x;
	// obs.Ry = ellipse.minor * h_y;	

	obs.Rx = cRadius[0];
	obs.Ry = cRadius[90];

	// // moving median of radii
	// vector<double> tRadius(360);
	// for(int i = 0; i < 360; i++){
	// 	int size = 9;
	// 	int middle = (size-1)/2;
	// 	vector<double> median(size);
	// 	for(int j = 0; j < size; j++){
	// 		int k = i + j - middle;
	// 		if(k < 0) k += 360;
	// 		if(k >= 360) k -= 360;
	// 		median[j] = cRadius[k];
	// 	}
	// 	std::sort(median.begin(),median.end());
	// 	tRadius[i] = median[middle];
	// }
	// cRadius = tRadius;

	// vector<double>::iterator maxRadius = std::max_element(cRadius.begin(), cRadius.end());
	// obs.r_max = *maxRadius;
	// obs.r_max_phi = std::distance(cRadius.begin(), maxRadius);

	// vector<double>::iterator minRadius = std::min_element(cRadius.begin(), cRadius.end());
	// obs.r_min = *minRadius;
	// obs.r_min_phi = std::distance(cRadius.begin(), minRadius);

	vector<double> tmp_ratio(360);
	for(int i = 0; i < 360; i++){
		int j = (i+90 < 360) ? i+90 : i-270;
		if(cRadius[j] > 0.0){
			tmp_ratio[i] = cRadius[i] / cRadius[j];
			obs.fixedAspectRatio(i) = tmp_ratio[i];
		} else {
			cout << "WARNING: Aspect-Ratio: Calculated Distance smaller than zero!" << endl;
		}
	}
}

Observables Eval::calculator(MatrixXcd DATA,int sampleindex){
	
	Observables obs = Observables(OBSERVABLES_DATA_POINTS_SIZE);

	aspectRatio(obs,sampleindex);

	obs.volume = data->meta.spacing[0] * data->meta.spacing[1] * densityCounter[sampleindex];

	double sum = 0;
	#pragma omp parallel for reduction(+:sum)
	for(int i = 0; i < data->meta.grid[0]; i++){
	    for(int j = 0; j < data->meta.grid[1]; j++){	    	    		
	      	sum += abs2(DATA(i,j));
	    }
	}

	obs.particle_count = sum;
	obs.particle_count *= data->meta.spacing[0] * data->meta.spacing[1];
	obs.density = obs.particle_count / obs.volume;

	

	// == Angular Density
	// vector<double> angularDensity;
	vector<double> phi;
	vector<double> polarDensity; // first entry is 
	vector<double> radius;
	vector<Coordinate<int32_t>> cartesianCoordinates;

	ArrayXd divisor2(obs.number.size());

	double r_index_factor = (obs.radialDensity.size() -1) / sqrt(data->meta.coord[0] * data->meta.coord[0] + data->meta.coord[1] * data->meta.coord[1]);
	// cout << "INDEX_FACTOR " << r_index_factor << endl;
	for(int i = 0; i < data->meta.grid[0]; i++){
	    for(int j = 0; j < data->meta.grid[1]; j++){
				int x_shift = i - data->meta.grid[0]/2;
				int y_shift = j - data->meta.grid[1]/2;
				phi.push_back( atan2(x_shift * data->meta.spacing[0] ,y_shift * data->meta.spacing[1]) * 180 / M_PI + 180);
				double r_tmp = sqrt(x_shift*x_shift * data->meta.spacing[0]*data->meta.spacing[0] + y_shift*y_shift *data->meta.spacing[1]*data->meta.spacing[1]);

				

				radius.push_back(r_tmp);
				double dens_tmp = abs2(DATA(i,j));
				polarDensity.push_back(dens_tmp);
				Coordinate<int32_t> tmpCoord =  Coordinate<int32_t>(i,j,0,data->meta.grid[0],data->meta.grid[1],1);
				cartesianCoordinates.push_back(tmpCoord);

				int index = r_index_factor * r_tmp;
				// cout << " radius " << r_tmp << " / " << data->meta.coord[0] << " index " << index << " / " << obs.radialDensity.size() << endl;
				// if(index >= obs.number.size()){
				// 	cout << "index " << index << " r_tmp " << r_tmp << " " << data->meta.spacing[0] << " " << endl;
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

	data->fft.Forward(DATA);

	DATA /= sqrt(DATA.cols() * DATA.rows());
	
	ArrayXd divisor(obs.number.size());
	divisor.setZero();
	
	vector<vector<double>> kspace;
	vector<double> Kmax(2);
	Kmax[0] = M_PI / data->meta.spacing[0];
	Kmax[1] = M_PI / data->meta.spacing[1];
	vector<double> deltaK(2);
	deltaK[0] = Kmax[0] / (data->meta.grid[0] / 2.0);
	deltaK[1] = Kmax[1] / (data->meta.grid[1] / 2.0);


	kspace.resize(2);
	for(int d = 0; d < 2; d++){
		// set k-space
		kspace[d].resize(data->meta.grid[d]);
		for(int i = 0; i <= data->meta.grid[d]/2; i++){
			// kspace[d][i] = (M_PI / rmax[d]) * (double)i;
			kspace[d][i] = deltaK[d] * (double)i;
		}
		for(int i = (data->meta.grid[d]/2)+1; i < data->meta.grid[d]; i++){
			// kspace[d][i] = -(M_PI / rmax[d]) * (double)(data->meta.grid[d] - i);
			kspace[d][i] = - deltaK[d] * (double)(data->meta.grid[d] - i);
		}
	}


	// DEFINITION OF KLENGTH FROM THE INTERNET! |||| ==>>     2*pi*i/(Nx*dx)
	// double kmax[2];

	// for(int i = 0; i < 2; i++){
	// 	kmax[i] = *max_element(kspace[i].begin(), kspace[i].end());
	// }
	
	double kwidth2[2];

	for(int i = 0; i < 2; i++)
		kwidth2[i] = (data->meta.grid[i] == 1) ? 0 : kspace[i][data->meta.grid[i]/2] * kspace[i][data->meta.grid[i]/2];
	
	double index_factor = (obs.number.size() - 1) / sqrt(kwidth2[0] + kwidth2[1]);

	for(int x = 0; x < data->meta.grid[0]; x++){
		for (int y = 0; y < data->meta.grid[1]; y++){
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

	// vector<double> x_tmp(densityCoordinates[0].size());
	// vector<double> y_tmp(densityCoordinates[0].size());
	int upperVal = data->meta.grid[0]*0.80;
	int lowerVal = data->meta.grid[0]*0.20;


	for(int i = 0; i < densityCoordinates[0].size(); i++){
		// x_tmp[i] = densityCoordinates[0][i].x();
		// y_tmp[i] = densityCoordinates[0][i].y();
		if(densityCoordinates[0][i].x() > upperVal || densityCoordinates[0][i].x() < lowerVal)
			return true;
		if(densityCoordinates[0][i].y() > upperVal || densityCoordinates[0][i].y() < lowerVal)
			return true;
	}
	return false;
	// auto x_minmax = std::minmax_element (x_tmp.begin(),x_tmp.end());
	// auto y_minmax = std::minmax_element (y_tmp.begin(),y_tmp.end());

	// edges[0] = (int)*x_minmax.first;
	// edges[1] = (int)*y_minmax.first;

	// edges[2] = (int)*x_minmax.second - (int)*x_minmax.first;
	// edges[3] = (int)*y_minmax.second - (int)*y_minmax.first;

	// return (edges[2] >= (data->meta.grid[0] * EDGE_RANGE_CHECK) || edges[3] >= (data->meta.grid[1] * EDGE_RANGE_CHECK));
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


	