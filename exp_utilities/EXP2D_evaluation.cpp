#include <EXP2D_evaluation.h>

#define OBSERVABLES_DATA_POINTS_SIZE opt.grid[1]*opt.grid[2]*opt.grid[3]
#define ANGULAR_AVERAGING_LENGTH 12
#define NUMBER_OF_VORTICES 100
#define VORTEX_SURROUND_DENSITY_RADIUS 5
#define EDGE_RANGE_CHECK 10

using namespace std;
using namespace Eigen;


Eval::Eval() {};

Eval::~Eval() {};

void Eval::saveData(vector<MatrixXcd> &wavefctVec,Options &external_opt,int external_snapshot_time,string external_runname){
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
	convertFromDimensionless();

}

void Eval::convertFromDimensionless(){
	opt.min_x *= opt.Ag;
	opt.min_y *= opt.Ag;
	opt.t_abs /= opt.OmegaG;
	opt.t_abs *= 1000.0; // conversion to ms
	opt.omega_x /= 2.0 * M_PI / opt.OmegaG;
	opt.omega_y /= 2.0 * M_PI / opt.OmegaG;
	opt.dispersion_x /= 2.0 * M_PI / opt.OmegaG;
	opt.dispersion_y /= 2.0 * M_PI / opt.OmegaG;
}

void Eval::saveData2DSlice(vector<ComplexGrid> &wavefctVec, Options & external_opt, int external_snapshot_time, string external_runname, int sliceNumber){
	runname = external_runname;
	opt = external_opt;
	snapshot_time = external_snapshot_time;
	PsiVec.resize(wavefctVec.size());
	#pragma omp parallel for
	for(int k = 0; k < wavefctVec.size(); k++){
		PsiVec[k] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
		for(int i = 0; i < opt.grid[1]; i++){
			for(int j = 0; j < opt.grid[2]; j++){
				PsiVec[k](0,i,j,0) = wavefctVec[k](0,i,j,sliceNumber);
			}
		}
	}
}

void Eval::saveData(MatrixXcd &wavefct,Options &external_opt,int external_snapshot_time,string external_runname){
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
	convertFromDimensionless();		
}

void Eval::saveDataFromEval(Options &external_opt,int &external_snapshot_time,string &external_runname,vector<Eval> &extEval){
	runname = external_runname;
	opt = external_opt;
	snapshot_time = external_snapshot_time;

	int numberOfSamples = extEval.size();


	totalResult = Observables(extEval[0].totalResult.number.size());
	contour.resize(numberOfSamples*opt.samplesize);
	pres.resize(numberOfSamples*opt.samplesize);
	
	for(int k = 0; k < numberOfSamples; k++){
		if(k == 0){
			totalResult = extEval[k].totalResult;
		}else{
			totalResult += extEval[k].totalResult;
		}
		for(int i = 0; i < opt.samplesize;i++){
			contour[i+k*opt.samplesize] = extEval[k].contour[i];
			pres[i+k*opt.samplesize].vlist = extEval[k].pres[i].vlist;
		}
	}
	totalResult /= numberOfSamples;
	CombinedEval();
	CombinedSpectrum();
}

void Eval::CombinedEval(){

	// // ONLY NEEDED UNTIL fixedAspectRatio is saved in BinaryFile::appendEval
	// Observables tmpResult = Observables(OBSERVABLES_DATA_POINTS_SIZE);
	// for(int sampleindex = 0; sampleindex < contour.size(); sampleindex++){
	// 	Observables obs = Observables(OBSERVABLES_DATA_POINTS_SIZE);
	// 	aspectRatio(obs,sampleindex);
	// 	tmpResult.fixedAspectRatio += obs.fixedAspectRatio;
	// }
	// tmpResult.fixedAspectRatio /= contour.size();
	// totalResult.fixedAspectRatio = tmpResult.fixedAspectRatio;
	// // ONLY NEEDED ONE TIME


	string dirname = "CombinedRunObservables";
    struct stat st;
    	if(stat(dirname.c_str(),&st) != 0){
        mkdir(dirname.c_str(),0755);
    }

	string filename = dirname + "/" + runname + "_Observables.dat";	
	
	struct stat buffer;
  	if(stat (filename.c_str(), &buffer) != 0){
  		ofstream datafile;
  		datafile.open(filename.c_str(), ios::out | ios::app);
  		datafile << std::left << std::setw(15) << "Timestep"
  						 << std::setw(15) << "X_max"
  						 << std::setw(15) << "Y_max"
  						 << std::setw(15) << "D_max"
  						 << std::setw(15) << "D_min"
  						 << std::setw(15) << "Rx"
						 << std::setw(15) << "Ry"
  						 << std::setw(15) << "D_max/D_min"
  						 << std::setw(15) << "D_max Angle"
  						 << std::setw(15) << "D_min Angle"
  						 << std::setw(15) << "Ratio"
  						 << std::setw(15) << "RatioAngle"
  						 << std::setw(15) << "N"
  						 << std::setw(15) << "V"
  						 << std::setw(15) << "N/V"
  						 << std::setw(15) << "E_kin"
  				 << endl;
  		datafile.close();
  	} 

  	ofstream datafile(filename.c_str(), std::ios_base::out | std::ios_base::app);
	// datafile.open;
	datafile << std::left << std::setw(15) << snapshot_time
					 << std::setw(15) << opt.min_x * opt.stateInformation[0]
					 << std::setw(15) << opt.min_y * opt.stateInformation[1]
 					 << std::setw(15) << totalResult.r_max
 					 << std::setw(15) << totalResult.r_min
 					 << std::setw(15) << totalResult.Rx
 					 << std::setw(15) << totalResult.Ry
 					 << std::setw(15) << totalResult.r_max / totalResult.r_min  
 					 << std::setw(15) << totalResult.r_max_phi
 					 << std::setw(15) << totalResult.r_min_phi
 					 << std::setw(15) << totalResult.aspectRatio 
 					 << std::setw(15) << totalResult.aspectRatioAngle 
					 << std::setw(15) << totalResult.particle_count
					 << std::setw(15) << totalResult.volume
					 << std::setw(15) << totalResult.density
					 << std::setw(15) << totalResult.Ekin
			 << endl;
	datafile.close();


	filename = dirname + "/" + runname + "_Observables.csv";	
	
  	if(stat (filename.c_str(), &buffer) != 0){
  		ofstream datafile1;
  		datafile1.open(filename.c_str(), ios::out | ios::app);
  		datafile1 << std::left << "," << "Timestep"
  						 << "," << "X_max"
  						 << "," << "Y_max"
  						 << "," << "D_max"
  						 << "," << "D_min"
  						 << "," << "Rx"
						 << "," << "Ry"
  						 << "," << "D_max/D_min"
  						 << "," << "D_max Angle"
  						 << "," << "D_min Angle"
  						 << "," << "Ratio"
  						 << "," << "RatioAngle"
  						 << "," << "N"
  						 << "," << "V"
  						 << "," << "N/V"
  						 << "," << "E_kin"
  				 << endl;
  		datafile1.close();
  	} 

  	ofstream datafile1(filename.c_str(), std::ios_base::out | std::ios_base::app);
	// datafile.open;
	datafile1 << std::left << "," << snapshot_time
					 << "," << opt.min_x * opt.stateInformation[0]
					 << "," << opt.min_y * opt.stateInformation[1]
 					 << "," << totalResult.r_max
 					 << "," << totalResult.r_min
 					 << "," << totalResult.Rx
					 << "," << totalResult.Ry
 					 << "," << totalResult.r_max / totalResult.r_min  
 					 << "," << totalResult.r_max_phi
 					 << "," << totalResult.r_min_phi
 					 << "," << totalResult.aspectRatio 
 					 << "," << totalResult.aspectRatioAngle 
					 << "," << totalResult.particle_count
					 << "," << totalResult.volume
					 << "," << totalResult.density
					 << "," << totalResult.Ekin
			 << endl;
	datafile1.close();

	filename = dirname + "/" + runname + "_Ratios.csv";	
	
  	if(stat (filename.c_str(), &buffer) != 0){
  		ofstream datafile2;
  		datafile2.open(filename.c_str(), ios::out | ios::app);
  		datafile2 << std::left << "," << "Timestep";
  		for(int i = 0; i < totalResult.fixedAspectRatio.size(); i++){
  			datafile2 << "," << std::left << i;
  		}
  		datafile2 << endl;
  		datafile2.close();
  	} 

  	ofstream datafile2(filename.c_str(), std::ios_base::out | std::ios_base::app);
	// datafile2.open;
	datafile2 << std::left << "," << snapshot_time;
  	for(int i = 0; i < totalResult.fixedAspectRatio.size(); i++){
  		datafile2 << "," << totalResult.fixedAspectRatio(i);
  	}
  	datafile2 << endl;
	datafile2.close();
}

void Eval::evaluateData(){
	// cout << "evaluateData" << endl;	
	// cout << "checkEdges call: " << endl;
	// checkEdges();

	pres.resize(PsiVec.size());
	densityCoordinates.clear();
	contour.resize(PsiVec.size());
	densityCounter.resize(PsiVec.size());

	densityLocationMap.resize(PsiVec.size());
	densityCoordinates.resize(PsiVec.size());
	phase.resize(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	zeros.resize(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	Contour tracker(opt);


	totalResult = Observables(OBSERVABLES_DATA_POINTS_SIZE);
	
	cout << endl << "Evaluating sample #: ";
	for(int k = 0; k < PsiVec.size(); k++){
		cout << k << " " ;
		getDensity(PsiVec[k],densityLocationMap[k],densityCoordinates[k],densityCounter[k]);
		cout << "-getDensity" << endl;
		contour[k] = tracker.trackContour(densityLocationMap[k]);
		cout << "-trackContour" << endl;
		totalResult += calculator(PsiVec[k],k);
		cout << "-calculator" << endl;
		getVortices(PsiVec[k],densityCoordinates[k],pres[k]);
		cout << "-getVortices" << endl;
		// getVortexDistance(pres[k]);
		// cout << "-getVortexDistance" << endl;
	}	
	totalResult /= PsiVec.size();

	cout << endl;

	string dirname = "runObservables";
    struct stat st;
    	if(stat(dirname.c_str(),&st) != 0){
        mkdir(dirname.c_str(),0755);
    }

	string filename = dirname + "/" + runname + "_Observables.dat";	
	
	struct stat buffer;   
  	if(stat (filename.c_str(), &buffer) != 0){
  		ofstream datafile;
  		datafile.open(filename.c_str(), ios::out | ios::app);
  		datafile << std::left << std::setw(15) << "Timestep"
  						 << std::setw(15) << "Time"
  						 << std::setw(15) << "X_max"
  						 << std::setw(15) << "Y_max"
  						 << std::setw(15) << "D_max"
  						 << std::setw(15) << "D_min"
  						 << std::setw(15) << "Rx"
						 << std::setw(15) << "Ry"
  						 << std::setw(15) << "D_max/D_min"
  						 << std::setw(15) << "D_max Angle"
  						 << std::setw(15) << "D_min Angle"
  						 << std::setw(15) << "Ratio"
  						 << std::setw(15) << "RatioAngle"
  						 << std::setw(15) << "N"
  						 << std::setw(15) << "V"
  						 << std::setw(15) << "N/V"
  						 << std::setw(15) << "E_kin"
  				 << endl;
  		datafile.close();
  	} 

  	ofstream datafile(filename.c_str(), std::ios_base::out | std::ios_base::app);
	// datafile.open;
	datafile << std::left << std::setw(15) << snapshot_time
					 << std::setw(15) << opt.min_x * opt.stateInformation[0]
					 << std::setw(15) << opt.min_y * opt.stateInformation[1]
					 << std::setw(15) << opt.t_abs
 					 << std::setw(15) << totalResult.r_max
 					 << std::setw(15) << totalResult.r_min
 					 << std::setw(15) << totalResult.Rx
 					 << std::setw(15) << totalResult.Ry
 					 << std::setw(15) << totalResult.r_max / totalResult.r_min  
 					 << std::setw(15) << totalResult.r_max_phi
 					 << std::setw(15) << totalResult.r_min_phi
 					 << std::setw(15) << totalResult.aspectRatio 
 					 << std::setw(15) << totalResult.aspectRatioAngle 
					 << std::setw(15) << totalResult.particle_count
					 << std::setw(15) << totalResult.volume
					 << std::setw(15) << totalResult.density
					 << std::setw(15) << totalResult.Ekin
			 << endl;
	datafile.close();

}

void Eval::evaluateDataITP(){

	// pres.vlist.clear();
	densityCoordinates.clear();
	contour.resize(PsiVec.size());
	densityCounter.resize(PsiVec.size());

	densityLocationMap.resize(PsiVec.size());
	densityCoordinates.resize(PsiVec.size());
	phase.resize(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	zeros.resize(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

	totalResult = Observables(OBSERVABLES_DATA_POINTS_SIZE);
		
	// for(int k = 0; k < PsiVec.size(); k++){
		getDensity(PsiVec[0],densityLocationMap[0],densityCoordinates[0],densityCounter[0]);
		totalResult = calculatorITP(PsiVec[0],0);
	// }
	// totalResult /= PsiVec.size();
}

void Eval::CombinedSpectrum(){
	string dirname = "CombinedRunPlots";
    struct stat st;
    	if(stat(dirname.c_str(),&st) != 0){
        mkdir(dirname.c_str(),0755);
    }
    
	std::string snapShotString = to_string(snapshot_time);
	std::stringstream ss;
	ss << std::setfill('0') << std::setw(5) << snapShotString;
	snapShotString = ss.str();
	
	runname = dirname + "/" + runname;

	string plotname = runname + "-Spectrum-" + snapShotString;
	string title = "Spectrum " + snapShotString; 
	plotSpectrum(plotname,title, totalResult);

}

void Eval::plotData(){
	string dirname = "runPlots_" + runname;
    struct stat st;
    	if(stat(dirname.c_str(),&st) != 0){
        mkdir(dirname.c_str(),0755);
    }
    
	std::string snapShotString = to_string(snapshot_time);
	std::stringstream ss;
	ss << std::setfill('0') << std::setw(5) << snapShotString;
	snapShotString = ss.str();
	
	runname = dirname + "/" + runname;

	// vector<double> Xexpanding(opt.grid[1]);
	// vector<double> Yexpanding(opt.grid[2]);
	// double b_x = opt.min_x * opt.stateInformation[0];
	// double b_y = opt.min_y * opt.stateInformation[1];
	// double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
	// double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2];
	// for(int i = 0; i < opt.grid[1]; i++){
	// 	Xexpanding[i] = -b_x + h_x * i;
	// }
	// for(int i = 0; i < opt.grid[2]; i++){
	// 	Yexpanding[i] = -b_y + h_y * i;
	// }

	// vector<double> ranges(2);
	// complex<double> tmp3 = complex<double>(opt.RTE_step * opt.n_it_RTE,0.0);
	// ranges[0] = opt.min_x * real(sqrt(complex<double>(1.0,0.0)+opt.exp_factor*opt.dispersion_x*opt.dispersion_x*tmp3*tmp3));
	// ranges[1] = opt.min_y * real(sqrt(complex<double>(1.0,0.0)+opt.exp_factor*opt.dispersion_y*opt.dispersion_y*tmp3*tmp3));
	

	string plotname = runname + "-Control-Plot-" + snapShotString + " " + to_string(opt.t_abs.real());
	string title = "Density " + snapShotString;
	plotDataToPngExpanding(plotname,title,PsiVec[0],opt);

	// if(opt.runmode.compare(1,1,"1") == 0){
	// 	title = "Density " + snapShotString;
	// 	plotname = runname + "-ExpandingFrame-" + snapShotString;
	// 	plotWithExpandingFrame(plotname,title,PsiVec[0],ranges,Xexpanding,Yexpanding,opt);
	// }

	plotname = runname + "-Spectrum-" + snapShotString;
	title = "Spectrum " + snapShotString; 
	plotSpectrum(plotname,title,totalResult);

	plotname = runname + "-Radial-Density-" + snapShotString;
	title = "Radial-Density " + snapShotString; 
	plotRadialDensity(plotname,title,totalResult);

	// if(pres[0].vlist.size() >= 0){
	// 	plotname = runname + "-PairDistance" + snapShotString;
	// 	title = "PairDistance" + snapShotString;
	// 	plotPairDistance(plotname,title,pres[0]);
	// }

	plotname = runname + "-Vortices-" + snapShotString;
	title = "Vortices " + snapShotString;
	plotVortexList(plotname,title,phase,pres[0],opt);	

	// plotname = runname + "-Density-" + snapShotString;
	// title = "Density " + snapShotString;
	// plotDataToPng(plotname,title,densityLocationMap[0],opt);

	// plotname = runname + "-Density-Axial-Distribution-Gradient-" + snapShotString;
	// title = "Density " + snapShotString;
	// plotVector(plotname,title,x_dist_grad,y_dist_grad,opt);

	plotname = runname + "-Angular-Dens-" + snapShotString;
	title = "Angular Density " + snapShotString;
	plotVector(plotname,title,totalResult.angularDensity,opt);	

	plotname = runname + "-Contour-" + snapShotString+ " " + to_string(opt.t_abs.real());
	title = "Contour " + snapShotString;
	plotContour(plotname,title,PsiVec[0],contour[0],opt);


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

void Eval::getVortices(ComplexGrid &data, vector<Coordinate<int32_t>> &densityCoordinates,PathResults &pres){
	
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

int Eval::get_phase_jump(const Coordinate<int32_t> &c, const Vector<int32_t> &v, const RealGrid &phase){
	if(phase.at(0,c + v) + M_PI < phase.at(0,c))	// Phase ueberschreitet 0/2pi von unten
		return 1;
	else if(phase.at(0,c) + M_PI < phase.at(0,c + v))	// Phase ueberschreitet 0/2pi von oben
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
	// vector< vector< vector<bool > > > checked(phase.width(), vector< vector<bool> >(phase.height(), vector<bool>(phase.depth(),false)));	// Welche felder schon ueberprueft wurden
	VortexData vortex;
					// Charakteristika eines gefundenen Vortex

	// cout << "findVortices_before iterator" << endl;



	for(vector<Coordinate<int32_t>>::const_iterator it = densityCoordinates.begin(); it != densityCoordinates.end(); ++it){

	// for (int z = 0; z < phase.depth(); z++)
	// {
		// for (int x = 0; x < phase.width(); x++)
	// 	{
			// for (int y = 0; y < phase.height(); y++)
	// 		{

				Coordinate<int32_t> c = *it; // phase.make_coord(x,y,z);
				Vector<int32_t> down = phase.make_vector(0,-1,0);
				Vector<int32_t> right = phase.make_vector(1,0,0);
				Vector<int32_t> up = phase.make_vector(0,1,0);
				Vector<int32_t> left = phase.make_vector(-1,0,0);
					
				int phase_winding = get_phase_jump(c, down, phase) + get_phase_jump(c+down, right, phase) + get_phase_jump(c+down+right, up, phase) + get_phase_jump(c+right, left, phase);
				// int mass_zeros = zeros.at(0,c) + zeros.at(0,c+down) + zeros.at(0,c+down+right) + zeros.at(0,c+right);
				// if(mass_zeros < 4 && mass_zeros >= 0)
				if(phase_winding != 0)
				{
					vortex.n = phase_winding;
					vortex.x = c + phase.make_vector(0.5, -0.5, 0);
					// cout << vortex.x << endl;
					vortex.points.clear();
					vortex.points.push_back(c);
					vortex.num_points = 1;
					vlist.push_back(vortex);					
				}
				/*if(get_vortex(phase.make_coord(x,y,z), phase, zeros, mass_zeros, checked, vortex)) // prueft auf Vortices
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

	 // cout << "findVortices before sorting" << endl;
	// list<VortexData> vlistCopy(vlist);
	// list<VortexData> vlistCopy1(vlist);
	vlist.sort([](VortexData &lhs, VortexData &rhs) {return lhs.surroundDens > rhs.surroundDens;});
	// vlist.sort([](VortexData &lhs, VortexData &rhs) {return lhs.zeroDensity < rhs.zeroDensity;});

	// for(list<VortexData>::iterator it = vlistCopy.begin(); it != vlistCopy.end(); ++it){
	// 	for(list<VortexData>::iterator et = vlistCopy1.begin(); et != vlistCopy.end(); ++et){
	// 		if(it->x == et->x)
	// 			vlist.push_back(*it);
	// 	}
	// }
	if(opt.initialRun == true){
		opt.vortexnumber = vlist.size();
		cout << "Evaluation found " << opt.vortexnumber << " Vortices." << endl;
	} else {
		if(vlist.size() > opt.vortexnumber){
			list<VortexData>::iterator it1 = vlist.begin();
			advance(it1,opt.vortexnumber);
			vlist.erase(it1,vlist.end());
		}
	}
}

inline double Eval::norm(Coordinate<double> &a, Coordinate<double> &b, double &h_x, double &h_y){
	return sqrt( (a.x() - b.x()) * (a.x() - b.x()) * h_x * h_x + (a.y() - b.y()) * (a.y() - b.y()) * h_y * h_y);
}

void Eval::getVortexDistance(PathResults &pres){

	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2];

	pres.histogram.resize(OBSERVABLES_DATA_POINTS_SIZE);
	pres.distance.resize(OBSERVABLES_DATA_POINTS_SIZE);
	for (list<VortexData>::iterator it = pres.vlist.begin(); it != pres.vlist.end(); ++it)
	{
		if(it->n != 0)
		{
			double shortest_distance = 65536.0;
			for(list<VortexData>::iterator oit = pres.vlist.begin(); oit != pres.vlist.end(); ++oit)
			{
				if((it != oit) && (oit->n!=0))
				{
					double distance = (it->x - oit->x).norm();
					double coordDistance = Eval::norm(it->x, oit->x,h_x,h_y);
					// cout << "Coordinate Distance: " << coordDistance << endl;
					// cout << "Grid Distance: " << distance << endl;
					if(distance < shortest_distance)
						shortest_distance = distance;
					pairDistanceHistogram(pres, distance, coordDistance);
				}
			}
			// if(shortest_distance != 65536.0)
			// {
			// 	it->pair_distance = shortest_distance;
			// 	//inc_pd_histogram(ares.pd_histogram_closest, shortest_distance);
			// 	av_pair_distance.average(it->pair_distance);
			// 	if(it->pair_distance > ares.max_pair_distance)
			// 		ares.max_pair_distance = it->pair_distance;
			// }
			// else
			// 	it->pair_distance = 0.0;
		}
		// else
		// 	it->pair_distance = 0.0;
	}
	
	// // mittlerer Abstand der Vortices
	// ares.pair_distance_all = av_pair_distance.av();
	// if(av_pair_distance.av() != 0)
	// 	ares.pair_distance_nonzero = av_pair_distance.av();
}

inline void Eval::pairDistanceHistogram(PathResults &pres, double &distance, double &coordDistance)
{
	double max_distance = sqrt(opt.grid[1]*opt.grid[1] + opt.grid[2]*opt.grid[2] + opt.grid[3]*opt.grid[3]) / 2.0;
	int index = (pres.histogram.size() - 1) * distance / max_distance;
	pres.histogram[index] += 0.5;
	pres.distance[index] = coordDistance;

}

int Eval::getVortexNumber(){
	return opt.vortexnumber;
}

void Eval::calc_fields(ComplexGrid &data, Options &opt){
	double LOWER_THRESHOLD = 0.0; // opt.N * 0.05 / (4. * opt.min_x * opt.stateInformation[0] * opt.min_y * opt.stateInformation[1]); //opt.N * 0.05 / data.width() / data.height() / data.depth();
	for(int x = 0; x < data.width(); x++)
	{
		for(int y = 0; y < data.height(); y++)
		{
			for(int z = 0; z < data.depth(); z++)
			{	
				phase.at(0,x,y,z) = arg(data(0,x,y,z));
				if(abs2(data(0,x,y,z)) <= LOWER_THRESHOLD)
					zeros.at(0,x,y,z) = 0.0;
				else
					zeros.at(0,x,y,z) = 1.0;
			}
		}
	}
	string name = "Test Zeros"+to_string(snapshot_time);
	plotDataToPng(name, name, zeros, opt);
}



void Eval::getDensity(ComplexGrid &data, RealGrid &densityLocationMap_local, vector<Coordinate<int32_t>> &densityCoordinates_local, int &densityCounter){
	 //  / 
	// double upper_threshold = 20.;
	// cout << "Threshold " << threshold << endl;

	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2];
	double threshold = opt.N * 0.001 / (4. * opt.min_x * opt.stateInformation[0] * opt.min_y * opt.stateInformation[1]);  //abs2(data(0,opt.grid[1]/2,opt.grid[2]/2,0))*0.9; 

	// RealGrid densityLocationMap = RealGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	densityLocationMap_local = RealGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);



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
	x_dist_grad[0] = x_dist_grad[opt.grid[1]-1] = y_dist_grad[0] = y_dist_grad[opt.grid[2]-1] = 0.0;

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

	// densityLocationMap = densityLocationMap_local;
}

void Eval::aspectRatio(Observables &obs, int &sampleindex){
	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2];
	double h[2];
	h[0] = h_x;
	h[1] = h_y;

		// Aspect-Ratio
	obs.r_max = 0;
	obs.r_min = (opt.grid[1] * h_x >= opt.grid[2] * h_y) ? opt.grid[1] * h_x : opt.grid[2] * h_y;
	vector<contourData> cData;
	for(c_set::iterator it = contour[sampleindex].begin(); it != contour[sampleindex].end(); ++it){
		contourData tmp;
		tmp.c = *it;
		int x_shift = tmp.c.x() - opt.grid[1]/2;
		int y_shift = tmp.c.y() - opt.grid[2]/2;
		tmp.phi = atan2(y_shift * h_y,x_shift * h_x) * 180 /M_PI + 180;
		tmp.r = sqrt(x_shift*x_shift * h_x*h_x + y_shift*y_shift *h_y*h_y);
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
		if(divisor_counter[i] == 0){
			divisor_counter[i] = 1;
		}

		cRadius[i] /= divisor_counter[i];
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

Observables Eval::calculator(ComplexGrid data,int sampleindex){
	
	Observables obs = Observables(OBSERVABLES_DATA_POINTS_SIZE);
	// R-Space
	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2];
	double h[2];
	double x_max = opt.stateInformation[0] * opt.min_x;
	double y_max = opt.stateInformation[1] * opt.min_y;
	double rmax[2];
	rmax[0] = x_max;
	rmax[1] = y_max;
	h[0] = h_x;
	h[1] = h_y; 
	// double raw_volume = h_x * opt.grid[1] * h_y * opt.grid[2];
	
	// double threshold = abs2(data(0,opt.grid[1]/2,opt.grid[2]/2,0))*0.9;

	// cout << "DensityCounter " << sampleindex << " : " << densityCounter[sampleindex] << endl;

	obs.volume = h_x * h_y * densityCounter[sampleindex];
	for(int i = 0; i < opt.grid[1]; i++){
	    for(int j = 0; j < opt.grid[2]; j++){	    	    		
	      	obs.particle_count += abs2(data(0,i,j,0));

	    }
	}
	obs.particle_count *= h_x * h_y;
	obs.density = obs.particle_count / obs.volume;

	aspectRatio(obs,sampleindex);

	// // Aspect-Ratio
	// obs.r_max = 0;
	// obs.r_min = (opt.grid[1] * h_x >= opt.grid[2] * h_y) ? opt.grid[1] * h_x : opt.grid[2] * h_y;
	// vector<contourData> cData;
	// for(c_set::iterator it = contour[sampleindex].begin(); it != contour[sampleindex].end(); ++it){
	// 	contourData tmp;
	// 	tmp.c = *it;
	// 	int x_shift = tmp.c.x() - opt.grid[1]/2;
	// 	int y_shift = tmp.c.y() - opt.grid[2]/2;
	// 	tmp.phi = atan2(y_shift * h_y,x_shift * h_x) * 180 /M_PI + 180;
	// 	tmp.r = sqrt(x_shift*x_shift * h_x*h_x + y_shift*y_shift *h_y*h_y);
	// 	cData.push_back(tmp);
	// 	// if(obs.r_max <= tmp.r){
	// 	// 	obs.r_max = tmp.r;
	// 	// 	obs.r_max_phi = tmp.phi;
	// 	// }
	// 	// if(obs.r_min >= tmp.r){
	// 	// 	obs.r_min = tmp.r;
	// 	// 	obs.r_min_phi = tmp.phi;
	// 	// }

	// }

	// std::sort(cData.begin(),cData.end(),[](const contourData &lhs, const contourData &rhs) -> bool {return (lhs.phi < rhs.phi);});

	// vector<double> cRadius(361);
	// vector<int> divisor_counter(361);
	// for(vector<contourData>::const_iterator it = cData.begin(); it != cData.end(); ++it){
	// 	int index = round(it->phi);
	// 	cRadius[index] += it->r;
	// 	divisor_counter[index]++;
	// }
	// // cRadius.erase(cRadius.begin());
	// // divisor_counter.erase(divisor_counter.begin());
	// cRadius[0] += cRadius[360]; cRadius.pop_back();
	// divisor_counter[0] += cRadius[360]; divisor_counter.pop_back();

	// for(int i = 0; i < 360; i++){
	// 	if(divisor_counter[i] == 0){
	// 		divisor_counter[i] = 1;
	// 	}

	// 	cRadius[i] /= divisor_counter[i];
	// }

	// // string name = "cRadius_" + to_string(sampleindex) + "_" + to_string(sampleindex);
	// // plotVector(name,cRadius,opt);

	// vector<double> cDistance(180);
	// for(int i = 0; i < 180; i++){
	// 	cDistance[i] = fabs(cRadius[i] + cRadius[i+180]);
	// }
	// vector<double>::iterator maxDistance = std::max_element(cDistance.begin(), cDistance.end());
	// obs.r_max = *maxDistance;
	// obs.r_max_phi = std::distance(cDistance.begin(), maxDistance);

	// vector<double>::iterator minDistance = std::min_element(cDistance.begin(), cDistance.end());
	// obs.r_min = *minDistance;
	// obs.r_min_phi = std::distance(cDistance.begin(), minDistance);

	// vector<double> tmp_ratio(90);
	// // double tmp_ratio = 0;
	// for(int i = 0; i < 89; i++){
	// 	if(cDistance[i+90] >= 0.0){
	// 		tmp_ratio[i] = cDistance[i] / cDistance[i+90];
	// 	} else {
	// 		cout << "WARNING: Aspect-Ratio: Calculated Distance smaller than zero!" << endl;
	// 	}
	// 	// double tmp1 = cDistance[i] / cDistance[i+90];
	// 	// double tmp2 = cDistance[i+1] / cDistance[i+91];
	// 	// tmp_ratio = (tmp1 > tmp2) ? tmp1 : tmp2;
	// }
	// obs.fixedAspectRatio = tmp_ratio[0];
	// vector<double>::iterator maxElement = std::max_element(tmp_ratio.begin(),tmp_ratio.end());
	// int maxAspectRatioIndex = std::distance(tmp_ratio.begin(), maxElement);
	// obs.aspectRatioAngle = maxAspectRatioIndex;
	// obs.aspectRatio = *maxElement; // zwischen 0 und 90 grad, also effektiv x und y richtung
	// // FIXME replace this with a check for the max and min values, save the corresponding angles and check if they change (= overall rotation in the gas!)

	// == Angular Density
	// vector<double> angularDensity;
	vector<double> phi;
	vector<double> polarDensity; // first entry is 
	vector<double> radius;
	vector<Coordinate<int32_t>> cartesianCoordinates;

	ArrayXd divisor2(obs.number.size());
	double r_index_factor = (obs.radialDensity.size() -1) / sqrt(x_max * x_max + y_max * y_max);
	for(int i = 0; i < opt.grid[1]; i++){
	    for(int j = 0; j < opt.grid[2]; j++){
				int x_shift = i - opt.grid[1]/2;
				int y_shift = j - opt.grid[2]/2;
				phi.push_back( atan2(x_shift * h_x ,y_shift * h_y) * 180 / M_PI + 180);
				double r_tmp = sqrt(x_shift*x_shift * h_x*h_x + y_shift*y_shift *h_y*h_y);
				radius.push_back(r_tmp);
				double dens_tmp = abs2(data(0,i,j,0));
				polarDensity.push_back(dens_tmp);
				cartesianCoordinates.push_back(data.make_coord(i,j,0));

				int index = r_index_factor * r_tmp;
				divisor2(index)++;
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
		for(int i = 0; i <= opt.grid[d+1]/2; i++){
		// for (int32_t i = 0; i < kspace[d].size()/2; i++){
			// kspace[d][i] = opt.klength[d]/**opt.stateInformation[0]*/*2.0*sin( M_PI*((double)i)/((double)opt.grid[d+1]) );
			// kspace[d][i] = opt.klength[d]*((double)i)/((double)(opt.grid[d+1]/2));
			// kspace[d][i] = opt.klength[d] * 2 * M_PI  * ((double)i) / ((double)(opt.grid[d+1]*opt.grid[d+1]*h[d]));
			kspace[d][i] = (M_PI / rmax[d]) * (double)i;
		}
		// for (int i=opt.grid[d+1]/2; i<opt.grid[d+1]; i++){
		for(int i = -(opt.grid[d+1]/2)+1; i < 0; i++){
		// for (int32_t i = kspace[d].size()/2; i < kspace[d].size(); i++){
			// kspace[d][i] = opt.klength[d]/**opt.stateInformation[1]*/*2.0*sin( M_PI*((double)(-opt.grid[d+1]+i))/((double)opt.grid[d+1]) );
			// kspace[d][i] = opt.klength[d]*((double)(opt.grid[d+1]-i))/((double)opt.grid[d+1]/2);
			// kspace[d][i] = opt.klength[d] * 2 * M_PI  * ((double)(-opt.grid[d+1]+i)) / ((double)(opt.grid[d+1]*opt.grid[d+1]*h[d]));
			kspace[d][i] = (M_PI / rmax[d]) * (double)i;
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

bool Eval::checkResizeCondition(){
	bool resize = false;
	double innerCircleDistance = opt.min_x + opt.min_y;
	if(totalResult.r_max >= (innerCircleDistance * 0.85)){
		resize = true;
	}
	return resize;
}


Observables Eval::calculatorITP(ComplexGrid data,int sampleindex){
	
	Observables obs = Observables(OBSERVABLES_DATA_POINTS_SIZE);
	// R-Space
	double h_x = 2. * opt.stateInformation[0] * opt.min_x / opt.grid[1];
	double h_y = 2. * opt.stateInformation[1] * opt.min_y / opt.grid[2];
	double h[2];
	h[0] = h_x;
	h[1] = h_y; 
	// double raw_volume = h_x * opt.grid[1] * h_y * opt.grid[2];
	
	// double threshold = abs2(data(0,opt.grid[1]/2,opt.grid[2]/2,0))*0.9;

	obs.volume = h_x * h_y * densityCounter[sampleindex];
	for(int i = 0; i < opt.grid[1]; i++){
	    for(int j = 0; j < opt.grid[2]; j++){	    	    		
	      	obs.particle_count += abs2(data(0,i,j,0));
	      	
	    }
	}
	obs.particle_count *= h_x * h_y;
	obs.density = obs.particle_count / obs.volume;

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
	
	