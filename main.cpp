/**************************************************************************
Title: Simulating the Expansion of Turbulent Bose-Einstein Condensates (2D) 
Author: Simon Sailer (This work is based on the work of Bartholomew Andrews who made this as his master thesis.)
Last Update: 22/07/13
**************************************************************************/
#define EIGEN_FFTW_DEFAULT

#include <boost/program_options.hpp>
#include <iostream>
#include <unistd.h>
#include <cstdlib>
// #include <cstring>
#include <string>
#include <cmath>
#include <complex>
#include <omp.h>
#include <sys/stat.h>
#include <dirent.h>

#include <complexgrid.h>
#include <bh3defaultgrid.h>
#include <averageclass.h>
#include <bh3observables.h>

#include <EXP2D_MatrixData.h>
#include <main.h>
#include <EXP2D_tools.h>
#include <EXP2D_itp.hpp>
#include <EXP2D_binaryfile.h>
#include <EXP2D_rk4.hpp>
#include <EXP2D_rte.hpp>
#include <EXP2D_evaluation.h>
#include <plot_with_mgl.h>
#include <EXP2D_startgrids.h>

// #include <typeinfo>

#define SUCCESS 0
#define ERROR_IN_COMMAND_LINE 1
#define ERROR_IN_CONFIG_FILE 2
#define ERROR_UNHANDLED_EXCEPTION 3
#define DEBUG_LOG 1

using namespace std;

int main( int argc, char** argv){	

try{

	InitMain initMain(argc,argv);	

	#if DEBUG_LOG
 		std::ofstream logstream("simulation.log");
 		redirecter redirectcout(logstream,std::cout);
 		// redirects cout to logstream, until termination of this program. If DEBUG_LOG 1 is set, use cerr for output to console.
 		// std::ofstream errorstream("error.log");
 		// redirecter redirectcerr(errorstream,std::cerr);
 	#endif

 	initMain.printInitVar();

	Options tmpOpt = initMain.getOptions();
	MainControl mC = initMain.getControl(); // controls the choise between expansion and trapped simulations, e.g. the lab setup, has to be set in console
	MainControl runMode = initMain.getRunMode(); // controls the chosen integration algorithm, has to be set in the cfg

	MatrixData* data = new MatrixData(initMain.getMeta());

	if(initMain.restart()){
		string runName = "ex";
		string filename = "LastGrid.h5";
		// Loading from existing HDF5
		binaryFile* dataFile = new binaryFile(filename,binaryFile::in);	
		vector<int> timeList = dataFile->getTimeList();
		dataFile->getSnapshot(runName,timeList[0],data,tmpOpt);
		delete dataFile;
		// temp fix to change behaviour of sim object;
		tmpOpt.initialRun = false;
		tmpOpt.n_it_RTE = initMain.getRunTime();
		tmpOpt.snapshots = initMain.getSnapShots();
	} else {
		// set MatrixData to specified initial conditions
		setGridToTF(data,initMain.getOptions());
		// save initial Grid
		string startGridName = initMain.getStartingGridName();
		binaryFile* startFile = new binaryFile(startGridName,binaryFile::out);
		startFile->appendSnapshot("StartGrid",0,data,initMain.getOptions());
		delete startFile;
	}

	if(runMode == RK4){
		// Object determining the integration method
		RungeKutta* dgl_algorithm;
		// Consequently choosing the DGL
		if(mC == ROT) dgl_algorithm = new RotatingTrap(tmpOpt);
		if(mC == EXP) dgl_algorithm = new Expansion(tmpOpt);
		// if(mC == TRAP)
		// 	dgl_algorithm = new Trap(tmpOpt);
		
		// create simulation object with all accumulated settings
		// and start the run
		RTE* sim = new RTE(data,dgl_algorithm,tmpOpt);
		sim->rteToTime("ex");

		delete sim, dgl_algorithm;
	}
	// if(runMode == SPLIT){
	// 	MatrixData* startGrid = new MatrixData(1,tmpOpt.grid[1],tmpOpt.grid[2],0,0,tmpOpt.min_x,tmpOpt.min_y);
	
	// 	cout << "EigenThreads: " << Eigen::nbThreads() << endl;
		
	// 	string startGridName = initMain.getStartingGridName(); // "StartGrid_2048x2048_N1000_alternatingVortices.h5";
	
	// 	MatrixData* data = new MatrixData(initMain.getMeta());
	
	// 	binaryFile* dataFile = new binaryFile(startGridName,binaryFile::in);
	// 	dataFile->getSnapshot("StartGrid",0,startGrid,tmpOpt);
	// 	delete dataFile;

	// 	for(int i = 0; i < data->meta.samplesize; i++){
	// 		data->wavefunction[i] = startGrid->wavefunction[0];
	// 	}
	// 	delete startGrid;
	
	// 	string runName = "split";
	// 	RTE* runExpanding = new RTE(data,initMain.getOptions());
	// 	cout << "splitToTime()" << endl;
	// 	runExpanding->splitToTime(runName);
	// 	delete runExpanding;
	// 	delete data;
	// }
	delete data;
}


catch(const std::exception& e){ 
  	std::cerr << "Unhandled Exception reached the top of main: " 
    	      << e.what() << ", application will now exit" << std::endl; 
	return ERROR_UNHANDLED_EXCEPTION; 
}
catch(expException& e){
	e.printString();
	std::cerr << " Terminating now." << endl;
	return ERROR_UNHANDLED_EXCEPTION;
}
catch (const std::string& errorMessage){ 
	std::cerr << errorMessage.c_str(); 
	std::cerr << " Terminating now." << endl; 
	return ERROR_UNHANDLED_EXCEPTION; 
}
cerr << "[END]" << endl; 
return SUCCESS; 	
}




