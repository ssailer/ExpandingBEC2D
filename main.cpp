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
	MainControl mC = initMain.getControl();

	if(!initMain.restart()){
		// runExpanding->noise();
		if(mC == RK4){

			MatrixData* startGrid = new MatrixData(1,tmpOpt.grid[1],tmpOpt.grid[2],0,0,tmpOpt.min_x,tmpOpt.min_y);

			setGridToTF(startGrid,initMain.getOptions());

			// int vnumber = 0;
			// addVorticesAlternating(startGrid,tmpOpt,vnumber);

			string startGridName = initMain.getStartingGridName();
			binaryFile* startFile = new binaryFile(startGridName,binaryFile::out);
			startFile->appendSnapshot("StartGrid",0,startGrid,tmpOpt);
			delete startFile;
	
			// cout << "EigenThreads: " << Eigen::nbThreads() << endl;
			
			MatrixData* data = new MatrixData(initMain.getMeta());
			
			binaryFile* dataFile = new binaryFile(startGridName,binaryFile::in);
			dataFile->getSnapshot("StartGrid",0,startGrid,tmpOpt);
			delete dataFile;

			for(int i = 0; i < data->meta.samplesize; i++){
				data->wavefunction[i] = startGrid->wavefunction[0];
			}
			delete startGrid;

			RotatingTrap* rt = new RotatingTrap(/*data,*/tmpOpt);
			
			string runName = "ex";
			RTE* runExpanding = new RTE(data,rt,tmpOpt);
			cout << "rteToTime()" << endl;
			runExpanding->rteToTime(runName);

			delete runExpanding;
			delete rt;
			delete data;

		}
		// if(mC == TRAP){

		// 	MatrixData* startGrid = new MatrixData(1,tmpOpt.grid[1],tmpOpt.grid[2],0,0,tmpOpt.min_x,tmpOpt.min_y);
	
		// 	cout << "EigenThreads: " << Eigen::nbThreads() << endl;
			
		// 	string startGridName = "StartGrid_2048_2048.h5"; // "StartGrid_2048x2048_N1000_alternatingVortices.h5";
			
		// 	MatrixData* data = new MatrixData(initMain.getMeta());
			
		// 	binaryFile* dataFile = new binaryFile(startGridName,binaryFile::in);
		// 	dataFile->getSnapshot("StartGrid",0,startGrid,tmpOpt);
		// 	delete dataFile;

		// 	for(int i = 0; i < data->meta.samplesize; i++){
		// 		data->wavefunction[i] = startGrid->wavefunction[0];
		// 	}
		// 	delete startGrid;
			
		// 	string runName = "trap";
		// 	RTE* runExpanding = new RTE(data,initMain.getOptions());
		// 	cout << "rteToTime()" << endl;
		// 	runExpanding->rteToTime(runName);

		// 	binaryFile* trapFile = new binaryFile(startGridName,binaryFile::out);
		// 	trapFile->appendSnapshot("StartGrid",0,data,tmpOpt);
		// 	delete trapFile;

		// 	delete runExpanding;
		// 	delete data;

		// }
		// if(mC == SPLIT){

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
		// if(mC == RK4_RESTART){

		// 	MatrixData* data = new MatrixData(initMain.getMeta());

		// 	string runName = "ex";
		// 	string filename = "LastGrid.h5";

		// 	binaryFile* dataFile = new binaryFile(filename,binaryFile::in);

		// 	vector<int> timeList = dataFile->getTimeList();
		// 	dataFile->getSnapshot(runName,timeList[0],data,tmpOpt);
		// 	delete dataFile;

		// 	tmpOpt = initMain.getOptions();

		// 	tmpOpt.initialRun = true;
		// 	tmpOpt.n_it_RTE = initMain.getRunTime();
		// 	tmpOpt.snapshots = initMain.getSnapShots();
	
		// 	data->meta.steps = 0;
		// 	data->meta.time = 0.0;
		// 	tmpOpt.t_abs = complex<double>(0.0,0.0);

	
		// 	RTE* runExpanding = new RTE(data,tmpOpt);

		// 	cout << "rteToTime()" << endl;
		// 	runExpanding->rteToTime(runName);

		// 	delete runExpanding;
		// 	delete data;
		// }



	}
	
	if(initMain.restart()){
		string runName = "ex";
		string filename = "LastGrid.h5";
		MatrixData* data = new MatrixData(initMain.getMeta());
		binaryFile* dataFile = new binaryFile(filename,binaryFile::in);
	
		vector<int> timeList = dataFile->getTimeList();
		dataFile->getSnapshot(runName,timeList[0],data,tmpOpt);
		delete dataFile;
		tmpOpt.initialRun = false;
		tmpOpt.n_it_RTE = initMain.getRunTime();
		tmpOpt.snapshots = initMain.getSnapShots();

		RotatingTrap* rt = new RotatingTrap(/*data,*/tmpOpt);

		RTE* runExpanding = new RTE(data,rt,tmpOpt);

		if(mC == RK4){
			cout << "rteToTime()" << endl;
			runExpanding->rteToTime(runName);
		}
		// if(mC == SPLIT){
		// 	cout << "splitToTime()" << endl;
		// 	runExpanding->splitToTime(runName);
		// }

		delete runExpanding;
		delete rt;
		delete data;
	}
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
	return SUCCESS; 
}
cerr << "[END]" << endl; 
return SUCCESS; 	
}




