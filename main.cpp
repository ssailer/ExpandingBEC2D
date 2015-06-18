/**************************************************************************
Title: Simulating the Expansion of Turbulent Bose-Einstein Condensates (2D) 
Author: Simon Sailer (This work is based on the work of Bartholomew Andrews who made this as his master thesis.)
Last Update: 22/07/13
**************************************************************************/
#define EIGEN_FFTW_DEFAULT


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

// #include <complexgrid.h>
#include <bh3defaultgrid.h>
// #include <averageclass.h>
// #include <bh3observables.h>

#include <EXP2D_MatrixData.h>
#include <main.h>
#include <EXP2D_tools.h>
// #include <EXP2D_itp.hpp>
#include <EXP2D_binaryfile.h>
#include <EXP2D_rk4.hpp>
// #include <EXP2D_rte.hpp>
#include <EXP2D_runner.hpp>
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

int simulation(InitMain &initMain){

	initMain.printInitVar();

	Options opt = initMain.getOptions();
	MainControl algo = initMain.getAlgorithm(); // controls the choise between expansion and trapped simulations, e.g. the lab setup, has to be set in console
	MainControl dgl = initMain.getDgl(); // controls the chosen integration algorithm, has to be set in the cfg
	MainControl start = initMain.getRestart();

	string filename;

	MatrixData* data = new MatrixData(initMain.getMeta());

	if(start == RESUME){
		// Options loadedOptions;		
		// Loading from existing HDF5
		if(dgl == ROT) filename = "rotdata.h5";
		if(dgl == EXP) filename = "expdata.h5";
		if(dgl == TRAP) filename = "trapdata.h5";
		binaryFile* dataFile = new binaryFile(filename,binaryFile::in);	
		vector<int> timeList = dataFile->getTimeList();
		dataFile->getLatestSnapshot("MatrixData",data,opt);
		delete dataFile;
		// temp fix to change behaviour of sim object;
		opt.initialRun = false;
	}
	if(start == RESTART){
		Options loadedOptions;
		// Loading from existing HDF5
		if(dgl == ROT) filename = "rotdata.h5";
		if(dgl == EXP) filename = "rotdata.h5";
		if(dgl == TRAP) filename = "rotdata.h5";		
		binaryFile* dataFile = new binaryFile(filename,binaryFile::in);	
		vector<int> timeList = dataFile->getTimeList();
		dataFile->getLatestSnapshot("MatrixData",data,loadedOptions);
		delete dataFile;
		// temp fix to change behaviour of sim object;
		opt.initialRun = true;
		data->meta.time = 0.0;
		// opt.n_it_RTE = initMain.getRunTime();
		// opt.snapshots = initMain.getSnapShots();
	}
	if(start == NEW){
		setGridToTF(data,initMain.getOptions());
		// addVorticesAlternating(data, opt, opt.vortexnumber);
	}

	if(algo == RK4){
		switch ( dgl ){
			case ROT : {
					Runner<RotatingTrap>* run = new Runner<RotatingTrap>(data,opt);
					run->runToTime("rot");
					delete run;
				}
				break;

			case EXP : {
					Runner<Expansion>* run = new Runner<Expansion>(data,opt);
					run->runToTime("exp");
					delete run;
				}
				break;

			// case TRAP : {
			// 		Runner<Trap>* run = new Runner<Trap>(data,opt);
			// 		run->runToTime("trap");
			// 		delete run;
			// 	}
			// 	break;

			default :
				cout << "No known runmode was recognized in main. Please revise." << endl;
				break;
		}
	}
	if(algo == SPLIT){
		switch ( dgl ){
			case ROT : {
					Runner<SplitRotStrang>* run = new Runner<SplitRotStrang>(data,opt);
					run->runToTime("rot");
					delete run;
				}
				break;

			case EXP : {
					Runner<SplitFree>* run = new Runner<SplitFree>(data,opt);
					run->runToTime("exp");
					delete run;
				}
				break;

			case TRAP : {
					Runner<SplitTrap>* run = new Runner<SplitTrap>(data,opt);
					run->runToTime("trap");
					delete run;
				}
				break;

			default :
				cout << "No known runmode was recognized in main. Please revise." << endl;
				break;
		}
	}  
	delete data;

	chdir("..");
	cerr << "[END]" << endl;
}

int main( int argc, char** argv){	

try{

	InitMain initMain(argc,argv);	

	#if DEBUG_LOG
 		std::ofstream logstream("simulation.log");
 		redirecter redirectcout(logstream,std::cout);
 	#endif

 	int size = initMain.getIterations();
 	for(int i = 0; i < size; ++i){
 		initMain.setIteration(i);
 		simulation(initMain);
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
	return ERROR_UNHANDLED_EXCEPTION; 
}
 
return SUCCESS; 	
}




