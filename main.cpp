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

	Options opt = initMain.getOptions();
	MainControl mC = initMain.getControl(); // controls the choise between expansion and trapped simulations, e.g. the lab setup, has to be set in console
	MainControl runMode = initMain.getRunMode(); // controls the chosen integration algorithm, has to be set in the cfg

	MatrixData* data = new MatrixData(initMain.getMeta());

	string runName = "run";

	// setenv("PYTHONPATH","../",1);
	// Py_Initialize();
 //    PyRun_SimpleString("import plot");

	if(initMain.restart()){
		
		Options loadedOptions;
		string filename = "rundata.h5";
		// Loading from existing HDF5
		MatrixData* loadedData = new MatrixData();
		binaryFile* dataFile = new binaryFile(filename,binaryFile::in);	
		vector<int> timeList = dataFile->getTimeList();
		dataFile->getLatestSnapshot("MatrixData",data,loadedOptions);
		delete dataFile;
		// data->checkedCopy(loadedData);
		// plotDataToPngEigen("checkedCopy", data->wavefunction[0],opt);
		delete loadedData;
		// temp fix to change behaviour of sim object;
		opt.initialRun = false;
		// opt.n_it_RTE = initMain.getRunTime();
		// opt.snapshots = initMain.getSnapShots();
	} else {
		// set MatrixData to specified initial conditions
		setGridToTF(data,initMain.getOptions());
		// addVorticesAlternating(data, opt, opt.vortexnumber);
		// save initial Grid
		// string startGridName = initMain.getStartingGridName();
		binaryFile* startFile = new binaryFile("rundata.h5",binaryFile::out);
		startFile->appendSnapshot("MatrixData",data,initMain.getOptions());
		delete startFile;
	}

	if(mC == RK4){
		switch ( runMode ){
			case ROT : {
					Runner<RotatingTrap>* run = new Runner<RotatingTrap>(data,opt);
					run->runToTime(runName);
					delete run;
				}
				break;

			case EXP : {
					// FIXME WATCH OUT time var reset for restarting, this can be very DANGEROUS!
					if(opt.initialRun == true)
						data->meta.time = 0.0;
					Runner<Expansion>* run = new Runner<Expansion>(data,opt);
					run->runToTime(runName);
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
	if(mC == SPLIT){
		switch ( runMode ){
			case ROT : {
					Runner<SplitRot>* run = new Runner<SplitRot>(data,opt);
					run->runToTime(runName);
					delete run;
				}
				break;

			case EXP : {
					Runner<SplitFree>* run = new Runner<SplitFree>(data,opt);
					run->runToTime(runName);
					delete run;
				}
				break;

			case TRAP : {
					Runner<SplitTrap>* run = new Runner<SplitTrap>(data,opt);
					run->runToTime(runName);
					delete run;
				}
				break;

			default :
				cout << "No known runmode was recognized in main. Please revise." << endl;
				break;
		}
	}
	if(mC == SPLITSTRANG){
		switch ( runMode ){
			case ROT : {
					Runner<SplitRotStrang>* run = new Runner<SplitRotStrang>(data,opt);
					run->runToTime(runName);
					delete run;
				}
				break;

			default :
				cout << "No known runmode was recognized in main. Please revise." << endl;
				break;
		}
	}    
	delete data;
	// Py_Finalize();
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




