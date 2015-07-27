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
#include <cstdio>
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
// #include <hydro.h>
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

int evaluate(InitMain &initMain){
	initMain.printInitVar();
	Options opt = initMain.getOptions();
	

	MainControl dgl = initMain.getDgl();
	string filename;
	if(dgl == ROT){
		filename = "rotdata.h5";
		if(remove("runObservables/ROT_Observables.dat")){ cerr << "deleted ROT_Observables.dat" << endl;} else { cerr << "could not delete ROT_Observables.dat" << endl;}
	}
	if(dgl == EXP){
		filename = "expdata.h5"; 
		if(remove("runObservables/EXP_Observables.dat")){ cerr << "deleted EXP_Observables.dat" << endl;} else { cerr << "could not delete EXP_Observables.dat" << endl;}
	}
	if(dgl == TRAP){
		filename = "trapdata.h5"; 
		if(remove("runObservables/TRAP_Observables.dat")){ cerr << "deleted TRAP_Observables.dat" << endl;} else { cerr << "could not delete TRAP_Observables.dat" << endl;}
	}

	binaryFile*dataFile = new binaryFile(filename,binaryFile::in);
	vector<int> timeList = dataFile->getTimeList();
	int size = timeList.size();

	// Eval* initEval = new Eval(*data,opt);
	// dataFile->getEval(timeList[size-1],*initEval,opt);
	// double maxTime = initEval->data.meta.time;
	// dataFile->getEval(timeList[0],*initEval,opt);	
	// cerr << endl << "hydro plotting" << endl;
	// initEval->opt.vortexnumber = opt.vortexnumber;
	// hydroSolver solver(initEval,maxTime);
	// solver.integrate();
	// solver.pyPlot();
	// delete initEval;
	MatrixData* data;

	#pragma omp parallel for ordered private(data) num_threads(4) schedule(static,1)
	for(int k = 0; k < size; ++k){
		// cerr << "loading: " << k << " / " << size-1;
		data = new MatrixData();
		#pragma omp critical
		{
			dataFile->getSnapshot("MatrixData", timeList[k], data, opt);
		}
		opt.isDimensionless = true;
		data->meta.Ag = opt.Ag;
		data->meta.OmegaG = opt.OmegaG;
		
	
		Eval eval(*data,opt);
		// cerr << " processing ";
		eval.process();
		#pragma omp ordered
		{	
			cerr << " " << data->meta.steps << " step " << data->meta.time << " time" << endl;
			eval.save();
		}

		delete data;
	}

	delete dataFile;



	chdir("..");
	cerr << "[END]" << endl;

}

int plotting(InitMain &initMain){
	initMain.printInitVar();
	Options opt = initMain.getOptions();
	

	MainControl dgl = initMain.getDgl();
	string filename;
	if(dgl == ROT){
		filename = "rotdata.h5";
	}
	if(dgl == EXP){
		filename = "expdata.h5"; 
	}
	if(dgl == TRAP){
		filename = "trapdata.h5"; 
	}

	binaryFile* dataFile = new binaryFile(filename,binaryFile::in);
	vector<int> timeList = dataFile->getTimeList();
	int size = timeList.size();

	// Eval* initEval = new Eval(*data,opt);
	// dataFile->getEval(timeList[size-1],*initEval,opt);
	// double maxTime = initEval->data.meta.time;
	// dataFile->getEval(timeList[0],*initEval,opt);	
	// cerr << endl << "hydro plotting" << endl;
	// initEval->opt.vortexnumber = opt.vortexnumber;
	// hydroSolver solver(initEval,maxTime);
	// solver.integrate();
	// solver.pyPlot();
	// delete initEval;
	Eval* eval;

	#pragma omp parallel for ordered private(eval) num_threads(4) schedule(static,1)
	for(int k = 0; k < size; ++k){
		// cerr << "loading: " << k << " / " << size-1;
		MatrixData* data = new MatrixData();
		eval = new Eval(*data,opt);
		#pragma omp critical
		{
			dataFile->getSnapshot("MatrixData", timeList[k], data, opt);
		}
		opt.isDimensionless = true;
		data->meta.Ag = opt.Ag;
		data->meta.OmegaG = opt.OmegaG;
		cerr << " " << data->meta.steps << " step " << data->meta.time << " time";
		eval = new Eval(*data,opt);
		eval->process();

		Plotter plotter(*eval,opt);
		cerr << "plotting" << endl;
		plotter.plotEval();
		delete data;
		delete eval;
	}

	delete dataFile;



	chdir("..");
	cerr << "[END]" << endl;

}

int hydro(InitMain &initMain){
	initMain.printInitVar();
	Options opt = initMain.getOptions();
	MatrixData* data = new MatrixData();

	MainControl dgl = initMain.getDgl();
	string filename;
	if(dgl == ROT) filename = "rotdata.h5";
	if(dgl == EXP) filename = "expdata.h5";
	if(dgl == TRAP) filename = "trapdata.h5";


	binaryFile* dataFile = new binaryFile(filename,binaryFile::in);
	vector<int> timeList = dataFile->getTimeList();
	int size = timeList.size();

	Eval* initEval = new Eval(*data,opt);
	dataFile->getEval(timeList[size-1],*initEval,opt);
	double maxTime = initEval->data.meta.time;
	dataFile->getEval(timeList[0],*initEval,opt);
	initEval->opt.vortexnumber = opt.vortexnumber;

	cerr << endl << "hydro plotting" << endl;
	hydroSolver solver(initEval,maxTime);
	solver.integrate();
	solver.integrate2();
	solver.pyPlot();
	delete initEval;

	chdir("..");
	cerr << "[END]" << endl;
}

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
			case ITP : {
					Runner<SplitITP>* run = new Runner<SplitITP>(data,opt);
					run->runToTime("itp");
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
 	MainControl mode = initMain.getRestart();
 	if(mode == PLOT){
 		for(int i = 0; i < size; ++i){
 			initMain.setIteration(i);
 			plotting(initMain);
 		}
 	} else if(mode == HYDRO){
 		for(int i = 0; i < size; ++i){
 			initMain.setIteration(i);
 			hydro(initMain);
 		}
 	} else if(mode == EVAL){
 		for(int i = 0; i < size; ++i){
 			initMain.setIteration(i);
 			evaluate(initMain);
 		}
 	} else {
 		for(int i = 0; i < size; ++i){
 			initMain.setIteration(i);
 			simulation(initMain);
 		}
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




