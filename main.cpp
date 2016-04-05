#define EIGEN_FFTW_DEFAULT

#include <iostream>
#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <complex>
#include <omp.h>
#include <sys/stat.h>
#include <dirent.h>

#include "main.h"

#include "matrixdata.h"
#include "tools.h"
#include "binaryfile.h"
#include "rk4.h"
#include "runner.h"
#include "evaluation.h"
#include "plot_with_mgl.h"
#include "startgrids.h"


#define SUCCESS 0
#define ERROR_IN_COMMAND_LINE 1
#define ERROR_IN_CONFIG_FILE 2
#define ERROR_UNHANDLED_EXCEPTION 3
#define DEBUG_LOG 0

using namespace std;

int evaluate(InitMain &initMain){
	// initMain.printInitVar();

	Options opt = initMain.getOptions();


	MainControl dgl = initMain.getDgl();
	string filename;
	string resultfilename;
	if(dgl == ROT){
		filename = "rotdata.h5";
		resultfilename = "rotobs.h5";
		if(!remove("runObservables/ROT_Observables.dat")){ cerr << "deleted ROT_Observables.dat" << endl << endl;} else { cerr << "could not delete ROT_Observables.dat" << endl << endl;}
	}
	else if(dgl == EXP){
		filename = "expdata.h5";
		resultfilename = "expobs.h5";
		if(!remove("runObservables/EXP_Observables.dat")){ cerr << "deleted EXP_Observables.dat" << endl << endl;} else { cerr << "could not delete EXP_Observables.dat" << endl << endl;}
	}
	else if(dgl == TRAP){
		filename = "trapdata.h5";
		resultfilename = "trapobs.h5";
		if(!remove("runObservables/TRAP_Observables.dat")){ cerr << "deleted TRAP_Observables.dat" << endl << endl;} else { cerr << "could not delete TRAP_Observables.dat" << endl << endl;}
	}
	else {
		cout << " No valid runmode, choose ROT, EXP or TRAP!" << endl;
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
	shared_ptr<MatrixData> data;

	// Initialize empty h5 file, if old file exists, it will be overwritten!
	binaryFile* initObsFile = new binaryFile(resultfilename,binaryFile::out);
	delete initObsFile;

	// #pragma omp parallel for ordered private(data) num_threads(1) schedule(static,1)
	for(int k = 0; k < size; ++k){
		data = make_shared<MatrixData>(new MatrixData());
		#pragma omp critical
		{
			dataFile->getSnapshot("MatrixData", timeList[k], data, opt);
		}
		opt.isDimensionless = true;
		data->meta.Ag = opt.Ag;
		data->meta.OmegaG = opt.OmegaG;


		auto eval = std::make_shared<Eval>(data,opt);
		eval->process();
		#pragma omp ordered
		{
			if(k == 0){
				 eval->data->meta.time = 0.0;
			}
			eval->save();
			binaryFile* obsFile = new binaryFile(resultfilename,binaryFile::append);
			obsFile->appendEval(eval,opt);
			delete obsFile;
		}
	}

	delete dataFile;

	chdir("..");
	cerr << "[END]" << endl;

}

int resave(InitMain &initMain){
	// initMain.printInitVar();
	Options opt = initMain.getOptions();


	MainControl dgl = initMain.getDgl();
	string filename;
	if(dgl == ROT){
		filename = "rotdata.h5";
		if(!remove("runObservables/ROT_Observables.dat")){ cerr << "deleted ROT_Observables.dat" << endl << endl;} else { cerr << "could not delete ROT_Observables.dat" << endl << endl;}
	}
	else if(dgl == EXP){
		filename = "expdata.h5";
		if(!remove("runObservables/EXP_Observables.dat")){ cerr << "deleted EXP_Observables.dat" << endl << endl;} else { cerr << "could not delete EXP_Observables.dat" << endl << endl;}
	}
	else if(dgl == TRAP){
		filename = "trapdata.h5";
		if(!remove("runObservables/TRAP_Observables.dat")){ cerr << "deleted TRAP_Observables.dat" << endl << endl;} else { cerr << "could not delete TRAP_Observables.dat" << endl << endl;}
	}
	else {
		cout << " No valid runmode, choose ROT, EXP or TRAP!" << endl;
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
	shared_ptr<MatrixData> data;

	string dataname = "exp2data.h5";
	binaryFile* bFile = new binaryFile(dataname,binaryFile::out);
	cerr << "opened new datafile" << endl;

	// #pragma omp parallel for ordered private(data) num_threads(2) schedule(static,1)
	for(int k = 0; k < size; ++k){
		// cerr << "loading: " << k << " / " << size-1;
		data = make_shared<MatrixData>(new MatrixData());
		// #pragma omp critical
		// {
			dataFile->getSnapshot("MatrixData", timeList[k], data, opt);
		// }
		opt.isDimensionless = true;
		data->meta.Ag = opt.Ag;
		data->meta.OmegaG = opt.OmegaG;

		bFile->appendSnapshot("MatrixData",data,opt);
			// if(opt.runmode != "EXP"){
		// bFile->appendEval(*eval, opt);
			// }

		// cerr << "saved step " << timeList[k] << endl;

		// Eval eval(*data,opt);
		// // cerr << " processing ";
		// eval.process();
		// #pragma omp ordered
		// {
		// 	if(k == 0){
		// 		 eval.data.meta.time = 0.0;
		// 	}
			// cerr << " " << data->meta.steps << " step " << eval.data.meta.time << " time" << endl;
		// 	eval.save();
		// }

	}

	delete dataFile;
	delete bFile;



	chdir("..");
	cerr << "[END]" << endl;

}


int plotting(InitMain &initMain){
	// initMain.printInitVar();
	Options opt = initMain.getOptions();


	MainControl dgl = initMain.getDgl();
	string filename;
	if(dgl == ROT){
		filename = "rotdata.h5";
	}
	else if(dgl == EXP){
		filename = "expdata.h5";
	}
	else if(dgl == TRAP){
		filename = "trapdata.h5";
	}
	else {
		cout << " No valid runmode, choose ROT, EXP or TRAP!" << endl;
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
	// shared_ptr<Eval> eval;

	// #pragma omp parallel for ordered private(eval) num_threads(1) schedule(static,1)
	for(int k = 0; k < size; ++k){
		// cerr << "loading: " << k << " / " << size-1;
		// MatrixData* data = new MatrixData();
		auto data = make_shared<MatrixData>(new MatrixData());
		// eval = new Eval(*data,opt);
		#pragma omp critical
		{
			dataFile->getSnapshot("MatrixData", timeList[k], data, opt);
		}
		opt.isDimensionless = true;
		data->meta.Ag = opt.Ag;
		data->meta.OmegaG = opt.OmegaG;
		cerr << " " << data->meta.steps << " step " << data->meta.time << " time";
		auto eval = make_shared<Eval>(data,opt);
		eval->process();

		Plotter plotter(eval,opt);
		cerr << "plotting" << endl;
		plotter.plotEval();
	}

	delete dataFile;



	chdir("..");
	cerr << "[END]" << endl;

}

int hydro(InitMain &initMain){
	// initMain.printInitVar();
	Options opt = initMain.getOptions();
	auto data = make_shared<MatrixData>(new MatrixData());
	// MatrixData* data = new MatrixData();

	MainControl dgl = initMain.getDgl();
	string filename;
	if(dgl == ROT) filename = "rotobs.h5";
	else if(dgl == EXP) filename = "expobs.h5";
	else if(dgl == TRAP) filename = "trapobs.h5";
	else cout << " No valid runmode, choose ROT, EXP or TRAP!" << endl;

	cout << "Filename: " << filename << endl;
	binaryFile* dataFile = new binaryFile(filename,binaryFile::in);
	vector<int> timeList = dataFile->getTimeList();
	int size = timeList.size();

	// Eval* initEval = new Eval(data,opt);
	auto initEval = make_shared<Eval>(data,opt);
	dataFile->getEval(timeList[size-1],*initEval,opt);
	double maxTime = initEval->data->meta.time;
	dataFile->getEval(timeList[0],*initEval,opt);
	initEval->opt.vortexnumber = opt.vortexnumber;

	cerr << endl << "hydro plotting" << endl;
	hydroSolver solver(initEval,maxTime);
	solver.integrate();
	// solver.integrate2();
	solver.pyPlot();
	// delete initEval;

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

	auto data = make_shared<MatrixData>(initMain.getMeta());

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
		data->convertToDimensionlessUnits();
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
		data->convertToDimensionlessUnits();
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
					// Runner<SplitRot>* run = new Runner<SplitRot>(data,opt);
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

	chdir("..");
	cerr << endl << "Ended application successfully, bye bye and thanks for the fish. :)" << endl;
}

int main( int argc, char** argv){

try{

	InitMain initMain(argc,argv);

	#if DEBUG_LOG
		cout << "DEBUG_LOG actived, cout will be found in simulation.log" << endl;
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
	std::cerr << " Terminating." << endl;
	return ERROR_UNHANDLED_EXCEPTION;
}
catch (const std::string& errorMessage){
	std::cerr << errorMessage.c_str();
	std::cerr << " Terminating." << endl;
	return ERROR_UNHANDLED_EXCEPTION;
}

return SUCCESS;
}
