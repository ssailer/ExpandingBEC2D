/**************************************************************************
Title: Simulating the Expansion of Turbulent Bose-Einstein Condensates (2D) 
Author: Simon Sailer (This work is based on the work of Bartholomew Andrews who made this as his master thesis.)
Last Update: 22/07/13
**************************************************************************/

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

#include <complexgrid.h>
#include <bh3defaultgrid.h>
#include <averageclass.h>
#include <bh3observables.h>

#include <main.h>
#include <EXP2D_tools.h>
#include <EXP2D_itp.hpp>
#include <EXP2D_rte.hpp>
#include <EXP2D_evaluation.h>
#include <plot_with_mgl.h>

// #include <typeinfo>

#define SUCCESS 0
#define ERROR_IN_COMMAND_LINE 1
#define ERROR_IN_CONFIG_FILE 2
#define ERROR_UNHANDLED_EXCEPTION 3
#define DEBUG_LOG 1

using namespace std;

int main( int argc, char** argv) 
{	
try{
 	std::streambuf *psbuf, *backup;
 	std::ofstream logstream;
 	backup = std::cout.rdbuf();     // back up cout's streambuf
 	// std::cout.rdbuf(backup);        // restore cout's original streambuf
	
	Options opt;
	
	read_cli_options(argc,argv,opt);
	read_config(argc,argv,opt);
	set_workingdirectory(opt);
	
	if(DEBUG_LOG == 1){
		logstream.open ("run.log");
		psbuf = logstream.rdbuf();        // get file's streambuf
		std::cout.rdbuf(psbuf);         // assign streambuf to cout
	}
	
	// Initialize the needed grid object 
	ComplexGrid* data = new ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	
	
	//////////// VORTICES ////////////////
	if(opt.runmode.compare(0,1,"0") == 0){

		printInitVar(opt); 
		ITP* itprun = new ITP(data,opt);
		
		double sigma_real[2];
		sigma_real[0] = opt.min_x/4;
		sigma_real[1] = opt.min_y/4;		
		data = set_grid_to_gaussian(data,opt,sigma_real[0],sigma_real[1]);
		
		// set the datafile identifier name and save the initial grid
		
		string runname = "INIT";
		plotDataToPng(runname,data,opt);
		
		
		
		//====> Imaginary Time Propagation (ITP)
		itprun->propagateToGroundState("ITP1");
		runname = "ITP1";
		plotDataToPng(runname,data,opt);
		
		// vector<MatrixXcd> tmpMatrix(1);
		// tmpMatrix[0] = MatrixXcd(opt.grid[1],opt.grid[2]);
		// for(int i = 0; i < opt.grid[1]; i++){
		// 	for(int j = 0; j < opt.grid[2]; j++){
		// 		tmpMatrix[0](i,j) = data->at(0,i,j,0);
		// 	}
		// }
		// binaryFile * ITP1 = new binaryFile("ITP1.h5",binaryFile::out);
		// ITP1->appendSnapshot("ITP1",0,tmpMatrix,opt);
		// delete ITP1;		
			
				// if the given value is true, add vortices to the startgrid
		if(opt.runmode.compare(3,1,"1") == 0){
			int sigma_grid[2];
			sigma_grid[0] = opt.grid[1]/8;
			sigma_grid[1] = opt.grid[2]/8;
			double r = (sigma_grid[0]+sigma_grid[1])/2.0; 
			
			// data = add_central_vortex(data,opt);	
			// data = add_circle_vortex(data,opt,r,4);
			// data = add_circle_vortex(data,opt,r/4.0,6);
			// data = add_circle_vortex(data,opt,r*2.0/4.0,12);
			// data = add_circle_vortex(data,opt,r*3.0/4.0,24);
			data = addVortices(data,opt);
			
			
			
			//====> Imaginary Time Propagation (ITP)
			itprun->formVortices("ITP2");
			runname = "ITP2";
			plotDataToPng(runname,data,opt);
			
			
			// vector<MatrixXcd> tmpMatrix1(1);
			// tmpMatrix1[0] = MatrixXcd(opt.grid[1],opt.grid[2]);
			// for(int i = 0; i < opt.grid[1]; i++){
			// 	for(int j = 0; j < opt.grid[2]; j++){
			// 		tmpMatrix1[0](i,j) = data->at(0,i,j,0);
			// 	}
			// }
			// binaryFile * ITP2 = new binaryFile("ITP2.h5",binaryFile::out);
			// ITP2->appendSnapshot("ITP2",0,tmpMatrix1,opt);
			// delete ITP2;
		}
	delete itprun;
	}
	////// END VORTICES //////////
	
	
	
	
	//====> Real Time Expansion (RTE)	
	RTE* rterun = new RTE(data,opt);
	string runname = "RT-Ex";
	

	
	if(opt.runmode.compare(0,1,"1") == 0){
		vector<string> snapShotFiles;
		string tmp;
		string fileNameListName = "runData/fileNameList.dat";
		ifstream fileNameList(fileNameListName);

		if(fileNameList.is_open()){
			while (getline (fileNameList,tmp)){
				snapShotFiles.push_back(tmp);
			}
		}

		string h5name = snapShotFiles.back();
		binaryFile *dataLoading = new binaryFile(h5name,binaryFile::in);
		int previousTimes = dataLoading->getTimeList().back();
		delete dataLoading;

		vector<int> snapshot_times(opt.snapshots);	
		for(int i = 0; i < opt.snapshots; i++){
			snapshot_times[i] = ((i+1) * opt.n_it_RTE / opt.snapshots) + previousTimes;

		}
		rterun->rteFromDataToTime(runname,snapshot_times,h5name);
		printInitVar(opt);
	}

	if(opt.runmode.compare(0,1,"0") == 0){
		vector<int> snapshot_times(opt.snapshots);	
		for(int i = 0; i < opt.snapshots; i++){
			snapshot_times[i] = (i+1) * opt.n_it_RTE / opt.snapshots;
		}
		ofstream runparameters;
		runparameters.open(("runparameters.txt")/*.c_str()*/, ios::out | ios::trunc);
		runparameters << opt;
		runparameters.close();

		complex<double> tmp;
		tmp = opt.omega_y;
		opt.omega_y = opt.omega_x;
		opt.omega_x = tmp;
		
		tmp = opt.omega_y;
		opt.dispersion_y = opt.dispersion_x;
		opt.dispersion_x = tmp;
	
		rterun->setOptions(opt);
		rterun->RunSetup();
		rterun->rteToTime(runname,snapshot_times);
	}


	delete rterun;
	delete data;
	
	cout << "Terminating successfully." << endl;

	if(DEBUG_LOG == 1){
		logstream.close();
	}
}  // exceptions catcher


catch(const std::exception& e) 
{ 
  	std::cerr << "Unhandled Exception reached the top of main: " 
    	      << e.what() << ", application will now exit" << std::endl; 
	return ERROR_UNHANDLED_EXCEPTION; 
}
catch(expException& e){
	e.printString();
	std::cerr << " Terminating now." << endl;
	return ERROR_UNHANDLED_EXCEPTION;
}
catch (const std::string& errorMessage) 
{ 
	std::cerr << errorMessage.c_str(); 
	std::cerr << " Terminating now." << endl; 
	return SUCCESS; 
// the code could be different depending on the exception message 
}
cout << "Run complete. Terminating successfully." << endl; 
return SUCCESS; 	
}




