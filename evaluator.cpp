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

#include <EXP2D_MatrixData.h>
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
#define DEBUG_LOG 0

using namespace std;

int main( int argc, char** argv) 
{	
try{

	StartUp startUp(argc,argv);

	if(DEBUG_LOG == 1){
 		std::ofstream logstream("run.log");
 		redirecter redirect(logstream,std::cout); // redirects cout to logstream, until termination of this program. If DEBUG_LOG 1 is set, use cerr for output to console.
 	}
	

		
	string runname = "RT-No-Ex";
	vector<string> snapShotFiles;
	string tmp;
	string fileNameListName = "runData/fileNameList.dat";
	ifstream fileNameList(fileNameListName);

	if(fileNameList.is_open()){
		while (getline (fileNameList,tmp)){
			snapShotFiles.push_back(tmp);
		}
	}

	cout << "Snapshot File Size: " << snapShotFiles.size() << endl;
	int counter = 0;		
	MatrixData* matrixData = new MatrixData(startUp.getMeta());
	for(int j = 0; j < snapShotFiles.size(); j++){
		
		counter++;
		string h5name = snapShotFiles[j];
		Options opt;
		Eval results;

		cout << counter << " Opening datafile " << h5name << flush;
		binaryFile data(h5name,binaryFile::in);	

		cout << " >> Reading Datafiles " << h5name << flush;
		vector<int> timeList = data.getTimeList();
		for(int i = 0; i < timeList.size(); i++){
			data.getSnapshot(runname,timeList[i],matrixData,opt);		
			results.saveData(matrixData->wavefunction,opt,timeList[i],runname);		
			cout << " >> Evaluating Datafiles "<< timeList[i] << flush;
			results.evaluateData();		
			cout << " >> Plotting Datafiles " << timeList[i] << endl;		
			results.plotData();
		}

	}	
	
	cout << "Terminating successfully." << endl;
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
	return ERROR_UNHANDLED_EXCEPTION; 
// the code could be different depending on the exception message 
}
return SUCCESS; 	
}



