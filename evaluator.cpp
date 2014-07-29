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
#define DEBUG_LOG 0

using namespace std;

int main( int argc, char** argv) 
{	
try{

 	std::streambuf *psbuf, *backup;
 	std::ofstream logstream;
 	backup = std::cout.rdbuf();     // back up cout's streambuf
 	// std::cout.rdbuf(backup);        // restore cout's original streambuf

	StartUp startUp(argc,argv);
	
	if(DEBUG_LOG == 1){
		logstream.open ("eval.log");
		psbuf = logstream.rdbuf();        // get file's streambuf
		std::cout.rdbuf(psbuf);         // assign streambuf to cout
	}
		
	string runname = "RT-Ex";
	vector<string> snapShotFiles;
	string tmp;
	string fileNameListName = "runData/fileNameList.dat";
	ifstream fileNameList(fileNameListName);

	if(fileNameList.is_open()){
		while (getline (fileNameList,tmp)){
			snapShotFiles.push_back(tmp);
		}
	}		

	#pragma omp parallel for
	for(int j = 0; j < snapShotFiles.size(); j++){
		
		string h5name = snapShotFiles[j];
		Options opt;
		Eval results;
		vector<MatrixXcd> wavefunction;	

		cout << "Opening Datafiles.." << h5name << endl;
		binaryFile data(h5name,binaryFile::in);	

		cout << "Reading Datafiles.. " << h5name << endl;
		vector<int> timeList = data.getTimeList();
		for(int i = 0; i < timeList.size(); i++){
			data.getSnapshot(runname,timeList[i],wavefunction,opt);		
			results.saveData(wavefunction,opt,timeList[i],runname);		
			cout << "Evaluating Datafiles.. "<< timeList[i] << endl;
			results.evaluateData();		
			cout << "Plotting Datafiles.. " << timeList[i] << endl;		
			results.plotData();
		}
	}	
	
	cout << "Terminating successfully." << endl;
	if(DEBUG_LOG == 1){
		logstream.close();
	}
// Everything finished here 
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




