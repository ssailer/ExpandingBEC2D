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
	if(DEBUG_LOG == 1){
 		std::ofstream logstream("evaluator.log");
 		redirecter redirect(logstream,std::cout); // redirects cout to logstream, until termination of this program. If DEBUG_LOG 1 is set, use cerr for output to console.
 	}
	int files = 1;
	vector<vector<Observables>> obs;	
	obs.resize(files);
	vector<Options> opt;
	MatrixData::MetaData meta;
	vector<vector<Eval>> results;
	results.resize(files);

	vector<int> timeList;

	for(int k = 0; k < files; k++){

		string runName = "Expanding-Set-"+to_string(k+1);
		string evalname = runName + "-Eval.h5";

		binaryFile* evalFile = new binaryFile(evalname,binaryFile::in);

		timeList = evalFile->getTimeList();

		results[k].resize(timeList.size());

		cout << "Loaded file " << evalname << " from runData/" << endl;

		for(int j = 0; j< timeList.size(); j++){
			Options tmpOpt;
			evalFile->getEval(timeList[j],tmpOpt,meta,results[k][j]);
			// obs[k-1].push_back(results[k].totalResult);
			if(k == 0){
				opt.push_back(tmpOpt);
			}
		}
		delete evalFile;
	}


	string finalRunName = "Expanding";
	Eval finalResult;
	for(int i = 0; i < timeList.size(); i++){
		cout << "Processing Time: " << timeList[i] << " .. " ;
		vector<Eval> tmpResults(files);
		for(int f = 0; f < files; f++){
			tmpResults[f] = results[f][i];
		}
		finalResult.saveDataFromEval(opt[i],timeList[i],finalRunName,tmpResults);
		cout << "\r" << flush;

	}
	cout << endl;
	cout << "Evaluation finished" << endl;

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




