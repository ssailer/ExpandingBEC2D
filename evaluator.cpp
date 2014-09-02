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
	
	vector<vector<Observables>> obs;	
	obs.resize(4);
	vector<Options> opt;
	MatrixData::MetaData meta;
	Eval* results = new Eval;

	vector<int> timeList;

	for(int k = 1; k <= 4; k++){

		string runName = "Expanding-Set-"+to_string(k);
		string evalname = runName + "-Eval.h5";

		binaryFile* evalFile = new binaryFile(evalname,binaryFile::in);

		timeList = evalFile->getTimeList();

		cout << "Loaded file " << evalname << " from runData/" << endl;

		for(int j = 0; j< timeList.size(); j++){
			Options tmpOpt;
			evalFile->getEval(timeList[j],tmpOpt,meta,*results);
			obs[k-1].push_back(results->totalResult);
			if(k == 1){
				opt.push_back(tmpOpt);
			}
		}
		delete evalFile;
	}

	delete results;

	vector<Observables> finalObs;
	for(int k = 0; k < obs[0].size(); k++){
		Observables tmpObs = obs[0][k] + obs[1][k] + obs[2][k] + obs[3][k];
		tmpObs /= 4;
		finalObs.push_back(tmpObs);
	}

	if(finalObs.size() != opt.size()){
		cout << "Obs and Opt are not the same size." << endl;
	}

	string filename = "Combined_Observables.dat";

	for(int i = 0; i < finalObs.size(); i++){
		
		struct stat buffer;   
	  	if(stat (filename.c_str(), &buffer) != 0){
	  		ofstream datafile;
	  		datafile.open(filename.c_str(), ios::out | ios::app);
	  		datafile << std::left << "Timestep"
	  						 << "," << "X_max"
	  						 << "," << "Y_max"
	  						 << "," << "D_max"
	  						 << "," << "D_min"
	  						 << "," << "D_Ratio"
	  						 << "," << "D_max_Angle"
	  						 << "," << "D_min_Angle"
	  						 << "," << "Ratio"
	  						 << "," << "RatioAngle"
	  						 << "," << "N"
	  						 << "," << "V"
	  						 << "," << "Density"
	  						 << "," << "E_kin"
	  				 << endl;
	  		datafile.close();
	  	} 
	
	  	ofstream datafile(filename.c_str(), std::ios_base::out | std::ios_base::app);
		// datafile.open;
		datafile << std::left << timeList[i]
						 << "," << opt[i].min_x * opt[i].stateInformation[0]
						 << "," << opt[i].min_y * opt[i].stateInformation[1]
	 					 << "," << finalObs[i].r_max
	 					 << "," << finalObs[i].r_min
	 					 << "," << finalObs[i].r_max / finalObs[i].r_min  
	 					 << "," << finalObs[i].r_max_phi
	 					 << "," << finalObs[i].r_min_phi
	 					 << "," << finalObs[i].aspectRatio 
	 					 << "," << finalObs[i].aspectRatioAngle 
						 << "," << finalObs[i].particle_count
						 << "," << finalObs[i].volume
						 << "," << finalObs[i].density
						 << "," << finalObs[i].Ekin
				 << endl;
		datafile.close();
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
	return ERROR_UNHANDLED_EXCEPTION; 
// the code could be different depending on the exception message 
}
return SUCCESS; 	
}




