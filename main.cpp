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
#include <dirent.h>

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
#include <EXP2D_startgrids.h>

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
	StartUp startUp(argc,argv);	

	#if DEBUG_LOG
 		std::ofstream logstream("run.log");
 		redirecter redirectcout(logstream,std::cout); // redirects cout to logstream, until termination of this program. If DEBUG_LOG 1 is set, use cerr for output to console.
 		// std::ofstream errorstream("error.log");
 		// redirecter redirectcerr(errorstream,std::cerr);
 	#endif

 	startUp.printInitVar();
	
	MatrixData* startGrid = new MatrixData(startUp.getMeta());
		
	setGridToGaussian(startGrid,startUp.getOptions());

	// cout << "value " << startGrid->wavefunction[0](1024,1024) << endl;

	ITP* groundStateITP = new ITP(startGrid->wavefunction[0],startUp.getOptions());
	string itpname = "ITP-Groundstate";
	groundStateITP->propagateToGroundState(itpname);
	startGrid->wavefunction[0] = groundStateITP->result();
	delete groundStateITP;

	string tmpRunMode = startUp.getRunMode();
	if(tmpRunMode.compare(3,1,"1") == 0){
		int vnumber = 0;
		addVorticesRegular(startGrid,startUp.getOptions(),vnumber);
		
		startUp.setVortexnumber(vnumber);
		cout << endl << "Set Vortices #: " << vnumber << endl;
	
		itpname = "ITP-Vortices";
		ITP* vorticesITP = new ITP(startGrid->wavefunction[0],startUp.getOptions());
		vorticesITP->formVortices(itpname);
		
		startGrid->wavefunction[0] = vorticesITP->result();
	
		delete vorticesITP;
	}

	// FIXME: To run RTE multiple times, go into RTE::RunSetup() and fix the expanding coordinates starting procedure. It has to be loaded from metaData, instead of calculating directly, not only the time.

	for( int k = 1; k <= 4; k++){
		MatrixData* data = new MatrixData(startUp.getMeta());
		for(int i = 0; i < data->meta.samplesize; i++){
			data->wavefunction[i] = startGrid->wavefunction[0];
		}	
		string runName = "Expanding-Set-"+to_string(k);
		RTE* runExpanding = new RTE(data,startUp.getOptions());
		runExpanding->noise();
		runExpanding->rteToTime(runName);
		delete runExpanding;
		delete data;
	}
	// }
	
	delete startGrid;	

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
cerr << "Run complete. Terminating successfully." << endl; 
return SUCCESS; 	
}




