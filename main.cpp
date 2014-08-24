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
 		redirecter redirect(logstream,std::cout); // redirects cout to logstream, until termination of this program. If DEBUG_LOG 1 is set, use cerr for output to console.
 	#endif

 	startUp.printInitVar();

 	// std::streambuf *filebuf, *coutbuf;

 	// coutbuf = std::cout.rdbuf();     // back up cout's streambuf
 	// std::cout.rdbuf(backup);        // restore cout's original streambuf
	
	// if(DEBUG_LOG == 1){
	// 	logstream.open ("run.log");
		// filebuf = logstream.rdbuf();        // get file's streambuf
		// std::cout.rdbuf(filebuf);         // assign streambuf to cout
	// }
	
	string runName = "RTE";


	MatrixData* data = new MatrixData(startUp.getMeta());

	if(startUp.newRun == false){
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
		binaryFile* loading = new binaryFile(h5name,binaryFile::in);
		vector<int> timeList = loading->getTimeList();
		Options opt = startUp.getOptions();
		loading->getSnapshot(runName,timeList.back(),data,opt);
		delete loading;

		RTE* run = new RTE(data,opt);
		run->rteToTime(runName);
		delete run;

	} else {
		
		setGridToGaussian(data,startUp.getOptions());

		ITP* itprun = new ITP(data->wavefunction[0],startUp.getOptions());
		string itpname = "ITP-Groundstate";
		itprun->propagateToGroundState(itpname);
		
		startUp.setVortexnumber(addVortices(data,startUp.getOptions()));

		itpname = "ITP-Vortices";
		itprun->formVortices(itpname);
			
		for(int i = 0; i < data->meta.samplesize; i++){
			data->wavefunction[i] = itprun->result();
		}
		delete itprun;
	
		RTE* run = new RTE(data,startUp.getOptions());
		
		run->noise();
		run->rteToTime(runName);
		delete run;
	}
	


	// cout << "Deleting objects." << endl;
	delete data;	

	cout << "Terminating successfully." << endl;
	// if(DEBUG_LOG == 1){
	// 	logstream.close();
	// 	std::cout.rdbuf(coutbuf);
	// }
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




