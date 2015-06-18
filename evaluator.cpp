/**************************************************************************
Title: Simulating the Expansion of Turbulent Bose-Einstein Condensates (2D) 
Author: Simon Sailer (This work is based on the work of Bartholomew Andrews who made this as his master thesis.)
Last Update: 22/07/13
**************************************************************************/

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
#define DEBUG_LOG 0

using namespace std;

int main( int argc, char** argv) 
{	
try{
	// if(DEBUG_LOG == 1){
 // 		std::ofstream logstream("evaluator.log");
 // 		redirecter redirect(logstream,std::cout); // redirects cout to logstream, until termination of this program. If DEBUG_LOG 1 is set, use cerr for output to console.
 // 	}
	// int files = 1;
	// vector<vector<Observables>> obs;	
	// obs.resize(files);
	// vector<Options> opt;
	// MatrixData::MetaData meta;
	// vector<vector<Eval>> results;
	// results.resize(files);

	// vector<int> timeList;
 	string filename = "rotdata.h5";
	MatrixData* data = new MatrixData();
	Options opt;
	

	binaryFile* dataFile = new binaryFile(filename,binaryFile::in);	
	vector<int> timeList = dataFile->getTimeList();
	int size = timeList.size();
	dataFile->getLatestSnapshot("MatrixData",data,opt);
	opt.isDimensionless = true;
	data->meta.Ag = opt.Ag;
	data->meta.OmegaG = opt.OmegaG;
	Eval eval(*data,opt);
	Observables obs;
	int counter = 0;

	dataFile->getEval(700000,eval,opt);
	obs = eval.totalResult;
	for(int k = 70; k < size; k++){
		dataFile->getEval(timeList[k],eval,opt);
		// cout << "timeList " << timeList[k] << endl;
		obs += eval.totalResult;
		counter++;
	}
	// dataFile->getLatestSnapshot("MatrixData",data,opt);
	obs /= (counter+1);
	eval.totalResult = obs;

		// spectrum
	vector<double> kval;
	vector<double> numberval;
	map<double,double> spectrum;
	pair<map<double,double>::iterator,bool> ret;    
    vector<double> tmpKval;
		
    for (int r = 0; r < eval.totalResult.number.size(); r++){
		if(eval.totalResult.k(r) != 0.0){
			if(eval.totalResult.number(r) != 0.0){
				ret = spectrum.insert(map<double,double>::value_type(eval.totalResult.k(r),eval.totalResult.number(r)));
				tmpKval.push_back(eval.totalResult.k(r));
				if(ret.second==false){
					cout << "Binning of spectrum failed, double value inserted." << endl;
				}

				// kval.push_back(k_int);
				// numberval.push_back(eval.totalResult.number(r));
			}
        }
	}

	// // estimate powerlaw
	// double gamma = 3.0;
	double k_max = 4.0;
	double k_min = 2.0;

	vector<double> klog;
	vector<double> nlog;
	for(int i = 0; i < eval.totalResult.k.size(); ++i){
		if(eval.totalResult.k(i) != 0.0 && eval.totalResult.number(i)){
			if(eval.totalResult.k(i) <= k_max && k_min <= eval.totalResult.k(i)){
				klog.push_back(log(eval.totalResult.k(i)));
				nlog.push_back(log(eval.totalResult.number(i)));
			}
		}
	}


	

	auto tmpMinMax = std::minmax_element(klog.begin(),klog.end());

	double linksX = *tmpMinMax.first;
	double rechtsX = *tmpMinMax.second;

   double SUMx = 0;     //sum of x values
   double SUMy = 0;     //sum of y values
   double SUMxy = 0;    //sum of x * y
   double SUMxx = 0;    //sum of x^2
   double SUMres = 0;   //sum of squared residue
   double res = 0;      //residue squared
   double slope = 0;    //slope of regression line
   double y_intercept = 0; //y intercept of regression line
   double SUM_Yres = 0; //sum of squared of the discrepancies
   double AVGy = 0;     //mean of y
   double AVGx = 0;     //mean of x
   double Yres = 0;     //squared of the discrepancies
   double Rsqr = 0;     //coefficient of determination
   int dataSize = nlog.size();
   //calculate various sums 
   for (int i = 0; i < dataSize; i++)
   {
      //sum of x
      SUMx = SUMx + klog[i];
      //sum of y
      SUMy = SUMy + nlog[i];
      //sum of squared x*y
      SUMxy = SUMxy + klog[i] * nlog[i];
      //sum of squared x
      SUMxx = SUMxx + klog[i] * klog[i];
   }

   //calculate the means of x and y
   AVGy = SUMy / dataSize;
   AVGx = SUMx / dataSize;

   //slope or a1
   slope = (dataSize * SUMxy - SUMx * SUMy) / (dataSize * SUMxx - SUMx*SUMx);

   //y intercept or a0
   y_intercept = AVGy - slope * AVGx;

   eval.totalResult.alpha = slope;


   eval.steigung = slope;
   double abschnitt = y_intercept;
   double linksY = slope * linksX + y_intercept;
   double rechtsY = slope * rechtsX + y_intercept;

   eval.punkte.push_back(exp(linksX));
   eval.punkte.push_back(exp(linksY));
   eval.punkte.push_back(exp(rechtsX));
   eval.punkte.push_back(exp(rechtsY));


   // printf("x mean(AVGx) = %0.5E\n", AVGx);

   // printf("y mean(AVGy) = %0.5E\n", AVGy);

   // printf ("\n");
   // printf ("The linear equation that best fits the given data:\n");
   // printf ("       y = %2.8lfx + %2.8f\n", slope, y_intercept);
   // printf ("------------------------------------------------------------\n");
   // printf ("   Original (x,y)   (y_i - y_avg)^2     (y_i - a_o - a_1*x_i)^2\n");
   // printf ("------------------------------------------------------------\n");

   // //calculate squared residues, their sum etc.
   for (int i = 0; i < dataSize; i++) 
   {
      //current (y_i - a0 - a1 * x_i)^2
      Yres = pow(nlog[i] - y_intercept - (slope * (klog[i])), 2);

      //sum of (y_i - a0 - a1 * x_i)^2
      SUM_Yres += Yres;

      //current residue squared (y_i - AVGy)^2
      res = pow(nlog[i] - AVGy, 2);

      //sum of squared residues
      SUMres += res;
      
      // printf ("   (%0.2f %0.2f)      %0.5E         %0.5E\n", 
      //  klog[i], nlog[i], res, Yres);
   }

   eval.fehler = sqrt(SUM_Yres / (dataSize - 2));

	Plotter plot(eval);
	plot.spectrum();
	plot.alphas();
	delete dataFile;
	delete data;

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




