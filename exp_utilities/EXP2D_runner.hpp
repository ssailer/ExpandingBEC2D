#ifndef EXP2D_RUNNER_H__
#define EXP2D_RUNNER_H__

#define EIGEN_VECTORIZE
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_NO_DEBUG


#include <iostream>
#include <complex>
#include <math.h>
#include <complexgrid.h>
#include <bh3binaryfile.h>
#include <gauss_random.h>
#include <vector>
#include <omp.h>
#include <string>
#include <iomanip>
#include <gauss_random.h>
#include <EXP2D_tools.h>
#include <EXP2D_evaluation.h>
#include <EXP2D_binaryfile.h>
#include <EXP2D_rk4.hpp>
#include <EXP2D_splitstep.hpp>
#include <EXP2D_constants.h>
#include <plot_with_mgl.h>
#include <EXP2D_MatrixData.h>
#include <EXP2D_plotter.hpp>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

typedef struct {
    int absoluteSteps;
    int lambdaSteps;
    int initialSteps;
} stepCounter;

template<class T>
class Runner{
public:
    // Runner(MatrixData* &d, RungeKutta* r, Options &opt);  
    // Runner(MatrixData* &d, SplitStep* s, Options &externaloptions);
    Runner(MatrixData* &d, Options &o);
    ~Runner();

    void setOptions(const Options &externaloptions);
    void RunSetup();
    
    // Propagatoren

    void runToTime(string runName);
    // void splitToTime(string runName);
   
    // StoragePointer to the wavefunction (MatrixData Object)
    MatrixData* pData;

    // Storage Variable for the runs

    vector<MatrixXcd> &wavefctVec;
    MatrixData::MetaData &meta;

    T* algorithm;

    // internal RunOptions, use setOptions(Options) to update from the outside
    Options opt;    

    // virtual void singleK(MatrixXcd &k, MatrixXcd &wavefctcp, int32_t &front, int32_t &end,int32_t &subx,int32_t & suby, int &t);

    vector<int> snapshot_times;
    // Coordinates
    vector<double> x_axis,y_axis;
    VectorXcd X,Y;
    MatrixXcd Xmatrix,Ymatrix;
    VectorXd Xexpanding, Yexpanding;

    int samplesize;

    // Plotting and progress functions 
    
    // void cli_plot(string name,int counter_state, int counter_max, double start,bool plot);
    void cli(string name, int index, double start);
    void plot(const string name);
    void noise();

    // inline double rotatingPotential(int &i, int &j, int &t);

    // void ComputeDeltaPsi(MatrixXcd &wavefct, MatrixXcd &wavefctcp, int &t,complex<double> delta_T);
    
    // void MSDBoundaries(MatrixXcd &U,MatrixXcd &Ut);
   

    // Variables
    complex<double> h_x, h_y;
    complex<double> pot_laplacian_x;
    complex<double> pot_laplacian_y;
    vector<double> ranges;
    // VectorXcd t_RTE;
    Matrix<std::complex<double>,Dynamic,Dynamic,ColMajor> wavefctcpX;
    Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor> wavefctcpY;
    MatrixXcd PotentialGrid,AbsorbingPotentialGrid;
    
    VectorXcd laplacian_coefficient_x,laplacian_coefficient_y,gradient_coefficient_x,gradient_coefficient_y;

   // little helper functions for stuff

   inline VectorXd x_expand(complex<double> &t){
    return X.real() * real(lambda_x(t));
   }

   inline VectorXd y_expand(complex<double> &t){
    return Y.real() * real(lambda_y(t));
   }

   inline complex<double> lambda_x(complex<double> &t){
    return sqrt(one+opt.exp_factor*opt.dispersion_x*opt.dispersion_x*t*t);
   }

   inline complex<double> lambda_x(const complex<double> &t){
    return sqrt(one+opt.exp_factor*opt.dispersion_x*opt.dispersion_x*t*t);
   }

   inline complex<double> lambda_x_dot(complex<double> &t){
    return (opt.exp_factor*opt.dispersion_x*opt.dispersion_x*t/sqrt(one+opt.exp_factor*opt.dispersion_x*opt.dispersion_x*t*t));
   }

   inline complex<double> lambda_y(complex<double> &t){
    return sqrt(one+opt.exp_factor*opt.dispersion_y*opt.dispersion_y*t*t);
   }

   inline complex<double> lambda_y_dot(complex<double> &t){
    return (opt.exp_factor*opt.dispersion_y*opt.dispersion_y*t/sqrt(one+opt.exp_factor*opt.dispersion_y*opt.dispersion_y*t*t));
   }
  
};


template<class T>
Runner<T>::Runner(MatrixData* &d, Options &o) : wavefctVec(d->wavefunction), meta(d->meta), pData(d)
{
	algorithm = new T(o);
	algorithm->assignMatrixData(pData);
	setOptions(o);
	RunSetup();
}

template<class T>
Runner<T>::~Runner(){
	delete algorithm;
}
// void Runner::setAlgorithm()

template<class T>
void Runner<T>::setOptions(const Options &o)
{
	opt = o;
}

template<class T>
void Runner<T>::RunSetup()
{
	int snapShotSize = opt.n_it_RTE / opt.snapshots;
	int nbTrueSnapShots =	opt.snapshots - meta.steps / snapShotSize; 

	snapshot_times.resize(nbTrueSnapShots);
	for(int k = 0; k < nbTrueSnapShots; k++){
		snapshot_times[k] = (k + 1) * snapShotSize + meta.steps;
	}

	double mu = sqrt(3.0  * opt.g * real(opt.omega_x) * real(opt.omega_y) * opt.N / 8.0);
    double Ry = sqrt(2.0 * mu / ( real(opt.omega_y)*real(opt.omega_y))) * opt.Ag;
    double Rx = sqrt(2.0 * mu / ( real(opt.omega_x)*real(opt.omega_x))) * opt.Ag;

    cout << "Initial Thomas Fermi Radii set to Rx = " << Rx << " and Ry = " << Ry << endl;
    double n0 = 2 * (opt.N / M_PI) * (1 / (Rx * Ry));
    cout << "n_0 = " << n0 << endl;


}

template<class T>
void Runner<T>::cli(string name, int index, double start)
{	
	// if(fmod((float)counter,(float)(counter_max/100))==0){
		int seconds, min, hour, total, expectedhour, expectedmin, expectedseconds;
		// double totalstate = 0;
		// double totalmaxpercent = (double)counter_max * (double)meta.samplesize / 100;
		// for(int i = 0; i < meta.samplesize; i++){
		// 	totalstate += stateOfLoops[i];
		// }
		// double totalPercent = totalstate/totalmaxpercent;

		// int overallStepState = keeperOfTime.absoluteSteps + totalstate / meta.samplesize;

		total = omp_get_wtime() - start;

		// overallStepState = (overallStepState == 0) ? 1 : overallStepState;

		// int remainingSeconds = (total * opt.n_it_RTE / pData->meta.steps + 1) - total;


		
		hour = total / 3600;
		min = (total / 60) % 60;
		seconds = total % 60;
		// expectedhour = (remainingSeconds / 3600);
		// expectedmin = (remainingSeconds / 60) % 60;
		// expectedseconds = remainingSeconds % 60;
		cout << "\r";
		cout << currentTime() <<  " " << name << " "
		 	 << std::setw(2) << std::setfill('0') << hour << ":"
			 << std::setw(2) << std::setfill('0') << min << ":"
			 << std::setw(2) << std::setfill('0') << seconds  << "    Steps: "
			 << std::setw(3) << std::setfill('0') << pData->meta.steps << " / " << opt.n_it_RTE
			 
			 // << " remaining runtime: "
			 // << std::setw(2) << std::setfill('0') << expectedhour << ":"
			 // << std::setw(2) << std::setfill('0') << expectedmin << ":"
			 // << std::setw(2) << std::setfill('0') << expectedseconds
		// }
		<< "    " << flush;
	// }
}

template<class T>
void Runner<T>::plot(const string name){

	if(opt.runmode == "EXP"){
		Xexpanding = x_expand(opt.t_abs);
		Yexpanding = y_expand(opt.t_abs);
		plotDataToPngEigenExpanding(name, wavefctVec[0],ranges,Xexpanding,Yexpanding,opt);
	}
	else {
		plotDataToPngEigen(name, wavefctVec[0],opt);
	}
}

template<class T>
void Runner<T>::noise(){
	for(int k = 0; k < wavefctVec.size(); k++){
		GaussRandom r (get_seed());
		double rvalue;
		for(int i = 0;i < wavefctVec[k].rows();i++){
			for(int j = 0; j < wavefctVec[k].cols();j++){
				rvalue = real(wavefctVec[k](i,j)) * 0.1;
				wavefctVec[k](i,j) += r.gauss_random(0.0,rvalue);
			}
		}
	}
}

template<class T>
void Runner<T>::runToTime(string runName)
{	
	// omp_set_num_threads(16);

	double start;  // starttime of the run

	if(opt.initialRun == true){

		Eval* initEval = new Eval(*pData,opt);
		initEval->process();

		Plotter* initPlot = new Plotter(*initEval,opt);
		initPlot->plotEval();
		delete initPlot;

		opt.vortexnumber = initEval->getVortexNumber();
		opt.initialRun = false;

		string evalname = "Evaluation.h5";
		binaryFile* evalFile = new binaryFile(evalname,binaryFile::out);
		evalFile->appendEval(meta.steps,opt,meta,*initEval);
		delete initEval;
		delete evalFile;
	}

	start = omp_get_wtime();

	cli(runName,pData->meta.steps,start);

	for(int j = 0; j < snapshot_times.size(); j++){

		while(pData->meta.steps < snapshot_times[j]){			
				
			algorithm->timeStep(opt.RTE_step);

			opt.t_abs += opt.RTE_step;
		}

		cli(runName,pData->meta.steps,start);

		// REMOVE
			if(opt.runmode == "EXP"){
				opt.stateInformation[0] = real(lambda_x(opt.t_abs)); // needed for expansion and the computing of the gradient etc.
				opt.stateInformation[1] = real(lambda_y(opt.t_abs));
				cout << opt.stateInformation[0] << "  " << opt.t_abs << " ";
				cout << pData->meta.initCoord[0] << " " << pData->meta.coord[0] << " " << pData->meta.time;
			}
			else {
				opt.stateInformation[0] = 1.0;
				opt.stateInformation[1] = 1.0;
			}
		// END REMOVE

		try{
			Eval* eval = new Eval(*pData,opt);
			eval->process();

			Plotter* plotter = new Plotter(*eval,opt);
			plotter->plotEval();
			delete plotter;

			string dataname = "LastGrid.h5";
			binaryFile* dataFile = new binaryFile(dataname,binaryFile::out);
			dataFile->appendSnapshot(runName,pData->meta.steps,pData,opt);
			delete dataFile;

			string evalname = "Evaluation.h5";
			binaryFile* evalFile = new binaryFile(evalname,binaryFile::append);
			evalFile->appendEval(pData->meta.steps,opt,pData->getMeta(),*eval);
			delete evalFile;

			delete eval;
		}
		catch(const std::exception& e) { 
			std::cerr 	<< "Unhandled Exception after dataFile.appendSnapshot() in rteToTime: " << std::endl; 
			throw e; 
		}
	}
}



#endif // EXP2D_RUNNER_H__