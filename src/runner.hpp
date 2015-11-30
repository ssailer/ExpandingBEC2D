#ifndef EXP2D_RUNNER_H__
#define EXP2D_RUNNER_H__

#define EIGEN_VECTORIZE
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_NO_DEBUG


#include <iostream>
#include <complex>
#include <math.h>
#include <vector>
#include <omp.h>
#include <string>
#include <iomanip>

#include <.archive/utilities/gauss_random.h>

#include <tools.h>
#include <evaluation.h>
#include <binaryfile.h>
#include <rk4.hpp>
#include <splitstep.hpp>
#include <constants.h>
#include <plot_with_mgl.h>
#include <matrixdata.h>
#include <plotter.hpp>
#include <hydro.h>
#include <eigen3/Eigen/Dense>

#define RESIZE 200

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
    Runner(shared_ptr<MatrixData> d, Options &o);
    ~Runner();

    void setOptions(const Options &externaloptions);
    void RunSetup();
    
    // Propagatoren

    void runToTime(string runName);
    // void splitToTime(string runName);
   
    // StoragePointer to the wavefunction (MatrixData Object)
    shared_ptr<MatrixData> pData;

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

   inline complex<double> lambda_y(const complex<double> &t){
    return sqrt(one+opt.exp_factor*opt.dispersion_y*opt.dispersion_y*t*t);
   }

   inline complex<double> lambda_y_dot(complex<double> &t){
    return (opt.exp_factor*opt.dispersion_y*opt.dispersion_y*t/sqrt(one+opt.exp_factor*opt.dispersion_y*opt.dispersion_y*t*t));
   }
  
};


template<class T>
Runner<T>::Runner(shared_ptr<MatrixData> d, Options &o) : wavefctVec(d->wavefunction), meta(d->meta), pData(d)
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
	snapshot_times.resize(opt.snapshots);
	for(int k = 0; k < opt.snapshots; k++){
		snapshot_times[k] = (k + 1) * opt.n_it_RTE + meta.steps;
	}

	// double mu = sqrt(3.0  * opt.g * real(opt.omega_x) * real(opt.omega_y) * opt.N / 8.0);
 //    double Ry = sqrt(2.0 * mu / ( real(opt.omega_y)*real(opt.omega_y))) * opt.Ag;
 //    double Rx = sqrt(2.0 * mu / ( real(opt.omega_x)*real(opt.omega_x))) * opt.Ag;

 //    cout << "Initial Thomas Fermi Radii set to Rx = " << Rx << " and Ry = " << Ry << endl;
 //    double n0 = 2 * (opt.N / M_PI) * (1 / (Rx * Ry));
 //    cout << "n_0 = " << n0 << endl;


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
		cerr << "\r";
		cerr << currentTime() <<  " " << name << " "
		 	 << std::setw(2) << std::setfill('0') << hour << ":"
			 << std::setw(2) << std::setfill('0') << min << ":"
			 << std::setw(2) << std::setfill('0') << seconds  << "    Steps: "
			 << std::setw(3) << std::setfill('0') << pData->meta.steps << " / " << snapshot_times.back()
			 
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
		Xexpanding = x_expand(complex<double>(pData->meta.time,0.0));
		Yexpanding = y_expand(complex<double>(pData->meta.time,0.0));
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
	double start;  // starttime of the run
	// int resizeByInt = 512;

	// Eval* initEval = new Eval(pData,opt);
	{
		auto initEval = std::make_shared<Eval>(pData,opt);
		initEval->process();
		initEval->save();
		vector<int> edges;
		bool should_i_resize = initEval->checkResizeCondition(edges);
	
		if(opt.initialRun == true){
	
			Plotter* initPlot = new Plotter(initEval,opt);
			initPlot->plotEval();
			delete initPlot;
	
			opt.vortexnumber = initEval->getVortexNumber();
			opt.initialRun = false;
	
			string filename = runName + "data.h5";
			string obsname  = runName + "obs.h5";
	
			binaryFile* bFile = new binaryFile(filename,binaryFile::out);
			bFile->appendSnapshot("MatrixData",pData,opt);
			binaryFile* obsFile = new binaryFile(obsname,binaryFile::out);
			obsFile->appendEval(initEval,opt);
			delete bFile;
			delete obsFile;
		}

		initEval.reset();

		if(opt.runmode == "EXP"){
			// hydroSolver solver(initEval);
			// solver.integrate();
			// solver.pyPlot();
	
			vector<int> edges;
			
			if(should_i_resize){
				pData->resizeBy(RESIZE);
				algorithm->setVariables();
				cout << "Resizing by " << RESIZE << endl;
			}
		}
	}

	start = omp_get_wtime();

	cli(runName,pData->meta.steps,start);

	for(int j = 0; j < snapshot_times.size(); j++){

		while(pData->meta.steps < snapshot_times[j]){			
				
			algorithm->timeStep(opt.RTE_step);

			// opt.t_abs += opt.RTE_step;

			cli(runName,pData->meta.steps,start);
		}

		// REMOVE
			if(opt.runmode == "EXP"){
				opt.stateInformation[0] = real(lambda_x(complex<double>(pData->meta.time,0.0))); // needed for expansion and the computing of the gradient etc.
				opt.stateInformation[1] = real(lambda_y(complex<double>(pData->meta.time,0.0)));
			}
			else {
				opt.stateInformation[0] = 1.0;
				opt.stateInformation[1] = 1.0;
			}
		
		// END REMOVE

		try{

			// if(opt.runmode != "EXP"){
				// Eval* eval = new Eval(pData,opt);
			auto eval = std::make_shared<Eval>(pData,opt);
			eval->process();
			eval->save();
			vector<int> edges;
			bool should_i_resize = eval->checkResizeCondition(edges);
			// }

			string dataname = runName + "data.h5";
			string obsname = runName + "obs.h5";
			binaryFile* bFile = new binaryFile(dataname,binaryFile::append);
			bFile->appendSnapshot("MatrixData",pData,opt);
			delete bFile;

			binaryFile* obsFile = new binaryFile(obsname,binaryFile::append);
			obsFile->appendEval(eval, opt);
			delete obsFile;

			// if(opt.runmode != "EXP"){
				Plotter* plotter = new Plotter(eval,opt);
				plotter->plotEval();
				delete plotter;
			// }
			// cout << "opt.RTE_step " << opt.RTE_step << endl;
			// cout << "spacing " << pData->meta.spacing[0] << endl;
			// cout << "grid " << pData->meta.grid[0] << endl;
			// cout << "size " << pData->meta.coord[0] << endl;

			eval.reset();

			if(opt.runmode == "EXP"){			

				if(should_i_resize){
					pData->resizeBy(RESIZE);
					algorithm->setVariables();
					cout << "Resizing by " << RESIZE << endl;
				}
				// hydroSolver solver;
				// solver.pyPlot();
			}

			// cout << "opt.RTE_step " << opt.RTE_step << endl;
			// cout << "spacing " << pData->meta.spacing[0] << endl;
			// cout << "grid " << pData->meta.grid[0] << endl;
			// cout << "size " << pData->meta.coord[0] << endl;


			// delete eval;


		}
		catch(const std::exception& e) { 
			std::cerr 	<< "Unhandled Exception after dataFile.appendSnapshot() in rteToTime: " << std::endl; 
			throw e; 
		}

		cli(runName,pData->meta.steps,start);
	}
}



#endif // EXP2D_RUNNER_H__