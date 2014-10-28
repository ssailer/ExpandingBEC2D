
#include <EXP2D_splitstep.hpp>
#include <omp.h>

#define EIGEN_DONT_VECTORIZE
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_NO_DEBUG

using namespace std;
using namespace Eigen;

void SplitStep::cli(string name,int &slowestthread, vector<int> threadinfo, vector<int> stateOfLoops, int counter_max, double start)
{	
	for(int i = 0;i < stateOfLoops.size();i++){
		slowestthread = (stateOfLoops[slowestthread] <= stateOfLoops[i]) ? slowestthread : threadinfo[i];
	}
	
	if(fmod((float)stateOfLoops[slowestthread],(float)(counter_max/10))==0){
		int seconds, min, hour, total, expectedhour, expectedmin, expectedseconds;
		double totalstate = 0;
		double totalmaxpercent = (double)counter_max * (double)meta.samplesize / 100;
		for(int i = 0; i < meta.samplesize; i++){
			totalstate += stateOfLoops[i];
		}
		double totalPercent = totalstate/totalmaxpercent;

		int overallStepState = keeperOfTime.absoluteSteps + totalstate / meta.samplesize;

		total = omp_get_wtime() - start;

		overallStepState = (overallStepState == 0) ? 1 : overallStepState;

		int remainingSeconds = (total * opt.n_it_RTE / overallStepState) - total;


		
		hour = total / 3600;
		min = (total / 60) % 60;
		seconds = total % 60;
		expectedhour = (remainingSeconds / 3600);
		expectedmin = (remainingSeconds / 60) % 60;
		expectedseconds = remainingSeconds % 60;
		cout << "\r";
		cout << currentTime() <<  " " << name << " "
		 	 << std::setw(2) << std::setfill('0') << hour << ":"
			 << std::setw(2) << std::setfill('0') << min << ":"
			 << std::setw(2) << std::setfill('0') << seconds  << "    "
			 << std::setw(3) << std::setfill('0') << (int)totalPercent << "% "
			 << " threads: " << stateOfLoops.size()
			 << " remaining runtime: "
			 << std::setw(2) << std::setfill('0') << expectedhour << ":"
			 << std::setw(2) << std::setfill('0') << expectedmin << ":"
			 << std::setw(2) << std::setfill('0') << expectedseconds
			// cout << " | Slowest Thread: " << std::setw(3) << std::setfill('0') << (float)(stateOfLoops[slowestthread])/(float)(counter_max/10) << "% ";
			// for(int k = 0; k < stateOfLoops.size(); k++){
			// cout << k << "_" << threadinfo[k] << ": " << std::setw(3) << std::setfill('0') << (float)stateOfLoops[k]/((float)counter_max/100) << "% ";
		// }
		<< "    " << flush;
	}
}

void SplitStep::plot(const string name){

	if(opt.runmode.compare(1,1,"1") == 0){
		complex<double> tmp = complex<double>(keeperOfTime.absoluteSteps,0.0) * t_RTE;
		Xexpanding = x_expand(tmp);
		Yexpanding = y_expand(tmp);
		plotDataToPngEigenExpanding(name, wavefctVec[0],ranges,Xexpanding,Yexpanding,opt);
	}
	if(opt.runmode.compare(1,1,"0") == 0){
		plotDataToPngEigen(name, wavefctVec[0],opt);
	}
}

SplitStep::SplitStep(MatrixData* &d,const Options &externaloptions) : wavefctVec(d->wavefunction), meta(d->meta), pData(d)
{	
	// Both essential Variables
	// pData = d;
  	setOptions(externaloptions);

  	// some constants used in computations to shorten stuff
	pi = M_PI;
 	zero=complex<double>(0,0);
 	half=complex<double>(0.5,0);
 	one=complex<double>(1,0);
 	two=complex<double>(2,0);
 	four=complex<double>(4,0);
 	six=complex<double>(6,0);
 	i_unit=complex<double>(0,1);

 	// setting up multithreading. Output to see what Eigen is doing.
	// omp_set_num_threads(omp_get_max_threads());
	cout << "Max Number of Threads in RTE: " << omp_get_max_threads() << endl;	
	// cout << "Eigenthreads: " << Eigen::nbThreads() << endl;

	// Using the setter function to initialize the stuff.
	RunSetup();

}

void SplitStep::setOptions(const Options &externaloptions){
	opt = externaloptions;
}


void SplitStep::RunSetup(){

	int snapShotSize = opt.n_it_RTE / opt.snapshots;
	int nbTrueSnapShots =	opt.snapshots - meta.steps / snapShotSize; 

	snapshot_times.resize(nbTrueSnapShots);
	for(int k = 0; k < nbTrueSnapShots; k++){
		snapshot_times[k] = (k + 1) * snapShotSize + meta.steps;
	}

	//Initialize and fill the Eigen Wavefunction Storage
	// wavefct = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);

	// the time-step sizes for Runge-Kutta integration for both schemes as complex valued variables
	t_RTE = complex<double>(opt.RTE_step,0.0);

	// Maximum x and y ranges of the grid after expanding for the full runtime.
	// Needed to compute the growing plots.
	ranges.resize(2);
	complex<double> tmp3 = complex<double>(opt.RTE_step * opt.n_it_RTE,0.0);
	ranges[0] = opt.min_x * real(lambda_x(tmp3));
	ranges[1] = opt.min_y * real(lambda_y(tmp3));

	if(ranges[0] > ranges[1])
		ranges[1] = ranges[0];
	else if(ranges[1] > ranges[0])
		ranges[0] = ranges[1];

	// Grid Spacing variables
	h_x = complex<double>((2.*opt.min_x/opt.grid[1]),0.0);
  	h_y = complex<double>((2.*opt.min_y/opt.grid[2]),0.0); 

  	// Coordinate vectors/arrays in different forms etc.
  	x_axis.resize(opt.grid[1]);
  	y_axis.resize(opt.grid[2]);
  	for(int i=0;i<opt.grid[1];i++){x_axis[i]=-opt.min_x+i*real(h_x);}
  	for(int j=0;j<opt.grid[2];j++){y_axis[j]=-opt.min_y+j*real(h_y);}

  	X = VectorXcd(opt.grid[1]); Y = VectorXcd(opt.grid[2]);
	for(int i = 0;i<opt.grid[1];i++){X(i) = complex<double>(x_axis[i],0.0);}
	for(int j = 0;j<opt.grid[2];j++){Y(j) = complex<double>(y_axis[j],0.0);}

	Xmatrix = MatrixXcd(meta.grid[0]-2,meta.grid[1]-2); Ymatrix = MatrixXcd(meta.grid[0]-2,meta.grid[1]-2);
	for( int i = 0; i < meta.grid[0]-2; i++){ Xmatrix.col(i) = X.segment(1,meta.grid[0]-2);	}
	for( int i = 0; i < meta.grid[0]-2; i++){ Ymatrix.row(i) = Y.segment(1,meta.grid[0]-2);	}

	Xexpanding = VectorXd::Zero(opt.grid[1]);
	Yexpanding = VectorXd::Zero(opt.grid[2]);

	// The laplacian and gradient coefficient needed for the RTE scheme.
	// These are precomputed here, to simplify the computations later
	int coefSize = 2 * opt.n_it_RTE + 1 - 2 * meta.steps;
   	laplacian_coefficient_x = VectorXcd::Zero(coefSize);
   	laplacian_coefficient_y = VectorXcd::Zero(coefSize);
   	gradient_coefficient_x = VectorXcd::Zero(coefSize);
   	gradient_coefficient_y = VectorXcd::Zero(coefSize);

   	complex<double> tmp;  	
   	for(int t = 0; t < coefSize; t++){
   	tmp = complex<double>(meta.time,0.0) + ( half * complex<double>(t,0.0) * t_RTE );   	
   	laplacian_coefficient_x(t) = i_unit / ( two * h_x * h_x * lambda_x(tmp) * lambda_x(tmp) );
   	laplacian_coefficient_y(t) = i_unit / ( two * h_y * h_y * lambda_y(tmp) * lambda_y(tmp) );
   	gradient_coefficient_x(t) = lambda_x_dot(tmp) / (two * h_x * lambda_x(tmp));
   	gradient_coefficient_y(t) = lambda_y_dot(tmp) / (two * h_y * lambda_y(tmp));
   	}

   	PotentialGrid = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
   	for(int i = 0; i< opt.grid[1]; i++){for(int j = 0; j < opt.grid[2]; j++){
	PotentialGrid(i,j) = complex<double>(opt.potFactor,0.0) * ( half * opt.omega_x * opt.omega_x * X(i) * X(i) +  half * opt.omega_y * opt.omega_y * Y(j) * Y(j) );}}
	

   	pot_laplacian_x = complex<double>(1.0,0.0) / (two * h_x * h_x);
	pot_laplacian_y = complex<double>(1.0,0.0) / (two * h_y * h_y);


}


void SplitStep::splitToTime(string runName){

	cout << "Starting SSFT run!" << endl;


	double start;  // starttime of the run
	int samplesize = wavefctVec.size();
	keeperOfTime.absoluteSteps = 0;
	keeperOfTime.lambdaSteps = 0;
	keeperOfTime.initialSteps = meta.steps;

	double timestepsize = opt.RTE_step;
	ComplexGrid kprop(opt.grid[0], opt.grid[1], opt.grid[2],opt.grid[3]);
	ComplexGrid rgrid(opt.grid[0], opt.grid[1], opt.grid[2],opt.grid[3]);
	ComplexGrid kgrid(opt.grid[0], opt.grid[1], opt.grid[2],opt.grid[3]);
	ComplexGrid potPlotGrid(opt.grid[0], opt.grid[1], opt.grid[2],opt.grid[3]);


	vector<vector<double>> kspace;	
	kspace.resize(2);
	for(int d = 0; d < 2; d++){
		kspace[d].resize(opt.grid[d+1]);
		for(int i = 0; i < opt.grid[d+1]/2; i++){
			kspace[d][i] = opt.klength[d]*2.0*sin( M_PI*((double)i)/((double)opt.grid[d+1]) );
			// kspace[d][i] = 2.0 * M_PI * i / (opt.grid[d+1]);
		}
		for(int i = opt.grid[d+1]/2; i < opt.grid[d+1]; i++){
			kspace[d][i] = opt.klength[d]*2.0*sin( M_PI*((double)(-opt.grid[d+1]+i))/((double)opt.grid[d+1]) );
			// kspace[d][i] = 2.0 * M_PI * (i - opt.grid[d+1]) / (opt.grid[d+1]);
		}
	}
	double beta = real(h_x) * real(h_x);
	for(int x = 0; x < kprop.width(); x++){
	    for(int y = 0; y < kprop.height(); y++){
	      	double k[2];
	      	// k[0] = opt.klength[1] * 2.0 * sin(M_PI * x / (double) opt.grid[1]);
	      	// k[1] = opt.klength[2] * 2.0 * sin(M_PI * y / (double) opt.grid[2]);
	      	// double T = - 0.5 * beta * (k[0] * k[0] + k[1] * k[1] ) * timestepsize;
	
	      	double T = - 0.5 * (kspace[0][x]*kspace[0][x] + kspace[1][y]*kspace[1][y]) * timestepsize / beta;	      
	      	
	      	kprop(0,x,y,0) = complex<double>(cos(T),sin(T)) / complex<double>((double)(opt.grid[1]*opt.grid[2]*opt.grid[3]),0.0);	    
	    }
	}
	plotDataToPng("RTE_Kprop","Control",kprop,opt);

	if(opt.initialRun == true){
		Eval* initialEval = new Eval;
		initialEval->saveData(wavefctVec,opt,meta.steps,runName);
		initialEval->evaluateData();
		initialEval->plotData();
		// Commenting out both lines below, to switch on behavior in evaluation
		// This basically counts every Vortex in each step, instead of capping at the initial value
		// opt.vortexnumber = initialEval->getVortexNumber();
		// opt.initialRun = false;

		string evalname = runName + "-Eval.h5";
		binaryFile* evalFile = new binaryFile(evalname,binaryFile::out);
		evalFile->appendEval(meta.steps,opt,meta,*initialEval);
		delete initialEval;
		delete evalFile;
	}
	
	start = omp_get_wtime();
	omp_set_num_threads(12);
	int previousTimes = meta.steps;
	for(int j = 0; j < snapshot_times.size(); j++){
		// some information about the computation status and stuff
		string stepname = runName + "-" + to_string(snapshot_times[j]);
		vector<int> stateOfLoops(samplesize);
		vector<int> threadinfo(samplesize);
		int slowestthread = 0;

		// omp_set_num_threads(12);
		// #pragma omp parallel for
		for(int i = 0; i < samplesize; i++){

			#pragma omp parallel for
			for(int x = 0; x < opt.grid[1];x++){
				for(int y = 0; y < opt.grid[2]; y++){
					rgrid(0,x,y,0) = wavefctVec[i](x,y);
				}
			}

			// list of which thread is working which iteration
			int lambdaSteps = keeperOfTime.lambdaSteps;
			threadinfo[i] = omp_get_thread_num();
			for(int m = previousTimes + 1; m <= snapshot_times[j]; m++){
		
				ComplexGrid::fft_unnormalized(rgrid, kgrid, true);
    
				#pragma omp parallel for
				for(int x = 0; x < kgrid.width(); x++){
					for(int y = 0; y < kgrid.height(); y++){
				   		kgrid(0,x,y,0) = kprop(0,x,y,0) * kgrid(0,x,y,0);
					}
				}

				// plotDataToPng("RTE_KGrid_"+to_string(m),"RTE_KGrid_"+to_string(m),kgrid,opt);
				
				ComplexGrid::fft_unnormalized(kgrid, rgrid, false);

				// plotDataToPng("RTE_RGrid_"+to_string(m),"RTE_RGrid_"+to_string(m),rgrid,opt); 

				// ComplexGrid potGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

				

				#pragma omp parallel for
				for(int x = 0; x < rgrid.width(); x++){
					for(int y = 0; y < rgrid.height(); y++){
				    	complex<double> value = rgrid(0,x,y,0);
				    	double V = - ( /*PotentialGrid(x,y).real()*/ /*rotatingPotential(x,y,m)*//* +*/ opt.g * abs2(value) ) * timestepsize;
				    	// potPlotGrid(0,x,y,0) = complex<double>(rotatingPotential(x,y,m) /*PotentialGrid(x,y).real()*/,0.0);
				    	// potGrid(0,x,y,0) = complex<double>(cos(V),sin(V));
				    	rgrid(0,x,y,0) = complex<double>(cos(V),sin(V)) * value;
					}
				}
				
				// plotDataToPng("RTE_PotGrid_"+to_string(m),"RTE_PotGrid_"+to_string(m),potPlotGrid,opt);

				// ComplexGrid::fft_unnormalized(rgrid, kgrid, true);
    
				// #pragma omp parallel for
				// for(int x = 0; x < kgrid.width(); x++){
				// 	for(int y = 0; y < kgrid.height(); y++){
				//    		kgrid(0,x,y,0) = kprop(0,x,y,0) * kgrid(0,x,y,0);
				// 	}
				// }

				// ComplexGrid::fft_unnormalized(kgrid, rgrid, false);					
		
				// progress to the cli from the slowest thread to always have an update. (otherwise progressbar would freeze until next snapshot computation starts)
   				stateOfLoops[i]= m - previousTimes;
   				if(omp_get_thread_num() == slowestthread){
   					int counter_max = snapshot_times[j] - previousTimes;
   					cli(stepname,slowestthread,threadinfo,stateOfLoops,counter_max,start);
   				}
	
			}

			#pragma omp parallel for
			for(int x = 0; x < opt.grid[1];x++){
				for(int y = 0; y < opt.grid[2]; y++){
					wavefctVec[i](x,y) = rgrid(0,x,y,0);
				}
			}
	
		}
		keeperOfTime.lambdaSteps += 2 * (snapshot_times[j] - previousTimes);
		keeperOfTime.absoluteSteps = snapshot_times[j] - keeperOfTime.initialSteps;	
		previousTimes = snapshot_times[j];

		complex<double> tmp = complex<double>(snapshot_times[j] * opt.RTE_step,0.0);
		opt.t_abs = tmp;  

		opt.stateInformation.resize(2);
		if(opt.runmode.compare(1,1,"1") == 0){
			opt.stateInformation[0] = real(lambda_x(tmp)); // needed for expansion and the computing of the gradient etc.
			opt.stateInformation[1] = real(lambda_y(tmp));
			cout << opt.stateInformation[0] << "  " << opt.stateInformation[1] << endl;
		}
		if(opt.runmode.compare(1,1,"0") == 0){
			opt.stateInformation[0] = 1.0;
			opt.stateInformation[1] = 1.0;
		}

		vector<double> coord(2);
		coord[0] = opt.min_x * opt.stateInformation[0];
		coord[1] = opt.min_y * opt.stateInformation[1];
		pData->update(real(tmp),snapshot_times[j],coord);

		// plot("3-"+to_string(snapshot_times[j]));
		
		try{
			// plotDataToPng("RTE_RGrid"+to_string(snapshot_times[j]),"Control"+to_string(snapshot_times[j]),rgrid,opt);
			Eval results;
	
			cout << " >> Evaluating Datafiles "<< snapshot_times[j] << " ";
			results.saveData(pData->wavefunction,opt,snapshot_times[j],runName);
			results.evaluateData();
			results.plotData();

			string dataname = runName + "-LastGrid.h5";
			binaryFile* dataFile = new binaryFile(dataname,binaryFile::out);
			dataFile->appendSnapshot(runName,snapshot_times[j],pData,opt);
			delete dataFile;

			string evalname = runName + "-Eval.h5";
			binaryFile* evalFile = new binaryFile(evalname,binaryFile::append);
			// evalFile->appendEval(snapshot_times[j],opt,pData->getMeta(),vec1Name,vec1Rank,vec1);
			evalFile->appendEval(snapshot_times[j],opt,pData->getMeta(),results);
			delete evalFile;

			cout << " ..Snapshot saved to runData/ ";

		}
		catch(const std::exception& e) { 
			std::cerr 	<< "Unhandled Exception after dataFile.appendSnapshot() in rteToTime: " << std::endl; 
			throw e; 
		}

	}
}