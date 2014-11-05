
#include <EXP2D_splitstep.hpp>
#include <omp.h>

#define EIGEN_DONT_VECTORIZE
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_NO_DEBUG

// #define SLICE_NUMBER 0

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

		plotDataToPng(name,"Control", wavefctVec[0],opt);

}

SplitStep::SplitStep(vector<ComplexGrid> &d,const MatrixData::MetaData &extMeta, const Options &externaloptions, int &extSLICE_NUMBER)
{	
	setOptions(externaloptions);
	SLICE_NUMBER = extSLICE_NUMBER;
	meta = extMeta;
	// Both essential Variables
	wavefctVec.resize(opt.samplesize);
	for(int i = 0; i < opt.samplesize; i++){
		wavefctVec[i] = ComplexGrid(opt.grid[0], opt.grid[1], opt.grid[2], opt.grid[3]);
		wavefctVec[i] = d[i];
	}
	
  	

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

	// Grid Spacing variables
	h_x = complex<double>((2.*opt.min_x/opt.grid[1]),0.0);
  	h_y = complex<double>((2.*opt.min_y/opt.grid[2]),0.0);
  	h_z = complex<double>((2.*opt.min_z/opt.grid[3]),0.0);

  	// Coordinate vectors/arrays in different forms etc.
  	x_axis.resize(opt.grid[1]);
  	y_axis.resize(opt.grid[2]);
  	z_axis.resize(opt.grid[3]);
  	for(int i=0;i<opt.grid[1];i++){x_axis[i]=-opt.min_x+i*real(h_x);}
  	for(int j=0;j<opt.grid[2];j++){y_axis[j]=-opt.min_y+j*real(h_y);}
  	for(int k=0;k<opt.grid[3];k++){z_axis[k]=-opt.min_z+k*real(h_z);}

  	X = VectorXcd(opt.grid[1]); Y = VectorXcd(opt.grid[2]); Z = VectorXcd(opt.grid[3]);
	for(int i = 0;i<opt.grid[1];i++){X(i) = complex<double>(x_axis[i],0.0);}
	for(int j = 0;j<opt.grid[2];j++){Y(j) = complex<double>(y_axis[j],0.0);}
	for(int k = 0;k<opt.grid[3];k++){Z(k) = complex<double>(z_axis[k],0.0);}


   	PotentialGrid = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
   	for(int i = 0; i< opt.grid[1]; i++){for(int j = 0; j < opt.grid[2]; j++){for(int k = 0; k < opt.grid[3]; k++){
	PotentialGrid(0,i,j,k) = complex<double>(opt.potFactor,0.0) * ( half * opt.omega_x * opt.omega_x * X(i) * X(i) +  half * opt.omega_y * opt.omega_y * Y(j) * Y(j) + half * opt.omega_z * opt.omega_z * Z(k) * Z(k) );}}}
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
	// ComplexGrid potPlotGrid(opt.grid[0], opt.grid[1], opt.grid[2],opt.grid[3]);




	// vector<vector<double>> kspace;	
	// kspace.resize(3);
	// for(int d = 0; d < 3; d++){
	// 	kspace[d].resize(opt.grid[d+1]);
	// 	for(int i = 0; i < opt.grid[d+1]/2; i++){
	// 		kspace[d][i] = opt.klength[d]*2.0*sin( M_PI*((double)i)/((double)opt.grid[d+1]) );
	// 		// kspace[d][i] = 2.0 * M_PI * i / (opt.grid[d+1]);
	// 	}
	// 	for(int i = opt.grid[d+1]/2; i < opt.grid[d+1]; i++){
	// 		kspace[d][i] = opt.klength[d]*2.0*sin( M_PI*((double)(-opt.grid[d+1]+i))/((double)opt.grid[d+1]) );
	// 		// kspace[d][i] = 2.0 * M_PI * (i - opt.grid[d+1]) / (opt.grid[d+1]);
	// 	}
	// }
	double beta = real(h_x) * real(h_x); // FIXME for all directions, not important if same size in all directions.

	#pragma omp parallel for
	for(int x = 0; x < opt.grid[1]; x++){
	    for(int y = 0; y < opt.grid[2]; y++){
	    	for(int z = 0; z < opt.grid[3]; z++){
	      		double k[3];
	      		k[0] = opt.klength[0] * 2.0 * sin(M_PI * x / (double) opt.grid[1]);
	      		k[1] = opt.klength[1] * 2.0 * sin(M_PI * y / (double) opt.grid[2]);
	      		k[2] = opt.klength[2] * 2.0 * sin(M_PI * z / (double) opt.grid[3]);
	      		double T = - 0.5 * (k[0] * k[0] + k[1] * k[1] + k[2] * k[2] ) * timestepsize; // / beta;
		
		      		// double T = - 0.5 * (kspace[0][x]*kspace[0][x] + kspace[1][y]*kspace[1][y] + kspace[2][z]*kspace[2][z]) * timestepsize / beta;	      
		      		
	      		kprop(0,x,y,z) = complex<double>(cos(T),sin(T)) / complex<double>((double)(opt.grid[1]*opt.grid[2]*opt.grid[3]),0.0);	    
	      	}
	    }
	}
	plotDataToPng("RTE_Kprop","Control",kprop,opt);

	if(opt.initialRun == true){
		Eval* initialEval = new Eval;
		initialEval->saveData2DSlice(wavefctVec,opt,meta.steps,runName,SLICE_NUMBER);
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
					for(int z = 0; z < opt.grid[3]; z++){
						rgrid(0,x,y,z) = wavefctVec[i](0,x,y,z);
					}
				}
			}

			// list of which thread is working which iteration
			int lambdaSteps = keeperOfTime.lambdaSteps;
			threadinfo[i] = omp_get_thread_num();
			for(int m = previousTimes + 1; m <= snapshot_times[j]; m++){
		
				ComplexGrid::fft_unnormalized(rgrid, kgrid, true);
    
				#pragma omp parallel for
				for(int x = 0; x < opt.grid[1]; x++){
					for(int y = 0; y < opt.grid[2]; y++){
						for(int z = 0; z < opt.grid[3]; z++){				
				   			kgrid(0,x,y,z) = kprop(0,x,y,z) * kgrid(0,x,y,z);
				   		}
					}
				}

				// plotDataToPng("RTE_KGrid_"+to_string(m),"RTE_KGrid_"+to_string(m),kgrid,opt);
				
				ComplexGrid::fft_unnormalized(kgrid, rgrid, false);

				// plotDataToPng("RTE_RGrid_"+to_string(m),"RTE_RGrid_"+to_string(m),rgrid,opt); 

				// ComplexGrid potGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

				

				#pragma omp parallel for
				for(int x = 0; x < opt.grid[1]; x++){
					for(int y = 0; y < opt.grid[2]; y++){
						for(int z = 0; z < opt.grid[3]; z++){	
				    		complex<double> value = rgrid(0,x,y,z);
				    		double V = - ( /*PotentialGrid(x,y).real()*/ /*rotatingPotential(x,y,m)*//* +*/ opt.g * abs2(value) ) * timestepsize;
				    		// potPlotGrid(0,x,y,0) = complex<double>(rotatingPotential(x,y,m) /*PotentialGrid(x,y).real()*/,0.0);
				    		// potGrid(0,x,y,0) = complex<double>(cos(V),sin(V));
				    		rgrid(0,x,y,z) = complex<double>(cos(V),sin(V)) * value;
				    	}
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
					for(int z = 0; z < opt.grid[3]; z++){
						wavefctVec[i](0,x,y,z) = rgrid(0,x,y,z);
					}
				}
			}
	
		}
		keeperOfTime.lambdaSteps += 2 * (snapshot_times[j] - previousTimes);
		keeperOfTime.absoluteSteps = snapshot_times[j] - keeperOfTime.initialSteps;	
		previousTimes = snapshot_times[j];

		complex<double> tmp = complex<double>(snapshot_times[j] * opt.RTE_step,0.0);
		opt.t_abs = tmp;  

		opt.stateInformation.resize(3);
		if(opt.runmode.compare(1,1,"1") == 0){
			cout << "Wrong runmode!! choose nonexpanding";
		}
		if(opt.runmode.compare(1,1,"0") == 0){
			opt.stateInformation[0] = 1.0;
			opt.stateInformation[1] = 1.0;
			opt.stateInformation[2] = 1.0;
		}

		// plot("3-"+to_string(snapshot_times[j]));
		
		try{
			// plotDataToPng("RTE_RGrid"+to_string(snapshot_times[j]),"Control"+to_string(snapshot_times[j]),rgrid,opt);
			Eval results;
			cout << " >> Evaluating Datafiles "<< snapshot_times[j] << " ";
			results.saveData2DSlice(wavefctVec,opt,snapshot_times[j],runName,SLICE_NUMBER);
			results.evaluateData();
			results.plotData();

			string dataname = runName + "-LastGrid.h5";
			binaryFile* dataFile = new binaryFile(dataname,binaryFile::out);
			dataFile->appendSnapshot(runName,snapshot_times[j],wavefctVec,meta,opt);
			delete dataFile;

			string evalname = runName + "-Eval.h5";
			binaryFile* evalFile = new binaryFile(evalname,binaryFile::append);
			// evalFile->appendEval(snapshot_times[j],opt,pData->getMeta(),vec1Name,vec1Rank,vec1);
			evalFile->appendEval(snapshot_times[j],opt,meta,results);
			delete evalFile;

			cout << " ..Snapshot saved to runData/ ";

		}
		catch(const std::exception& e) { 
			std::cerr 	<< "Unhandled Exception after dataFile.appendSnapshot() in rteToTime: " << std::endl; 
			throw e; 
		}

	}
}