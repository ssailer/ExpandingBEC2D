
// 

#include <EXP2D_rte.hpp>
#include <omp.h>

#define EIGEN_VECTORIZE
#define EIGEN_PARALLELIZE
#define EIGEN_NO_DEBUG

using namespace std;
using namespace Eigen;

RTE::RTE(MatrixData* &d,const Options &externaloptions) : wavefctVec(d->wavefunction), meta(d->meta), pData(d)
{	
	opt = externaloptions;
	RunSetup();
}

void RTE::setOptions(const Options &externaloptions)
{
	opt = externaloptions;
}

void RTE::RunSetup()
{
	int snapShotSize = opt.n_it_RTE / opt.snapshots;
	int nbTrueSnapShots =	opt.snapshots - meta.steps / snapShotSize; 

	snapshot_times.resize(nbTrueSnapShots);
	for(int k = 0; k < nbTrueSnapShots; k++){
		snapshot_times[k] = (k + 1) * snapShotSize + meta.steps;
	}

 //   	ranges.resize(2);
	// // complex<double> tmp3 = complex<double>(opt.RTE_step * opt.n_it_RTE,0.0);
	// ranges[0] = opt.min_x * real(lambda_x(absTime));
	// ranges[1] = opt.min_y * real(lambda_y(absTime));
	// double arr = real(absTime) * 1000.0 / opt.OmegaG;
	// cout << "Max ExpFactor: " << real(lambda_x(absTime)) << "  " << real(lambda_y(absTime)) << " with a time of " << std::setprecision(15) << arr << " ms" << " " << t_RTE(0) / opt.OmegaG << " " <<  t_RTE(l-1) / opt.OmegaG<< endl;

	// if(ranges[0] > ranges[1])
	// 	ranges[1] = ranges[0];
	// else if(ranges[1] > ranges[0])
	// 	ranges[0] = ranges[1];
	// cout << "Max Ranges: " << ranges[0] << "  " << ranges[1] << endl;


 //   	double x0 = meta.grid[0] * 0.20 / 2;
 //   	double y0 = meta.grid[1] * 0.20 / 2;
 //   	double rx = opt.min_x - x0;
	// double ry = opt.min_y - y0;
	// double strength = 5;


	double mu = sqrt(3.0  * opt.g * real(opt.omega_x) * real(opt.omega_y) * opt.N / 8.0);
    double Ry = sqrt(2.0 * mu / ( real(opt.omega_y)*real(opt.omega_y))) * opt.Ag;
    double Rx = sqrt(2.0 * mu / ( real(opt.omega_x)*real(opt.omega_x))) * opt.Ag;

    cout << "Initial Thomas Fermi Radii set to Rx = " << Rx << " and Ry = " << Ry << endl;
    double n0 = 2 * (opt.N / M_PI) * (1 / (Rx * Ry));
    cout << "n_0 = " << n0 << endl;


}

void RTE::cli(string name, int index, double start)
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

void RTE::plot(const string name){

	if(opt.runmode.compare(1,1,"1") == 0){
		// complex<double> tmp = complex<double>(keeperOfTime.absoluteSteps,0.0) * t_RTE;
		Xexpanding = x_expand(opt.t_abs);
		Yexpanding = y_expand(opt.t_abs);
		plotDataToPngEigenExpanding(name, wavefctVec[0],ranges,Xexpanding,Yexpanding,opt);
	}
	if(opt.runmode.compare(1,1,"0") == 0){
		plotDataToPngEigen(name, wavefctVec[0],opt);
	}
}
void RTE::noise(){
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

void RTE::rteToTime(string runName)
{	
	omp_set_num_threads(16);

	double start;  // starttime of the run

	RK4 *rk4 = new RK4(pData, opt);


	if(opt.initialRun == true){
		Eval* initialEval = new Eval;
		initialEval->saveData(wavefctVec,opt,meta.steps,runName);
		initialEval->evaluateData();
		initialEval->plotData();
		opt.vortexnumber = initialEval->getVortexNumber();
		opt.initialRun = false;

		string evalname = runName + "-Eval.h5";
		binaryFile* evalFile = new binaryFile(evalname,binaryFile::out);
		evalFile->appendEval(meta.steps,opt,meta,*initialEval);
		delete initialEval;
		delete evalFile;
	}
	
	start = omp_get_wtime();

	for(int j = 0; j < snapshot_times.size(); j++){

		while(pData->meta.steps < snapshot_times[j]){

			cli(runName,j,start);
				
			rk4->timeStep(opt.RTE_step);

			opt.t_abs += opt.RTE_step;


		}			

		// REMOVE
			if(opt.runmode.compare(1,1,"1") == 0){
				opt.stateInformation[0] = real(lambda_x(opt.t_abs)); // needed for expansion and the computing of the gradient etc.
				opt.stateInformation[1] = real(lambda_y(opt.t_abs));
				cout << opt.stateInformation[0] << "  " << opt.t_abs << " ";
				cout << pData->meta.initCoord[0] << " " << pData->meta.coord[0] << " " << pData->meta.time;
			}
			if(opt.runmode.compare(1,1,"0") == 0){
				opt.stateInformation[0] = 1.0;
				opt.stateInformation[1] = 1.0;
			}
		// END REMOVE


		try{

			

			// One Result file
			// results
			// Eval results;
			// results.saveData(pData->wavefunction,opt,snapshot_times[j],runName);
			// results.evaluateData();
			// results.plotData();

			// // one save file
			// // save snapshot
			// // save eval
			// string dataname = runName + "-LastGrid.h5";
			// binaryFile* dataFile = new binaryFile(dataname,binaryFile::out);
			// dataFile->appendSnapshot(runName,snapshot_times[j],pData,opt);
			// delete dataFile;

			// string evalname = runName + "-Eval.h5";
			// binaryFile* evalFile = new binaryFile(evalname,binaryFile::append);
			// evalFile->appendEval(snapshot_times[j],opt,pData->getMeta(),results);
			// delete evalFile;

			// BinaryFile* bin = new BinaryFile(BinaryFile::append);
			// Evaluation* eval = new Evaluation(*pData,opt);
			// eval->process();
			// bin->appendSnapshot(*pData,*eval);
			// // bin->appendWavefunction(pData);
			// // bin->appendEvaluation(eval);
			// delete bin, eval;


		}
		catch(const std::exception& e) { 
			std::cerr 	<< "Unhandled Exception after dataFile.appendSnapshot() in rteToTime: " << std::endl; 
			throw e; 
		}
	}
	delete rk4;
}







// lower order O(h^2) for reference.
// void RTE::singleK(MatrixXcd &k, MatrixXcd &wavefctcp, int32_t &front, int32_t &end,int32_t &subx,int32_t & suby, int &t){
// 	k.block(front,1,end,suby).noalias() =  (wavefctcp.block(front-1,1,end,suby)	- two * wavefctcp.block(front  ,1,end,suby) + wavefctcp.block(front+1,1,end,suby)) * laplacian_coefficient_x(t)
// 										 + (wavefctcp.block(front  ,0,end,suby) - two * wavefctcp.block(front  ,1,end,suby) + wavefctcp.block(front  ,2,end,suby)) * laplacian_coefficient_y(t);

// 	k.block(front,1,end,suby).array() += (wavefctcp.block(front+1,1,end,suby).array() - wavefctcp.block(front-1,1,end,suby).array()) * Xmatrix.block(front,1,end,suby).array() * gradient_coefficient_x(t)
// 									   + (wavefctcp.block(front  ,2,end,suby).array() - wavefctcp.block(front  ,0,end,suby).array()) * Ymatrix.block(front,1,end,suby).array() * gradient_coefficient_y(t);

// 	k.block(front,1,end,suby).array() -= i_unit * ( PotentialGrid.block(front,1,end,suby).array() + (complex<double>(opt.g,0.0) * ( wavefctcp.block(front,1,end,suby).conjugate().array() * wavefctcp.block(front,1,end,suby).array() ))) * wavefctcp.block(front,1,end,suby).array();
// }

// void Expansion::singleK(MatrixXcd &k, MatrixXcd &wavefctcp, int32_t &front, int32_t &end,int32_t &subx,int32_t & suby, int &t){
// 	k.block(front,2,end,suby).noalias() =  (-wavefctcp.block(front-2,2,end,suby) + sixteen * wavefctcp.block(front-1,2,end,suby) - thirty * wavefctcp.block(front  ,2,end,suby) + sixteen * wavefctcp.block(front+1,2,end,suby) - wavefctcp.block(front+2,2,end,suby)) * laplacian_coefficient_x(t)
// 										 + (-wavefctcp.block(front  ,0,end,suby) + sixteen * wavefctcp.block(front  ,1,end,suby) - thirty * wavefctcp.block(front  ,2,end,suby) + sixteen * wavefctcp.block(front  ,3,end,suby) - wavefctcp.block(front  ,4,end,suby)) * laplacian_coefficient_y(t);

// 	k.block(front,2,end,suby).array() += (-wavefctcp.block(front+2,2,end,suby).array() + eight * wavefctcp.block(front+1,2,end,suby).array() - eight * wavefctcp.block(front-1,2,end,suby).array() + wavefctcp.block(front-2,2,end,suby).array()) * Xmatrix.block(front,2,end,suby).array() * gradient_coefficient_x(t)
// 									   + (-wavefctcp.block(front  ,4,end,suby).array() + eight * wavefctcp.block(front  ,3,end,suby).array() - eight * wavefctcp.block(front  ,1,end,suby).array() + wavefctcp.block(front  ,0,end,suby).array()) * Ymatrix.block(front,2,end,suby).array() * gradient_coefficient_y(t);

// 	k.block(front,2,end,suby).array() -= i_unit * ( /*PotentialGrid.block(front,2,end,suby).array() +*/ (complex<double>(opt.g,0.0) * ( wavefctcp.block(front,2,end,suby).conjugate().array() * wavefctcp.block(front,2,end,suby).array() ))) * wavefctcp.block(front,2,end,suby).array();
// }

// void Trap::singleK(MatrixXcd &k, MatrixXcd &wavefctcp, int32_t &front, int32_t &end,int32_t &subx,int32_t & suby, int &t){
// 	k.block(front,2,end,suby).noalias() =  (-wavefctcp.block(front-2,2,end,suby) + sixteen * wavefctcp.block(front-1,2,end,suby) - thirty * wavefctcp.block(front  ,2,end,suby) + sixteen * wavefctcp.block(front+1,2,end,suby) - wavefctcp.block(front+2,2,end,suby)) * laplacian_coefficient_x(t)
// 										 + (-wavefctcp.block(front  ,0,end,suby) + sixteen * wavefctcp.block(front  ,1,end,suby) - thirty * wavefctcp.block(front  ,2,end,suby) + sixteen * wavefctcp.block(front  ,3,end,suby) - wavefctcp.block(front  ,4,end,suby)) * laplacian_coefficient_y(t);

// 	k.block(front,2,end,suby).array() += (-wavefctcp.block(front+2,2,end,suby).array() + eight * wavefctcp.block(front+1,2,end,suby).array() - eight * wavefctcp.block(front-1,2,end,suby).array() + wavefctcp.block(front-2,2,end,suby).array()) * Xmatrix.block(front,2,end,suby).array() * gradient_coefficient_x(t)
// 									   + (-wavefctcp.block(front  ,4,end,suby).array() + eight * wavefctcp.block(front  ,3,end,suby).array() - eight * wavefctcp.block(front  ,1,end,suby).array() + wavefctcp.block(front  ,0,end,suby).array()) * Ymatrix.block(front,2,end,suby).array() * gradient_coefficient_y(t);

// 	k.block(front,2,end,suby).array() -= i_unit * ( PotentialGrid.block(front,2,end,suby).array() + (complex<double>(opt.g,0.0) * ( wavefctcp.block(front,2,end,suby).conjugate().array() * wavefctcp.block(front,2,end,suby).array() ))) * wavefctcp.block(front,2,end,suby).array();
// }



void RTE::splitToTime(string runName){

	cout << "Starting SSFT run!" << endl;


	double start;  // starttime of the run
	int samplesize = wavefctVec.size();
	keeperOfTime.absoluteSteps = 0;
	keeperOfTime.lambdaSteps = 0;
	keeperOfTime.initialSteps = meta.steps;

	double timestepsize = opt.RTE_step;
	ComplexGrid kprop(opt.grid[0], meta.grid[0], meta.grid[1],opt.grid[3]);
	ComplexGrid rgrid(opt.grid[0], meta.grid[0], meta.grid[1],opt.grid[3]);
	ComplexGrid kgrid(opt.grid[0], meta.grid[0], meta.grid[1],opt.grid[3]);
	ComplexGrid potPlotGrid(opt.grid[0], meta.grid[0], meta.grid[1],opt.grid[3]);


	vector<vector<double>> kspace(2);
	for(int d = 0; d < 2; d++){
		kspace[d].resize(opt.grid[d+1]);
		for(int i = 0; i <= opt.grid[d+1]/2; i++){
			kspace[d][i] = (M_PI / opt.min_x) * (double)i;
		}

		for(int i = (opt.grid[d+1]/2)+1; i < opt.grid[d+1]; i++){
			kspace[d][i] = -(M_PI / opt.min_x) * (double)(opt.grid[d+1] - i);
		}
	}

	#pragma omp parallel for
	for(int x = 0; x < meta.grid[0]; x++){
	    for(int y = 0; y < meta.grid[1]; y++){		
		    double T = - 0.5 * (kspace[0][x]*kspace[0][x] + kspace[1][y]*kspace[1][y]) * timestepsize; // / beta;	      
		      		
	      	kprop(0,x,y,0) = complex<double>(cos(T),sin(T)) / complex<double>((double)(meta.grid[0]*meta.grid[1]),0.0);	    	      	
	    }
	}
	// plotDataToPng("RTE_Kprop","Control",kprop,opt);

	if(opt.initialRun == true){
		Eval* initialEval = new Eval;
		initialEval->saveData(wavefctVec,opt,meta.steps,runName);
		initialEval->evaluateData();
		initialEval->plotData();
		// Commenting out both lines below, to switch on behavior in evaluation
		// This basically counts every Vortex in each step, instead of capping at the initial value
		opt.vortexnumber = initialEval->getVortexNumber();
		opt.initialRun = false;

		string evalname = runName + "-Eval.h5";
		binaryFile* evalFile = new binaryFile(evalname,binaryFile::out);
		evalFile->appendEval(meta.steps,opt,meta,*initialEval);
		delete initialEval;
		delete evalFile;
	}
	
	start = omp_get_wtime();
	omp_set_num_threads(16);
	int previousTimes = meta.steps;
	for(int j = 0; j < snapshot_times.size(); j++){
		// some information about the computation status and stuff
		string stepname = runName + "-" + to_string(snapshot_times[j]);
		vector<int> stateOfLoops(samplesize);
		vector<int> threadinfo(samplesize);
		int slowestthread = 0;

		for(int i = 0; i < samplesize; i++){

			#pragma omp parallel for
			for(int x = 0; x < meta.grid[0];x++){
				for(int y = 0; y < meta.grid[1]; y++){
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

				// ComplexGrid potGrid(opt.grid[0],meta.grid[0],meta.grid[1],opt.grid[3]);

				

				#pragma omp parallel for
				for(int x = 0; x < rgrid.width(); x++){
					for(int y = 0; y < rgrid.height(); y++){
				    	complex<double> value = rgrid(0,x,y,0);
				    	double V = - ( PotentialGrid(x,y).real() /*rotatingPotential(x,y,m)*/ + opt.g * abs2(value) ) * timestepsize;
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
   					// cli(stepname,slowestthread,threadinfo,stateOfLoops,counter_max,start);
   				}
	
			}

			#pragma omp parallel for
			for(int x = 0; x < meta.grid[0];x++){
				for(int y = 0; y < meta.grid[1]; y++){
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

inline double RTE::rotatingPotential(int &i, int &j, int &t){
	double potential;
	if(t <= 3000){
		double alpha = 2 * M_PI / 500;
		double x = X(i).real() * cos(alpha * t) + Y(j).real() * sin(alpha * t);
		double y = - X(i).real() * sin(alpha * t) + Y(j).real() * cos(alpha * t);
		potential = 0.5 * opt.omega_x.real() * opt.omega_x.real() * x * x + 0.5 * opt.omega_y.real() * opt.omega_y.real() * y * y;
	} else {
		potential = 0.5 * opt.omega_x.real() * opt.omega_x.real() * X(i).real() * X(i).real() + 0.5 * opt.omega_x.real() * opt.omega_x.real() * Y(j).real() * Y(j).real();
	}
	return potential;
}



