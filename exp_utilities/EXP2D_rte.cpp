#define EIGEN_VECTORIZE
#define EIGEN_NO_DEBUG

#include <EXP2D_rte.hpp>
#include <omp.h>

using namespace std;
using namespace Eigen;

RTE::RTE()
{
  	// some constants used in computations to shorten stuff
	pi = M_PI;
 	zero=complex<double>(0,0);
 	half=complex<double>(0.5,0);
 	one=complex<double>(1,0);
 	two=complex<double>(2,0);
 	four=complex<double>(4,0);
 	six=complex<double>(6,0);
 	i_unit=complex<double>(0,1);

 	// Use this to control the program flow: first char determines if the program is loading from a dataset or using ITP to generate the necessary datafile
    // second char determines if expanding coordinates are used or not
    // third char determines if potential is switch on for the differential equation
	opt.runmode[0] = 0;
	opt.runmode[1] = 1;
	opt.runmode[2] = 0;
	opt.runmode[3] = 1;


}

RTE::RTE(ComplexGrid* &c,Options &externaloptions)
{	
	// Both essential Variables
	pPsi = c;
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
	omp_set_num_threads(omp_get_max_threads());
	cout << "Max Number of Threads: " << omp_get_max_threads() << endl;	
	// cout << "Eigenthreads: " << Eigen::nbThreads() << endl;

	// Using the setter function to initialize the stuff.
	RunSetup();

}

RTE::~RTE(){

}

void RTE::setOptions(Options &externaloptions){
	opt = externaloptions;

}

void RTE::RunSetup(){

	//Initialize and fill the Eigen Wavefunction Storage
	wavefct = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);

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

	Xexpanding = VectorXd::Zero(opt.grid[1]);
	Yexpanding = VectorXd::Zero(opt.grid[2]);

	// The laplacian and gradient coefficient needed for the RTE scheme.
	// These are precomputed here, to simplify the computations later

   	laplacian_coefficient_x = VectorXcd::Zero(2 * opt.n_it_RTE + 1);
   	laplacian_coefficient_y = VectorXcd::Zero(2 * opt.n_it_RTE + 1);
   	gradient_coefficient_x = VectorXcd::Zero(2 * opt.n_it_RTE + 1);
   	gradient_coefficient_y = VectorXcd::Zero(2 * opt.n_it_RTE + 1);

   	complex<double> tmp;  	
   	for(int t = 0; t < ( 2 * opt.n_it_RTE); t++){
   	tmp = half * complex<double>(t+1,0.0) * t_RTE;   	
   	laplacian_coefficient_x(t) = i_unit / ( two * h_x * h_x * lambda_x(tmp) * lambda_x(tmp) );
   	laplacian_coefficient_y(t) = i_unit / ( two * h_y * h_y * lambda_y(tmp) * lambda_y(tmp) );
   	gradient_coefficient_x(t) = lambda_x_dot(tmp) / (two * h_x * lambda_x(tmp));
   	gradient_coefficient_y(t) = lambda_y_dot(tmp) / (two * h_y * lambda_y(tmp));
   	}

   	PotentialGrid = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
   	for(int i = 0; i< opt.grid[1]; i++){for(int j = 0; j < opt.grid[2]; j++){
	PotentialGrid(i,j) = half * opt.omega_x * opt.omega_x * X(i) * X(i) +  half * opt.omega_y * opt.omega_y * Y(j) * Y(j);}}

   	pot_laplacian_x = complex<double>(1.0,0.0) / (two * h_x * h_x);
	pot_laplacian_y = complex<double>(1.0,0.0) / (two * h_y * h_y);


}

void RTE::cli(string name,int &slowestthread, vector<int> threadinfo, vector<int> stateOfLoops, int counter_max, double start)
{	
	for(int i = 0;i < stateOfLoops.size();i++){
		slowestthread = (stateOfLoops[slowestthread] <= stateOfLoops[i]) ? slowestthread : threadinfo[i];
	}
	
	if(fmod((float)stateOfLoops[slowestthread],(float)(counter_max/10))==0){
		int seconds;
		int min;
		int hour;
		int total;
		double totalstate = 0;
		double totalmaxpercent = (double)counter_max * (double)opt.samplesize / 100;
		for(int i = 0; i < opt.samplesize; i++){
			totalstate += stateOfLoops[i];
		}

		total = omp_get_wtime() - start;
		hour = total / 3600;
		min = (total / 60) % 60;
		seconds = total % 60;
		cout << "\r" << flush;
		cout << currentTime() <<  " " << name << " "
		 	 << std::setw(2) << std::setfill('0') << hour << ":"
			 << std::setw(2) << std::setfill('0') << min << ":"
			 << std::setw(2) << std::setfill('0') << seconds  << "    "
			 << std::setw(3) << std::setfill('0') << std::setprecision(2) << (totalstate/totalmaxpercent) << "% | "
			 << "Number of Threads: " << stateOfLoops.size() ;
			// cout << " | Slowest Thread: " << std::setw(3) << std::setfill('0') << (float)(stateOfLoops[slowestthread])/(float)(counter_max/10) << "% ";
			// for(int k = 0; k < stateOfLoops.size(); k++){
			// cout << k << "_" << threadinfo[k] << ": " << std::setw(3) << std::setfill('0') << (float)stateOfLoops[k]/((float)counter_max/100) << "% ";
		// }
		cout << "                          ";
	}
}

void RTE::plot(string name,int counter_state, int counter_max){
	wavefct = wavefctVec[0];

	if(opt.runmode.compare(1,1,"1") == 0){
		complex<double> tmp = complex<double>(keeperOfTime.absoluteSteps,0.0) * t_RTE;
		Xexpanding = x_expand(tmp);
		Yexpanding = y_expand(tmp);
		plotDataToPngEigenExpanding(name, wavefct,ranges,Xexpanding,Yexpanding,opt);
	}
	if(opt.runmode.compare(1,1,"0") == 0){
		plotDataToPngEigen(name, wavefct,opt);
	}
}


void RTE::CopyComplexGridToEigen(){
	for(int i = 0; i < opt.grid[1]; i++){for(int j = 0; j < opt.grid[2]; j++){ wavefct(i,j) = pPsi->at(0,i,j,0);}}
}

void RTE::CopyEigenToComplexGrid(){
	for(int i = 0; i < opt.grid[1]; i++){for(int j = 0; j < opt.grid[2]; j++){ pPsi->at(0,i,j,0) = wavefct(i,j);}}
}

void RTE::toEigenAndNoise(ComplexGrid g,MatrixXcd &wavefct){
	noiseTheGrid(g);
	for(int i = 0; i < opt.grid[1]; i++){for(int j = 0; j < opt.grid[2]; j++){ wavefct(i,j) = g(0,i,j,0);}}
}

void RTE::rteToTime(string runname, vector<int> snapshot_times)
{
	double start;  // starttime of the run
	// int t = 0;		// counter for the expanding lambdavectors with coefficients
	// int step_counter = 0;
	keeperOfTime.absoluteSteps = 0;
	keeperOfTime.lambdaSteps = 0;

	wavefctVec.resize(opt.samplesize);

	ComplexGrid Psi(*pPsi);

	// ComplexGrid Psi(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	// for(int i = 0; i < opt.grid[1];i++){for(int j = 0; j < opt.grid[2]; j++){Psi(0,i,j,0) = pPsi->at(0,i,j,0);}}

	vector<MatrixXcd> wavefctcp(opt.samplesize);	
	vector<MatrixXcd> k0(opt.samplesize);
	vector<MatrixXcd> k1(opt.samplesize);
	vector<MatrixXcd> k2(opt.samplesize);
	vector<MatrixXcd> k3(opt.samplesize);

	// CopyComplexGridToEigen();

	for(int i = 0; i < opt.samplesize;i++){
	wavefctVec[i] = MatrixXcd(opt.grid[1],opt.grid[2]);
	toEigenAndNoise(Psi,wavefctVec[i]);
	wavefctcp[i] = MatrixXcd(opt.grid[1],opt.grid[2]);
	k0[i] = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	k1[i] = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	k2[i] = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	k3[i] = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);	
	}


	binaryFile *dataFile = new binaryFile("00000.h5",binaryFile::out);
	dataFile->appendSnapshot(runname,0,wavefctVec,opt);
	delete dataFile;
	
	start = omp_get_wtime();


	//start loop here
	Eigen::initParallel();
	int previousTimes = 0;
	for(int j = 0; j < snapshot_times.size(); j++){
		// some information about the computation status and stuff
		string stepname = runname + "-" + to_string(snapshot_times[j]);
		vector<int> stateOfLoops(opt.samplesize);
		vector<int> threadinfo(opt.samplesize);
		int slowestthread = 0;

		#pragma omp parallel for
		for(int i = 0; i < opt.samplesize; i++){
			// list of which thread is working which iteration
			int lambdaSteps = keeperOfTime.lambdaSteps;
			threadinfo[i] = omp_get_thread_num();
			for(int m = previousTimes + 1; m <= snapshot_times[j]; m++){
		
				wavefctcp[i] = wavefctVec[i];
		
				// boundary conditions -- Dirichlet
		
				wavefctVec[i].row(0) = VectorXcd::Zero(opt.grid[1]);
				wavefctVec[i].row(opt.grid[1]-1) = VectorXcd::Zero(opt.grid[1]);
				wavefctVec[i].col(0) = VectorXcd::Zero(opt.grid[2]);
				wavefctVec[i].col(opt.grid[2]-1) = VectorXcd::Zero(opt.grid[2]);
		
				// boundary conditions end
		
				RTE_compute_k(k0[i],wavefctcp[i],lambdaSteps);
				wavefctcp[i] = wavefctVec[i] + half * t_RTE * k0[i];

				lambdaSteps++;
				RTE_compute_k(k1[i],wavefctcp[i],lambdaSteps);
				wavefctcp[i] = wavefctVec[i] + half * t_RTE * k1[i];
		
				RTE_compute_k(k2[i],wavefctcp[i],lambdaSteps);		
				wavefctcp[i] = wavefctVec[i] + t_RTE * k2[i];
		
				lambdaSteps++;
				RTE_compute_k(k3[i],wavefctcp[i],lambdaSteps);
		
				wavefctVec[i] += (t_RTE/six) * ( k0[i] + two * k1[i] + two * k2[i] + k3[i]);
		
				// // Neumann Boundaries
		
				// wavefctVec[i].col(0).real() = wavefctVec[i].col(1).real();
				// wavefctVec[i].col(opt.grid[2]-1).real() = wavefctVec[i].col(opt.grid[2]-2).real();
				// wavefctVec[i].row(0).real() = wavefctVec[i].row(0).real();
				// wavefctVec[i].row(opt.grid[1]-1).real() = wavefctVec[i].row(opt.grid[1]-2).real();
		
				// // Boundaries
		
				// progress to the cli from the slowest thread to always have an update. (otherwise progressbar would freeze until next snapshot computation starts)
   				stateOfLoops[i]= m - previousTimes;
   				if(omp_get_thread_num() == slowestthread){
   					int counter_max = snapshot_times[j] - previousTimes;
   					cli(stepname,slowestthread,threadinfo,stateOfLoops,counter_max,start);
   				}
	
			}
			// if(omp_get_thread_num() == 0){
				
			// }
	
		}
		keeperOfTime.lambdaSteps = 2 * snapshot_times[j];
		keeperOfTime.absoluteSteps = snapshot_times[j];	
		previousTimes = snapshot_times[j];

		complex<double> tmp = complex<double>(keeperOfTime.absoluteSteps * opt.RTE_step,0.0);  

		opt.stateInformation.resize(2);
		if(opt.runmode.compare(1,1,"1") == 0){
			opt.stateInformation[0] = real(lambda_x(tmp)); // needed for expansion and the computing of the gradient etc.
			opt.stateInformation[1] = real(lambda_y(tmp));
		}
		if(opt.runmode.compare(1,1,"0") == 0){
			opt.stateInformation[0] = 1.0;
			opt.stateInformation[1] = 1.0;
		}

		try{
			std::string h5name = to_string(snapshot_times[j]);
			std::stringstream ss;
			ss << std::setfill('0') << std::setw(5) << h5name;
			h5name = ss.str() + ".h5";

			binaryFile dataFile(h5name,binaryFile::out);
			dataFile.appendSnapshot(runname,snapshot_times[j],wavefctVec,opt);
			// dataFile.close();
			cout << endl << currentTime() << " Snapshot saved to runData.h5" << endl;

		}
		catch(const std::exception& e) { 
			std::cerr 	<< "Unhandled Exception after dataFile.appendSnapshot() in rteToTime: " 
					    << e.what() << ", application will now exit" << std::endl; 
			throw; 
		}

	}

// update the ComplexGrid* DATA object outside of this.
// if(opt.samplesize == 1){
// 	CopyEigenToComplexGrid();
// }
	

}

void RTE::rteFromDataToTime(string runname, vector<int> snapshot_times)
{
	double start;  // starttime of the run
	// int t = 0;		// counter for the expanding lambdavectors with coefficients
	// int step_counter = 0;
	keeperOfTime.absoluteSteps = 0;
	keeperOfTime.lambdaSteps = 0;

	binaryFile *dataLoading = new binaryFile("runData/00000.h5",binaryFile::in);
	dataLoading->getSnapshot(runname,0,wavefctVec,opt);
	delete dataLoading;

	// wavefctVec.resize(opt.samplesize);

	// ComplexGrid Psi(*pPsi);

	// ComplexGrid Psi(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	// for(int i = 0; i < opt.grid[1];i++){for(int j = 0; j < opt.grid[2]; j++){Psi(0,i,j,0) = pPsi->at(0,i,j,0);}}

	vector<MatrixXcd> wavefctcp(opt.samplesize);	
	vector<MatrixXcd> k0(opt.samplesize);
	vector<MatrixXcd> k1(opt.samplesize);
	vector<MatrixXcd> k2(opt.samplesize);
	vector<MatrixXcd> k3(opt.samplesize);

	// CopyComplexGridToEigen();

	for(int i = 0; i < opt.samplesize;i++){
	// wavefctVec[i] = MatrixXcd(opt.grid[1],opt.grid[2]);
	// toEigenAndNoise(Psi,wavefctVec[i]);
	wavefctcp[i] = MatrixXcd(opt.grid[1],opt.grid[2]);
	k0[i] = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	k1[i] = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	k2[i] = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	k3[i] = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);	
	}

	start = omp_get_wtime();


	//start loop here
	Eigen::initParallel();
	int previousTimes = 0;
	for(int j = 0; j < snapshot_times.size(); j++){
		// some information about the computation status and stuff
		string stepname = runname + "-" + to_string(snapshot_times[j]);
		vector<int> stateOfLoops(opt.samplesize);
		vector<int> threadinfo(opt.samplesize);
		int slowestthread = 0;

		#pragma omp parallel for
		for(int i = 0; i < opt.samplesize; i++){
			// list of which thread is working which iteration
			int lambdaSteps = keeperOfTime.lambdaSteps;
			threadinfo[i] = omp_get_thread_num();
			for(int m = previousTimes + 1; m <= snapshot_times[j]; m++){
		
				wavefctcp[i] = wavefctVec[i];
		
				// boundary conditions -- Dirichlet
		
				wavefctVec[i].row(0) = VectorXcd::Zero(opt.grid[1]);
				wavefctVec[i].row(opt.grid[1]-1) = VectorXcd::Zero(opt.grid[1]);
				wavefctVec[i].col(0) = VectorXcd::Zero(opt.grid[2]);
				wavefctVec[i].col(opt.grid[2]-1) = VectorXcd::Zero(opt.grid[2]);
		
				// boundary conditions end
		
				RTE_compute_k(k0[i],wavefctcp[i],lambdaSteps);
				wavefctcp[i] = wavefctVec[i] + half * t_RTE * k0[i];

				lambdaSteps++;
				RTE_compute_k(k1[i],wavefctcp[i],lambdaSteps);
				wavefctcp[i] = wavefctVec[i] + half * t_RTE * k1[i];
		
				RTE_compute_k(k2[i],wavefctcp[i],lambdaSteps);		
				wavefctcp[i] = wavefctVec[i] + t_RTE * k2[i];
		
				lambdaSteps++;
				RTE_compute_k(k3[i],wavefctcp[i],lambdaSteps);
		
				wavefctVec[i] += (t_RTE/six) * ( k0[i] + two * k1[i] + two * k2[i] + k3[i]);
		
				// // Neumann Boundaries
		
				// wavefctVec[i].col(0).real() = wavefctVec[i].col(1).real();
				// wavefctVec[i].col(opt.grid[2]-1).real() = wavefctVec[i].col(opt.grid[2]-2).real();
				// wavefctVec[i].row(0).real() = wavefctVec[i].row(0).real();
				// wavefctVec[i].row(opt.grid[1]-1).real() = wavefctVec[i].row(opt.grid[1]-2).real();
		
				// // Boundaries
		
				// progress to the cli from the slowest thread to always have an update. (otherwise progressbar would freeze until next snapshot computation starts)
   				stateOfLoops[i]= m - previousTimes;
   				if(omp_get_thread_num() == slowestthread){
   					int counter_max = snapshot_times[j] - previousTimes;
   					cli(stepname,slowestthread,threadinfo,stateOfLoops,counter_max,start);
   				}
	
			}
			// if(omp_get_thread_num() == 0){
				
			// }
	
		}
		keeperOfTime.lambdaSteps = 2 * snapshot_times[j];
		keeperOfTime.absoluteSteps = snapshot_times[j];	
		previousTimes = snapshot_times[j];

		complex<double> tmp = complex<double>(keeperOfTime.absoluteSteps * opt.RTE_step,0.0);  

		opt.stateInformation.resize(2);
		if(opt.runmode.compare(1,1,"1") == 0){
			opt.stateInformation[0] = real(lambda_x(tmp)); // needed for expansion and the computing of the gradient etc.
			opt.stateInformation[1] = real(lambda_y(tmp));
		}
		if(opt.runmode.compare(1,1,"0") == 0){
			opt.stateInformation[0] = 1.0;
			opt.stateInformation[1] = 1.0;
		}

		try{
			std::string h5name = to_string(snapshot_times[j]);
			std::stringstream ss;
			ss << std::setfill('0') << std::setw(5) << h5name;
			h5name = ss.str() + ".h5";
			
			binaryFile dataFile(h5name,binaryFile::out);
			dataFile.appendSnapshot(runname,snapshot_times[j],wavefctVec,opt);
			// dataFile.close();
			cout << endl << currentTime() << " Snapshot saved to runData.h5" << endl;

		}
		catch(const std::exception& e) { 
			std::cerr 	<< "Unhandled Exception after dataFile.appendSnapshot() in rteToTime: " 
					    << e.what() << ", application will now exit" << std::endl; 
			throw; 
		}

	}

// update the ComplexGrid* DATA object outside of this.
// if(opt.samplesize == 1){
// 	CopyEigenToComplexGrid();
// }
	

}

inline void RTE::RTE_compute_k(MatrixXcd &k,MatrixXcd &wavefctcp,int &t)
	{
	Matrix<std::complex<double>,Dynamic,Dynamic,ColMajor> wavefctcpX = Matrix<std::complex<double>,Dynamic,Dynamic,ColMajor>::Zero(opt.grid[1],opt.grid[2]);
	Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor> wavefctcpY = Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor>::Zero(opt.grid[1],opt.grid[2]);
	k = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);                                                                                                                                            

	if(opt.runmode.compare(1,1,"1") == 0){
	//laplacian
	for(int j = 1;j<opt.grid[2]-1;j++){
	for(int i = 1;i<opt.grid[1]-1;i++){
	wavefctcpX(i,j) = wavefctcp(i-1,j) - two * wavefctcp(i,j) + wavefctcp(i+1,j);
	wavefctcpY(i,j) = wavefctcp(i,j-1) - two * wavefctcp(i,j) + wavefctcp(i,j+1);
	}}
	k.noalias() +=   wavefctcpX * laplacian_coefficient_x(t) + wavefctcpY * laplacian_coefficient_y(t);

	// gradient
	for(int j = 1;j<opt.grid[2]-1;j++){
	for(int i = 1;i<opt.grid[1]-1;i++){
	wavefctcpX(i,j) = wavefctcp(i+1,j) - wavefctcp(i-1,j);
	wavefctcpY(i,j) = wavefctcp(i,j+1) - wavefctcp(i,j-1);
	}}

	for(int i = 0;i<opt.grid[1];i++){ wavefctcpY.row(i).array() *= Y.array(); }
	for(int j = 0;j<opt.grid[2];j++){ wavefctcpX.col(j).array() *= X.array(); }
	k.noalias() += wavefctcpX * gradient_coefficient_x(t) + wavefctcpY * gradient_coefficient_y(t);
	}

	if(opt.runmode.compare(1,1,"0") == 0){
	//laplacian
	for(int j = 1;j<opt.grid[2]-1;j++){
	for(int i = 1;i<opt.grid[1]-1;i++){
	wavefctcpX(i,j) = wavefctcp(i-1,j) - two * wavefctcp(i,j) + wavefctcp(i+1,j);
	wavefctcpY(i,j) = wavefctcp(i,j-1) - two * wavefctcp(i,j) + wavefctcp(i,j+1);
	}}
	k.noalias() +=   wavefctcpX * pot_laplacian_x * i_unit + wavefctcpY * pot_laplacian_x * i_unit;
	}

	if(opt.runmode.compare(2,1,"0") == 0){
	//interaction
	k.array() -= complex<double>(0.0,opt.g) * ( wavefctcp.conjugate().array() * wavefctcp.array() ) * wavefctcp.array();
	}

	if(opt.runmode.compare(2,1,"1") == 0){
		//Potential + Interaction
	k.array() -= ( i_unit * PotentialGrid.array() + complex<double>(0.0,opt.g) * ( wavefctcp.conjugate().array() * wavefctcp.array() )) * wavefctcp.array();
	}

	}

// inline void RTE::RTE_compute_k_pot(MatrixXcd &k,MatrixXcd &wavefctcp,int &t)
// 	{
// 	Matrix<std::complex<double>,Dynamic,Dynamic,ColMajor> wavefctcpX = Matrix<std::complex<double>,Dynamic,Dynamic,ColMajor>::Zero(opt.grid[1],opt.grid[2]);
// 	Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor> wavefctcpY = Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor>::Zero(opt.grid[1],opt.grid[2]);
// 	k = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);


// 	//laplacian
// 	for(int j = 1;j<opt.grid[2]-1;j++){
// 	for(int i = 1;i<opt.grid[1]-1;i++){
// 	wavefctcpX(i,j) = wavefctcp(i-1,j) - two * wavefctcp(i,j) + wavefctcp(i+1,j);
// 	wavefctcpY(i,j) = wavefctcp(i,j-1) - two * wavefctcp(i,j) + wavefctcp(i,j+1);
// 	}}
// 	k.noalias() +=   wavefctcpX * pot_laplacian_x * i_unit + wavefctcpY * pot_laplacian_x * i_unit;

// 	//Potential + Interaction
// 	k.array() -= ( i_unit * PotentialGrid.array() + complex<double>(0.0,opt.g) * ( wavefctcp.conjugate().array() * wavefctcp.array() )) * wavefctcp.array();

// 	}

