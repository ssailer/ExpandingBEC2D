#define EIGEN_VECTORIZE
#define EIGEN_NO_DEBUG

#include <EXP2D_itp.hpp>
#include <omp.h>

#define VORTICES_BUILD_TIME 1000

using namespace std;
using namespace Eigen;

ITP::ITP()
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
}

ITP::ITP(ComplexGrid* &c,Options &externaloptions)
{	
	// Both essential Variables
	pPsi = c;
  	opt = externaloptions;


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
	cout << "Eigenthreads: " << Eigen::nbThreads() << endl;

	// Using the setter function to initialize the stuff.
	RunSetup();

}

ITP::~ITP(){

}

void ITP::setOptions(Options &externaloptions){
	opt = externaloptions;

}

void ITP::RunSetup(){

	//Initialize and fill the Eigen Wavefunction Storage
	wavefct = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);

	// the time-step sizes for Runge-Kutta integration for both schemes as complex valued variables
 	t_ITP = complex<double>(opt.ITP_step,0.0);

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

   	// The laplacian and potential coefficients needed for the ITP scheme
   	// Precomputing

   	PotentialGrid = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
   	for(int i = 0; i< opt.grid[1]; i++){for(int j = 0; j < opt.grid[2]; j++){
	PotentialGrid(i,j) = half * opt.omega_x * opt.omega_x * X(i) * X(i) +  half * opt.omega_y * opt.omega_y * Y(j) * Y(j);}}

	itp_laplacian_x = complex<double>(1.0,0.0) / (two * h_x * h_x);
	itp_laplacian_y = complex<double>(1.0,0.0) / (two * h_y * h_y);

	opt.stateInformation.resize(2);
	opt.stateInformation[0] = 1.0;
	opt.stateInformation[1] = 1.0;

}

inline void ITP::rescale(MatrixXcd &wavefct)
{	
	double Integral = 0.;  
	for(int i=0;i<opt.grid[1]-1;i++){
    	for(int j=0;j<opt.grid[2]-1;j++){
      		Integral += real(h_x)*real(h_y)*(abs2(wavefct(i,j))+abs2(wavefct(i+1,j))+abs2(wavefct(i,j+1))+abs2(wavefct(i+1,j+1)))/real(four);      
    	}
    }
    
	opt.scale_factor = opt.N/Integral;	
	// cout << "Integral : " << Integral << " scalefactor: " << opt.scale_factor << " " << sqrt(opt.scale_factor) << endl;
	// wavefct.array() *= sqrt(opt.scale_factor);
	wavefct.array() *= sqrt(opt.scale_factor);
}

void ITP::cli(string name,int counter_state, int counter_max, double start)
{
	if(counter_state%(counter_max/100)==0)
		{
			int seconds;
			int min;
			int hour;
			int total;

			total = omp_get_wtime() - start;
			hour = total / 3600;
			min = (total / 60) % 60;
			seconds = total % 60;

			cout << "  " << name << " with " << VORTICES_BUILD_TIME << " Steps: "
				 << std::setw(2) << std::setfill('0') << hour << ":"
				 << std::setw(2) << std::setfill('0') << min << ":"
				 << std::setw(2) << std::setfill('0') << seconds  << "    "
				 << std::setw(3) << std::setfill('0') << (counter_state/(counter_max/100)) << "%\r" << flush;
		}
	if(counter_state == counter_max)
	{
		cout << endl;
	}
}

void ITP::CopyComplexGridToEigen(){
	for(int i = 0; i < opt.grid[1]; i++){for(int j = 0; j < opt.grid[2]; j++){ wavefct(i,j) = pPsi->at(0,i,j,0);}}
}

void ITP::CopyEigenToComplexGrid(){
	for(int i = 0; i < opt.grid[1]; i++){for(int j = 0; j < opt.grid[2]; j++){ pPsi->at(0,i,j,0) = wavefct(i,j);}}
}

void ITP::formVortices(string runname){
	double start;  // starttime of the run

	CopyComplexGridToEigen();

	MatrixXcd wavefctcp = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k0 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k1 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k2 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k3 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);

	start = omp_get_wtime();

	//start loop here
	Eigen::initParallel();

	for(int m = 1; m <= VORTICES_BUILD_TIME; m++){

		wavefct.row(0) = VectorXcd::Zero(opt.grid[1]);
		wavefct.row(opt.grid[1]-1) = VectorXcd::Zero(opt.grid[1]);
		wavefct.col(0) = VectorXcd::Zero(opt.grid[2]);
		wavefct.col(opt.grid[2]-1) = VectorXcd::Zero(opt.grid[2]);


		wavefctcp = wavefct;

		ITP_compute_k(k0,wavefctcp);

		wavefctcp = wavefct + half * t_ITP * k0;
		ITP_compute_k(k1,wavefctcp);

		wavefctcp = wavefct + half * t_ITP * k1;
		ITP_compute_k(k2,wavefctcp);

		wavefctcp = wavefct + t_ITP * k2;
		ITP_compute_k(k3,wavefctcp);

		wavefct += (t_ITP/six) * ( k0 + two * k1 + two * k2 + k3);



		rescale(wavefct);

		cli(runname,m,VORTICES_BUILD_TIME,start);	
	}
	cout << endl;
	// update the ComplexGrid* DATA object outside of this.
	CopyEigenToComplexGrid();
}
void ITP::propagateToGroundState(string runname)
{
	double start;  // starttime of the run
	bool finished = false;
	int counter_finished = 0;
	int state = 0;
	int old_Ekin = 0;

	// load external Data into wavefct
	CopyComplexGridToEigen();

	Eval breakCondition;
	MatrixXcd wavefctcp = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k0 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k1 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k2 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k3 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);

	start = omp_get_wtime();

	//start loop here
	Eigen::initParallel();

	// do{
	// 	rescale(wavefct);
	// 	if(opt.scale_factor == 1){
	// 		finished = true;
	// 	}
	// } while (finished == false);

	// finished = false;

	// for(int m = 1; opt.scale_factor < 0.99 && opt.scale_factor > 1.01; m++){
	do {
		for(int m = 0; m < 100; m++){

			wavefct.row(0) = VectorXcd::Zero(opt.grid[1]);
			wavefct.row(opt.grid[1]-1) = VectorXcd::Zero(opt.grid[1]);
			wavefct.col(0) = VectorXcd::Zero(opt.grid[2]);
			wavefct.col(opt.grid[2]-1) = VectorXcd::Zero(opt.grid[2]);

			rescale(wavefct);

			wavefctcp = wavefct;
	
			ITP_compute_k(k0,wavefctcp);
	
			wavefctcp = wavefct + half * t_ITP * k0;
			ITP_compute_k(k1,wavefctcp);
	
			wavefctcp = wavefct + half * t_ITP * k1;
			ITP_compute_k(k2,wavefctcp);
	
			wavefctcp = wavefct + t_ITP * k2;
			ITP_compute_k(k3,wavefctcp);
	
			wavefct += (t_ITP/six) * ( k0 + two * k1 + two * k2 + k3);			
	
			state++;	
		}
		

		breakCondition.saveData(wavefct,opt,state,runname);
		breakCondition.evaluateDataITP();

		cli_groundState(runname,start,state,breakCondition.totalResult);
		int difference = breakCondition.totalResult.Ekin - old_Ekin;
		if(difference == 0){
		// if(opt.scale_factor == 0){
			counter_finished++;
		}else{
			counter_finished = 0;
		}
		old_Ekin = breakCondition.totalResult.Ekin;
		if(counter_finished >= 3){
			finished = true;
		}		
		// cli_plot(runname,m,runtime,start,plot);
	} while (finished == false);

	cout << endl;



	// update the ComplexGrid* DATA object outside of this.
	CopyEigenToComplexGrid();
}

void ITP::cli_groundState(string name, double start,int state,Observables totalResult){	

			int seconds;
			int min;
			int hour;
			int total;

			total = omp_get_wtime() - start;
			hour = total / 3600;
			min = (total / 60) % 60;
			seconds = total % 60;	

		cout << "  " << name << ":"
			<< std::setw(5) << std::setfill(' ') << state
			<< " Kinetic Energy: " << std::setw(5) << std::setfill(' ') << totalResult.Ekin << " "
			<< " Particles: " << std::setw(12) << std::setfill(' ') << totalResult.particle_count << " "
			<< " Volume: " << std::setw(12) << std::setfill(' ') << totalResult.volume << " "
			<< " Density: " << std::setw(12) << std::setfill(' ') << totalResult.density << " "
			<< " Scaling Factor: " << std::setw(12) << std::setfill(' ') << opt.scale_factor << " "
			<< std::setw(2) << std::setfill('0') << hour << ":"
			<< std::setw(2) << std::setfill('0') << min << ":"
			<< std::setw(2) << std::setfill('0') << seconds  << "\r" << flush;


}

inline void ITP::ITP_compute_k(MatrixXcd &k,MatrixXcd &wavefctcp){
	Matrix<std::complex<double>,Dynamic,Dynamic,ColMajor> wavefctcpX = Matrix<std::complex<double>,Dynamic,Dynamic,ColMajor>::Zero(opt.grid[1],opt.grid[2]);
	Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor> wavefctcpY = Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor>::Zero(opt.grid[1],opt.grid[2]);
	
	k = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);

	// laplacian
	for(int j = 1;j<opt.grid[2]-1;j++){
	for(int i = 1;i<opt.grid[1]-1;i++){
	wavefctcpX(i,j) = wavefctcp(i-1,j) - two * wavefctcp(i,j) + wavefctcp(i+1,j);
	wavefctcpY(i,j) = wavefctcp(i,j-1) - two * wavefctcp(i,j) + wavefctcp(i,j+1);
	}}
	k.noalias() += wavefctcpX * itp_laplacian_x + wavefctcpY * itp_laplacian_y;

	// interaction + potential
	k.array() -= (PotentialGrid.array() + complex<double>(opt.g,0.0) * ( wavefctcp.conjugate().array() * wavefctcp.array() )) * wavefctcp.array();

	}


