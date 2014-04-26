#define EIGEN_VECTORIZE
#define EIGEN_NO_DEBUG

#include <EXP2D_itp.hpp>
#include <omp.h>

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

}

inline void ITP::rescale(MatrixXcd &wavefct)
{	
	Integral= complex<double>(0,0);  
	for(int i=0;i<opt.grid[1]-1;i++){
    for(int j=0;j<opt.grid[2]-1;j++)
    {
      Integral += h_x*h_y*(norm(wavefct(i,j))+norm(wavefct(i+1,j))+norm(wavefct(i,j+1))+norm(wavefct(i+1,j+1)))/four;      
    }}
	opt.scale_factor=complex<double>(opt.N,0)/Integral;	
	wavefct.array() *= sqrt(opt.scale_factor);
}

void ITP::cli_plot(string name,int counter_state, int counter_max, double start,bool plot)
{
	if(counter_state%(counter_max/100)==0)
		{
			int seconds;
			int min;
			int hour;
			int total;

			if(plot == true)
				{
					opt.name = name; //+ "-" + std::to_string(counter_state/(counter_max/100));
					// plotdatatopng(pPsi,opt);
					plotdatatopngEigen(wavefct,opt);

					// // kvalue analysis
					// CopyEigenToComplexGrid();
					// ComplexGrid::fft(*pPsi,*pK,true);
					// opt.name = "kvalues -" + std::to_string(counter_state/(counter_max/100));
					// plotdatatopng(pK,opt);


				}

			total = omp_get_wtime() - start;
			hour = total / 3600;
			min = total / 60;
			seconds = total % 60;

			cout << "  " << name << " reached: "
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

void ITP::itpToTime(string runname, bool plot)
{
	double start;  // starttime of the run
	bool finished = false;
	int counter_finished = 0;
	int state = 0;
	double oldabsolute = 0.0;

	// load external Data into wavefct
	CopyComplexGridToEigen();

	MatrixXcd wavefctcp = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k0 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k1 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k2 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k3 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);

	start = omp_get_wtime();

	//start loop here
	Eigen::initParallel();
	// for(int m = 1; opt.scale_factor < 0.99 && opt.scale_factor > 1.01; m++){
	do {

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

		state++;
		counter_finished += cli_itp(runname, start,state,oldabsolute);

		if(counter_finished >= 1){
			finished = true;
		}
		// cli_plot(runname,m,runtime,start,plot);
	} while (finished == false);

	cout << endl;



	// update the ComplexGrid* DATA object outside of this.
	CopyEigenToComplexGrid();
}

int ITP::cli_itp(string name, double start,int state, double &oldabsolute){	
	int counter;
	double absolute = abs(opt.scale_factor.real() - oldabsolute);
		if (absolute  <= 1){
		counter = 1;
	}else{counter = 0;}

			int seconds;
			int min;
			int hour;
			int total;

			total = omp_get_wtime() - start;
			hour = total / 3600;
			min = total / 60;
			seconds = total % 60;	

		cout << "  " << name << " Runstep: " << state << "  Scalefactor: "<< std::setw(5) << opt.scale_factor.real() << "  " << "Absolute: " << std::setw(5) << absolute << "  "
			<< std::setw(2) << std::setfill('0') << hour << ":"
			<< std::setw(2) << std::setfill('0') << min << ":"
			<< std::setw(2) << std::setfill('0') << seconds  << "\r" << flush;


			oldabsolute = opt.scale_factor.real();
	return counter;

}

inline void ITP::ITP_compute_k(MatrixXcd &k,MatrixXcd &wavefctcp)
	{
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


