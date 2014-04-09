#define EIGEN_VECTORIZE
#define EIGEN_NO_DEBUG

#include <EXP2D_tools.h>
#include <2dexpan.h>
#include <omp.h>
#include <plot_with_mgl.h>


using namespace std;
using namespace Eigen;

EXP2D::EXP2D()
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

EXP2D::EXP2D(ComplexGrid* &c,Options &externaloptions)
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

EXP2D::~EXP2D(){}

void EXP2D::setOptions(Options &externaloptions){
	opt = externaloptions;

}

void EXP2D::RunSetup(){

	//Initialize and fill the Eigen Wavefunction Storage
	wavefct = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);

	// the time-step sizes for Runge-Kutta integration for both schemes as complex valued variables
	t_RTE = complex<double>(opt.RTE_step,0.0);
 	t_ITP = complex<double>(opt.ITP_step,0.0);

	// Maximum x and y ranges of the grid after expanding for the full runtime.
	// Needed to compute the growing plots.
	ranges.resize(2);
	complex<double> tmp3 = complex<double>(opt.RTE_step * opt.n_it_RTE,0.0);
	ranges[0] = opt.min_x * real(lambda_x(tmp3));
	ranges[1] = opt.min_x * real(lambda_y(tmp3));

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

   	laplacian_coefficient_x = VectorXcd::Zero(2 * opt.n_it_RTE);
   	laplacian_coefficient_y = VectorXcd::Zero(2 * opt.n_it_RTE);
   	gradient_coefficient_x = VectorXcd::Zero(2 * opt.n_it_RTE);
   	gradient_coefficient_y = VectorXcd::Zero(2 * opt.n_it_RTE);

   	complex<double> tmp;  	
   	for(int t = 0; t< ( 2* opt.n_it_RTE); t++){
   	tmp = half * complex<double>(t,0.0) * t_RTE;   	
   	laplacian_coefficient_x(t) = i_unit / ( two * h_x * h_x * lambda_x(tmp) * lambda_x(tmp) );
   	laplacian_coefficient_y(t) = i_unit / ( two * h_y * h_y * lambda_y(tmp) * lambda_y(tmp) );
   	gradient_coefficient_x(t) = lambda_x_dot(tmp) / (two * h_x * lambda_x(tmp));
   	gradient_coefficient_y(t) = lambda_y_dot(tmp) / (two * h_y * lambda_y(tmp));
   	}

   	// The laplacian and potential coefficients needed for the ITP scheme
   	// Precomputing

   	PotentialGrid = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
   	for(int i = 0; i< opt.grid[1]; i++){for(int j = 0; j < opt.grid[2]; j++){
	PotentialGrid(i,j) = half * opt.omega_x * opt.omega_x * X(i) * X(i) +  half * opt.omega_y * opt.omega_y * Y(j) * Y(j);}}

	itp_laplacian_x = complex<double>(1.0,0.0) / (two * h_x * h_x);
	itp_laplacian_y = complex<double>(1.0,0.0) / (two * h_y * h_y);

}



void EXP2D::rescale(MatrixXcd &wavefct)
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



void EXP2D::cli_plot(string name,int counter_state, int counter_max, double start,bool plot)
{
	int seconds;
	int min;
	int hour;
	int total;

	if(counter_state%(counter_max/100)==0)
		{
			if(plot == true)
				{
					opt.name = name + "-" + std::to_string(counter_state/(counter_max/100));
					// plotdatatopng(pPsi,opt);
					plotdatatopngEigen(wavefct,opt);
				}

			total = omp_get_wtime() - start;
			hour = total / 3600;
			min = total / 60;
			seconds = total % 60;

			cout << std::setw(2) << std::setfill('0') << hour << ":"
				 << std::setw(2) << std::setfill('0') << min << ":"
				 << std::setw(2) << std::setfill('0') << seconds  << "    "
				 << std::setw(3) << std::setfill('0') << (counter_state/(counter_max/100)) << "%\r" << flush;
		}
	if(counter_state == counter_max)
	{
		cout << endl;
	}
}

void EXP2D::cli_plot_expanding(vector<double> &ranges,string name,int counter_state, int counter_max, double start,bool plot)
{
	int seconds;
	int min;
	int hour;
	int total;

	if(counter_state%(counter_max/100)==0)
		{
			if(plot == true)
				{
					opt.name = name + "-" + std::to_string(counter_state/(counter_max/100));
					// plotdatatopng(pPsi,opt);
					plotdatatopngEigenExpanding(wavefct,ranges,Xexpanding,Yexpanding,opt);
				}

			total = omp_get_wtime() - start;
			hour = total / 3600;
			min = total / 60;
			seconds = total % 60;

			cout << std::setw(2) << std::setfill('0') << hour << ":"
				 << std::setw(2) << std::setfill('0') << min << ":"
				 << std::setw(2) << std::setfill('0') << seconds  << "    "
				 << std::setw(3) << std::setfill('0') << (counter_state/(counter_max/100)) << "%\r" << flush;
		}
	if(counter_state == counter_max)
	{
		cout << endl;
	}
}

void EXP2D::CopyComplexGridToEigen(){
	for(int i = 0; i < opt.grid[1]; i++){for(int j = 0; j < opt.grid[2]; j++){ wavefct(i,j) = pPsi->at(0,i,j,0);}}
}

void EXP2D::CopyEigenToComplexGrid(){
	for(int i = 0; i < opt.grid[1]; i++){for(int j = 0; j < opt.grid[2]; j++){ pPsi->at(0,i,j,0) = wavefct(i,j);}}
}

void EXP2D::itpToTime(string runname, int runtime, bool plot)
{
	double start;  // starttime of the run

	// load external Data into wavefct
	CopyComplexGridToEigen();

	MatrixXcd wavefctcp = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k0 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k1 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k2 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k3 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);

	start = omp_get_wtime();

	cout << " " << runname << endl;	
	//start loop here
	Eigen::initParallel();
	for(int m = 1; m <= runtime; m++){

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

		cli_plot(runname,m,runtime,start,plot);
	}

	// update the ComplexGrid* DATA object outside of this.
	CopyEigenToComplexGrid();


	cout << "\n";
}

inline void EXP2D::ITP_compute_k(MatrixXcd &k,MatrixXcd &wavefctcp)
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





void EXP2D::rteToTime(string runname, int runtime, bool plot)
{
	double start;  // starttime of the run
	int t = 0;		// counter for the expanding lambdavectors with coefficients

	CopyComplexGridToEigen();

	MatrixXcd wavefctcp(opt.grid[1],opt.grid[2]);
	MatrixXcd k0 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k1 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k2 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);
	MatrixXcd k3 = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);	
	
	start = omp_get_wtime();

	cout << " " << runname << endl;	
	//start loop here
	Eigen::initParallel();

	if(opt.startgrid[2] == false){
	for(int m = 1; m <= runtime; m++){

		wavefctcp = wavefct;

		//boundary conditions -- Dirichlet

		// wavefct.row(0) = VectorXcd::Zero(opt.grid[1]);
		// wavefct.row(opt.grid[1]-1) = VectorXcd::Zero(opt.grid[1]);
		// wavefct.col(0) = VectorXcd::Zero(opt.grid[2]);
		// wavefct.col(opt.grid[2]-1) = VectorXcd::Zero(opt.grid[2]);

		//boundary conditions end

		RTE_compute_k(k0,wavefctcp,t);
		wavefctcp = wavefct + half * t_RTE * k0;

		t += 1;
		RTE_compute_k(k1,wavefctcp,t);
		wavefctcp = wavefct + half * t_RTE * k1;

		RTE_compute_k(k2,wavefctcp,t);		
		wavefctcp = wavefct + t_RTE * k2;

		t += 1;
		RTE_compute_k(k3,wavefctcp,t);

		wavefct += (t_RTE/six) * ( k0 + two * k1 + two * k2 + k3);

		complex<double> tmp2 = complex<double>(m,0.0) * t_RTE;
		Xexpanding = x_expand(tmp2);
   		Yexpanding = y_expand(tmp2);

		cli_plot_expanding(ranges,runname,m,runtime,start,plot);

	}}

	if(opt.startgrid[2] == true){ // with potential and no expansion
	for(int m = 1; m <= runtime; m++){

		wavefctcp = wavefct;

		//boundary conditions -- Dirichlet

		// wavefct.row(0) = VectorXcd::Zero(opt.grid[1]);
		// wavefct.row(opt.grid[1]-1) = VectorXcd::Zero(opt.grid[1]);
		// wavefct.col(0) = VectorXcd::Zero(opt.grid[2]);
		// wavefct.col(opt.grid[2]-1) = VectorXcd::Zero(opt.grid[2]);

		//boundary conditions end

		RTE_compute_k_pot(k0,wavefctcp,t);
		wavefctcp = wavefct + half * t_RTE * k0;

		t += 1;
		RTE_compute_k_pot(k1,wavefctcp,t);
		wavefctcp = wavefct + half * t_RTE * k1;

		RTE_compute_k_pot(k2,wavefctcp,t);		
		wavefctcp = wavefct + t_RTE * k2;

		t += 1;
		RTE_compute_k_pot(k3,wavefctcp,t);

		wavefct += (t_RTE/six) * ( k0 + two * k1 + two * k2 + k3);

		cli_plot(runname,m,runtime,start,plot);
	}}

// update the ComplexGrid* DATA object outside of this.
CopyEigenToComplexGrid();

cout << "\n";

}

inline void EXP2D::RTE_compute_k(MatrixXcd &k,MatrixXcd &wavefctcp,int &t)
	{
	Matrix<std::complex<double>,Dynamic,Dynamic,ColMajor> wavefctcpX = Matrix<std::complex<double>,Dynamic,Dynamic,ColMajor>::Zero(opt.grid[1],opt.grid[2]);
	Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor> wavefctcpY = Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor>::Zero(opt.grid[1],opt.grid[2]);
	k = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);


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

	//interaction
	k.array() -= complex<double>(0.0,opt.g) * ( wavefctcp.conjugate().array() * wavefctcp.array() ) * wavefctcp.array();

	}

inline void EXP2D::RTE_compute_k_pot(MatrixXcd &k,MatrixXcd &wavefctcp,int &t)
	{
	Matrix<std::complex<double>,Dynamic,Dynamic,ColMajor> wavefctcpX = Matrix<std::complex<double>,Dynamic,Dynamic,ColMajor>::Zero(opt.grid[1],opt.grid[2]);
	Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor> wavefctcpY = Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor>::Zero(opt.grid[1],opt.grid[2]);
	k = MatrixXcd::Zero(opt.grid[1],opt.grid[2]);


	//laplacian
	for(int j = 1;j<opt.grid[2]-1;j++){
	for(int i = 1;i<opt.grid[1]-1;i++){
	wavefctcpX(i,j) = wavefctcp(i-1,j) - two * wavefctcp(i,j) + wavefctcp(i+1,j);
	wavefctcpY(i,j) = wavefctcp(i,j-1) - two * wavefctcp(i,j) + wavefctcp(i,j+1);
	}}
	k.noalias() +=   wavefctcpX * itp_laplacian_x * i_unit + wavefctcpY * itp_laplacian_x * i_unit;

	//Potential + Interaction
	k.array() -= ( i_unit * PotentialGrid.array() + complex<double>(0.0,opt.g) * ( wavefctcp.conjugate().array() * wavefctcp.array() )) * wavefctcp.array();

	}

