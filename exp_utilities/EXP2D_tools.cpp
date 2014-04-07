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
	h_x = 0;
	h_y = 0;
	Integral = 0;
	Integral_aux = 0;
	pi = M_PI;
	zero=complex<double>(0,0),half=complex<double>(0.5,0),one=complex<double>(1,0),two=complex<double>(2,0),four=complex<double>(4,0),six=complex<double>(6,0),i_unit=complex<double>(0,1);
}

EXP2D::EXP2D(ComplexGrid* &c,Options &opt)
{	
	omp_set_num_threads(omp_get_max_threads());
	cout << "Max Number of Threads: " << omp_get_max_threads() << endl;	
	cout << "Eigenthreads: " << Eigen::nbThreads() << endl;

  	pPsi = c;  
  	h_x = complex<double>((2.*opt.min_x/opt.grid[1]),0.0);
  	h_y = complex<double>((2.*opt.min_y/opt.grid[2]),0.0); 
  	x_axis.resize(opt.grid[1]);
  	y_axis.resize(opt.grid[2]);
  	for(int i=0;i<opt.grid[1];i++){x_axis[i]=-opt.min_x+i*real(h_x);}
  	for(int j=0;j<opt.grid[2];j++){y_axis[j]=-opt.min_y+j*real(h_y);}
  	pi = M_PI;
  	dispersion_x = opt.dispersion_x;
  	dispersion_y = opt.dispersion_y;
  	exp_factor = opt.exp_factor;
  	g = opt.g;
 	zero=complex<double>(0,0),half=complex<double>(0.5,0),one=complex<double>(1,0),two=complex<double>(2,0),four=complex<double>(4,0),six=complex<double>(6,0),i_unit=complex<double>(0,1);
 	t_RTE = complex<double>(opt.RTE_step,0);

 	//compute lambda_squareds

   laplacian_coefficient_x = vector<complex<double>>(2 * opt.n_it_RTE);
   laplacian_coefficient_y = vector<complex<double>>(2 * opt.n_it_RTE);
   gradient_coefficient_x = vector<complex<double>>(2 * opt.n_it_RTE);
   gradient_coefficient_y = vector<complex<double>>(2 * opt.n_it_RTE);


  complex<double> tmp;
   for(int t = 0; t< ( 2* opt.n_it_RTE); t++)
   {
   	tmp = half * complex<double>(t,0.0) * t_RTE;
   	laplacian_coefficient_x[t] = i_unit / ( two * h_x * h_x * lambda_x(tmp) * lambda_x(tmp) );
   	laplacian_coefficient_y[t] = i_unit / ( two * h_y * h_y * lambda_y(tmp) * lambda_y(tmp) );
   	gradient_coefficient_x[t] = lambda_x_dot(tmp) / (two * h_x * lambda_x(tmp));
   	gradient_coefficient_y[t] = lambda_y_dot(tmp) / (two * h_y * lambda_y(tmp));   	
   }

}

EXP2D::~EXP2D(){}



void EXP2D::rescale(ComplexGrid* & pPsi, Options &opt)
{	
	Integral=integral(pPsi,opt);
	opt.scale_factor=complex<double>(opt.N,0)/Integral;
	
	for(int i=0;i<opt.grid[1];i++)
	{
		for(int j=0;j<opt.grid[2];j++)
		{
		  pPsi->at(0,i,j,0)*=sqrt(opt.scale_factor);
		}
	}
}



void EXP2D::cli_plot(MatrixXcd& mPsi, Options &opt, string name,int counter_state, int counter_max, double start,bool plot)
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
					plotdatatopngEigen(mPsi,opt);
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

void EXP2D::run_status(string name,int counter_state, int counter_max, double start)
{
	int seconds;
	int min;
	int hour;
	int total;

	if(counter_state%(counter_max/100)==0)
		{
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


void EXP2D::Neumann(ComplexGrid &k,ComplexGrid &PsiCopy,Options &opt){
      //Fixed derivative (Neumann) boundary conditions for the ITP


    for(int i=0;i<opt.grid[1];i++) 
    { 
		k(0,i,0,0)=itp_potential(PsiCopy,i,0,opt);
		k(0,i,opt.grid[1]-1,0)= itp_potential(PsiCopy,i,opt.grid[1]-1,opt);
    }

    for(int j=0;j<opt.grid[2];j++)
    { 
		k(0,0,j,0)=itp_potential(PsiCopy,0,j,opt);
		k(0,opt.grid[2]-1,j,0)= itp_potential(PsiCopy,opt.grid[2]-1,j,opt);
    }
}


void EXP2D::computeK_ITP(ComplexGrid* &pPsi, vector<ComplexGrid> &k,Options &opt,complex<double> &t_ITP){ 
	
	ComplexGrid PsiCopy(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	for(int i = 0; i < opt.grid[1]; i++){for(int j = 0; j < opt.grid[2]; j++){ PsiCopy(0,i,j,0) = pPsi->at(0,i,j,0);}}

	
	

	for(int d=0; d<4; d++){

  // The k's have to be computed differently, this is decided by int d
		switch (d){
			case 1 : case 2 :
				#pragma omp parallel for
      			for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ PsiCopy(0,i,j,0)=pPsi->at(0,i,j,0)+half*t_ITP*k[d-1](0,i,j,0) ;}}
      			break;
      		case 3 :
      			#pragma omp parallel for
    			for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ PsiCopy(0,i,j,0)=pPsi->at(0,i,j,0)+t_ITP*k[d-1](0,i,j,0) ;}}	
    			break;
		}

    // Compute all k[1]..k[4] for the EXP2D
		#pragma omp parallel for 
    	for(int i=1;i<opt.grid[1]-1;i++) {for(int j=1;j<opt.grid[2]-1;j++) { k[d](0,i,j,0) = itp_kinetic(PsiCopy,i,j) + itp_potential(PsiCopy,i,j,opt) ;}}
    
    	Neumann(k[d],PsiCopy,opt);

	}
  
}



void EXP2D::ITP(ComplexGrid* & pPsi, Options &opt)
{
	vector<ComplexGrid> k(4);
	complex<double> t_ITP(opt.ITP_step,0); //Timetep size for ITP (the equations already assume imaginary time so a real t_ITP should be used)
	
	for ( int d=0;d<4;d++){k[d] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);}
	
	// Compute all the k's for the EXP2D	

	computeK_ITP(pPsi,k,opt,t_ITP);
	
	// Use the k's of EXP2D to timestep Psi	

	#pragma omp parallel for
	for(int i=0;i<opt.grid[1];i++){
	  for(int j=0;j<opt.grid[2];j++){ 	
	    pPsi->at(0,i,j,0)+=(t/six)*(k[0](0,i,j,0)+two*k[1](0,i,j,0)+two*k[2](0,i,j,0)+k[3](0,i,j,0));
	   }
	}

	// rescale Psi to conserve number of particles
	
	rescale(pPsi,opt); 
}

void EXP2D::itpToTime(Options &opt,bool plot)
{
	int counter_ITP = 0;
	double start;

	start = omp_get_wtime();
	
	cout << " " << opt.name << endl;
	string tmp = opt.name;
	for(int k=0;k<opt.n_it_ITP;k++)
	{ 
		ITP(pPsi,opt);		

  		counter_ITP += 1;

  		run_status(tmp,counter_ITP,opt.n_it_ITP,start);

  			
	}
	cout << "\n";	
}




void EXP2D::rteToTime(Options &opt, bool plot)
{
	int grid_x = opt.grid[1]; // Shorts for the Gridsize
	int grid_y = opt.grid[2];
	double start;  // starttime of the run
	int t = 0;		// counter for the expanding lambdavectors with coefficients

	// rte coefficients list
	// gridx gridy
	// t_RTE

	MatrixXcd wavefct(grid_x,grid_y);
	MatrixXcd wavefctcp(grid_x,grid_y);
	MatrixXcd k0 = MatrixXcd::Zero(grid_x,grid_y);
	MatrixXcd k1 = MatrixXcd::Zero(grid_x,grid_y);
	MatrixXcd k2 = MatrixXcd::Zero(grid_x,grid_y);
	MatrixXcd k3 = MatrixXcd::Zero(grid_x,grid_y);

	for(int i = 0; i < grid_x; i++){for(int j = 0; j < grid_y; j++){ wavefct(i,j) = pPsi->at(0,i,j,0);}}
	
	VectorXcd X(opt.grid[1]), Y(opt.grid[2]);
	for(int i = 0;i<grid_x;i++){X(i) = complex<double>(x_axis[i],0.0);}
	for(int j = 0;j<grid_y;j++){Y(j) = complex<double>(y_axis[j],0.0);}

	start = omp_get_wtime();

	cout << " " << opt.name << endl;	
	//start loop here
	Eigen::initParallel();
	for(int m = 1; m <= opt.n_it_RTE; m++){

		wavefctcp = wavefct;

		//boundary conditions -- Dirichlet

		wavefct.row(0) = VectorXcd::Zero(grid_x);
		wavefct.row(grid_x-1) = VectorXcd::Zero(grid_x);
		wavefct.col(0) = VectorXcd::Zero(grid_y);
		wavefct.col(grid_y-1) = VectorXcd::Zero(grid_y);

		//boundary conditions end

		RTE_compute_k(k0,wavefctcp,X,Y,t,grid_x,grid_y);
		wavefctcp = wavefct + half * t_RTE * k0;

		t += 1;
		RTE_compute_k(k1,wavefctcp,X,Y,t,grid_x,grid_y);
		wavefctcp = wavefct + half * t_RTE * k1;

		RTE_compute_k(k2,wavefctcp,X,Y,t,grid_x,grid_y);		
		wavefctcp = wavefct + t_RTE * k2;

		t += 1;
		RTE_compute_k(k3,wavefctcp,X,Y,t,grid_x,grid_y);

		wavefct += (t_RTE/six) * ( k0 + two * k1 + two * k2 + k3);

		cli_plot(wavefct,opt,"RTE",m,opt.n_it_RTE,start,plot);

	}

cout << "\n";

}

inline void EXP2D::RTE_compute_k(MatrixXcd &k,MatrixXcd &wavefctcp, VectorXcd &X,VectorXcd &Y,int &t,int &grid_x,int &grid_y)
	{
	Matrix<std::complex<double>,Dynamic,Dynamic,ColMajor> wavefctcpX = Matrix<std::complex<double>,Dynamic,Dynamic,ColMajor>::Zero(grid_x,grid_y);
	Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor> wavefctcpY = Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor>::Zero(grid_x,grid_y);
	k = MatrixXcd::Zero(grid_x,grid_y);


	//laplacian

	// #pragma omp parallel for
	for(int j = 1;j<grid_y-1;j++){
	for(int i = 1;i<grid_x-1;i++){
	wavefctcpX(i,j) = wavefctcp(i-1,j) - two * wavefctcp(i,j) + wavefctcp(i+1,j);
	wavefctcpY(i,j) = wavefctcp(i,j-1) - two * wavefctcp(i,j) + wavefctcp(i,j+1);
	}}
	k.noalias() +=   wavefctcpX * laplacian_coefficient_x[t] + wavefctcpY * laplacian_coefficient_y[t];
	//laplacian end

	// gradient

	// #pragma omp parallel for
	for(int j = 1;j<grid_y-1;j++){
	for(int i = 1;i<grid_x-1;i++){
	wavefctcpX(i,j) = wavefctcp(i+1,j) - wavefctcp(i-1,j);
	wavefctcpY(i,j) = wavefctcp(i,j+1) - wavefctcp(i,j-1);
	}}

	for(int i = 0;i<grid_x;i++){ wavefctcpY.row(i).array() *= Y.array(); }
	for(int j = 0;j<grid_y;j++){ wavefctcpX.col(j).array() *= X.array(); }

	k.noalias() += wavefctcpX * gradient_coefficient_x[t] + wavefctcpY * gradient_coefficient_y[t];
	// // gradient end

	//interaction
	k.array() -= complex<double>(0.0,g) * ( wavefctcp.conjugate().array() * wavefctcp.array() ) * wavefctcp.array();

	//interaction end

	}

