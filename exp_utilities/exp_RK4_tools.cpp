#include <exp_RK4_tools.h>
#include <2dexpan.h>
#include <omp.h>
#include <plot_with_mgl.h>


using namespace std;

RK4::RK4()
{
	h_x = 0;
	h_y = 0;
	Integral = 0;
	Integral_aux = 0;
	pi = M_PI;
	zero=complex<double>(0,0),half=complex<double>(0.5,0),one=complex<double>(1,0),two=complex<double>(2,0),four=complex<double>(4,0),six=complex<double>(6,0),i_unit=complex<double>(0,1);
}

RK4::RK4(ComplexGrid* &c,Options &opt)
{	
	opt.threads = omp_get_max_threads();
	cout << "Max Number of Threads: " << opt.threads << endl;
	omp_set_num_threads(opt.threads);
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

   lambda_x_squared = vector<complex<double>>(2 * opt.n_it_RTE);
   lambda_y_squared = vector<complex<double>>(2 * opt.n_it_RTE);
   lambda_dot_y = vector<complex<double>>(2 * opt.n_it_RTE);
   lambda_dot_x = vector<complex<double>>(2 * opt.n_it_RTE);


   for(int t = 0; t< ( 2* opt.n_it_RTE); t++)
   {
   	tmp = half * complex<double>(t,0.0) * t_RTE;
   	lambda_x_squared[t] = l_x(tmp) * l_x(tmp);
   	lambda_y_squared[t] = l_y(tmp) * l_y(tmp);
   	lambda_dot_x[t] = l_x_dot(tmp) / l_x(tmp);
   	lambda_dot_y[t] = l_y_dot(tmp) / l_y(tmp);
   	
  //  	if(t%((2*opt.n_it_RTE)/100)==0)
		// {
		// cout << lambda_x_squared[t] << endl;;
		// }


   }




// // X.resize(opt.grid[1],opt.grid[2]);
// // Y.resize(opt.grid[1],opt.grid[2]);  
// // L.resize(opt.grid[1],opt.grid[2]);  
// Eigen::MatrixXcd G(opt.grid[1],opt.grid[2]);
// Eigen::MatrixXcd X(opt.grid[1],opt.grid[2]); 
// Eigen::MatrixXcd Y(opt.grid[1],opt.grid[2]); 
// Eigen::MatrixXcd L(opt.grid[1],opt.grid[2]);   
// // X.reserve(opt.grid[1]);
// // Y.reserve(opt.grid[2]);
// // L.reserve(3 * (opt.grid[1]-2));

// // G.reserve(2 * (opt.grid[1]-2));
//  	for(int i = 0; i<opt.grid[1];i++)
// 	{
// 		X(i,i) = complex<double>(x_axis[i],0.0);
// 	}
// 	// X.makeCompressed();

// 	cout << " X: \n" << X << endl;

// 	for(int j = 0; j<opt.grid[2];j++)
// 	{
// 		Y(j,j) = complex<double>(y_axis[j],0.0);
// 	}
// 	// Y.makeCompressed();

// 	cout << " Y: \n"<< Y << endl;

// 	for(int i = 1;i<opt.grid[1]-1;i++){
// 	L(i,i-1) = complex<double>(1.0,0.0);
// 	L(i,i) = complex<double>(-2.0,0.0);
// 	L(i,i+1) = complex<double>(1.0,0.0);
// 	}
// 	// L.makeCompressed();

// 	cout<< " L: \n" << L << endl;


// 	for(int i = 1;i<opt.grid[1]-1;i++){
// 	G(i,i-1) = complex<double>(-1.0,0.0);
// 	G(i,i+1) = complex<double>(1.0,0.0);
// 	}
// 	// G.makeCompressed();
// 	cout << " G: \n"<< G << endl;

}

RK4::~RK4(){}

// double RK4::gauss(double x,double y){return (exp(-x*x-y*y));}

complex<double> RK4::interaction(complex<double> a,Options &opt)
{return (opt.g*norm(a));} //Interaction term in the GPE Hamiltonian 

complex<double> RK4::integral(ComplexGrid* & pPsi,Options &opt)
{	
	Integral_aux = complex<double>(0,0);	
	for(int i=0;i<opt.grid[1]-1;i++)
	{
		for(int j=0;j<opt.grid[2]-1;j++)
		{
			Integral_aux+=h_x*lambda_x(opt)*h_y*lambda_y(opt)*(norm(pPsi->at(0,i,j,0))+norm(pPsi->at(0,i+1,j,0))+norm(pPsi->at(0,i,j+1,0))+norm(pPsi->at(0,i+1,j+1,0)))/four;
			
		}
	}
	return Integral_aux;
}

void RK4::rescale(ComplexGrid* & pPsi, Options &opt)
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

double RK4::phase_save(ComplexGrid* &pPsi,int a,int b) //Definition of phase 
{
 //        if(arg(pPsi->at(0,a,b,0))<0){ return 2*pi+arg(pPsi->at(0,a,b,0)); } //arg uses the atan2 function, so the same applies
	// else{ 
		return arg(pPsi->at(0,a,b,0)); 
	// }
}

void RK4::save_2D(ComplexGrid* &pPsi,Options &opt) //Function to save the data to a file (with compatible blocks for gnuplot)
{	
	openDataFiles_obdm(opt.name,1,opt.times); //Open file with name (1,name)
	
        for(int i=0;i<opt.grid[1];i++)
        {
	        for(int j=0;j<opt.grid[2];j++)
	        {
			double x;
			double y;
			x = x_axis[i];
			y = y_axis[j];
		    save_obdm(x*real(lambda_x(opt)),y*real(lambda_y(opt)),norm(pPsi->at(0,i,j,0)),phase_save(pPsi,i,j));
	        }
                blank_line();
        }
        blank_line();
	closeDataFiles_obdm();
}

complex<double> RK4::itp_kinetic(ComplexGrid &PsiCopy,int i, int j)
{
    return half*((PsiCopy(0,i+1,j,0)-(two*PsiCopy(0,i,j,0))+PsiCopy(0,i-1,j,0))/(h_x*h_x))+half*((PsiCopy(0,i,j+1,0)-(two*PsiCopy(0,i,j,0))+PsiCopy(0,i,j-1,0))/(h_x*h_x)); 
}

complex<double> RK4::itp_potential(ComplexGrid & PsiCopy,int i, int j,Options & opt)
{ 
    complex<double> xvalue = complex<double>(x_axis[i],0);
    complex<double> yvalue = complex<double>(y_axis[j],0);
  
    return -(half*opt.omega_x*opt.omega_x*xvalue*xvalue+half*opt.omega_y*opt.omega_y*yvalue*yvalue+complex<double>(opt.g,0)*norm(PsiCopy(0,i,j,0)))*PsiCopy(0,i,j,0);
}




void RK4::TimeStepRK4(ComplexGrid* &pPsi,vector<ComplexGrid> &k,Options &opt,complex<double> &t)
{	

	#pragma omp parallel for
	for(int i=0;i<opt.grid[1];i++){
	  for(int j=0;j<opt.grid[2];j++){ 	
	    pPsi->at(0,i,j,0)+=(t/six)*(k[0](0,i,j,0)+two*k[1](0,i,j,0)+two*k[2](0,i,j,0)+k[3](0,i,j,0));
	   }
	}
}	


void RK4::cli_plot(Eigen::MatrixXcd& mPsi, Options &opt, string name,int counter_state, int counter_max, double start,bool plot)
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


void RK4::Neumann(ComplexGrid &k,ComplexGrid &PsiCopy,Options &opt){
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


void RK4::computeK_ITP(ComplexGrid* &pPsi, vector<ComplexGrid> &k,Options &opt,complex<double> &t_ITP){ 
	
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

    // Compute all k[1]..k[4] for the RK4
		#pragma omp parallel for 
    	for(int i=1;i<opt.grid[1]-1;i++) {for(int j=1;j<opt.grid[2]-1;j++) { k[d](0,i,j,0) = itp_kinetic(PsiCopy,i,j) + itp_potential(PsiCopy,i,j,opt) ;}}
    
    	Neumann(k[d],PsiCopy,opt);

	}
  
}



void RK4::ITP(ComplexGrid* & pPsi, Options &opt)
{
	vector<ComplexGrid> k(4);
	complex<double> t_ITP(opt.ITP_step,0); //Timetep size for ITP (the equations already assume imaginary time so a real t_ITP should be used)
	
	for ( int d=0;d<4;d++){k[d] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);}
	
	// Compute all the k's for the RK4	

	computeK_ITP(pPsi,k,opt,t_ITP);
	
	// Use the k's of RK4 to timestep Psi	

	TimeStepRK4(pPsi,k,opt,t_ITP);
	
	// rescale Psi to conserve number of particles
	
	rescale(pPsi,opt); 
}

void RK4::itpToTime(Options &opt,bool plot)
{
	// int counter_ITP = 0;
	// double start;

	// start = omp_get_wtime();
	
	cout << " " << opt.name << endl;
	string tmp = opt.name;
	for(int k=0;k<opt.n_it_ITP;k++)
	{ 
		ITP(pPsi,opt);		

  		// counter_ITP += 1;

  		// cli_plot(opt,tmp,counter_ITP,opt.n_it_ITP,start,plot);
  			
	}
	cout << "\n";	
}

// complex<double> RK4::function_RTE(ComplexGrid &wavefct,int i, int j, Options &opt)
// {
// 	complex<double> tmp;

// 	tmp = rte_kinetic(wavefct,i,j,opt);

// 	if(opt.g != 0.0)
// 	{tmp -= (rte_interaction(wavefct,i,j,opt) * wavefct(0,i,j,0));}

// 	if(opt.exp_factor != 0.0)
// 	{tmp += rte_expandingframe(wavefct,i,j,opt);}

// 	if(opt.startgrid[2] == true)
// 	{tmp -= rte_potential(i,j,opt) * wavefct(0,i,j,0);}

// 	return tmp;
// }



void RK4::Dirichlet(ComplexGrid* &pPsi,Options &opt){

        //Fixed (Dirichlet) boundary conditions for the RTE
    for(int l=0;l<opt.grid[1];l++){ 
	  pPsi->at(0,l,0,0).real(0.0);
	  pPsi->at(0,l,opt.grid[2]-1,0).real(0.0);
	}
    for(int m=0;m<opt.grid[2];m++){
	  pPsi->at(0,0,m,0).real(0.0);
	  pPsi->at(0,opt.grid[1]-1,m,0).real(0.0);
	}
}

// void RK4::computeK_RTE(ComplexGrid* &pPsi, vector<ComplexGrid> &k,Options &opt,complex<double> &t_RTE){

// 	int grid_x = opt.grid[1];
// 	int grid_y = opt.grid[2];

// 	ComplexGrid PsiCopy(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

// 	for(int i = 0; i < grid_x; i++){for(int j = 0; j < grid_y; j++){ PsiCopy(0,i,j,0) = pPsi->at(0,i,j,0);}}

// 	Dirichlet(pPsi,opt);
	
// 	complex<double> timestep[3];
// 	timestep[0] = half * t_RTE;
// 	timestep[1] = half * t_RTE;
// 	timestep[2] = t_RTE;

// 	#pragma omp parallel for
// 	for(int i=1;i<grid_x-1;i++){for(int j=1;j<grid_y-1;j++){ k[0](0,i,j,0)=function_RTE(PsiCopy,i,j,opt);}}

// 	opt.t_abs += half * t_RTE;
      
//     for(int d=0;d<3;d++)
// 	{  

// 		if(d == 2){opt.t_abs += half * t_RTE;}

// 		#pragma omp parallel for
// 		for(int i=0;i<grid_x;i++){for(int j=0;j<grid_y;j++){ PsiCopy(0,i,j,0) = pPsi->at(0,i,j,0) + timestep[d]*k[d](0,i,j,0) ;}}

// 		#pragma omp parallel for
// 		for(int i=1;i<grid_x-1;i++){for(int j=1;j<grid_y-1;j++){ k[d+1](0,i,j,0)=function_RTE(PsiCopy,i,j,opt);}}

// 	}
// }

// void RK4::functionEigen_RTE(Eigen::MatrixXcd& k,Eigen::MatrixXcd& mPsiCopy,complex<double>& t)
// {	
// 	cout << "1" << endl;
// 	k.noalias() += complex<double>(0.0,0.5) * ((L * mPsiCopy) / ( l_x(t) * l_x(t) ) + ( mPsiCopy * LT ) / (l_y(t) * l_y(t) ) );
// 	cout << "2" << endl;
// 	k.noalias() += ( (l_x_dot(t) / l_x(t)) * X * ( G * mPsiCopy) );
// 	cout << "3" << endl;
// 	k.noalias() += ( (l_y_dot(t) / l_y(t)) * ( mPsiCopy * GT ) * Y );
// 	cout << "4" << endl;
// 	k.array() += - complex<double>(0.0, g ) * (mPsiCopy.adjoint().array() * mPsiCopy.array()) * mPsiCopy.array();
// 	cout << "5" << endl;
// }

// void RK4::computeWithEigen_RTE(Options &opt, bool plot)
// {	
// 	int grid_x = opt.grid[1];
// 	int grid_y = opt.grid[2];
// 	complex<double> t(0.0,0.0);
// 	double start;
// 	complex<double> t_RTE(opt.RTE_step,0); 


// 	Eigen::MatrixXcd mPsi(grid_x,grid_y);
// 	Eigen::MatrixXcd mPsiCopy(grid_x,grid_y);
// 	Eigen::MatrixXcd k0(grid_x,grid_y);
// 	Eigen::MatrixXcd k1(grid_x,grid_y);
// 	Eigen::MatrixXcd k2(grid_x,grid_y);
// 	Eigen::MatrixXcd k3(grid_x,grid_y);
// 	Eigen::MatrixXcd L(grid_x,grid_y);
// 	Eigen::MatrixXcd G(grid_x,grid_y);
// 	Eigen::MatrixXcd X(grid_x,grid_y);
// 	Eigen::MatrixXcd Y(grid_x,grid_y);

// 	Eigen::SparseMatrix<std::complex<double>,1,std::ptrdiff_t >  Lsparse, Gsparse, Xsparse, Ysparse;
// 	Lsparse.resize(grid_x,grid_y);
// 	Gsparse.resize(grid_x,grid_y);
// 	Xsparse.resize(grid_x,grid_y);
// 	Ysparse.resize(grid_x,grid_y);

// 	Lsparse.reserve(3 * (grid_x-2));
// 	Xsparse.reserve(opt.grid[1]);
// 	Ysparse.reserve(opt.grid[2]);
// 	Gsparse.reserve(2 * (opt.grid[1]-2));



// 		for(int i = 1;i<opt.grid[1]-1;i++){
// 	Lsparse.insert(i,i-1) = complex<double>(1.0,0.0);
// 	Lsparse.insert(i,i) = complex<double>(-2.0,0.0);
// 	Lsparse.insert(i,i+1) = complex<double>(1.0,0.0);
// 	}
// 	Lsparse.makeCompressed();

// 	// cout<< " Lsparse: \n" << Lsparse << endl;

	
// 	for(int i = 0; i < grid_x; i++){for(int j = 0; j < grid_y; j++){ mPsiCopy(i,j) = pPsi->at(0,i,j,0);}}
// 	for(int i = 0; i < grid_x; i++){for(int j = 0; j < grid_y; j++){ mPsi(i,j) = pPsi->at(0,i,j,0);}}

// 	// 	Eigen::MatrixXcd result;
// 	// result = Lsparse * mPsiCopy;

// 	// cout << "result: " << result << endl;

// // Eigen::MatrixXcd G(opt.grid[1],opt.grid[2]);
// // Eigen::MatrixXcd X(opt.grid[1],opt.grid[2]); 
// // Eigen::MatrixXcd Y(opt.grid[1],opt.grid[2]); 
// // Eigen::MatrixXcd L(opt.grid[1],opt.grid[2]);   

//  	for(int i = 0; i<opt.grid[1];i++)
// 	{
// 		Xsparse.insert(i,i) = complex<double>(x_axis[i],0.0);
// 	}
// 	Xsparse.makeCompressed();

// 	// cout << " X: \n" << X << endl;

// 	for(int j = 0; j<opt.grid[2];j++)
// 	{
// 		Ysparse.insert(j,j) = complex<double>(y_axis[j],0.0);
// 	}
// 	Ysparse.makeCompressed();

// 	// cout << " Y: \n"<< Y << endl;

// 	// for(int i = 1;i<opt.grid[1]-1;i++){
// 	// L(i,i-1) = complex<double>(1.0,0.0);
// 	// L(i,i) = complex<double>(-2.0,0.0);
// 	// L(i,i+1) = complex<double>(1.0,0.0);
// 	// }
// 	// L.makeCompressed();

// 	// cout<< " L: \n" << L << endl;
// 	Gsparse.resize(grid_x,grid_y);
// 	Gsparse.reserve(2 * (opt.grid[1]-2));
// 	for(int i = 1;i<opt.grid[1]-1;i++){
// 	Gsparse.insert(i,i-1) = complex<double>(-1.0,0.0);
// 	Gsparse.insert(i,i+1) = complex<double>(1.0,0.0);
// 	}
// 	Gsparse.makeCompressed();
// 	// cout << " G: \n"<< G << endl;

// 	Eigen::SparseMatrix<std::complex<double>,0,std::ptrdiff_t > LsparseT = Eigen::SparseMatrix<std::complex<double>,0,std::ptrdiff_t >(Lsparse.transpose());
// 	Eigen::SparseMatrix<std::complex<double>,0,std::ptrdiff_t > GsparseT = Eigen::SparseMatrix<std::complex<double>,0,std::ptrdiff_t >(Gsparse.transpose());

// 	cout << "All Matrix Sizes: " << endl
// 		<< "Lsparse: " << Lsparse.size() << endl
// 		<< "LsparseT: " << LsparseT.size() << endl
// 		<< "GsparseT: " << GsparseT.size() << endl
// 		<< "Gsparse: " << Gsparse.size() << endl
// 		<< "XSparse: " << Xsparse.size() << endl
// 		<< "YSparse: " << Ysparse.size() << endl
// 		<< "mPsiCopy: " << mPsiCopy.size() << endl;

// 	start = omp_get_wtime();

// 	cout << " " << opt.name << endl;

// for(int m = 1; m <= opt.n_it_RTE; m++)
// {
// 	// cout << "before function" << endl;
// 	// functionEigen_RTE(k0,mPsiCopy,t);
// 		// cout << "L " << L.size() << " mPsiCopy:  " << mPsiCopy.size() << " Ltrans:  " << L.transpose().size()<< endl;
// 	k0.noalias() += complex<double>(0.0,0.5) * ((Lsparse * mPsiCopy) / ( l_x(t) * l_x(t) ) + ( mPsiCopy * LsparseT ) / (l_y(t) * l_y(t) ) );
// 	// cout << "2" << endl;
// 	k0.noalias() += ( (l_x_dot(t) / l_x(t)) * /*Xsparse **/ ( Gsparse * mPsiCopy) );
// 	// cout << "3" << endl;
// 	k0.noalias() += ( (l_y_dot(t) / l_y(t)) * ( mPsiCopy * GsparseT ) /** Ysparse*/ );
// 	// cout << "4" << endl;
// 	k0.array() += - complex<double>(0.0, g ) * (mPsiCopy.adjoint().array() * mPsiCopy.array()) * mPsiCopy.array();
// 	// cout << "5" << endl;
// 	// cout << "after function" << endl;

// 	t += half * t_RTE;
// 	mPsiCopy = mPsi + half * t_RTE * k0;
// 	// functionEigen_RTE(k1,mPsiCopy,t);
// 	k1.noalias() += complex<double>(0.0,0.5) * ((Lsparse * mPsiCopy) / ( l_x(t) * l_x(t) ) + ( mPsiCopy * LsparseT ) / (l_y(t) * l_y(t) ) );
// 	k1.noalias() += ( (l_x_dot(t) / l_x(t)) * /*Xsparse **/ ( Gsparse * mPsiCopy) );
// 	k1.noalias() += ( (l_y_dot(t) / l_y(t)) * ( mPsiCopy * GsparseT ) /** Ysparse*/ );
// 	k1.array() += - complex<double>(0.0, g ) * (mPsiCopy.adjoint().array() * mPsiCopy.array()) * mPsiCopy.array();

// 	mPsiCopy = mPsi + half * t_RTE * k1;
// 	// functionEigen_RTE(k2,mPsiCopy,t);
// 	k2.noalias() += complex<double>(0.0,0.5) * ((Lsparse * mPsiCopy) / ( l_x(t) * l_x(t) ) + ( mPsiCopy * LsparseT ) / (l_y(t) * l_y(t) ) );
// 	k2.noalias() += ( (l_x_dot(t) / l_x(t)) * /*Xsparse **/ ( Gsparse * mPsiCopy) );
// 	k2.noalias() += ( (l_y_dot(t) / l_y(t)) * ( mPsiCopy * GsparseT ) /** Ysparse*/ );
// 	k2.array() += - complex<double>(0.0, g ) * (mPsiCopy.adjoint().array() * mPsiCopy.array()) * mPsiCopy.array();

// 	t += half * t_RTE;
// 	mPsiCopy = mPsi + t_RTE * k2;
// 	// functionEigen_RTE(k3,mPsiCopy,t);
// 	k3.noalias() += complex<double>(0.0,0.5) * ((Lsparse * mPsiCopy) / ( l_x(t) * l_x(t) ) + ( mPsiCopy * LsparseT ) / (l_y(t) * l_y(t) ) );
// 	k3.noalias() += ( (l_x_dot(t) / l_x(t)) * /*Xsparse **/ ( Gsparse * mPsiCopy) );
// 	k3.noalias() += ( (l_y_dot(t) / l_y(t)) * ( mPsiCopy * GsparseT ) /** Ysparse*/ );
// 	k3.array() += - complex<double>(0.0, g ) * (mPsiCopy.adjoint().array() * mPsiCopy.array()) * mPsiCopy.array();

// 	mPsi += (t_RTE/six) * ( k0 + two * k1 + two * k2 + k3);

// 	mPsiCopy = mPsi;

// 	// opt.min_x = x_expand(opt.grid[1]-1,opt);
//   	// opt.min_y = y_expand(opt.grid[2]-1,opt);

//   	cli_plot(mPsi,opt,"RTE",m,opt.n_it_RTE,start,plot);

// }
// 	// for(int i = 0; i < grid_x; i++){for(int j = 0; j < grid_y; j++){ pPsi->at(0,i,j,0) = mPsi(i,j) ;}}
// 	// cli_plot(opt,"RTE",100,opt.n_it_RTE,start,plot);

// cout << "\n";


	
// }




void RK4::functionEigen2_RTE(Options &opt, bool plot)
{
	int grid_x = opt.grid[1];
	int grid_y = opt.grid[2];
	double start;

	Eigen::MatrixXcd wavefct(grid_x,grid_y);
	Eigen::MatrixXcd wavefctcp(grid_x,grid_y);
	Eigen::MatrixXcd k0 = Eigen::MatrixXcd::Zero(grid_x,grid_y);
	Eigen::MatrixXcd k1 = Eigen::MatrixXcd::Zero(grid_x,grid_y);
	Eigen::MatrixXcd k2 = Eigen::MatrixXcd::Zero(grid_x,grid_y);
	Eigen::MatrixXcd k3 = Eigen::MatrixXcd::Zero(grid_x,grid_y);
	for(int i = 0; i < grid_x; i++){for(int j = 0; j < grid_y; j++){ wavefct(i,j) = pPsi->at(0,i,j,0);}}
	
	Eigen::VectorXcd X(opt.grid[1]), Y(opt.grid[2]);
	for(int i = 0;i<grid_x;i++){
		X(i) = complex<double>(x_axis[i],0.0);
	}
	for(int j = 0;j<grid_y;j++){
		Y(j) = complex<double>(y_axis[j],0.0);
	}

	start = omp_get_wtime();

	cout << " " << opt.name << endl;

	int t = 0;
	//start loop here
	for(int m = 1; m <= opt.n_it_RTE; m++){

		wavefctcp = wavefct;

		//boundary conditions -- Dirichlet

		wavefct.row(0) = Eigen::VectorXcd::Zero(grid_x);
		wavefct.row(grid_x-1) = Eigen::VectorXcd::Zero(grid_x);
		wavefct.col(0) = Eigen::VectorXcd::Zero(grid_y);
		wavefct.col(grid_y-1) = Eigen::VectorXcd::Zero(grid_y);

		//boundary conditions end


		// if(m%((opt.n_it_RTE)/100)==0)
		// {	
		// 	cout << "  " << m << " wavefct " << wavefct.block<1,4>(grid_x/2,grid_y/2) << endl;
		// 	cout << "  " << m << " norm " << wavefct.norm() << endl << endl;
		// }

		k0 = compute_the_eigen_k(wavefctcp,X,Y,t,grid_x,grid_y);
		wavefctcp = wavefct + half * t_RTE * k0;

		t += 1;
		k1 = compute_the_eigen_k(wavefctcp,X,Y,t,grid_x,grid_y);
		wavefctcp = wavefct + half * t_RTE * k1;

		k2 = compute_the_eigen_k(wavefctcp,X,Y,t,grid_x,grid_y);		
		wavefctcp = wavefct + t_RTE * k2;

		t += 1;
		k3 = compute_the_eigen_k(wavefctcp,X,Y,t,grid_x,grid_y);

		wavefct += (t_RTE/six) * ( k0 + two * k1 + two * k2 + k3);

		// if(m%((opt.n_it_RTE)/100)==0)
		// {
		// 	cout << "1 " << k0.block<1,4>(grid_x/2,grid_y/2) << endl
		// 		 << "2 " << k1.block<1,4>(grid_x/2,grid_y/2) << endl
		// 		 << "3 " << k2.block<1,4>(grid_x/2,grid_y/2) << endl
		// 		 << "4 " << k3.block<1,4>(grid_x/2,grid_y/2) << endl << endl
		// 		 << t_RTE/six << endl;
		// }

		cli_plot(wavefct,opt,"RTE",m,opt.n_it_RTE,start,plot);

	}

cout << "\n";

}

Eigen::MatrixXcd RK4::compute_the_eigen_k(Eigen::MatrixXcd &wavefctcp, Eigen::VectorXcd &X,Eigen::VectorXcd &Y,int &t,int &grid_x,int &grid_y)
	{
	Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> wavefctcpX = Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>::Zero(grid_x,grid_y);
	Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> wavefctcpY = Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::Zero(grid_x,grid_y);
	Eigen::MatrixXcd k = Eigen::MatrixXcd::Zero(grid_x,grid_y);


	//laplacian
	// wavefctcpX = complex<double>(0.0,1.0) * wavefctcp / (two * h_x * h_x /** lambda_x_squared[t]*/);
	// wavefctcpY = complex<double>(0.0,1.0) * wavefctcp / (two * h_y * h_y /** lambda_y_squared[t]*/);

	// #pragma omp parallel for
	for(int j = 1;j<grid_y-1;j++){
	for(int i = 1;i<grid_x-1;i++){
	wavefctcpX(i,j) = wavefctcp(i-1,j) - two * wavefctcp(i,j) + wavefctcp(i+1,j);
	wavefctcpY(i,j) = wavefctcp(i,j-1) - two * wavefctcp(i,j) + wavefctcp(i,j+1);
	}}
	k += i_unit * ( wavefctcpX / (two * h_x * h_x /** lambda_x_squared[t]*/) +  wavefctcpY / (two * h_y * h_y /** lambda_y_squared[t]*/)  );
	// k += i_unit * wavefctcpX / (two * h_x * h_x /** lambda_x_squared[t]*/);
	// #pragma omp parallel for
	// for(int i = 1;i<grid_x-1;i++){
	// for(int j = 1;j<grid_y-1;j++){
	// wavefctcpY(i,j) = wavefctcp(i,j-1) - two * wavefctcp(i,j) + wavefctcp(i,j+1);
	// }}
	// k += i_unit * wavefctcpY / (two * h_y * h_y /** lambda_y_squared[t]*/);
	//laplacian end

	// gradient

	// wavefctcpX = wavefctcp * lambda_dot_x[t] / (two * h_x );
	// wavefctcpY = wavefctcp * lambda_dot_y[t] / (two * h_y );

	// #pragma omp parallel for
	// for(int j = 0;j<grid_y;j++){ wavefctcpX.col(j).array() *= X.array(); }
	// #pragma omp parallel for
	// for(int i = 0;i<grid_x;i++){ wavefctcpY.row(i).array() *= Y.array(); }

	// #pragma omp parallel for
	// for(int j = 0;j<grid_y;j++){
	// for(int i = 1;i<grid_x-1;i++){
	// k(i,j) += wavefctcpX(i+1,j) - wavefctcpX(i-1,j);
	// }}
	// #pragma omp parallel for
	// for(int i = 0;i<grid_x;i++){
	// for(int j = 1;j<grid_y-1;j++){
	// k(i,j) += wavefctcpY(i,j+1) - wavefctcpY(i,j-1);
	// }}
	// // gradient end

	//interaction
	k.array() -= complex<double>(0.0,g) * ( wavefctcp.conjugate().array() * wavefctcp.array() ) * wavefctcp.array();

	//interaction end

	return k;
	}


// Propagation Wrapper Functions

// void RK4::rteToTime(Options &opt, bool plot)
// {
// 	double start;
// 	complex<double> t_RTE(opt.RTE_step,0); //Time-step size for RTE


// 	vector<ComplexGrid> k(4);	// kvectors for RK4
// 	for(int d=0;d<4;d++){
// 		k[d] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
// 		for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ k[d](0,i,j,0) = zero; }}
// 	}

// 	start = omp_get_wtime();

// 	cout << " " << opt.name << endl;

// 	for(int m = 1; m <= opt.n_it_RTE; m++)
// 	{
// 		computeK_RTE(pPsi,k,opt,t_RTE);

// 		TimeStepRK4(pPsi,k,opt,t_RTE);

//   		opt.min_x = x_expand(opt.grid[1]-1,opt);
//   		opt.min_y = y_expand(opt.grid[2]-1,opt);

//   		// cli_plot(opt,"RTE",m,opt.n_it_RTE,start,plot);

// 	}

// 	cout << "\n";
// }


