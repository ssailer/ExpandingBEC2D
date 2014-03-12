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

/*
RK4::RK4(Options &opt)
{
  	real(h_x) = 2.*opt.min_x/opt.grid[1]; 
  	real(h_y) = 2.*opt.min_y/opt.grid[2]; 

    x_axis.resize(opt.grid[1]);
  	y_axis.resize(opt.grid[2]);
  	for(int i=0;i<opt.grid[1];i++){x_axis[i]=-opt.min_x+i*real(h_x);} //Initialisation of the x-axis from -min_x to +min_x
  	for(int j=0;j<opt.grid[2];j++){y_axis[j]=-opt.min_y+j*real(h_y); 

    Integral=0;
    Integral_aux=0;

    pi = M_PI;
    zero=complex<double>(0,0),half=complex<double>(0.5,0),one=complex<double>(1,0),two=complex<double>(2,0),four=complex<double>(4,0),six=complex<double>(6,0),i_unit=complex<double>(0,1);

    ComplexGrid* pPsi = new ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

    for(int i=0;i<opt.grid[1];i++) //Initialise the wavefunction as gaussian by default with default_gridsize
 	{
 		for(int j=0;j<opt.grid[2];j++)
 		{
 						
			double xfactor;
 			xfactor = gauss(x_axis[i],y_axis[j]); 			
 			complex<double> factor (xfactor,0); 			
 			pPsi->at(0,i,j,0) = factor;
 			pPhase->at(0,i,j,0) = (0,0); //!!!! Check if this should be double, and not complexdouble !!!!
 			 
 		}	
 	}
    
}
*/

RK4::RK4(ComplexGrid* &c,Options &opt)
{	
	opt.threads = omp_get_max_threads();
	cout << "Max Number of Threads: " << opt.threads << endl;
	omp_set_num_threads(opt.threads);
  	pPsi = c;  
  	h_x.real(2.*opt.min_x/opt.grid[1]);
  	h_y.real(2.*opt.min_y/opt.grid[2]); 
  	x_axis.resize(opt.grid[1]);
  	y_axis.resize(opt.grid[2]);
  	for(int i=0;i<opt.grid[1];i++){x_axis[i]=-opt.min_x+i*real(h_x);}
  	for(int j=0;j<opt.grid[2];j++){y_axis[j]=-opt.min_y+j*real(h_y);}
  	pi = M_PI;
 	zero=complex<double>(0,0),half=complex<double>(0.5,0),one=complex<double>(1,0),two=complex<double>(2,0),four=complex<double>(4,0),six=complex<double>(6,0),i_unit=complex<double>(0,1);
}

RK4::~RK4(){};

// double RK4::gauss(double x,double y){return (exp(-x*x-y*y));}

complex<double> RK4::interaction(complex<double> a,Options &opt)
{return (opt.g*norm(a));} //Interaction term in the GPE Hamiltonian   

complex<double> RK4::grad_x(complex<double> a, complex<double> b)
{return ((a-b)/(two*h_x));} //Central-difference x-grad approximation

complex<double> RK4::grad_y(complex<double> a, complex<double> b)
{return ((a-b)/(two*h_y));} //Central-difference y-grad approximation

complex<double> RK4::lambda_x(Options &opt)
{return sqrt(one+opt.exp_factor*opt.omega_x*opt.omega_x*opt.t_abs*opt.t_abs);}

complex<double> RK4::lambda_x_dot(Options &opt)
{return (opt.exp_factor*opt.omega_x*opt.omega_x*opt.t_abs/sqrt(one+opt.exp_factor*opt.omega_x*opt.omega_x*opt.t_abs*opt.t_abs));}

complex<double> RK4::lambda_y(Options &opt)
{return sqrt(one+opt.exp_factor*opt.omega_y*opt.omega_y*opt.t_abs*opt.t_abs);}

complex<double> RK4::lambda_y_dot(Options &opt)
{return (opt.exp_factor*opt.omega_y*opt.omega_y*opt.t_abs/sqrt(one+opt.exp_factor*opt.omega_y*opt.omega_y*opt.t_abs*opt.t_abs));}

complex<double> RK4::x_expand(complex<double> a,Options &opt)
{return ((-complex<double>(opt.grid[1],0)/two+a)*h_x*lambda_x(opt));}

complex<double> RK4::y_expand(complex<double> a,Options &opt)
{return ((-complex<double>(opt.grid[2],0)/two+a)*h_y*lambda_y(opt));}

complex<double> RK4::integral(ComplexGrid* & pPsi,Options &opt)
{	
	Integral_aux=(0,0);	
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
	// cout << "Particle Number: " << opt.N << endl;
 	// cout  << "opt.scale_factor " << opt.scale_factor<< "   Integral : " <<  Integral << endl;
	
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
        if(arg(pPsi->at(0,a,b,0))<0){ return 2*pi+arg(pPsi->at(0,a,b,0)); } //arg uses the atan2 function, so the same applies
	else{ return arg(pPsi->at(0,a,b,0)); }
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
// 			cout << "x = " << x << "  |  y = " << y << endl;
		    save_obdm(x*real(lambda_x(opt)),y*real(lambda_y(opt)),norm(pPsi->at(0,i,j,0)),phase_save(pPsi,i,j));
// 			cout << "lambda: " << real(lambda_x(opt.t_abs,opt)) << " |  return of lambda_x  "  << sqrt(one+opt.exp_factor*opt.omega_x*opt.omega_x*opt.t_abs*opt.t_abs) << "  Psi = " << pPsi->at(0,i,j,0) << "  |  NORM: " << norm(pPsi->at(0,i,j,0)) <<  "      ";
// 			cout << x*real(lambda_x(opt.t_abs,opt)) << "  " << y*real(lambda_y(opt.t_abs,opt))  << "  " << norm(pPsi->at(0,i,j,0)) << "  " << phase_save(pPsi,i,j) << endl << endl;
	        }
                blank_line();
        }
        blank_line();
	closeDataFiles_obdm();
}


complex<double> RK4::T(ComplexGrid &PsiCopy,int i, int j)
{
	return half*((PsiCopy(0,i+1,j,0)-(two*PsiCopy(0,i,j,0))+PsiCopy(0,i-1,j,0))/(h_x*h_x))+half*((PsiCopy(0,i,j+1,0)-(two*PsiCopy(0,i,j,0))+PsiCopy(0,i,j-1,0))/(h_x*h_x)); 
}

complex<double> RK4::V(ComplexGrid & PsiCopy,int i, int j,Options & opt)
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

void RK4::Neumann(ComplexGrid &k,ComplexGrid &PsiCopy,Options &opt){
      //Fixed derivative (Neumann) boundary conditions for the ITP


    for(int i=0;i<opt.grid[1];i++) 
    { 
		k(0,i,0,0)=V(PsiCopy,i,0,opt);
		k(0,i,opt.grid[1]-1,0)=V(PsiCopy,i,opt.grid[1]-1,opt);
    }

    for(int j=0;j<opt.grid[2];j++)
    { 
		k(0,0,j,0)=V(PsiCopy,0,j,opt);
		k(0,opt.grid[2]-1,j,0)=V(PsiCopy,opt.grid[2]-1,j,opt);
    }
}


void RK4::computeK_ITP(ComplexGrid* &pPsi, vector<ComplexGrid> &k,Options &opt,complex<double> &t_ITP){ 
	
	ComplexGrid PsiCopy = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	PsiCopy = *pPsi;
	

	for(int d=1; d<=4; d++){

  // The k's have to be computed differently, this is decided by int d
    if (d == 2 || d == 3){
      	for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ PsiCopy(0,i,j,0)=pPsi->at(0,i,j,0)+half*t_ITP*k[d-2](0,i,j,0) ;}}
    }
    else if (d == 4){
    	for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ PsiCopy(0,i,j,0)=pPsi->at(0,i,j,0)+t_ITP*k[d-2](0,i,j,0) ;}}	
	}
   
    // Compute all k[1]..k[4] for the RK4
	#pragma omp parallel for 
    for(int i=1;i<opt.grid[1]-1;i++) {
		for(int j=1;j<opt.grid[2]-1;j++) {
	   		k[d-1](0,i,j,0) = T(PsiCopy,i,j)+V(PsiCopy,i,j,opt); 
	 	}
	}
	
    
    Neumann(k[d-1],PsiCopy,opt);

	}
  
}



void RK4::ITP(ComplexGrid* & pPsi, Options &opt)
{
	vector<ComplexGrid> k(4);
	complex<double> t_ITP(opt.ITP_step,0); //Timetep size for ITP (the equations already assume imaginary time so a real t_ITP should be used)
	
	for ( int i=0;i<4;i++){k[i] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);}
	
	// Compute all the k's for the RK4	

	computeK_ITP(pPsi,k,opt,t_ITP);
	
	// Use the k's of RK4 to timestep Psi	

	TimeStepRK4(pPsi,k,opt,t_ITP);
	
	// rescale Psi to conserve number of particles
	
	rescale(pPsi,opt); 
}

void RK4::itpToTime(Options &opt)
{
	int counter_ITP = 0;
	double start, end;

	start = omp_get_wtime();
	

	for(int k=0;k<opt.n_it_ITP;k++)
	{ 
		ITP(pPsi,opt);

		// if(k==vortex_start){add_vortex();} //Add vortex at ~80% of the ITP

		// if(k>0 && k%opt.n_save_ITP==0)
		// 	{
		// 		save_2D(pPsi,opt);
		// 	}

  		counter_ITP+=1;

		if(counter_ITP%(opt.n_it_ITP/100)==0)
			{	
				end = omp_get_wtime();
				cout << "    " << opt.name << " " << (counter_ITP/(opt.n_it_ITP/100)) << "%   " << end - start << "s        \r" << flush;
				opt.times = counter_ITP;
			}
	
	}
	cout << "\n";	
}


complex<double> RK4::function_RTE(ComplexGrid &PsiCopy,int i, int j,Options &opt) //Function used for the RTE Runge-Kutta evolution (expanding frame)
{
	return 
 ((i_unit*((PsiCopy(0,i+1,j,0)-(two*PsiCopy(0,i,j,0))+PsiCopy(0,i-1,j,0))/(h_x*h_x)))/(two*lambda_x(opt)*lambda_x(opt)))
+((i_unit*((PsiCopy(0,i,j+1,0)-(two*PsiCopy(0,i,j,0))+PsiCopy(0,i,j-1,0))/(h_y*h_y)))/(two*lambda_y(opt)*lambda_y(opt)))
-(i_unit*(complex<double>(opt.g,0)*norm((PsiCopy(0,i,j,0))))*PsiCopy(0,i,j,0))
+((lambda_x_dot(opt)*real(x_expand(i,opt))*((PsiCopy(0,i+1,j,0)-PsiCopy(0,i-1,j,0))/(two*h_x)))/lambda_x(opt))
+((lambda_y_dot(opt)*real(y_expand(j,opt))*((PsiCopy(0,i,j+1,0)-PsiCopy(0,i,j-1,0))/(two*h_y)))/lambda_y(opt));
}

void RK4::Dirichlet(ComplexGrid* &pPsi,Options &opt){

        //Fixed (Dirichlet) boundary conditions for the RTE
    for(int l=0;l<opt.grid[1];l++){ 
	  pPsi->at(0,l,0,0)=zero;
	  pPsi->at(0,l,opt.grid[2]-1,0)=zero;
	}
    for(int m=0;m<opt.grid[2];m++){
	  pPsi->at(0,0,m,0)=zero;
	  pPsi->at(0,opt.grid[1]-1,m,0)=zero;
	}
}

void RK4::computeK_RTE(ComplexGrid* &pPsi, vector<ComplexGrid> &k,Options &opt,complex<double> &t_RTE){

	ComplexGrid PsiCopy = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

	Dirichlet(pPsi,opt);
	PsiCopy = *pPsi;
	
      
    for(int d=1;d<=4;d++)
	{  
		if(d==2){
			opt.t_abs += half*t_RTE;
			for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ PsiCopy(0,i,j,0)=pPsi->at(0,i,j,0)+half*t_RTE*k[d-1](0,i,j,0) ;}}		
		}else if(d ==3){	
			for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ PsiCopy(0,i,j,0)=pPsi->at(0,i,j,0)+half*t_RTE*k[d-1](0,i,j,0) ;}}
		} else if (d==4){
			opt.t_abs += half*t_RTE;
			for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ PsiCopy(0,i,j,0)=pPsi->at(0,i,j,0)+t_RTE*k[d-1](0,i,j,0) ;}}
		}

		#pragma omp parallel for
		for(int i=1;i<opt.grid[1]-1;i++)
			{
			for(int j=1;j<opt.grid[2]-1;j++)
			{
				k[d-1](0,i,j,0)=function_RTE(PsiCopy,i,j,opt) ;
			}
		}
	}

}

void RK4::RTE(ComplexGrid* &pPsi,Options &opt)
{	
	complex<double> t_RTE(opt.RTE_step,0); //Time-step size for RTE
	vector<ComplexGrid> k(4);
	
	for(int i=0;i<4;i++){k[i] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);}
	
	computeK_RTE(pPsi,k,opt,t_RTE);

	TimeStepRK4(pPsi,k,opt,t_RTE);

	// rescale(pPsi,opt);

}

// Propagation Wrapper Functions



void RK4::rteToTime(Options &opt)
{
	int counter_RTE = 0;
	double start, end;

	start = omp_get_wtime();
	

	for(int k=0;k<opt.n_it_RTE;k++)
	{
		RTE(pPsi,opt);
		opt.name = "RTE" + std::to_string(k);
		plotdatatopng(pPsi,opt);


		// if(k>0 && k%opt.n_save_RTE==0)
		// 	{
		// 		save_2D(pPsi,opt);
		// 	}

  		counter_RTE+=1;

		if(counter_RTE%(opt.n_it_RTE/100)==0)
			{
				end = omp_get_wtime();
				cout << "    " << "RTE" << " " << (counter_RTE/(opt.n_it_RTE/100)) << "%   " << end - start << "s        \r" << flush;
				opt.times = counter_RTE;
			}
	}

	cout << "\n";
}	