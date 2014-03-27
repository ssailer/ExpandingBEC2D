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
 	zero=complex<double>(0,0),half=complex<double>(0.5,0),one=complex<double>(1,0),two=complex<double>(2,0),four=complex<double>(4,0),six=complex<double>(6,0),i_unit=complex<double>(0,1);
}

RK4::~RK4(){};

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


void RK4::cli_plot(Options &opt, string name,int counter_state, int counter_max, double start,bool plot)
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
					plotdatatopng(pPsi,opt);
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
	int counter_ITP = 0;
	double start, end;

	start = omp_get_wtime();
	
	cout << " " << opt.name << endl;
	string tmp = opt.name;
	for(int k=0;k<opt.n_it_ITP;k++)
	{ 
		ITP(pPsi,opt);		

  		counter_ITP += 1;

  		cli_plot(opt,tmp,counter_ITP,opt.n_it_ITP,start,plot);
  			
	}
	cout << "\n";	
}

complex<double> RK4::function_RTE(ComplexGrid &wavefct,int i, int j, Options &opt)
{
	complex<double> tmp;

	tmp = rte_kinetic(wavefct,i,j,opt);

	if(opt.g != 0.0)
	{tmp -= (rte_interaction(wavefct,i,j,opt) * wavefct(0,i,j,0));}

	if(opt.exp_factor != 0.0)
	{tmp += rte_expandingframe(wavefct,i,j,opt);}

	if(opt.startgrid[2] == true)
	{tmp -= rte_potential(i,j,opt) * wavefct(0,i,j,0);}

	return tmp;
}



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

void RK4::DirichletK(ComplexGrid &pPsi,Options &opt){

        //Fixed (Dirichlet) boundary conditions for the RTE
    for(int l=0;l<opt.grid[1];l++){ 
	  pPsi.at(0,l,0,0)=zero;
	  pPsi.at(0,l,opt.grid[2]-1,0)=zero;
	}
    for(int m=0;m<opt.grid[2];m++){
	  pPsi.at(0,0,m,0)=zero;
	  pPsi.at(0,opt.grid[1]-1,m,0)=zero;
	}
}

void RK4::NeumannRTE(ComplexGrid &k,ComplexGrid &wavefct,Options &opt){
      //Fixed derivative (Neumann) boundary conditions for the RTE


    for(int i=0;i<opt.grid[1];i++) 
    { 
		k(0,i,0,0)  = - rte_interaction(wavefct,i,0,opt) * wavefct(0,i,0,0);
		k(0,i,opt.grid[1]-1,0) = - rte_interaction(wavefct,opt.grid[1]-1,0,opt) * wavefct(0,opt.grid[1]-1,0,0);
    }

    for(int j=0;j<opt.grid[2];j++)
    { 
		k(0,0,j,0)= - rte_interaction(wavefct,0,j,opt) * wavefct(0,0,j,0);
		k(0,opt.grid[2]-1,j,0)= - rte_interaction(wavefct,0,opt.grid[2]-1,opt) * wavefct(0,0,opt.grid[2]-1,0);
    }
}

void RK4::computeK_RTE(ComplexGrid* &pPsi, vector<ComplexGrid> &k,Options &opt,complex<double> &t_RTE){

	ComplexGrid PsiCopy(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	for(int i = 0; i < opt.grid[1]; i++){for(int j = 0; j < opt.grid[2]; j++){ PsiCopy(0,i,j,0) = pPsi->at(0,i,j,0);}}

	Dirichlet(pPsi,opt);
	// PsiCopy = *pPsi;	
      
    for(int d=0;d<4;d++)
	{  
		// NeumannRTE(k[d-1],PsiCopy,opt);

		switch ( d ){
			case 1:
				opt.t_abs += half*t_RTE;
				#pragma omp parallel for
				for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ PsiCopy(0,i,j,0)=pPsi->at(0,i,j,0)+half*t_RTE*k[d-1](0,i,j,0) ;}}				
				break;
			case 2:
				#pragma omp parallel for	
				for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ PsiCopy(0,i,j,0)=pPsi->at(0,i,j,0)+half*t_RTE*k[d-1](0,i,j,0) ;}}
				break;
			case 3:
				opt.t_abs += half*t_RTE;
				#pragma omp parallel for
				for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ PsiCopy(0,i,j,0)=pPsi->at(0,i,j,0)+t_RTE*k[d-1](0,i,j,0) ;}}
				break;
		}

		#pragma omp parallel for
			for(int i=1;i<opt.grid[1]-1;i++){for(int j=1;j<opt.grid[2]-1;j++){ k[d](0,i,j,0)=function_RTE(PsiCopy,i,j,opt);}}
		// DirichletK(k[d],opt);
	}
}

void RK4::computeK_RTE_v2(ComplexGrid* &pPsi, vector<ComplexGrid> &k,Options &opt,complex<double> &t_RTE){

	ComplexGrid PsiCopy(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	for(int i = 0; i < opt.grid[1]; i++){for(int j = 0; j < opt.grid[2]; j++){ PsiCopy(0,i,j,0) = pPsi->at(0,i,j,0);}}

	Dirichlet(pPsi,opt);
	// PsiCopy = *pPsi;	
      
    for(int d=0;d<4;d++)
	{  
		// NeumannRTE(k[d-1],PsiCopy,opt);

		switch ( d ){
			case 1:
				opt.t_abs += half*t_RTE;
				#pragma omp parallel for
				for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ PsiCopy(0,i,j,0)=pPsi->at(0,i,j,0)+half*t_RTE*k[d-1](0,i,j,0) ;}}				
				break;
			case 2:
				#pragma omp parallel for	
				for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ PsiCopy(0,i,j,0)=pPsi->at(0,i,j,0)+half*t_RTE*k[d-1](0,i,j,0) ;}}
				break;
			case 3:
				opt.t_abs += half*t_RTE;
				#pragma omp parallel for
				for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ PsiCopy(0,i,j,0)=pPsi->at(0,i,j,0)+t_RTE*k[d-1](0,i,j,0) ;}}
				break;
		}

		#pragma omp parallel for
			for(int i=1;i<opt.grid[1]-1;i++){for(int j=1;j<opt.grid[2]-1;j++){ k[d](0,i,j,0)=function_RTE(PsiCopy,i,j,opt);}}
		// DirichletK(k[d],opt);
	}
}

void RK4::RTE(ComplexGrid* &pPsi,Options &opt)
{	
	complex<double> t_RTE(opt.RTE_step,0); //Time-step size for RTE
	vector<ComplexGrid> k(4);
	
	for(int d=0;d<4;d++){k[d] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

		for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ k[d](0,i,j,0) = zero; }}

	}
	
	computeK_RTE(pPsi,k,opt,t_RTE);

	TimeStepRK4(pPsi,k,opt,t_RTE);

}

// Propagation Wrapper Functions

void RK4::rteToTime(Options &opt, bool plot)
{
	int counter_RTE = 0;
	double start;

	start = omp_get_wtime();

	cout << " " << opt.name << endl;

	for(int k=0;k<opt.n_it_RTE;k++)
	{
		RTE(pPsi,opt);

  		counter_RTE+=1;
  		opt.min_x = x_expand(opt.grid[1]-1,opt);
  		opt.min_y = y_expand(opt.grid[2]-1,opt);

  		cli_plot(opt,"RTE",counter_RTE,opt.n_it_RTE,start,plot);

	}

	cout << "\n";
}


