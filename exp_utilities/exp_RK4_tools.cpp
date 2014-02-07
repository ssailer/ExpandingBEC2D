#include <exp_RK4_tools.h>
#include <2dexpan.h>


using namespace std;



RK4::RK4()
{
    h_x =0;
    h_y=0;
    x_axis;
    y_axis;
    Integral=0;
    Integral_aux=0;
    pi = M_PI; //acos(-1.0L);
    zero=complex<double>(0,0),half=complex<double>(0.5,0),one=complex<double>(1,0),two=complex<double>(2,0),four=complex<double>(4,0),six=complex<double>(6,0),i_unit=complex<double>(0,1);

    
}

RK4::RK4(ComplexGrid* &c,Options &opt)
{
  pPsi = c;  
  real(h_x) = 2.*opt.min_x/opt.grid[1]; 
  real(h_y) = 2.*opt.min_y/opt.grid[2]; 
  x_axis.resize(opt.grid[1]);
  y_axis.resize(opt.grid[2]);
  for(int i=0;i<opt.grid[1];i++){x_axis[i]=-opt.min_x+i*real(h_x); /*cout << "x_axis["<<i<<"] = "<< x_axis[i] << endl;*/} //Initialisation of the x-axis from -min_x to +min_x
  for(int j=0;j<opt.grid[2];j++){y_axis[j]=-opt.min_y+j*real(h_y); /*cout << "y_axis["<<j<<"] = "<< y_axis[j] << endl;*/} //Initialisation of the y-axis from -min_y to +min_y
  pi = M_PI; //acos(-1.0L);
 zero=complex<double>(0,0),half=complex<double>(0.5,0),one=complex<double>(1,0),two=complex<double>(2,0),four=complex<double>(4,0),six=complex<double>(6,0),i_unit=complex<double>(0,1);
}
RK4::~RK4(){};


complex<double> RK4::interaction(complex<double> a,Options &opt)
{return (opt.g*norm(a));} //Interaction term in the GPE Hamiltonian   

complex<double> RK4::grad_x(complex<double> a, complex<double> b)
{return ((a-b)/(two*h_x));} //Central-difference x-grad approximation

complex<double> RK4::grad_y(complex<double> a, complex<double> b)
{return ((a-b)/(two*h_y));} //Central-difference y-grad approximation

complex<double> RK4::lambda_x(complex<double> t, Options &opt)
{return sqrt(one+opt.exp_factor*opt.omega_x*opt.omega_x*t*t);}

complex<double> RK4::lambda_x_dot(complex<double> t, Options &opt)
{return (opt.exp_factor*opt.omega_x*opt.omega_x*t/sqrt(one+opt.exp_factor*opt.omega_x*opt.omega_x*t*t));}

complex<double> RK4::lambda_y(complex<double> t, Options &opt)
{return sqrt(one+opt.exp_factor*opt.omega_y*opt.omega_y*t*t);}

complex<double> RK4::lambda_y_dot(complex<double> t, Options &opt)
{return (opt.exp_factor*opt.omega_y*opt.omega_y*t/sqrt(one+opt.exp_factor*opt.omega_y*opt.omega_y*t*t));}

complex<double> RK4::x_expand(complex<double> a, complex<double> t, Options &opt)
{return ((-complex<double>(opt.grid[1],0)/two+a)*h_x*lambda_x(t,opt));}

complex<double> RK4::y_expand(complex<double> a, complex<double> t, Options &opt)
{return ((-complex<double>(opt.grid[2],0)/two+a)*h_y*lambda_y(t,opt));}

complex<double> RK4::integral(ComplexGrid* & pPsi,Options &opt)
{	
	Integral_aux=(0,0);	
	for(int i=0;i<opt.grid[1]-1;i++)
	{
		for(int j=0;j<opt.grid[2]-1;j++)
		{
			Integral_aux+=h_x*lambda_x(opt.t_abs,opt)*h_y*lambda_y(opt.t_abs,opt)*(norm(pPsi->at(0,i,j,0))+norm(pPsi->at(0,i+1,j,0))+norm(pPsi->at(0,i,j+1,0))+norm(pPsi->at(0,i+1,j+1,0)))/four;
			
		}
	}
	return Integral_aux;
}

void RK4::rescale(ComplexGrid* & pPsi,ComplexGrid* & pPsiCopy, Options &opt)
{	
	Integral=integral(pPsi,opt);
	opt.scale_factor=opt.N/Integral;
// 	cout  << "opt.scale_factor " << opt.scale_factor<< "   Integral : " <<  Integral << endl;
	
	for(int i=0;i<opt.grid[1];i++)
	{
		for(int j=0;j<opt.grid[2];j++)
		{
		  pPsi->at(0,i,j,0)*=sqrt(opt.scale_factor);
		  pPsiCopy->at(0,i,j,0)*=sqrt(opt.scale_factor);
		}
	}
}

double RK4::vortex(int a, int b, int x, int y) //Vortex with phase [0,2*pi)          
{
        if(atan2(b-y,a-x)<0){ return 2*pi+atan2(b-y,a-x); } //atan2 is defined from [-pi,pi) so it needs to be changed to [0,2*pi)
	else{ return atan2(b-y,a-x); }        
}

double RK4::phase_save(ComplexGrid* & pPsi,int a,int b) //Definition of phase 
{
        if(arg(pPsi->at(0,a,b,0))<0){ return 2*pi+arg(pPsi->at(0,a,b,0)); } //arg uses the atan2 function, so the same applies
	else{ return arg(pPsi->at(0,a,b,0)); }
}

// void RK4::add_vortex() // computes the phasefield by adding all vortices together and saves the final psi to use for timeevolution
// {
//    	for(i=0;i<n_x;i++)
// 	{
//      		for(j=0;j<n_y;j++)
//       		{
// 		        phase[i][j]=phase_save(i,j)+vortex(i,j,x_1,y_1)
// 			  /*Ring 1*/
// +vortex(i,j,x_2,y_2)+vortex(i,j,x_3,y_3)+vortex(i,j,x_4,y_4)+vortex(i,j,x_5,y_5)+vortex(i,j,x_6,y_6)+vortex(i,j,x_7,y_7)
// 			  /*Ring 2*/
// +vortex(i,j,x_8,y_8)+vortex(i,j,x_9,y_9)+vortex(i,j,x_10,y_10)+vortex(i,j,x_11,y_11)+vortex(i,j,x_12,y_12)+vortex(i,j,x_13,y_13)+vortex(i,j,x_14,y_14)+vortex(i,j,x_15,y_15)+vortex(i,j,x_16,y_16)+vortex(i,j,x_17,y_17)+vortex(i,j,x_18,y_18)+vortex(i,j,x_19,y_19)
// 			  /*Ring 3*/
// +vortex(i,j,x_20,y_20)+vortex(i,j,x_21,y_21)+vortex(i,j,x_22,y_22)+vortex(i,j,x_23,y_23)+vortex(i,j,x_24,y_24)+vortex(i,j,x_25,y_25)+vortex(i,j,x_26,y_26)+vortex(i,j,x_27,y_27)+vortex(i,j,x_28,y_28)+vortex(i,j,x_29,y_29)+vortex(i,j,x_30,y_30)+vortex(i,j,x_31,y_31)+vortex(i,j,x_32,y_32)+vortex(i,j,x_33,y_33)+vortex(i,j,x_34,y_34)+vortex(i,j,x_35,y_35)+vortex(i,j,x_36,y_36)+vortex(i,j,x_37,y_37) 
// 			  /*Ring 4*/
// 			  +vortex(i,j,x_38,y_38)+vortex(i,j,x_39,y_39)+vortex(i,j,x_40,y_40)+vortex(i,j,x_41,y_41)+vortex(i,j,x_42,y_42)+vortex(i,j,x_43,y_43)+vortex(i,j,x_44,y_44)+vortex(i,j,x_45,y_45)+vortex(i,j,x_46,y_46)+vortex(i,j,x_47,y_47)+vortex(i,j,x_48,y_48)+vortex(i,j,x_49,y_49)+vortex(i,j,x_50,y_50)+vortex(i,j,x_51,y_51)+vortex(i,j,x_52,y_52)+vortex(i,j,x_53,y_53)+vortex(i,j,x_54,y_54)+vortex(i,j,x_55,y_55)+vortex(i,j,x_56,y_56)+vortex(i,j,x_57,y_57)+vortex(i,j,x_58,y_58)+vortex(i,j,x_59,y_59)+vortex(i,j,x_60,y_60)+vortex(i,j,x_61,y_61)
// ;
// 			psi[i][j]=polar(abs(psi_copy[i][j]),phase[i][j]); // compute psi by using the initial psi^2 and adding the phase 
// 			psi_copy[i][j]=psi[i][j];
//       		}
//    	}
// }

void RK4::save_2D(ComplexGrid* & pPsi,Options &opt) //Function to save the data to a file (with compatible blocks for gnuplot)
{	
	openDataFiles_obdm(1,opt.name); //Open file with name (1,name)
	
        for(int i=0;i<opt.grid[1];i++)
        {
	        for(int j=0;j<opt.grid[2];j++)
	        {
			double x;
			double y;
			x = x_axis[i];
			y = y_axis[j];
// 			cout << "x = " << x << "  |  y = " << y << endl;
		        save_obdm(x*real(lambda_x(opt.t_abs,opt)),y*real(lambda_y(opt.t_abs,opt)),norm(pPsi->at(0,i,j,0)),phase_save(pPsi,i,j));
// 			cout << "lambda: " << real(lambda_x(opt.t_abs,opt)) << " |  return of lambda_x  "  << sqrt(one+opt.exp_factor*opt.omega_x*opt.omega_x*opt.t_abs*opt.t_abs) << "  Psi = " << pPsi->at(0,i,j,0) << "  |  NORM: " << norm(pPsi->at(0,i,j,0)) <<  "      ";
// 			cout << x*real(lambda_x(opt.t_abs,opt)) << "  " << y*real(lambda_y(opt.t_abs,opt))  << "  " << norm(pPsi->at(0,i,j,0)) << "  " << phase_save(pPsi,i,j) << endl << endl;
	        }
                blank_line();
        }
        blank_line();
	closeDataFiles_obdm();
}


complex<double> RK4::T(ComplexGrid* & pPsiCopy,int i, int j){return half*((pPsiCopy->at(0,i+1,j,0)-(two*pPsiCopy->at(0,i,j,0))+pPsiCopy->at(0,i-1,j,0))/(h_x*h_x))+half*((pPsiCopy->at(0,i,j+1,0)-(two*pPsiCopy->at(0,i,j,0))+pPsiCopy->at(0,i,j-1,0))/(h_x*h_x)); }
complex<double> RK4::V(ComplexGrid* & pPsiCopy,int i, int j,Options & opt){ 
  complex<double> xvalue;
  complex<double> yvalue;
  xvalue = (x_axis[i],0);
  yvalue = (y_axis[j],0);
  
  return -(half*opt.omega_x*opt.omega_x*xvalue*xvalue+half*opt.omega_y*opt.omega_y*yvalue*yvalue+opt.g*norm(pPsiCopy->at(0,i,j,0)))*pPsiCopy->at(0,i,j,0);} 

void RK4::computeK(ComplexGrid* & pPsiCopy,ComplexGrid* & pPsi, vector<ComplexGrid*> & k,Options & opt,complex<double> & t_ITP, int d){ 
	
  // The k's have to be computed differently, this is decided by int d
    if (d == 1){
      
      pPsiCopy = pPsi;
    }
    else if ( d== 2){
      for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ pPsiCopy->at(0,i,j,0)=pPsi->at(0,i,j,0)+half*t_ITP*k[d-2]->at(0,i,j,0) ;}}    
    }
    else if (d == 3){
      for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ pPsiCopy->at(0,i,j,0)=pPsi->at(0,i,j,0)+half*t_ITP*k[d-2]->at(0,i,j,0) ;}}
    }
    else if (d == 4){
     for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ pPsiCopy->at(0,i,j,0)=pPsi->at(0,i,j,0)+t_ITP*k[d-2]->at(0,i,j,0) ;}}	
   }
   
      // Compute all k[1]..k[4] for the RK4
      for(int i=1;i<opt.grid[1]-1;i++) 
      {
	 for(int j=1;j<opt.grid[2]-1;j++)
	 {
// 	   cout << "PsiCopy = " << pPsiCopy->at(0,i,j,0) << endl;
	   k[d-1]->at(0,i,j,0) = T(pPsiCopy,i,j)+V(pPsiCopy,i,j,opt); 
// 	   cout << "k[" << d-1 << "] = " << k[d-1]->at(0,i,j,0) << endl;
	 }	    
      }
      
      //Fixed derivative (Neumann) boundary conditions for the ITP
      for(int i=0;i<opt.grid[1];i++) 
      { 
	k[d-1]->at(0,i,0,0)=V(pPsiCopy,i,0,opt);
	k[d-1]->at(0,i,opt.grid[1]-1,0)=V(pPsiCopy,i,opt.grid[1]-1,opt);
      }
      for(int j=0;j<opt.grid[2];j++)
      { 
	k[d-1]->at(0,0,j,0)=V(pPsiCopy,0,j,opt);
	k[d-1]->at(0,opt.grid[2]-1,j,0)=V(pPsiCopy,opt.grid[2]-1,j,opt);
      }
   
  
}


void RK4::ITP(ComplexGrid* & pPsi, Options &opt)
{
	ComplexGrid* pPsiCopy;
	vector<ComplexGrid*> k(4);
	complex<double> t_ITP(opt.ITP_step,0); //Timetep size for ITP (the equations already assume imaginary time so a real t_ITP should be used)
	pPsiCopy = new ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	for ( int i=0;i<4;i++){k[i] = new ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);}
	
	// Compute all the k's for the RK4
	
	for(int d=1; d<=4; d++){
	computeK(pPsiCopy,pPsi,k,opt,t_ITP,d);
	}
	
	// Use RK4 to timestep Psi
	
	for(int i=0;i<opt.grid[1];i++){
	  for(int j=0;j<opt.grid[2];j++){
	    pPsi->at(0,i,j,0)+=(t_ITP)*(k[0]->at(0,i,j,0)+two*k[1]->at(0,i,j,0)+two*k[2]->at(0,i,j,0)+k[3]->at(0,i,j,0));
// 	    cout << "final Psi of computeK: " << pPsi->at(0,i,j,0) << endl;
	   }
	}
	
	
	rescale(pPsi,pPsiCopy,opt);
	
}


complex<double> RK4::function_RTE(ComplexGrid* & pPsiCopy,int i, int j, complex<double> t,Options &opt) //Function used for the RTE Runge-Kutta evolution (expanding frame)
{return ((i_unit*((pPsiCopy->at(0,i+1,j,0)-(two*pPsiCopy->at(0,i,j,0))+pPsiCopy->at(0,i-1,j,0))/(h_x*h_x)))/(two*lambda_x(t,opt)*lambda_x(t,opt)))+((i_unit*((pPsiCopy->at(0,i,j+1,0)-(two*pPsiCopy->at(0,i,j,0))+pPsiCopy->at(0,i,j-1,0))/(h_x*h_x)))/(two*lambda_y(t,opt)*lambda_y(t,opt)))
-(i_unit*(opt.g*norm((pPsiCopy->at(0,i,j,0))))*pPsiCopy->at(0,i,j,0))+((lambda_x_dot(t,opt)*real(x_expand(i,opt.t_abs,opt))*((pPsiCopy->at(0,i+1,j,0)-pPsiCopy->at(0,i-1,j,0))/(two*h_x)))/lambda_x(t,opt))+((lambda_y_dot(t,opt)*real(y_expand(j,opt.t_abs,opt))*((pPsiCopy->at(0,i,j+1,0)-pPsiCopy->at(0,i,j-1,0))/(two*h_y)))/lambda_y(t,opt)) ;}

void RK4::RTE(ComplexGrid* & pPsi,Options &opt)
{	
	ComplexGrid* pPsiCopy;
	complex<double> t_RTE(opt.RTE_step,0); //Time-step size for RTE
	pPsiCopy = new ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	vector<ComplexGrid*> k(4);
	
	for(int i=0;i<4;i++){
	  k[i] = new ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	  
	}
	
        //Fixed (Dirichlet) boundary conditions for the RTE
        for(int l=0;l<opt.grid[1];l++){ 
	  pPsi->at(0,l,0,0)=zero;
	  pPsiCopy->at(0,l,0,0)=zero;
	  pPsi->at(0,l,opt.grid[2]-1,0)=zero;
	  pPsiCopy->at(0,l,opt.grid[2]-1,0)=zero;
	  
	}
        for(int m=0;m<opt.grid[2];m++){
	  pPsi->at(0,0,m,0)=zero;
	  pPsiCopy->at(0,0,m,0)=zero;
	  pPsi->at(0,opt.grid[1]-1,m,0)=zero;
	  pPsiCopy->at(0,opt.grid[1]-1,m,0)=zero;
	  
	}
        
		//k1		
	for(int i=1;i<opt.grid[1]-1;i++){for(int j=1;j<opt.grid[2]-1;j++){ k[0]->at(0,i,j,0)=function_RTE(pPsiCopy,i,j,opt.t_abs,opt) ;}}
		//k2
	for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ pPsiCopy->at(0,i,j,0)=pPsi->at(0,i,j,0)+half*t_RTE*k[1]->at(0,i,j,0) ;}}		
	
	for(int i=1;i<opt.grid[1]-1;i++){for(int j=1;j<opt.grid[2]-1;j++){ k[1]->at(0,i,j,0)=function_RTE(pPsiCopy,i,j,opt.t_abs+half*t_RTE,opt) ;}}
		//k3		
	for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ pPsiCopy->at(0,i,j,0)=pPsi->at(0,i,j,0)+half*t_RTE*k[2]->at(0,i,j,0) ;}}
	
	for(int i=1;i<opt.grid[1]-1;i++){for(int j=1;j<opt.grid[2]-1;j++){ k[2]->at(0,i,j,0)=function_RTE(pPsiCopy,i,j,opt.t_abs+half*t_RTE,opt) ;}}
		//k4		
	for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ pPsiCopy->at(0,i,j,0)=pPsi->at(0,i,j,0)+t_RTE*k[3]->at(0,i,j,0) ;}}
	
	for(int i=1;i<opt.grid[1]-1;i++){for(int j=1;j<opt.grid[2]-1;j++){ k[3]->at(0,i,j,0)=function_RTE(pPsiCopy,i,j,opt.t_abs+t_RTE,opt) ;}}
	
	
		//Time-step psi		
	for(int i=0;i<opt.grid[1];i++)
	{		
		for(int j=0;j<opt.grid[2];j++)		
		{
			pPsi->at(0,i,j,0)+=(t_RTE/six)*(k[0]->at(0,i,j,0)+two*k[1]->at(0,i,j,0)+two*k[2]->at(0,i,j,0)+k[3]->at(0,i,j,0));
// 			pPsiCopy->at(0,i,j,0)=pPsi->at(0,i,j,0);
		}
	}
	delete pPsiCopy;
	for(int i=0;i<4;i++){
	  delete k[i];
	  
	}
	
	}


