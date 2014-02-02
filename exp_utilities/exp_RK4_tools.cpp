#include <exp_RK4_tools.h>



using namespace std;

RK4::RK4()
{
    h_x =0;
    h_y=0;
    x_axis[0]=0;
    y_axis[0]=0;
    Integral=0;
    Integral_aux=0;
}

RK4::RK4(ComplexGrid* &c,Options &opt)
{
  h_x = (((2.*opt.min_x)/opt.grid[1]),0);
  h_y = (((2.*opt.min_y)/opt.grid[2]),0);
  x_axis[opt.grid[1]],y_axis[opt.grid[2]];
  for(int i=0;i<opt.grid[1];i++){x_axis[i]=-opt.min_x+i*real(h_x);} //Initialisation of the x-axis from -min_x to +min_x
  for(int j=0;j<opt.grid[2];j++){y_axis[j]=-opt.min_y+j*real(h_y);} //Initialisation of the y-axis from -min_y to +min_y
}

double RK4::gauss(double x,double y)
{return (exp(-x*x-y*y));} //A simple Gaussian

// complex<double> RK4::laplacian_x(complex<double> a, complex<double> b, complex<double> c)
// {return ((a-(two*b)+c)/(h_x*h_x));} //Central-difference x-laplacian approximation

// complex<double> RK4::laplacian_y(complex<double> a, complex<double> b, complex<double> c)
// {return ((a-(two*b)+c)/(h_y*h_y));} //Central-difference y-laplacian approximation

// complex<double> RK4::potential(double x,double y, Options &opt)
// {return (half*opt.omega_x*opt.omega_x*complex<double>(x,0)*complex<double>(x,0)+half*opt.omega_y*opt.omega_y*complex<double>(y,0)*complex<double>(y,0));} //SHO potential

/*complex<double> RK4::interaction(complex<double> a,Options &opt)
{return (opt.g*norm(a));} //Interaction term in the GPE Hamiltonian  */ 

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

complex<double> RK4::rescale(ComplexGrid* & pPsi,ComplexGrid* & pPsiCopy, Options &opt)
{	
	Integral=integral(pPsi,opt);
	opt.scale_factor=opt.N/Integral;
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
        for(int i=0;i<opt.grid[1];i++)
        {
	        for(int j=0;j<opt.grid[2];j++)
	        {
		        save_obdm(x_axis[i]*real(lambda_x(opt.t_abs,opt)),y_axis[j]*real(lambda_y(opt.t_abs,opt)),norm(pPsi->at(0,i,j,0)),phase_save(pPsi,i,j));
	        }
                blank_line();
        }
        blank_line();
}
/*
complex<double> RK4::function_ITP_BC(int i,int j) //Function used for the ITP fixed-derivative (Neumann) boundary conditions
{return -(potential(x_axis[i],y_axis[j])+interaction(psi_copy[i][j]))*psi_copy[i][j];}

complex<double> RK4::function_ITP(int i,int j) //Function used for the ITP Runge-Kutta evolution 
{return half*laplacian_x(psi_copy[i+1][j],psi_copy[i][j],psi_copy[i-1][j])+half*laplacian_y(psi_copy[i][j+1],psi_copy[i][j],psi_copy[i][j-1])-(potential(x_axis[i],y_axis[j])+interaction(psi_copy[i][j]))*psi_copy[i][j];} 

complex<double> RK4::function_RTE(int i, int j, complex<double> t) //Function used for the RTE Runge-Kutta evolution (expanding frame)
{return ((i_unit*laplacian_x(psi_copy[i+1][j],psi_copy[i][j],psi_copy[i-1][j]))/(two*lambda_x(t)*lambda_x(t)))+((i_unit*laplacian_y(psi_copy[i][j+1],psi_copy[i][j],psi_copy[i][j-1]))/(two*lambda_y(t)*lambda_y(t)))-(i_unit*interaction(psi_copy[i][j])*psi_copy[i][j])+((lambda_x_dot(t)*real(x_expand(i,t_abs))*grad_x(psi_copy[i+1][j],psi_copy[i-1][j]))/lambda_x(t))+((lambda_y_dot(t)*real(y_expand(j,t_abs))*grad_y(psi_copy[i][j+1],psi_copy[i][j-1]))/lambda_y(t)) ;}*/

complex<double> RK4::T(ComplexGrid* & pPsiCopy,int i, int j){return half*((pPsiCopy->at(0,i+1,j,0)-(two*pPsiCopy->at(0,i,j,0))+pPsiCopy->at(0,i-1,j,0))/(h_x*h_x))+half*((pPsiCopy->at(0,i,j+1,0)-(two*pPsiCopy->at(0,i,j,0))+pPsiCopy->at(0,i,j-1,0))/(h_x*h_x)); }
complex<double> RK4::V(ComplexGrid* & pPsiCopy,int i, int j,Options & opt){ return -(half*opt.omega_x*opt.omega_x*complex<double>(x_axis[i],0)*complex<double>(x_axis[i],0)+half*opt.omega_y*opt.omega_y*complex<double>(y_axis[j],0)*complex<double>(y_axis[j],0)+opt.g*norm(pPsiCopy->at(0,i,j,0)))*pPsiCopy->at(0,i,j,0);} 

void RK4::computeK(ComplexGrid* & pPsiCopy,ComplexGrid* & pPsi, ComplexGrid** k,Options & opt,complex<double> & t_ITP, int d){ 

    if (d == 1){
      pPsiCopy = pPsi;
    }
    else if ( d== 2){
      for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ pPsiCopy->at(0,i,j,0)=pPsi->at(0,i,j,0)+half*t_ITP*k[d-1]->at(0,i,j,0) ;}}    
    }
    else if (d == 3){
      for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ pPsiCopy->at(0,i,j,0)=pPsi->at(0,i,j,0)+half*t_ITP*k[d-1]->at(0,i,j,0) ;}}
    }
    else if (d == 4){
     for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ pPsiCopy->at(0,i,j,0)=pPsi->at(0,i,j,0)+t_ITP*k[d-1]->at(0,i,j,0) ;}}	
   }
   
   
      for(int i=1;i<opt.grid[1]-1;i++) // Compute all k[1]..k[4] for the RK4
      {
	 for(int j=1;j<opt.grid[2]-1;j++)
	 {
	   k[d]->at(0,i,j,0) = T(pPsiCopy,i,j)+V(pPsiCopy,i,j,opt);  
	 }	    
      }
      for(int i=0;i<opt.grid[1];i++) //Fixed derivative (Neumann) boundary conditions for the ITP
      { 
	k[d]->at(0,i,0,0)=V(pPsiCopy,i,0,opt);
	k[d]->at(0,i,opt.grid[1]-1,0)=V(pPsiCopy,i,opt.grid[1]-1,opt);
      }
      for(int j=0;j<opt.grid[2];j++)
      { 
	k[d]->at(0,0,j,0)=V(pPsiCopy,0,j,opt);
	k[d]->at(0,opt.grid[2]-1,j,0)=V(pPsiCopy,opt.grid[2]-1,j,opt);
      }
   
  
}
void RK4::ITP(ComplexGrid* & pPsi, Options &opt, complex<double> & t_ITP)
{
	ComplexGrid* pPsiCopy;
	ComplexGrid* k[4];
	pPsiCopy = new ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	k[4] = new ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	
	// Compute all the k's for the RK4
	
	for(int d=1; d<=4; d++){
	computeK(pPsiCopy,pPsi,k,opt,t_ITP,d);
	 }
	
	// Use RK4 to timestep Psi
	
	for(int i=0;i<opt.grid[1];i++){
	  for(int j=0;j<opt.grid[2];j++){
	    pPsi->at(0,i,j,0)+=(t_ITP)*(k[1]->at(0,i,j,0)+two*k[2]->at(0,i,j,0)+two*k[3]->at(0,i,j,0)+k[4]->at(0,i,j,0));
	  }
	}
}
/*
void RK4::RTE(const int & n_x, const int & n_y, complex<double> & t_RTE)
{
        //Fixed (Dirichlet) boundary conditions for the RTE
        for(l=0;l<n_x;l++){ psi[l][0]=zero; psi_copy[l][0]=zero; psi[l][n_y-1]=zero; psi_copy[l][n_y-1]=zero; }
        for(m=0;m<n_y;m++){ psi[0][m]=zero; psi_copy[0][m]=zero; psi[n_x-1][m]=zero; psi_copy[n_x-1][m]=zero; }
		//k1		
	for(i=1;i<n_x-1;i++){for(j=1;j<n_y-1;j++){ k1[i][j]=function_RTE(i,j,t_abs) ;}}
		//k2
	for(i=0;i<n_x;i++){for(j=0;j<n_y;j++){ psi_copy[i][j]=psi[i][j]+half*t_RTE*k1[i][j] ;}}		
	
	for(i=1;i<n_x-1;i++){for(j=1;j<n_y-1;j++){ k2[i][j]=function_RTE(i,j,t_abs+half*t_RTE) ;}}
		//k3		
	for(i=0;i<n_x;i++){for(j=0;j<n_y;j++){ psi_copy[i][j]=psi[i][j]+half*t_RTE*k2[i][j] ;}}
	
	for(i=1;i<n_x-1;i++){for(j=1;j<n_y-1;j++){ k3[i][j]=function_RTE(i,j,t_abs+half*t_RTE) ;}}
		//k4		
	for(i=0;i<n_x;i++){for(j=0;j<n_y;j++){ psi_copy[i][j]=psi[i][j]+t_RTE*k3[i][j] ;}}
	
	for(i=1;i<n_x-1;i++){for(j=1;j<n_y-1;j++){ k4[i][j]=function_RTE(i,j,t_abs+t_RTE) ;}}
		//Time-step psi		
	for(i=0;i<n_x;i++)
	{		
		for(j=0;j<n_y;j++)		
		{
			psi[i][j]+=(t_RTE/six)*(k1[i][j]+two*k2[i][j]+two*k3[i][j]+k4[i][j]);
			psi_copy[i][j]=psi[i][j];
		}
	}
	}*/

const double RK4::pi = M_PI; //acos(-1.0L);
const complex<double> RK4::zero=(0,0),RK4::half=(0.5,0),RK4::one=(1,0),RK4::two=(2,0),RK4::four=(4,0),RK4::six=(6,0),RK4::i_unit=(0,1);
