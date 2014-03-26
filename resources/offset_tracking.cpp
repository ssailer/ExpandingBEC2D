/**************************************************************************
Title: Simulating the Expansion of Turbulent Bose-Einstein Condensates (2D) 
Author: Bartholomew Andrews
Last Update: 03/07/13
Website: www.bartholomewandrews.com
**************************************************************************/

#include <iostream>
#include <cmath>
#include <complex>
#include "2dexpan.h"

using namespace std;

//*****Parameter Initialisation*****

const double ITP_step=0.000001; //Time-step for the ITP (0.000001)
const double RTE_step=0.00001; //Time-step for the RTE (0.00001)
const int n_x=200,n_y=200; //Lattice size (200x200)  
const complex<double> N(1000,0); //Particle number (1000)
const complex<double> g(1,0); //Interaction constant (1)
const complex<double> omega_x(100,0),omega_y(150,0); //Trap frequency (100,150)
const double min_x=2,min_y=2; //Symmetric axis ranges start from -min_x and -min_y (2,2)

const int vortex_start=8000; //ITP iterations before the phase disturbances are added (8000<n_it_ITP) 
//const int n_save_ITP=1000; //Save ITP after every n_save_ITP iterations (initial state is auto saved)
const int n_it_ITP=10000; //Number of iterations for ITP (10000)
const int n_save_RTE=500; //Save RTE after every n_save_RTE iterations - intial state is auto saved (500)
const int n_it_RTE=100001; //Number of iterations for RTE (100001)
const int name=3; //Name of output (must be an integer)

//*****Varible Declarations*****

int i,j,k,l,m; //Declare the loop variables
int counter_ITP=0, counter_RTE=0; //Initialise the variables for the percentage loading

const complex<double> zero(0,0),half(0.5,0),one(1,0),two(2,0),four(4,0),six(6,0),i_unit(0,1); //Some useful numbers (to shorten code)
const long double pi=acos(-1.0L); //Precise definition of pi
		
complex<double> h_x(((2.*min_x)/n_x),0),h_y(((2.*min_y)/n_y),0); //Lattice constants (chosen for symmetric initial axis ranges)
complex<double> psi[n_x][n_y],psi_copy[n_x][n_y]; //Wavefunction (and copy for the RK4 algorithm)
complex<double> k1[n_x][n_y],k2[n_x][n_y],k3[n_x][n_y],k4[n_x][n_y]; //Runge-Kutta coefficients
double x_axis[n_x],y_axis[n_y];
double phase[n_x][n_y]; //Sum of all phases (used in the add_vortex() function) 

complex<double> scale_factor(0,0); //Scale factor
complex<double> Integral(0,0);
complex<double> Integral_aux(0,0); //Result of the integral function

//*****Vortex Initial Coordinates*****

//Spacing (16% across and 4% up relative to the lattice grid size)
const int across=16,up=4;
int half_across=across/2;

//Central Vortex (at the origin)
int x_1=n_x/2;
int y_1=n_y/2; 

//Ring 1 
int x_2=(100-across)*n_x/200,x_3=(100+across)*n_x/200,x_4=(100+half_across)*n_x/200,x_5=(100-half_across)*n_x/200,x_6=(100-half_across)*n_x/200,x_7=(100+half_across)*n_x/200;
int y_2=n_y/2,y_3=n_y/2,y_4=(100+up)*n_y/200,y_5=(100+up)*n_y/200,y_6=(100-up)*n_y/200,y_7=(100-up)*n_y/200; 

//Ring 2
int x_8=(100-2*across)*n_x/200,x_9=(100-across-half_across)*n_x/200,x_10=(100-across)*n_x/200,x_11=n_x/2,x_12=(100+across)*n_x/200,x_13=(100+across+half_across)*n_x/200,x_14=(100+2*across)*n_x/200,x_15=(100+across+half_across)*n_x/200,x_16=(100+across)*n_x/200,x_17=n_x/2,x_18=(100-across)*n_x/200,x_19=(100-across-half_across)*n_x/200;
int y_8=n_y/2,y_9=(100+up)*n_y/200,y_10=(100+2*up)*n_y/200,y_11=(100+2*up)*n_y/200,y_12=(100+2*up)*n_x/200,y_13=(100+up)*n_y/200,y_14=n_y/2,y_15=(100-up)*n_y/200,y_16=(100-2*up)*n_y/200,y_17=(100-2*up)*n_y/200,y_18=(100-2*up)*n_y/200,y_19=(100-up)*n_y/200;

//*****Function Definitions*****

double gauss(double x,double y){return (exp(-x*x-y*y));} //A simple Gaussian

complex<double> laplacian_x(complex<double> a, complex<double> b, complex<double> c)
{return ((a-(two*b)+c)/(h_x*h_x));} //Central-difference x-laplacian approximation

complex<double> laplacian_y(complex<double> a, complex<double> b, complex<double> c)
{return ((a-(two*b)+c)/(h_y*h_y));} //Central-difference y-laplacian approximation

complex<double> potential(double x,double y)
{return (half*omega_x*omega_x*complex<double>(x,0)*complex<double>(x,0)+half*omega_y*omega_y*complex<double>(y,0)*complex<double>(y,0));} //SHO potential

complex<double> interaction(complex<double> a)
{return (g*norm(a));} //Interaction term in the GPE Hamiltonian   

complex<double> grad_x(complex<double> a, complex<double> b)
{return ((a-b)/(two*h_x));} //Central-difference x-grad approximation

complex<double> grad_y(complex<double> a, complex<double> b)
{return ((a-b)/(two*h_y));} //Central-difference y-grad approximation

complex<double> function_ITP_BC(int i,int j) //Function used for the ITP fixed-derivative (Neumann) boundary conditions
{return -(potential(x_axis[i],y_axis[j])+interaction(psi_copy[i][j]))*psi_copy[i][j];}

complex<double> function_ITP(int i,int j) //Function used for the ITP Runge-Kutta evolution 
{return half*laplacian_x(psi_copy[i+1][j],psi_copy[i][j],psi_copy[i-1][j])+half*laplacian_y(psi_copy[i][j+1],psi_copy[i][j],psi_copy[i][j-1])-(potential(x_axis[i],y_axis[j])+interaction(psi_copy[i][j]))*psi_copy[i][j];} 

complex<double> function_RTE(int i, int j) //Function used for the RTE Runge-Kutta evolution
{return i_unit*half*laplacian_x(psi_copy[i+1][j],psi_copy[i][j],psi_copy[i-1][j])+i_unit*half*laplacian_y(psi_copy[i][j+1],psi_copy[i][j],psi_copy[i][j-1])-i_unit*interaction(psi_copy[i][j])*psi_copy[i][j]-i_unit*potential(x_axis[i],y_axis[j])*psi_copy[i][j];}

complex<double> integral()
{
	Integral_aux=0;	
	for(i=0;i<n_x-1;i++)
	{
		for(j=0;j<n_y-1;j++)
		{
			Integral_aux+=h_x*h_y*(norm(psi[i][j])+norm(psi[i+1][j])+norm(psi[i][j+1])+norm(psi[i+1][j+1]))/four;
		}
	}
	return Integral_aux;
}

complex<double> rescale()
{	
	Integral=integral();
	scale_factor=N/Integral;
	for(i=0;i<n_x;i++)
	{
		for(j=0;j<n_y;j++)
		{
		        psi[i][j]*=sqrt(scale_factor);
		        psi_copy[i][j]*=sqrt(scale_factor);
		}
	}
}

double vortex(int a, int b, int x, int y) //Vortex with phase [0,2*pi)          
{
        if(atan2(b-y,a-x)<0){ return 2*pi+atan2(b-y,a-x); } 
	else{  return atan2(b-y,a-x); }        
}

double phase_save(int a,int b) //Definition of phase 
{
        if(arg(psi[a][b])<0){ return 2*pi+arg(psi[a][b]); }
	else{ return arg(psi[a][b]); }
}

void add_vortex()
{
   	for(i=0;i<n_x;i++)
	{
     		for(j=0;j<n_y;j++)
      		{
		        phase[i][j]=phase_save(i,j)+vortex(i,j,x_2,y_2);
			psi[i][j]=polar(abs(psi_copy[i][j]),phase[i][j]);
			psi_copy[i][j]=psi[i][j];
      		}
   	}
}

void track()
{ 
        for(i=5;i<(n_x-5);i++)
        {
	        for(j=5;j<(n_y-5);j++)
	        {
		        if(norm(psi[i][j])<300 && (/*norm(psi[i][j+1])>300 &&*/ norm(psi[i][j-1])>300 && norm(psi[i-1][j-1])>300 && norm(psi[i-1][j])>300 && norm(psi[i-1][j+1])>300 && norm(psi[i+1][j-1])>300 && norm(psi[i+1][j])>300 && norm(psi[i+1][j+1])>300)){save_track(x_axis[i],y_axis[j],k*RTE_step,norm(psi[i][j]));break;}

		        else if(norm(psi[i][j])<300 && (norm(psi[i][j+1])>300 /*&& norm(psi[i][j-1])>300*/ && norm(psi[i-1][j-1])>300 && norm(psi[i-1][j])>300 && norm(psi[i-1][j+1])>300 && norm(psi[i+1][j-1])>300 && norm(psi[i+1][j])>300 && norm(psi[i+1][j+1])>300)){save_track(x_axis[i],y_axis[j],k*RTE_step,norm(psi[i][j]));break;}

		        else if(norm(psi[i][j])<300 && (norm(psi[i][j+1])>300 && norm(psi[i][j-1])>300 /*&& norm(psi[i-1][j-1])>300*/ && norm(psi[i-1][j])>300 && norm(psi[i-1][j+1])>300 && norm(psi[i+1][j-1])>300 && norm(psi[i+1][j])>300 && norm(psi[i+1][j+1])>300)){save_track(x_axis[i],y_axis[j],k*RTE_step,norm(psi[i][j]));break;}

		        else if(norm(psi[i][j])<300 && (norm(psi[i][j+1])>300 && norm(psi[i][j-1])>300 && norm(psi[i-1][j-1])>300 /*&& norm(psi[i-1][j])>300*/ && norm(psi[i-1][j+1])>300 && norm(psi[i+1][j-1])>300 && norm(psi[i+1][j])>300 && norm(psi[i+1][j+1])>300)){save_track(x_axis[i],y_axis[j],k*RTE_step,norm(psi[i][j]));break;}

		        else if(norm(psi[i][j])<300 && (norm(psi[i][j+1])>300 && norm(psi[i][j-1])>300 && norm(psi[i-1][j-1])>300 && norm(psi[i-1][j])>300 /*&& norm(psi[i-1][j+1])>300*/ && norm(psi[i+1][j-1])>300 && norm(psi[i+1][j])>300 && norm(psi[i+1][j+1])>300)){save_track(x_axis[i],y_axis[j],k*RTE_step,norm(psi[i][j]));break;}

		        else if(norm(psi[i][j])<300 && (norm(psi[i][j+1])>300 && norm(psi[i][j-1])>300 && norm(psi[i-1][j-1])>300 && norm(psi[i-1][j])>300 && norm(psi[i-1][j+1])>300 /*&& norm(psi[i+1][j-1])>300*/ && norm(psi[i+1][j])>300 && norm(psi[i+1][j+1])>300)){save_track(x_axis[i],y_axis[j],k*RTE_step,norm(psi[i][j]));break;}

		        else if(norm(psi[i][j])<300 && (norm(psi[i][j+1])>300 && norm(psi[i][j-1])>300 && norm(psi[i-1][j-1])>300 && norm(psi[i-1][j])>300 && norm(psi[i-1][j+1])>300 && norm(psi[i+1][j-1])>300 /*&& norm(psi[i+1][j])>300*/ && norm(psi[i+1][j+1])>300)){save_track(x_axis[i],y_axis[j],k*RTE_step,norm(psi[i][j]));break;}

		        else if(norm(psi[i][j])<300 && (norm(psi[i][j+1])>300 && norm(psi[i][j-1])>300 && norm(psi[i-1][j-1])>300 && norm(psi[i-1][j])>300 && norm(psi[i-1][j+1])>300 && norm(psi[i+1][j-1])>300 && norm(psi[i+1][j])>300 /*&& norm(psi[i+1][j+1])>300*/)){save_track(x_axis[i],y_axis[j],k*RTE_step,norm(psi[i][j]));break;}

	        }      
		blank_line();
        }       
        blank_line();
}

void save_2D() //Function to save the data to a file (with compatible blocks for gnuplot)
{
        for(i=0;i<n_x;i++)
        {
	        for(j=0;j<n_y;j++)
	        {
		        save_obdm(x_axis[i],y_axis[j],norm(psi[i][j]),phase_save(i,j));
	        }
                blank_line();
        }
        blank_line();
}

//>>>>>Main Program<<<<< 

int main()
{
	openDataFiles_obdm(1,name); //Open file with name (1,name)

	for(i=0;i<n_x;i++){x_axis[i]=-min_x+i*real(h_x);} //Initialisation of the x-axis from -min_x to +min_x
	for(j=0;j<n_y;j++){y_axis[j]=-min_y+j*real(h_y);} //Initialisation of the y-axis from -min_y to +min_y

	for(i=0;i<n_x;i++) //Initialise the wavefunction
	{
		for(j=0;j<n_y;j++)
		{
		        psi[i][j]=complex<double>(gauss(x_axis[i],y_axis[j]),0);
			psi_copy[i][j]=psi[i][j];
			phase[i][j]=0;
		}	
	}
	
	//====> Imaginary Time Propagation (ITP)

	complex<double> t_ITP(ITP_step,0); //Time-step size for ITP (the equations already assume imaginary time so a real t_ITP should be used)				

	for(k=0;k<n_it_ITP;k++)	//Time-evolution given by the 4th-order Runge-Kutta iteration			
	{
		//k1		
		for(i=1;i<n_x-1;i++){for(j=1;j<n_y-1;j++){ k1[i][j]=function_ITP(i,j) ;}}
		
		for(i=0;i<n_x;i++) //Fixed-derivative (Neumann) boundary conditions for the ITP
		{ 
		        k1[i][0]=function_ITP_BC(i,0); 
		        k1[i][n_x-1]=function_ITP_BC(i,n_x-1);
		}

		for(j=0;j<n_y;j++)
		{ 
		        k1[0][j]=function_ITP_BC(0,j); 
		        k1[n_y-1][j]=function_ITP_BC(n_y-1,j);
		}
		
		//k2
		for(i=0;i<n_x;i++){for(j=0;j<n_y;j++){ psi_copy[i][j]=psi[i][j]+half*t_ITP*k1[i][j] ;}}		
		
		for(i=1;i<n_x-1;i++){for(j=1;j<n_y-1;j++){ k2[i][j]=function_ITP(i,j) ;}}
		
		for(i=0;i<n_x;i++) //Fixed-derivative (Neumann) boundary conditions for the ITP
		{ 
		        k2[i][0]=function_ITP_BC(i,0); 
		        k2[i][n_x-1]=function_ITP_BC(i,n_x-1);
		}

		for(j=0;j<n_y;j++)
		{ 
		        k2[0][j]=function_ITP_BC(0,j); 
		        k2[n_y-1][j]=function_ITP_BC(n_y-1,j);
		}
		
		//k3		
		for(i=0;i<n_x;i++){for(j=0;j<n_y;j++){ psi_copy[i][j]=psi[i][j]+half*t_ITP*k2[i][j] ;}}
		
		for(i=1;i<n_x-1;i++){for(j=1;j<n_y-1;j++){ k3[i][j]=function_ITP(i,j) ;}}
		
		for(i=0;i<n_x;i++) //Fixed-derivative (Neumann) boundary conditions for the ITP
		{ 
		        k3[i][0]=function_ITP_BC(i,0); 
		        k3[i][n_x-1]=function_ITP_BC(i,n_x-1);
		}

		for(j=0;j<n_y;j++)
		{ 
		        k3[0][j]=function_ITP_BC(0,j); 
		        k3[n_y-1][j]=function_ITP_BC(n_y-1,j);
		}
		
		//k4		
		for(i=0;i<n_x;i++){for(j=0;j<n_y;j++){ psi_copy[i][j]=psi[i][j]+t_ITP*k3[i][j] ;}}
		
		for(i=1;i<n_x-1;i++){for(j=1;j<n_y-1;j++){ k4[i][j]=function_ITP(i,j) ;}}
		
		for(i=0;i<n_x;i++) //Fixed-derivative (Neumann) boundary conditions for the ITP
		{ 
		        k4[i][0]=function_ITP_BC(i,0); 
		        k4[i][n_x-1]=function_ITP_BC(i,n_x-1);
		}
	       
		for(j=0;j<n_y;j++)
		{ 
		        k4[0][j]=function_ITP_BC(0,j); 
		        k4[n_y-1][j]=function_ITP_BC(n_y-1,j);
		}
		
		//Time-step psi		
		for(i=0;i<n_x;i++)
		{		
			for(j=0;j<n_y;j++)		
			{
				psi[i][j]+=(t_ITP/six)*(k1[i][j]+two*k2[i][j]+two*k3[i][j]+k4[i][j]);
				psi_copy[i][j]=psi[i][j];
			}
		}
		
		rescale();

		if(k==vortex_start){add_vortex();} //Add the vortex near the end of the imaginary time propagation 
		
	       	//if(k>0 && k%n_save_ITP==0){save_2D();} //New block every multiple of n_save_ITP
  
		counter_ITP+=1; //Loading counter for ITP
   		if(counter_ITP%(n_it_ITP/100)==0){cout<<"ITP "<<(counter_ITP/(n_it_ITP/100))<<"%"<<endl;}
	}
	
	//====> Real Time Expansion (RTE)

	complex<double> t_RTE(RTE_step,0); //Time-step size for RTE				

	for(k=0;k<n_it_RTE;k++)	//Time-evolution given by the 4th-order Runge-Kutta iteration			
	{

	        //Fixed (Dirichlet) boundary conditions for the RTE

                for(l=0;l<n_x;l++){ psi[l][0]=zero; psi_copy[l][0]=zero; psi[l][n_y-1]=zero; psi_copy[l][n_y-1]=zero; }
	        for(m=0;m<n_y;m++){ psi[0][m]=zero; psi_copy[0][m]=zero; psi[n_x-1][m]=zero; psi_copy[n_x-1][m]=zero; }

		//k1		
		for(i=1;i<n_x-1;i++){for(j=1;j<n_y-1;j++){ k1[i][j]=function_RTE(i,j) ;}}

		//k2
		for(i=0;i<n_x;i++){for(j=0;j<n_y;j++){ psi_copy[i][j]=psi[i][j]+half*t_RTE*k1[i][j] ;}}		
		
		for(i=1;i<n_x-1;i++){for(j=1;j<n_y-1;j++){ k2[i][j]=function_RTE(i,j) ;}}

		//k3		
		for(i=0;i<n_x;i++){for(j=0;j<n_y;j++){ psi_copy[i][j]=psi[i][j]+half*t_RTE*k2[i][j] ;}}
		
		for(i=1;i<n_x-1;i++){for(j=1;j<n_y-1;j++){ k3[i][j]=function_RTE(i,j) ;}}

		//k4		
		for(i=0;i<n_x;i++){for(j=0;j<n_y;j++){ psi_copy[i][j]=psi[i][j]+t_RTE*k3[i][j] ;}}
		
		for(i=1;i<n_x-1;i++){for(j=1;j<n_y-1;j++){ k4[i][j]=function_RTE(i,j) ;}}

		//Time-step psi		
		for(i=0;i<n_x;i++)
		{		
			for(j=0;j<n_y;j++)		
			{
				psi[i][j]+=(t_RTE/six)*(k1[i][j]+two*k2[i][j]+two*k3[i][j]+k4[i][j]);
				psi_copy[i][j]=psi[i][j];
			}
		}
		
		//if(real(Integral) < 3495 || real(Integral) > 3505){rescale();}
		//rescale();

		if(k%n_save_RTE==0){track();}

		//if(k>=0 && k%n_save_RTE==0){save_2D();} //New block every multiple of n_save_RTE
  
		counter_RTE+=1; //Loading counter for RTE
   		if(counter_RTE%(n_it_RTE/100)==0){cout<<"RTE "<<(counter_RTE/(n_it_RTE/100))<<"%"<<endl;}
	}

	closeDataFiles_obdm(); //Close the file
	return 0;
}
