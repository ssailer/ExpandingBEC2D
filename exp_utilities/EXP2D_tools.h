#ifndef EXP2D_TOOLS_H__
#define EXP2D_TOOLS_H__

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <complex>
#include <math.h>
#include <complexgrid.h>
#include <bh3binaryfile.h>
#include <vector>
#include <omp.h>
#include <string>
#include <iomanip>

using namespace std;
using namespace Eigen;

typedef struct {
        // From bh3binaryfile
    double N; // Number of particles
    int32_t grid[4];  // gridsize
        // my own
    complex<double> omega_x,omega_y; // Frequency of the harmonic trap
    complex<double> dispersion_x, dispersion_y; // dispersion relation for the expandion frame
    double min_x,min_y; // Coordinate boundaries
    complex<double> scale_factor; //Scale factor
    complex<double> t_abs; //Absolute time // remove from opt! put into the function, don't need it here
    complex<double> exp_factor; //Expansion factor
    double g; // coupling constant
    double ITP_step, RTE_step; // stepsize for the timeiteration
    int n_it_ITP; // number of timesteps
    int n_it_ITP1; // number of timesteps
    int n_it_ITP2; // number of timesteps
    int n_it_RTE; // number of timesteps
    int n_save_RTE; // times, when to save the process // don't need it here anymore
    int n_save_ITP; // replace with snapshot_times     // don't need it here anymore
    string name; // naming of the datafile      // think about that naming system remove it from here
    string config; // name of the config file 
    string workingdirectory;   // remove it from here, only needed in the program itself
    string workingfile;
    bool startgrid[3];
    //Vortex Positions and winding Number
    int Q;
    bool RTE_only;
    
} Options;


class EXP2D
{
  public:
    EXP2D();
    EXP2D(ComplexGrid* &c,Options &opt);    
    ~EXP2D();
    
    // Propagatoren
    void itpToTime(Options &opt,bool plot);
    void rteToTime(Options &opt, bool plot);    
   
    // StoragePointer for the wavefunction
    ComplexGrid* pPsi;
    
    // Coordinates
    vector<double> x_axis,y_axis; 

  private:

  void ITP(ComplexGrid* & pPsi,Options &opt);
  inline void RTE_compute_k(MatrixXcd &k,MatrixXcd &wavefctcp, VectorXcd &X,VectorXcd &Y,int &t,int &grid_x,int &grid_y);

    // double gauss(double x,double y); //A simple Gaussian
   
    // Scaling of Wavefunction after every timestep in ITP
    void rescale(ComplexGrid* & pPsi, Options &opt);
    
    // Hilfsfunktionen fuer ITP
    void computeK_ITP(ComplexGrid* &pPsi, vector<ComplexGrid> &k,Options &opt,complex<double> &t_ITP);
    void Neumann(ComplexGrid &k,ComplexGrid &PsiCopy,Options &opt);    
   
    // Hilfsfunktionen 
    void cli_plot(MatrixXcd& mPsi, Options &opt,string name,int counter_state, int counter_max, double start,bool plot);  
    void run_status(string name,int counter_state, int counter_max, double start);
    
    // Hilfsvariablen
    complex<double> h_x, h_y;
    complex<double> Integral;
    complex<double> Integral_aux;
    
    
    // some used constants
    
  double pi;
  complex<double>  zero,half,one,two,four,six,i_unit;
  complex<double> exp_factor, dispersion_x, dispersion_y;
  double g;
  complex<double> t_RTE;

  vector<complex<double>> laplacian_coefficient_x,laplacian_coefficient_y,gradient_coefficient_x,gradient_coefficient_y;  




    // inline functions for stuff

    inline complex<double> interaction(complex<double> a,Options &opt)
{return (opt.g*norm(a));} //Interaction term in the GPE Hamiltonian 

inline complex<double> integral(ComplexGrid* & pPsi,Options &opt)
{ 
  Integral_aux = complex<double>(0,0);  
  for(int i=0;i<opt.grid[1]-1;i++)
  {
    for(int j=0;j<opt.grid[2]-1;j++)
    {
      Integral_aux+=h_x*h_y*(norm(pPsi->at(0,i,j,0))+norm(pPsi->at(0,i+1,j,0))+norm(pPsi->at(0,i,j+1,0))+norm(pPsi->at(0,i+1,j+1,0)))/four;
      
    }
  }
  return Integral_aux;
}

    inline complex<double> itp_kinetic(ComplexGrid &PsiCopy,int i, int j)
{
    return half*((PsiCopy(0,i+1,j,0)-(two*PsiCopy(0,i,j,0))+PsiCopy(0,i-1,j,0))/(h_x*h_x))+half*((PsiCopy(0,i,j+1,0)-(two*PsiCopy(0,i,j,0))+PsiCopy(0,i,j-1,0))/(h_x*h_x)); 
}

inline complex<double> itp_potential(ComplexGrid & PsiCopy,int i, int j,Options & opt)
{ 
    complex<double> xvalue = complex<double>(x_axis[i],0);
    complex<double> yvalue = complex<double>(y_axis[j],0);
  
    return -(half*opt.omega_x*opt.omega_x*xvalue*xvalue+half*opt.omega_y*opt.omega_y*yvalue*yvalue+complex<double>(opt.g,0)*norm(PsiCopy(0,i,j,0)))*PsiCopy(0,i,j,0);
}

   

   // inline double x_expand(int i, Options &opt)
   // {
   //  return x_axis[i] * real(lambda_x(opt));
   // }

   // inline double y_expand(int j, Options &opt)
   // {
   //  return y_axis[j] * real(lambda_y(opt));
   // }


  inline complex<double> lambda_x(complex<double>& t)
   {
    return sqrt(one+exp_factor*dispersion_x*dispersion_x*t*t);
   }

   inline complex<double> lambda_x_dot(complex<double>& t)
   {
   return (exp_factor*dispersion_x*dispersion_x*t/sqrt(one+exp_factor*dispersion_x*dispersion_x*t*t));
   }

   inline complex<double> lambda_y(complex<double>& t)
   {
   return sqrt(one+exp_factor*dispersion_y*dispersion_y*t*t);
   }

   inline complex<double> lambda_y_dot(complex<double>& t)
   {
    return (exp_factor*dispersion_y*dispersion_y*t/sqrt(one+exp_factor*dispersion_y*dispersion_y*t*t));
   }



  
};

#endif // EXP2D_TOOLS_H__