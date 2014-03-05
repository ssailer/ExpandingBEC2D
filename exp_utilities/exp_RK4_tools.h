#ifndef EXP_RK4_TOOLS_H__
#define EXP_RK4_TOOLS_H__

#include <iostream>
#include <complex>
#include <math.h>
#include <complexgrid.h>
#include <bh3binaryfile.h>
#include <vector>
#include <omp.h>
#include <cstring>


using namespace std;

// Inherit PathOptions from bh3binaryfile.h with additional Options for RK4 and the Potential
typedef struct : PathOptions {
  complex<double> omega_x,omega_y; // Frequency of the harmonic trap
  double min_x,min_y; // Coordinate boundaries
  complex<double> scale_factor; //Scale factor
  complex<double> t_abs; //Absolute time 
  complex<double> exp_factor; //Expansion factor
  double g; // coupling constant
  double ITP_step, RTE_step; // stepsize for the timeiteration
  int n_it_ITP; // number of timesteps
  int n_it_ITP1; // number of timesteps
  int n_it_ITP2; // number of timesteps
  int n_it_RTE; // number of timesteps
  int n_save_RTE; // times, when to save the process 
  int n_save_ITP; // replace with snapshot_times
  int times; // naming of the datafile - time of the snapshot
  std::string name; // naming of the datafile
  std::string config; // name of the config file
  bool startgrid[2];
  int threads;
  //Vortex Positions and winding Number
  int Q;
	
} Options;



class RK4
{
  public:
    RK4();
    // RK4(Options &opt);
    RK4(ComplexGrid* &c,Options &opt);    
    ~RK4();
    
    // Propagatoren
    void itpToTime(Options &opt);
    void rteToTime(Options &opt);
    void ITP(ComplexGrid* & pPsi,Options &opt);
    void RTE(ComplexGrid* & pPsi,Options &opt);
    
    
    // Hilfsfunktionen  
  
    double phase_save(ComplexGrid* & pPsi,int a,int b);
    
    // save the Grid to file
    void save_2D(ComplexGrid* & pPsi,Options &opt);
   
    // StorageObjects for the wavefunction and its phase
    ComplexGrid* pPsi;
    ComplexGrid* pPhase;
    
    // Coordinates
    vector<double> x_axis,y_axis;
    
    
  
         

  private:

    // double gauss(double x,double y); //A simple Gaussian
   
    // Scaling of Wavefunction after every timestep in ITP and RTE
    void rescale(ComplexGrid* & pPsi, Options &opt);
    
    // Hilfsfunktionen fuer ITP
    void computeK_ITP(ComplexGrid* &pPsi, vector<ComplexGrid> &k,Options &opt,complex<double> &t_ITP);
    void computeK_RTE(ComplexGrid* &pPsi, vector<ComplexGrid> &k,Options &opt,complex<double> &t_RTE);
    void TimeStepRK4(ComplexGrid* &pPsi,vector<ComplexGrid> &k,Options &opt,complex<double> &t);
    void Neumann(ComplexGrid &k,ComplexGrid &PsiCopy,Options &opt);
    void Dirichlet(ComplexGrid* &pPsi,Options &opt);

    complex<double> T(ComplexGrid & pPsiCopy,int i, int j);
    complex<double> V(ComplexGrid & pPsicopy,int i, int j,Options &opt);   
    
    // Hilfsfunktionen fuer RTE
    complex<double> function_RTE(ComplexGrid & pPsiCopy,int i, int j,Options &opt);
    complex<double> interaction(complex<double> a,Options &opt);
    complex<double> grad_x(complex<double> a, complex<double> b);
    complex<double> grad_y(complex<double> a, complex<double> b);
    complex<double> lambda_x(Options &opt);
    complex<double> lambda_x_dot(Options &opt);
    complex<double> lambda_y(Options &opt);
    complex<double> lambda_y_dot(Options &opt);
    complex<double> x_expand(complex<double> a,Options &opt);
    complex<double> y_expand(complex<double> a,Options &opt);
    
    // Hilfsvariablen
    complex<double> h_x, h_y;
    complex<double> integral(ComplexGrid* & pPsi,Options &opt);
    complex<double> Integral;
    complex<double> Integral_aux;
    
    
    // some used constants
    
    double pi; //acos(-1.0L);
    complex<double>  zero,half,one,two,four,six,i_unit;
      

  
};

#endif // EXP_RK4_TOOLS_H__