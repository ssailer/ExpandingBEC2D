#ifndef EXP_RK4_TOOLS_H__
#define EXP_RK4_TOOLS_H__

#include <iostream>
#include <complex>
#include <math.h>
#include <complexgrid.h>
#include <bh3binaryfile.h>
#include <vector>


using namespace std;

// Inherit PathOptions from bh3binaryfile.h with additional Options for RK4 and the Potential
typedef struct : PathOptions {
  complex<double> omega_x,omega_y; // Frequency of the harmonic trap
  double min_x,min_y; // Coordinate boundaries
  complex<double> scale_factor; //Scale factor
  complex<double> t_abs; //Absolute time 
  complex<double> exp_factor; //Expansion factor
  double g;
  double ITP_step, RTE_step;
  int name;
	
} Options;



class RK4
{
  public:
    RK4();
    RK4(ComplexGrid* &c,Options &opt);    
    ~RK4();
    
    // Propagatoren
    void ITP(ComplexGrid* & pPsi,Options &opt);
    void RTE(ComplexGrid* & pPsi,Options &opt);
    
    
    // Hilfsfunktionen

  
  
    double vortex(int a, int b, int x, int y);
    double phase_save(ComplexGrid* & pPsi,int a,int b);
//     void add_vortex(ComplexGrid* & pPsi,Options &opt);
    
    // save the Grid to file
    void save_2D(ComplexGrid* & pPsi,Options &opt);
   
    // StorageObjects for the wavefunction and its phase
    ComplexGrid* pPsi;
    ComplexGrid* pPhase;
    
    // Coordinates
    vector<double> x_axis,y_axis;
    
    
  
         

  private:
    
   
    // Scaling of Wavefunction after every timestep in ITP and RTE
    void rescale(ComplexGrid* & pPsi,ComplexGrid* & pPsiCopy, Options &opt);
    
    // Hilfsfunktionen fuer ITP
    void computeK(ComplexGrid* & pPsiCopy,ComplexGrid* & pPsi, vector<ComplexGrid*> & k,Options & opt,complex<double> & t_ITP, int d);
    complex<double> T(ComplexGrid* & pPsiCopy,int i, int j);
    complex<double> V(ComplexGrid* & pPsicopy,int i, int j,Options &opt);   
    
    // Hilfsfunktionen fuer RTE
    complex<double> function_RTE(ComplexGrid* & pPsiCopy,int i, int j, complex<double> t,Options &opt);
    complex<double> interaction(complex<double> a,Options &opt);
    complex<double> grad_x(complex<double> a, complex<double> b);
    complex<double> grad_y(complex<double> a, complex<double> b);
    complex<double> lambda_x(complex<double> t, Options &opt);
    complex<double> lambda_x_dot(complex<double> t, Options &opt);
    complex<double> lambda_y(complex<double> t, Options &opt);
    complex<double> lambda_y_dot(complex<double> t, Options &opt);
    complex<double> x_expand(complex<double> a, complex<double> t, Options &opt);
    complex<double> y_expand(complex<double> a, complex<double> t, Options &opt);
    
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