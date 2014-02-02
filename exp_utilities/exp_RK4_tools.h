#ifndef EXP_RK4_TOOLS_H__
#define EXP_RK4_TOOLS_H__

#include <iostream>
#include <complex>
#include <math.h>
#include <complexgrid.h>
#include <bh3binaryfile.h>
#include <2dexpan.h>





using namespace std;

// Inherit PathOptions from bh3binaryfile.h with additional Options for RK4 and the Potential
typedef struct : PathOptions {
  complex<double> omega_x,omega_y; // Frequency of the harmonic trap
  double min_x,min_y; // Coordinate boundaries
  complex<double> scale_factor; //Scale factor
  complex<double> t_abs; //Absolute time 
  complex<double> exp_factor; //Expansion factor
  double g;
	
} Options;

class RK4
{
  public:
    RK4();
    RK4(ComplexGrid* &c,Options &opt);
    
    ~RK4();
    
    void ITP(ComplexGrid* & pPsi, Options &opt, complex<double> & t_ITP);
    void RTE(const int & n_x, const int & n_y, complex<double> & t_RTE);
    
    
    // Hilfsfunktionen
//     complex<double> laplacian_x(complex<double> a, complex<double> b, complex<double> c);
//     complex<double> laplacian_y(complex<double> a, complex<double> b, complex<double> c);
//     complex<double> potential(double x,double y, Options &opt);
//     complex<double> interaction(complex<double> a);
    complex<double> grad_x(complex<double> a, complex<double> b);
    complex<double> grad_y(complex<double> a, complex<double> b);
    complex<double> lambda_x(complex<double> t, Options &opt);
    complex<double> lambda_x_dot(complex<double> t, Options &opt);
    complex<double> lambda_y(complex<double> t, Options &opt);
    complex<double> lambda_y_dot(complex<double> t, Options &opt);
    complex<double> x_expand(complex<double> a, complex<double> t, Options &opt);
    complex<double> y_expand(complex<double> a, complex<double> t, Options &opt);
    complex<double> integral(ComplexGrid* & pPsi,Options &opt);
    complex<double> rescale(ComplexGrid* & pPsi,ComplexGrid* & pPsiCopy, Options &opt);
    double vortex(int a, int b, int x, int y);
    double phase_save(ComplexGrid* & pPsi,int a,int b);
//    void add_vortex(ComplexGrid* & pPsi,Options &opt);
    void save_2D(ComplexGrid* & pPsi,Options &opt);
//     complex<double> function_RTE(int i, int j, complex<double> t);
//     complex<double> function_ITP(int i,int j);
//     complex<double> function_ITP_BC(int i,int j);
    double gauss(double x,double y);

         

  private:
    void computeK(ComplexGrid* & pPsiCopy,ComplexGrid* & pPsi, ComplexGrid** k,Options & opt,complex<double> & t_ITP, int d);
    complex<double> T(ComplexGrid* & pPsiCopy,int i, int j);
    complex<double> V(ComplexGrid* & pPsicopy,int i, int j,Options &opt);    
    complex<double> h_x, h_y;
    double x_axis[],y_axis[];
    complex<double> Integral;
    complex<double> Integral_aux;
    
    // some used constants
    
    static const double pi; //acos(-1.0L);
    static const complex<double>  zero,half,one,two,four,six,i_unit;
    
    

    // const complex<double> N(1000,0); //Particle number (1000)
    // const complex<double> g(15,0); //Interaction constant (15)
    // const complex<double> omega_x(100,0),omega_y(150,0); //Trap frequency (100,150)
    // const double min_x=4,min_y=4; //Symmetric axis ranges start from -min_x and -min_y (6,6)
    
   

  
};

#endif // EXP_RK4_TOOLS_H__