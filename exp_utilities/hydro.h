/*
 Solver for second order Ordinary Differential Equations

 x"(t) = f(t,x,x')    equation
 x(ti) = xi           initial value 1
 x'(ti)= vi           initial value 2

 the second order ODE is solved as a set of two first order ODEs
 x'(t) = f1(t,x)
 x"(t) = f2(t,x,x')

 Methods (select one by a key)
 key = 0; simple Euler
 key = 1; modified Euler (predictor-corrector)
 key = 2; 4-th order Runge-Kutta

 Alex Godunov: Last revision March 2007
*/

#ifndef HYDRO_H__
#define HYDRO_H__

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <cstdlib>
#include <EXP2D_tools.h>
#include <EXP2D_MatrixData.h>
#include <python2.7/Python.h>
#include <EXP2D_evaluation.h>
// #include "gnuplot-iostream.h"
// #include <../exp_utilities/plot_with_mgl.h>

using namespace std;

typedef struct {
    double n0;
    double a;
    double sigma_x;
    double sigma_y;
    double omega;
    double alpha;
    double alpha_x;
    double alpha_y;
    double phi;
    double ratio;
} hydroParams;

class hydroSolver {
public:
  hydroSolver(Eval* &e, double &maxTime) : eval(e), hbar(1.0545718e-22), m(86.9091835 *  1.660538921e-27), ti(0), tmax(maxTime) {
    // g = eval->opt.g * (hbar * hbar / (m * m)) * (4.0) * eval->opt.N / M_PI;
    g = eval->opt.g * (hbar * hbar ) / (m * m );
    // g = eval->opt.g * (hbar * hbar / (m * m)) * (15.0/8.0) * eval->opt.N / M_PI;
  }

  hydroSolver() : hbar(1.0545718e-22), m(86.9091835 *  1.660538921e-27) {}
  double rk4_2nd(double, double, double, double,double&, double&);
  void rk4_1st(double ti, double xi, double tf, const hydroParams& params, double& xf, double (hydroSolver::*func)(double, const hydroParams&));
  double f1(double, double, double);
  double f2(double, double, double);

  double ode_sigma_x(double x,const hydroParams& params);
  double ode_sigma_y(double x,const hydroParams& params);
  double ode_phi(double x,const hydroParams& params);
  double ode_n0(double x,const hydroParams& params);
  double ode_alpha_x(double x,const hydroParams& params);
  double ode_alpha_y(double x,const hydroParams& params);
  double ode_alpha(double x,const hydroParams& params);
  double ode_a(double x,const hydroParams& params);
  double ode_omega(double x,const hydroParams& params);
  void calc_phi(const hydroParams& params, double& result);
  void calc_ratio(const hydroParams& params, double& result);

  void printParams(const hydroParams& params); // debugging

  void integrate();
  void integrate2();
  void integrate3();
  void pyPlot();
  Eval* eval;

  double beta;
  double zeta;
  double* PchangingValue;
  double* PchangingValueQ;
  double hbar;
  double m;
  double g;

  bool xxx;

  double ti;
  double tmax;
};

#endif // HYDRO_H__
