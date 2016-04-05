#ifndef HYDRO_H__
#define HYDRO_H__

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <cstdlib>

#include <python2.7/Python.h>

#include "tools.h"
#include "matrixdata.h"
#include "evaluation.h"

using namespace std;

typedef struct {
    double n0;
    double mu;
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
  hydroSolver(shared_ptr<Eval> e, double &maxTime) : eval(e), hbar(1.0545718e-22), m(86.9091835 *  1.660538921e-27), ti(0), tmax(maxTime) {
    // g = eval->opt.g * (hbar * hbar / (m * m)) * (4.0) * eval->opt.N / M_PI;
    // double kappa = eval->opt.vortexnumber / 10.0;
    // g = sqrt(1 + kappa) * eval->opt.g * (hbar * hbar ) / (m * m );
    // g = eval->opt.g * (hbar * hbar ) / (m * m);
    g = eval->opt.g * (hbar * hbar ) / (m );

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
  double ode_mu(double x,const hydroParams& params);
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
  shared_ptr<Eval> eval;

  double beta;
  double zeta;
  double* PchangingValue;
  double* PchangingValueQ;
  double hbar;
  double m;
  double g;
  double omega_average;

  bool xxx;

  double ti;
  double tmax;
};

#endif // HYDRO_H__
