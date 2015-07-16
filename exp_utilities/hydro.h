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
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <cstdlib>
#include <EXP2D_tools.h>
#include <EXP2D_MatrixData.h>
#include <python2.7/Python.h>
// #include "gnuplot-iostream.h"
// #include <../exp_utilities/plot_with_mgl.h>

using namespace std;

class hydroSolver {
public:
  hydroSolver(Eval* &e, double &maxTime) : eval(e), hbar(1.054e-22), m(87 * 1.66e-27), ti(0), tmax(maxTime) {
    g = eval->opt.g * (hbar * hbar / (m * m)) * (4.0) * eval->opt.N / M_PI;
    // g = eval->opt.g * (hbar * hbar / (m * m)) * (15.0/8.0) * eval->opt.N / M_PI;
  }

  hydroSolver() : hbar(1.054e-22), m(87 * 1.66e-27) {}
  double rk4_2nd(double, double, double, double,double&, double&);
  double f1(double, double, double);
  double f2(double, double, double);

  void integrate();
  void integrate2();
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



void hydroSolver::integrate()
{   
    PchangingValue = new double;
    PchangingValueQ = new double;
    bool xxx = true;
    // cout << "g = " << g << endl;
    cerr << eval->opt.vortexnumber << endl;
    double Nv = eval->opt.vortexnumber;
    // int Nv = 0;
    double r[2] = {eval->totalResult.Rx,eval->totalResult.Ry};
    double v[2] = {0.0,0.0};
    
    double xi, vi, tf, xf, vf, dt;
    double energy;
    vector<double> T, X, Y, Xdot, Ydot;

    string name = "runObservables/hydro.dat";
    ofstream file;
    file.open (name);
    file.precision(20);
    file.setf(ios::fixed | ios::showpoint);

    beta = 4 * hbar * hbar * Nv * Nv / (m * m);
    // zeta = 4 * hbar * Nv / m;

    dt = 1.0e-6;             // step size for integration

    X.push_back(r[0]);
    Y.push_back(r[1]);

/* integration of ODE */
    while (ti <= tmax)
    {   
        tf = ti + dt;
        T.push_back(tf);

        xxx = true;
        xi = r[0];
        vi = v[0];
        *PchangingValue = r[1];
        *PchangingValueQ = v[1];        
        rk4_2nd(ti,xi,vi,tf,xf,vf);
        r[0] = xf;
        v[0] = vf;
        X.push_back(xf);
        Xdot.push_back(vf);

        xxx = false;
        xi = r[1];
        vi = v[1];
        *PchangingValue = r[0];
        *PchangingValueQ = v[0];         
        rk4_2nd(ti,xi,vi,tf,xf,vf);
        r[1] = xf;
        v[1] = vf;
        Y.push_back(xf);
        Ydot.push_back(vf);

        ti = tf;
    }

  ti = 0.0;
  for(int i = 0; i < X.size(); ++i){
    file << setw(12) << ti << "," << setw(12) << X[i] << setw(12) << "," << Y[i]   << endl;
    ti += dt;
  }

  delete PchangingValue, PchangingValueQ;
}

void hydroSolver::integrate2()
{   
    ti = 0.0;
    PchangingValue = new double;
    PchangingValueQ = new double;
    bool xxx = true;
    double Nv = 0;
    double r[2] = {eval->totalResult.Rx,eval->totalResult.Ry};
    double v[2] = {0.0,0.0};
    
    double xi, vi, tf, xf, vf, dt;
    double energy;
    vector<double> T, X, Y, Xdot, Ydot;

    string name = "runObservables/hydro2.dat";
    ofstream file;
    file.open (name);
    file.precision(20);
    file.setf(ios::fixed | ios::showpoint);

    beta = 4 * hbar * hbar * Nv * Nv / (m * m);
    zeta = 4 * hbar * Nv / m;

    dt = 1.0e-6;             // step size for integration

    X.push_back(r[0]);
    Y.push_back(r[1]);

/* integration of ODE */
    while (ti <= tmax)
    {   
        tf = ti + dt;
        T.push_back(tf);

        xxx = true;
        xi = r[0];
        vi = v[0];
        *PchangingValue = r[1];
        *PchangingValueQ = v[1];        
        rk4_2nd(ti,xi,vi,tf,xf,vf);
        r[0] = xf;
        v[0] = vf;
        X.push_back(xf);
        Xdot.push_back(vf);

        xxx = false;
        xi = r[1];
        vi = v[1];
        *PchangingValue = r[0];
        *PchangingValueQ = v[0];         
        rk4_2nd(ti,xi,vi,tf,xf,vf);
        r[1] = xf;
        v[1] = vf;
        Y.push_back(xf);
        Ydot.push_back(vf);

        ti = tf;

    }

  ti = 0.0;
  for(int i = 0; i < X.size(); ++i){
    file << setw(12) << ti << "," << setw(12) << X[i] << setw(12) << "," << Y[i]   << endl;
    ti += dt;
  }

  delete PchangingValue, PchangingValueQ;
}

void hydroSolver::pyPlot(){    
    FILE * f = popen( "python ../plot.py", "r" );
    if ( f == 0 ) {
        fprintf( stderr, "Could not execute pyPlot()\n" );
    } else {
      const int BUFSIZE = 1000;
      char buf[ BUFSIZE ];
      while( fgets( buf, BUFSIZE,  f ) ) {
          fprintf( stdout, "%s", buf  );
      }
    }

    pclose( f );
    
}

/*
  Definition of the x'(t) = f1(t,x,x') = x' by the definition
*/
    double hydroSolver::f1(double t, double x, double v)
{
    double d1x;
    d1x = v;
    return d1x;
}
/*
 *  Definition of the x"(t) = f2(t,x,x')
*/
    double hydroSolver::f2(double t, double x, double v)
{
    double d2x;
    double secondTerm;
    d2x = g / (x * x * *PchangingValue) + beta * x / ((x * x + *PchangingValue * *PchangingValue ) * (x * x + *PchangingValue * *PchangingValue));
    // secondTerm = /*zeta * */( x * *PchangingValue *  *PchangingValueQ - *PchangingValue * *PchangingValue * v ) / ((x * x + *PchangingValue * *PchangingValue ) * (x * x + *PchangingValue * *PchangingValue));
    // cerr << "first " << d2x << " second " << secondTerm  << " * " << zeta << endl;
    // if(xxx == false){
    // // secondTerm = zeta * ( *PchangingValue * *PchangingValue * v - x * *PchangingValue * *PchangingValueQ) / ((x * x + *PchangingValue * *PchangingValue ) * (x * x + *PchangingValue * *PchangingValue));
    //   result = d2x * cos(t) 
    // } else {
    
    // }
    // secondTerm = x * *PchangingValueQ - *PchangingValue * v;
    // cerr << "secondTerm: " << secondTerm << endl;
    // double result;

    return d2x;
}

double hydroSolver::rk4_2nd(double ti, double xi, double vi, double tf, double& xf, double& vf)
{
      double h,t,k1x,k2x,k3x,k4x,k1v,k2v,k3v,k4v;

      h = tf-ti;
      t = ti;

      k1x = h*f1(t,xi,vi);
      k1v = h*f2(t,xi,vi);

      k2x = h*f1(t+h/2.0,xi+k1x/2.0,vi+k1v/2.0);
      k2v = h*f2(t+h/2.0,xi+k1x/2.0,vi+k1v/2.0);

      k3x = h*f1(t+h/2.0,xi+k2x/2.0,vi+k2v/2.0);
      k3v = h*f2(t+h/2.0,xi+k2x/2.0,vi+k2v/2.0);

      k4x = h*f1(t+h,xi+k3x,vi+k3v);
      k4v = h*f2(t+h,xi+k3x,vi+k3v);

      xf = xi + (k1x + 2.0*(k2x+k3x) + k4x)/6.0;
      vf = vi + (k1v + 2.0*(k2v+k3v) + k4v)/6.0;

      return 0.0;
}
