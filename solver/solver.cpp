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
#include "gnuplot-iostream.h"
// #include <../exp_utilities/plot_with_mgl.h>

using namespace std;

/* function prototypes */
double f1(double, double, double);
double f2(double, double, double);
double euler2d(double(*)(double, double, double),
               double(*)(double, double, double),
               double, double, double, double,
               double&, double&);
double euler2m(double(*)(double, double, double),
               double(*)(double, double, double),
               double, double, double, double,
               double&, double&);
double rk4_2nd(double(*)(double, double, double),
               double(*)(double, double, double),
               double, double, double, double,
               double&, double&);



double beta; 
double* PchangingValue;
const double hbar = 1.054e-22;
const double m = 87 * 1.66e-27;
const double N = 2.0e5;
const double g = (hbar * hbar / (m * m)) * 0.145 * (4.0) * 1.8 * N / M_PI;

int main( int argc, char** argv)
{   
    PchangingValue = new double;
    cout << "g = " << g << endl;
    int Nv = 50;
    double r[2] = {30.0,40.0};
    double v[2] = {0.0,0.0};
    if(argc >= 4){
      Nv = atof(argv[3]);
      r[0] = atof(argv[1]);
      r[1] = atof(argv[2]);
      cout << "Got input: Nv =" << Nv << endl;
      cout << "Got input: Rx =" << r[0] << endl;
      cout << "Got input: Ry =" << r[1] << endl;
    }else{
      cout << "Got no input, using default Nv = " << Nv << endl;
    }
    // double r[2] = {30.0,40.0}; // 2.116e-09;
    
    // alpha = 3.12625609723e-13;
    // beta = 2.1305244952e-18 * Nv * Nv;
    
    double ti, xi, vi, tf, xf, vf, dt, tmax;
    double energy;
    int key;
    const string method[3] = {"simple Euler","modified Euler","4th order Runge-Kutta"};
    vector<double> T, X, Y, Xdot, Ydot;




/* output: file and formats */
    ofstream file;
    string name = "ode_Rx_Ry.dat";
    file.open (name);
    file.precision(20);
    file.setf(ios::fixed | ios::showpoint);
    cout.precision(20);
    cout.setf(ios::fixed | ios::showpoint);

/* initial information */
    key =  2;             // select a method (key = 0, 1, 2)
    ti = 0.0;             // initial value for variable
               // initial value for function x(t)
    beta = 4 * hbar * hbar * Nv * Nv / (m * m);
    cout << "beta " << beta << endl;
            // initial
    dt = 5.0e-8;             // step size for integration
    tmax = 50.0e-3;          // integrate from ti till tmax

    cout << "xi = " << xi << endl;
    X.push_back(r[0]);
    Y.push_back(r[1]);

/* end of initial information */

    file << setw(30) << method[key] << endl;
    file << setw(12) << "t" << "," << setw(12) << "Rx"<< "," << setw(12) << "Ry" << endl;
    


/* integration of ODE */
    while (ti <= tmax)
    {   
        tf = ti + dt;
        T.push_back(tf);

        xi = r[0];
        vi = v[0];
        *PchangingValue = r[1];        
        rk4_2nd(f1,f2,ti,xi,vi,tf,xf,vf);
        r[0] = xf;
        v[0] = vf;
        X.push_back(xf);
        Xdot.push_back(vf);

        xi = r[1];
        vi = v[1];
        *PchangingValue = r[0];        
        rk4_2nd(f1,f2,ti,xi,vi,tf,xf,vf);
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


  // std::vector<std::pair<double, double> > xy_pts;
  // for(int i = 0; i < T.size(); i++){
  //   xy_pts.push_back(std::make_pair(T[i],X[i]));
  // }
  // int number = (int)Nv;
  // string name = "aspect_Nv_" + to_string(number) + "RX.png";
  // cout << "Now plotting" << endl;
  // Gnuplot gp;
  // gp << "set term pngcairo\n";
  // gp << "set output \"" + name + "\" \n";
  // gp << "plot '-' with lines\n";
  // gp.send1d(xy_pts);

  // for(int i = 0; i < T.size(); i++){
  //   xy_pts.push_back(std::make_pair(T[i],Y[i]));
  // }
  // name = "aspect_Nv_" + to_string(number) + "RY.png";
  // cout << "Now plotting" << endl;
  // gp << "set term pngcairo\n";
  // gp << "set output \"" + name + "\" \n";
  // gp << "plot '-' with lines\n";
  // gp.send1d(xy_pts);


    // system ("pause");
    // cout << "Done" << endl;
    return 0;
}

/*
  Definition of the x'(t) = f1(t,x,x') = x' by the definition
*/
    double f1(double t, double x, double v)
{
    double d1x;
    d1x = v;
    return d1x;
}
/*
 *  Definition of the x"(t) = f2(t,x,x')
*/
    double f2(double t, double x, double v)
{
    double d2x;
    d2x = g / (x * x * *PchangingValue) + beta * x / ((x * x + *PchangingValue * *PchangingValue ) * (x * x + *PchangingValue * *PchangingValue));
    return d2x;
}

double rk4_2nd(double(*d1x)(double, double, double),
               double(*d2x)(double, double, double),
               double ti, double xi, double vi, double tf,
               double& xf, double& vf)
{
      double h,t,k1x,k2x,k3x,k4x,k1v,k2v,k3v,k4v;

      h = tf-ti;
      t = ti;

      k1x = h*d1x(t,xi,vi);
      k1v = h*d2x(t,xi,vi);

      k2x = h*d1x(t+h/2.0,xi+k1x/2.0,vi+k1v/2.0);
      k2v = h*d2x(t+h/2.0,xi+k1x/2.0,vi+k1v/2.0);

      k3x = h*d1x(t+h/2.0,xi+k2x/2.0,vi+k2v/2.0);
      k3v = h*d2x(t+h/2.0,xi+k2x/2.0,vi+k2v/2.0);

      k4x = h*d1x(t+h,xi+k3x,vi+k3v);
      k4v = h*d2x(t+h,xi+k3x,vi+k3v);

      xf = xi + (k1x + 2.0*(k2x+k3x) + k4x)/6.0;
      vf = vi + (k1v + 2.0*(k2v+k3v) + k4v)/6.0;

      return 0.0;
}
