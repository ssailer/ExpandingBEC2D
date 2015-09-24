#include "hydro.h"

void hydroSolver::integrate3()
{   
    PchangingValue = new double;
    PchangingValueQ = new double;
    // bool xxx = true;
    // cout << "g = " << g << endl;
    // cerr << eval->opt.vortexnumber << endl;
    double Nv = eval->opt.vortexnumber;

    // int Nv = 0;
    double r[2] = {eval->totalResult.Rx,eval->totalResult.Ry};
    double v[2] = {0.0,0.0};

    // double Nv = ( eval->opt.omega_w.real() * 2.0 * M_PI ) / ( hbar / ( m * r[0] * r[1]));
    cerr << "VortexNumber: " << Nv << endl;
    
    double xi, vi, tf, xf, vf, dt;
    double energy;
    vector<double> T, X, Y, Xdot, Ydot;

    string name = "runObservables/hydro.dat";
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

        double tmp0,tmp1;
        


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

        // tmp0 = v[0];
        // tmp1 = v[1];
        // v[1] = r[0] * tmp0 / r[1];

        // v[0] = r[1] * tmp1 / r[0];

        ti = tf;
    }

  ti = 0.0;
  for(int i = 0; i < X.size(); ++i){
    file << setw(12) << ti << "," << setw(12) << X[i] << setw(12) << "," << Y[i] << "," << Xdot[i] << "," << Ydot[i] << "," << X[i] * Ydot[i] - Y[i] * Xdot[i]  << endl;
    ti += dt;
  }

  delete PchangingValue, PchangingValueQ;
}

void hydroSolver::integrate2()
{   
    ti = 0.0;
    PchangingValue = new double;
    PchangingValueQ = new double;
    // bool xxx = true;
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
    file << setw(12) << ti << "," << setw(12) << X[i] << setw(12) << "," << Y[i] << "," << Xdot[i] << "," << Ydot[i] << "," << X[i] * Ydot[i] - Y[i] * Xdot[i]  << endl;
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

    f = popen( "python ../ellipse.py", "r" );
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
    // d1x = 0.0;
    // d1x = (x / *PchangingValue) * *PchangingValueQ;
    return d1x;
}
/*
 *  Definition of the x"(t) = f2(t,x,x')
*/
    double hydroSolver::f2(double t, double x, double v)
{
    double d2x;
    double secondTerm = 0;
    double factor0 = *PchangingValue / x;
    double delta = ((x * x + *PchangingValue * *PchangingValue ) * (x * x + *PchangingValue * *PchangingValue));
    d2x = g / (x * x * *PchangingValue) + beta * x / delta;
    // d2x = g / (*PchangingValue * *PchangingValue * *PchangingValue) + beta * x / delta + v * v / x - v * *PchangingValueQ / *PchangingValue;
    if(xxx == true){
    secondTerm = zeta * factor0 * ( *PchangingValue * v - x * *PchangingValueQ) / delta;
         // d2x = g / ( *PchangingValue * *PchangingValue * *PchangingValue) + beta * (x / delta) + (*PchangingValueQ / *PchangingValue) * v - (v * v) / x;
    } else {
    secondTerm = + zeta * factor0 * (*PchangingValue * v - x * *PchangingValueQ) / delta;
    //   // d2x = g / ( *PchangingValue * *PchangingValue * *PchangingValue) + beta * (x / delta) - (*PchangingValueQ / *PchangingValue) * v + (v * v) / x;
    }
    // std::cerr << std::setprecision (24) << secondTerm << std::endl;
    // CHECK FOR DIFFERENT MODEL PDF 

    return d2x + secondTerm;
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











double hydroSolver::ode_n0(double x,const hydroParams& params)
{
    double result;
    result = - x * (params.alpha_x + params.alpha_y);
    return result;
}

double hydroSolver::ode_a(double x,const hydroParams& params)
{
    double result;
    result = - ( 2.0 * ( params.alpha - params.omega)) / (params.sigma_x * params.sigma_x)
             - ( 2.0 * ( params.alpha + params.omega)) / (params.sigma_y * params.sigma_y)
             - x * ( params.alpha_x + params.alpha_y);
    return result;
}

double hydroSolver::ode_sigma_x(double x,const hydroParams& params)
{
    double result;
    result = x * params.alpha_x + 0.5 * x * x * x * params.a * (params.alpha + params.omega);
    return result;
}

double hydroSolver::ode_sigma_y(double x,const hydroParams& params)
{
    double result;
    result = x * params.alpha_y + 0.5 * x * x * x * params.a * (params.alpha - params.omega);
    return result;
}

double hydroSolver::ode_omega(double x,const hydroParams& params)
{
    double result;
    result = - x * (params.alpha_x + params.alpha_y);
    return result;
}

double hydroSolver::ode_alpha(double x,const hydroParams& params)
{
    double result;
    result = - x * (params.alpha_x + params.alpha_y) + 2.0 * params.n0 * g * params.a;
    return result;
}

double hydroSolver::ode_alpha_x(double x,const hydroParams& params)
{
    double result;
    result = - x * x - params.alpha * params.alpha + params.omega + params.omega + 2.0 * params.n0 * g / (params.sigma_x * params.sigma_x);
    return result;
}

double hydroSolver::ode_alpha_y(double x,const hydroParams& params)
{
    double result;
    result = - x * x - params.alpha * params.alpha + params.omega + params.omega + 2.0 * params.n0 * g  / (params.sigma_y * params.sigma_y);
    return result;
}

void hydroSolver::calc_phi(const hydroParams& params, double& result)
{
    // xxx = true;
    // double tmp = - params.a * params.sigma_x * params.sigma_x * params.sigma_y * params.sigma_y / (params.sigma_x * params.sigma_x - params.sigma_y * params.sigma_y);
    // if(atan(tmp) < 0)
    //     result = atan(tmp) * 180 / M_PI;
    // else
    //     result = atan(tmp) * 180 / M_PI;
    // result /= 2.0;
    double tmp1 = - params.a * params.sigma_x * params.sigma_x * params.sigma_y * params.sigma_y;
    double tmp2 = (params.sigma_x * params.sigma_x - params.sigma_y * params.sigma_y);
    // double at = atan2(tmp1,tmp2);
    double at = atan(tmp1/tmp2);
    // double at = atan(tmp1/tmp2);
    at *= ( 180 / M_PI ) / 2.0;
    // if(at < 0.0){
    //     at += 90;
    // }
    // if(at < 0.0){
    //     at += 180 ;
    // //         if(at > 90){
    // //             cerr << "what?";
    // //         }
    // }
    result = at;
}

void hydroSolver::calc_ratio(const hydroParams& params, double& result)
{
    double tmp1 = (params.sigma_x * params.sigma_x + params.sigma_y * params.sigma_y);
    double tmp2 = sqrt((params.sigma_x * params.sigma_x - params.sigma_y * params.sigma_y) * (params.sigma_x * params.sigma_x - params.sigma_y * params.sigma_y)
         + params.a * params.a * params.sigma_x * params.sigma_x * params.sigma_x * params.sigma_x * params.sigma_y * params.sigma_y * params.sigma_y * params.sigma_y);
    double avar = tmp1 - tmp2;
    double bvar = tmp1 + tmp2;
    result = sqrt(bvar/avar);
}


void hydroSolver::rk4_1st(double ti, double xi, double tf, const hydroParams& params, double& xf, double (hydroSolver::*func)(double, const hydroParams&))
{
      double h,t,k1x,k2x,k3x,k4x;

      h = tf-ti;
      t = ti;

      k1x = h * (this->*func)(xi,params);

      k2x = h * (this->*func)(xi+k1x/2.0,params);

      k3x = h * (this->*func)(xi+k2x/2.0,params);

      k4x = h * (this->*func)(xi+k3x,params);

      xf = xi + (k1x + 2.0*(k2x+k3x) + k4x)/6.0;
}

void hydroSolver::printParams(const hydroParams& params)
{
    cerr << "hydroParams" << endl
         << "sigma_x " << params.sigma_x << endl
         << "sigma_y " << params.sigma_y << endl
         << "phi " << params.phi << endl
         << "n0 " << params.n0 << endl
         << "alpha_x " << params.alpha_x << endl
         << "alpha_y " << params.alpha_y << endl
         << "alpha " << params.alpha << endl
         << "a " << params.a << endl
         << "omega " << params.omega << endl;
    std::cin.ignore();
}

void hydroSolver::integrate()
{   


    double tf, dt;

    ti = 0.0;
    dt = 1.0e-5;

    hydroParams varnew;
    hydroParams varold;

    varold.sigma_x = eval->totalResult.r_max;
    varold.sigma_y = eval->totalResult.r_min;
    varold.phi = eval->totalResult.r_max_phi;
    // varold.n0 = 2.0 * (/*eval->opt.N*/ eval->totalResult.particle_count / M_PI) / (varold.sigma_x * varold.sigma_y);
    varold.n0 = eval->totalResult.n0;
    varold.alpha_x = 0.0;
    varold.alpha_y = 0.0;
    varold.a = 0.0;

    // varold.alpha = 0.0;
    varold.alpha = hbar * eval->opt.vortexnumber / (m * varold.sigma_x * varold.sigma_y);
    // varold.alpha = 2.0 * M_PI * real(eval->opt.omega_w);

    // cerr << " alpha " << varold.alpha << " vs omega " << eval->opt.omega_w * 2.0 * M_PI << endl;
    
    varold.omega = 0.0; // 2.0 * M_PI * real(eval->opt.omega_w);
    // varold.omega = hbar * eval->opt.vortexnumber / (m * varold.sigma_x * varold.sigma_y);

    varold.ratio = varold.sigma_x / varold.sigma_y;

    // printParams(varold);


    string name = "runObservables/hydro.dat";
    ofstream file;
    file.open (name);
    file.precision(20);
    file.setf(ios::fixed | ios::showpoint);

    file << setw(12) << "ti      " << "," 
         << setw(12) << "sigma_x " << ","
         << setw(12) << "sigma_y " << ","
         << setw(12) << "majorAxis" << ","
         << setw(12) << "minorAxis" << ","
         << setw(12) << "phi     " << ","
         << setw(12) << "n0      " << ","
         << setw(12) << "alpha_x " << ","
         << setw(12) << "alpha_y " << ","
         << setw(12) << "alpha   " << ","
         << setw(12) << "a       " << ","
         << setw(12) << "omega   " << ","
         << setw(12) << "ratio   "
        << endl;

    file << setw(12) << ti << "," 
         << setw(12) << varold.sigma_x << ","
         << setw(12) << varold.sigma_y << ","
         << setw(12) << varold.sigma_x << ","
         << setw(12) << varold.sigma_y << ","
         << setw(12) << varold.phi     << ","
         << setw(12) << varold.n0      << ","
         << setw(12) << varold.alpha_x << ","
         << setw(12) << varold.alpha_y << ","
         << setw(12) << varold.alpha   << ","
         << setw(12) << varold.a       << ","
         << setw(12) << varold.omega   << ","
         << setw(12) << varold.ratio
         << endl;
/* integration of ODE */
    while (ti <= tmax)
    {   
        tf = ti + dt;

        rk4_1st(ti,varold.n0,tf,varold,varnew.n0,&hydroSolver::ode_n0);
        rk4_1st(ti,varold.a,tf,varold,varnew.a,&hydroSolver::ode_a);
        rk4_1st(ti,varold.sigma_x,tf,varold,varnew.sigma_x,&hydroSolver::ode_sigma_x);
        rk4_1st(ti,varold.sigma_y,tf,varold,varnew.sigma_y,&hydroSolver::ode_sigma_y);
        rk4_1st(ti,varold.alpha_x,tf,varold,varnew.alpha_x,&hydroSolver::ode_alpha_x);
        rk4_1st(ti,varold.alpha_y,tf,varold,varnew.alpha_y,&hydroSolver::ode_alpha_y);
        rk4_1st(ti,varold.alpha,tf,varold,varnew.alpha,&hydroSolver::ode_alpha);
        rk4_1st(ti,varold.omega,tf,varold,varnew.omega,&hydroSolver::ode_omega);
        calc_phi(varnew,varnew.phi);

        calc_ratio(varnew,varnew.ratio);



        ti = tf;
        varold = varnew;

        // printParams(varold);
        double majorAxis, minorAxis;
        (varold.sigma_x >= varold.sigma_y) ? (majorAxis = varold.sigma_x, minorAxis = varold.sigma_y) : (majorAxis = varold.sigma_y, minorAxis = varold.sigma_x);


        file << setw(12) << ti << "," 
             << setw(12) << varold.sigma_x << ","
             << setw(12) << varold.sigma_y << ","
             << setw(12) << majorAxis << ","
             << setw(12) << minorAxis << ","
             << setw(12) << eval->totalResult.r_max_phi + varold.phi     << ","
             << setw(12) << varold.n0      << ","
             << setw(12) << varold.alpha_x << ","
             << setw(12) << varold.alpha_y << ","
             << setw(12) << varold.alpha   << ","
             << setw(12) << varold.a       << ","
             << setw(12) << varold.omega   << ","
             << setw(12) << varold.ratio
        << endl;     

    }
}
