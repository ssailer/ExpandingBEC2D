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

    void setOptions(Options &externaloptions);
    void RunSetup();
    
    // Propagatoren
    void itpToTime(string runname, int runtime, bool plot);
    void rteToTime(string runname, int runtime, bool plot);    
   
    // StoragePointer for the wavefunction
    ComplexGrid* pPsi;
    // Storage Variable for the runs
    MatrixXcd wavefct;

    void CopyComplexGridToEigen();
    void CopyEigenToComplexGrid();
    
    // Coordinates
    vector<double> x_axis,y_axis;
    VectorXcd X,Y;
    VectorXd Xexpanding, Yexpanding;

    // internal RunOptions, use setOptions(Options) to update from the outside
    Options opt; 


  private:

    //
    void RTE_compute_k(MatrixXcd &k,MatrixXcd &wavefctcp,int &t);
    void RTE_compute_k_pot(MatrixXcd &k,MatrixXcd &wavefctcp,int &t);
    void ITP_compute_k(MatrixXcd &k,MatrixXcd &wavefctcp);
   
    // Scaling of Wavefunction after every timestep in ITP
    void rescale(MatrixXcd &wavefct);   
   
    // Plotting and progress functions 
    void cli_plot(string name,int counter_state, int counter_max, double start,bool plot);
    void cli_plot_expanding(vector<double> &ranges,string name,int counter_state, int counter_max, double start,bool plot);
    
    // Variables
    complex<double> h_x, h_y;
    complex<double> Integral;
    complex<double> itp_laplacian_x;
    complex<double> itp_laplacian_y;
    vector<double> ranges;
    Matrix<std::complex<double>,Dynamic,Dynamic,ColMajor> wavefctcpX;
    Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor> wavefctcpY;
    MatrixXcd PotentialGrid;
    VectorXcd laplacian_coefficient_x,laplacian_coefficient_y,gradient_coefficient_x,gradient_coefficient_y;
    
    
    // some used constants
    
    double pi;
    complex<double>  zero,half,one,two,four,six,i_unit;
    complex<double> t_RTE;
    complex<double> t_ITP;

    // little helper functions for stuff

   inline VectorXd x_expand(complex<double> &t){
    return X.real() * real(lambda_x(t));
   }

   inline VectorXd y_expand(complex<double> &t){
    return Y.real() * real(lambda_y(t));
   }

   inline complex<double> lambda_x(complex<double> &t){
    return sqrt(one+opt.exp_factor*opt.dispersion_x*opt.dispersion_x*t*t);
   }

   inline complex<double> lambda_x_dot(complex<double> &t){
    return (opt.exp_factor*opt.dispersion_x*opt.dispersion_x*t/sqrt(one+opt.exp_factor*opt.dispersion_x*opt.dispersion_x*t*t));
   }

   inline complex<double> lambda_y(complex<double> &t){
    return sqrt(one+opt.exp_factor*opt.dispersion_y*opt.dispersion_y*t*t);
   }

   inline complex<double> lambda_y_dot(complex<double> &t){
    return (opt.exp_factor*opt.dispersion_y*opt.dispersion_y*t/sqrt(one+opt.exp_factor*opt.dispersion_y*opt.dispersion_y*t*t));
   }
  
};

#endif // EXP2D_TOOLS_H__