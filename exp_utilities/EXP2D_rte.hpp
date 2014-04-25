#ifndef EXP2D_RTE_H__
#define EXP2D_RTE_H__

#include <iostream>
#include <complex>
#include <math.h>
#include <complexgrid.h>
#include <bh3binaryfile.h>
#include <gauss_random.h>
#include <vector>
#include <omp.h>
#include <string>
#include <iomanip>
#include <EXP2D_tools.h>
#include <EXP2D_observables.h>
#include <plot_with_mgl.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

typedef struct {
        int absoluteSteps;
        int lambdaSteps;
    } StepCounter;

class RTE
{
  public:
    RTE();
    RTE(ComplexGrid* &c,Options &opt);    
    ~RTE();

    void setOptions(Options &externaloptions);
    void RunSetup();
    
    // Propagatoren

    void rteToTime(string runname, vector<int> snapshot_times,Averages* &eval);    
   
    // StoragePointer for the wavefunction
    ComplexGrid* pPsi;

    // Storage Variable for the runs
    MatrixXcd wavefct;
    vector<MatrixXcd> wavefctVec;

    void CopyComplexGridToEigen();
    void CopyEigenToComplexGrid();
    
    // Coordinates
    vector<double> x_axis,y_axis;
    VectorXcd X,Y;
    VectorXd Xexpanding, Yexpanding;

    int samplesize;

    // Plotting and progress functions 
    
    void cli_plot(string name,int counter_state, int counter_max, double start,bool plot);
    void cli(string name,int &slowestthread, vector<int> threadinfo, vector<int> stateOfLoops, int counter_max, double start);
    void plot(string name,int counter_state, int counter_max);
    

    // internal RunOptions, use setOptions(Options) to update from the outside
    Options opt;

    

    


  private:

    //
    inline void RTE_compute_k(MatrixXcd &k,MatrixXcd &wavefctcp,int &t);
    // inline void RTE_compute_k_pot(MatrixXcd &k,MatrixXcd &wavefctcp,int &t);
    void ToEigenAndNoise(ComplexGrid g,MatrixXcd &wavefct);


    // Variables
    complex<double> h_x, h_y;
    complex<double> pot_laplacian_x;
    complex<double> pot_laplacian_y;
    vector<double> ranges;
    Matrix<std::complex<double>,Dynamic,Dynamic,ColMajor> wavefctcpX;
    Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor> wavefctcpY;
    MatrixXcd PotentialGrid;
    VectorXcd laplacian_coefficient_x,laplacian_coefficient_y,gradient_coefficient_x,gradient_coefficient_y;

    StepCounter KeeperOfTime;
    
    
    // some used constants
    
    double pi;
    complex<double>  zero,half,one,two,four,six,i_unit;
    complex<double> t_RTE;

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

#endif // EXP2D_RTE_H__