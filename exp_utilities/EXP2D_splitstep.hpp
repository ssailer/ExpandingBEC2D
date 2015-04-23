#ifndef EXP2D_SPLITSTEP_H__
#define EXP2D_SPLITSTEP_H__

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
#include <gauss_random.h>
#include <EXP2D_tools.h>
#include <EXP2D_evaluation.h>
#include <EXP2D_binaryfile.h>
#include <plot_with_mgl.h>
#include <EXP2D_MatrixData.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

class SplitStep
{
  public:
    void assignMatrixData(MatrixData* &d);
    void setVariables();
    void timeStep(double delta_t);


    SplitStep(vector<ComplexGrid> &d,const MatrixData::MetaData &extMeta, const Options &externaloptions, int &extSLICE_NUMBER);

    void setOptions(const Options &externaloptions);
    void RunSetup();
    
    // Propagatoren

    void splitToTime(string runName);
    // void rteFromDataToTime(string runname, vector<int> snapshot_times, string h5name);    
   
    // StoragePointer for the wavefunction
    // MatrixData* pData;
    vector<ComplexGrid> wavefctVec;

    
    // Coordinates
    vector<double> x_axis,y_axis, z_axis;
    VectorXcd X,Y,Z;


    int samplesize;

    // Plotting and progress functions 
    
    // void cli_plot(string name,int counter_state, int counter_max, double start,bool plot);
    void cli(string name,int &slowestthread, vector<int> threadinfo, vector<int> stateOfLoops, int counter_max, double start);
    void plot(const string name);
    
    MatrixData* w;

    // internal RunOptions, use setOptions(Options) to update from the outside
    Options opt;
    vector<int> snapshot_times;

  private:

    vector<vector<double>> kspace;
    MatrixXcd kprop;
    MatrixXcd Vgrid;  

    MatrixData::MetaData meta;

    // inline void RTE_compute_k_pot(MatrixXcd &k,MatrixXcd &wavefctcp,int &t);
    inline double rotatingPotential(int &i, int &j, int &t);
   

    // Variables
    complex<double> h_x, h_y,h_z;
    ComplexGrid PotentialGrid;

    int SLICE_NUMBER;
    
    
    // some used constants
    
    double pi;
    complex<double>  zero,half,one,two,four,six,i_unit;
    complex<double> t_RTE;

    // little helper functions for stuff
};


#endif // EXP2D_SPLITSTEP_H__