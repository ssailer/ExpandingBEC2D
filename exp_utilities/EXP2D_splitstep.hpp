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
#include <EXP2D_constants.h>
#include <plot_with_mgl.h>
#include <EXP2D_MatrixData.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

class SplitStep
{
  public:
    SplitStep(Options &o);
    void assignMatrixData(MatrixData* &d);
    void setVariables();
    virtual void timeStep(double delta_t) = 0;


    // SplitStep(vector<ComplexGrid> &d,const MatrixData::MetaData &extMeta, const Options &externaloptions, int &extSLICE_NUMBER);


    void setOptions(const Options &externaloptions);
    void RunSetup();
    
    // Propagatoren

    void splitToTime(string runName);
    // void rteFromDataToTime(string runname, vector<int> snapshot_times, string h5name);    
   
    // StoragePointer for the wavefunction
    // MatrixData* pData;
    vector<ComplexGrid> wavefctVec;

    int samplesize;

    // Plotting and progress functions 
    
    // void cli_plot(string name,int counter_state, int counter_max, double start,bool plot);
    void cli(string name,int &slowestthread, vector<int> threadinfo, vector<int> stateOfLoops, int counter_max, double start);
    void plot(const string name);
    
    MatrixData* w;

    // internal RunOptions, use setOptions(Options) to update from the outside
    Options opt;
    vector<int> snapshot_times;

  protected:

    vector<vector<double>> kspace;
    vector<double> x_axis,y_axis;
    MatrixXcd kprop, kprop_x, kprop_y, Vgrid, PotentialGrid;

    MatrixData::MetaData meta;

    // inline void RTE_compute_k_pot(MatrixXcd &k,MatrixXcd &wavefctcp,int &t);
    inline double rotatingPotential(int &i, int &j, int &t);
   

    // Variables
    complex<double> h_x, h_y,h_z;

    int SLICE_NUMBER;
    
    
    // some used constants
    
    double pi;
    complex<double> t_RTE;

    // little helper functions for stuff
};

class SplitRot : public SplitStep {
public:
    SplitRot(Options &o) : SplitStep(o) {}
    virtual void timeStep(double delta_t);
};

class SplitTrap : public SplitStep {
public:
    SplitTrap(Options &o) : SplitStep(o) {}
    virtual void timeStep(double delta_t);
};

class SplitFree : public SplitStep {
public:
    SplitFree(Options &o) : SplitStep(o) {}
    virtual void timeStep(double delta_t);
};




#endif // EXP2D_SPLITSTEP_H__