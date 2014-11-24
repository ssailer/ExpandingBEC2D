#ifndef EXP2D_ITP_H__
#define EXP2D_ITP_H__

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
#include <EXP2D_evaluation.h>
#include <plot_with_mgl.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

class ITP
{
public:
    ITP();
    ITP(MatrixXcd &wavedata,const Options &opt);
    ~ITP();
    MatrixXcd result();

    void setOptions(Options &externaloptions);
    void RunSetup();

    void propagateToGroundState(string runname);
    void formVortices(string runname);
    void findVortices(string runname);

    // StoragePointer for the wavefunction
    ComplexGrid* pPsi;
    // Storage Variable for the runs
    MatrixXcd wavefct;

        // Coordinates
    vector<double> x_axis,y_axis;
        VectorXcd X,Y;

    void CopyComplexGridToEigen();
    void CopyEigenToComplexGrid();

        // Plotting and progress functions 
    void cli(string name,int counter_state, int counter_max, double start);
    void cli_groundState(string name, double start,int state,Observables totalResult);

        // internal RunOptions, use setOptions(Options) to update from the outside
    Options opt;


private:
    void plot(const string name);

    inline void ITP_compute_k(MatrixXcd &k,MatrixXcd &wavefctcp);
    void ITP_compute_k_parallel(MatrixXcd &k, MatrixXcd &wavefctcp);

    void ComputeDeltaPsi(MatrixXcd &wavefct, MatrixXcd &wavefctcp);
    void singleK(MatrixXcd &k, MatrixXcd &wavefctcp, int32_t &front, int32_t &end,int32_t &subx,int32_t & suby);

    inline double rotatingPotential(int i, int j, int angle);

    inline void rescale(MatrixXcd &wavefct); 

    // Variables
    complex<double> h_x, h_y;
    complex<double> itp_laplacian_x;
    complex<double> itp_laplacian_y;

    Matrix<std::complex<double>,Dynamic,Dynamic,ColMajor> wavefctcpX;
    Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor> wavefctcpY;
    MatrixXcd PotentialGrid,wavefctcp,k0,k1,k2,k3;

    double pi, scaleFactor;
    complex<double>  zero,half,one,two,four,six,i_unit;
    complex<double> t_ITP;

};

#endif // EXP2D_ITP_H__