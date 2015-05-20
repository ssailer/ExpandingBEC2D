#ifndef EXP2D_RK4_H__
#define EXP2D_RK4_H__

#include <iostream>
#include <complex>
#include <omp.h>
#include <EXP2D_MatrixData.h>
#include <EXP2D_constants.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

class RungeKutta
{
public:
	RungeKutta(/*MatrixData* &d,*/ Options &extOpt);
    void assignMatrixData(MatrixData* &d);
    void setVariables();
	void timeStep(double delta_T);
// private:
	MatrixData* w;
	MatrixXcd wavefctcp, k0, k1, k2, k3;
	Options opt;
	double absTime;

    vector<int32_t> frontx, endx;
    int32_t threads, subx, suby, partx;

	// FIXME REMOVE EVENTUALLY
	vector<double> x_axis,y_axis;
    VectorXcd X,Y;
    MatrixXcd Xmatrix,Ymatrix,PotentialGrid;
    // END FIXME

	vector<complex<double>> laplacian_coefficient_x,laplacian_coefficient_y,gradient_coefficient_x,gradient_coefficient_y;

	virtual void computeCoefficients(double &delta_T) = 0;
	virtual void singleK(MatrixXcd &k,int32_t &front, int32_t &end,int32_t &subx,int32_t &suby, int &t) = 0;

	void MSDBoundaries(MatrixXcd &U,MatrixXcd &Ut);

    inline complex<double> lambda_x(complex<double> absTime){
    	return sqrt(one+opt.exp_factor*opt.dispersion_x*opt.dispersion_x*absTime*absTime);
    }

   	inline complex<double> lambda_x_dot(complex<double> absTime){
    	return (opt.exp_factor*opt.dispersion_x*opt.dispersion_x*absTime/sqrt(one+opt.exp_factor*opt.dispersion_x*opt.dispersion_x*absTime*absTime));
   	}

    inline complex<double> lambda_y(complex<double> absTime){
   		return sqrt(one+opt.exp_factor*opt.dispersion_y*opt.dispersion_y*absTime*absTime);
    }

    inline complex<double> lambda_y_dot(complex<double> absTime){
    	return (opt.exp_factor*opt.dispersion_y*opt.dispersion_y*absTime/sqrt(one+opt.exp_factor*opt.dispersion_y*opt.dispersion_y*absTime*absTime));
    }

};

class Expansion : public RungeKutta {
public:
    Expansion(Options &opt) : RungeKutta(opt) {}
private:
    virtual void singleK(MatrixXcd &k, int32_t &front, int32_t &end,int32_t &subx,int32_t &suby, int &t);
    virtual void computeCoefficients(double &delta_T);
};

class RotatingTrap : public RungeKutta {
public:
    RotatingTrap(Options &opt) : RungeKutta(opt) {}
private:
    virtual void singleK(MatrixXcd &k, int32_t &front, int32_t &end,int32_t &subx,int32_t &suby, int &t);
    virtual void computeCoefficients(double &delta_T);
};

#endif // EXP2D_RK4_H__