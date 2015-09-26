#ifndef EXP2D_LMFITTER_H__
#define EXP2D_LMFITTER_H__

#include <iostream>
// #include <sys/stat.h>
// #include <complex>
// #include <cmath>
// #include <numeric>
// #include <stack>
// #include <algorithm>
// #include <complexgrid.h>
// #include <realgrid.h>
// #include <bh3binaryfile.h>
// #include <coordinate.h>
// #include <vector>
// #include <unordered_set>
// #include <omp.h>
// #include <string>
// #include <sstream>
// #include <EXP2D_Contour.h>
// #include <plot_with_mgl.h>
// #include <EXP2D_tools.h>
#include <EXP2D_MatrixData.h>
// #include <EXP2D_observables.h>

// #include <vector>
#include <eigen3/Eigen/Dense>
// #include <eigen3/unsupported/Eigen/NonLinearOptimization>
// #include <eigen3/unsupported/Eigen/NumericalDiff>



#include <plot_with_mgl.h>

#include <dlib/optimization.h>


// #define PARAMETER_SIZE 4

// #include <gsl/gsl_sf_zeta.h>

// #include "TF2.h"
// #include "TH2.h"
// #include "TMath.h"
// #include "TCanvas.h"

using namespace std;
using namespace Eigen;
// using namespace dlib;

typedef dlib::matrix<double,2,1> input_vector;
typedef dlib::matrix<double,4,1> parameter_vector;

class lmfitter
{

public:
	lmfitter(MatrixXd dens, MatrixData::MetaData& m) : density(dens) , meta(m) {};





	// int nbins = density.cols();

	// Double_t g2(Double_t *x, Double_t *par);
	// Double_t fun2(Double_t *x, Double_t *par);
	vector<double> optimize();

private:
	// int x = density.cols()/2;
	// int y = density.rows()/2;
	MatrixXd density;
	MatrixData::MetaData meta;

	parameter_vector set_initial_parameters();
	void plotQuerschnitte(const string& str, const parameter_vector& params_tf);


	static double tf_dist ( const input_vector& input, const parameter_vector& params );
	static double gauss_dist ( const input_vector& input, const parameter_vector& params );

	static double tf_residual ( const std::pair<input_vector, double>& data, const parameter_vector& params );
	static double gauss_residual ( const std::pair<input_vector, double>& data, const parameter_vector& params );
	
	static parameter_vector tf_residual_derivative ( const std::pair<input_vector, double>& data, const parameter_vector& params );
	static parameter_vector gauss_residual_derivative ( const std::pair<input_vector, double>& data, const parameter_vector& params );

	

};

#endif