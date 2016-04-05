#ifndef EXP2D_LMFITTER_H__
#define EXP2D_LMFITTER_H__

#include <iostream>
#include <sstream>

#include "matrixdata.h"
#include "tools.h"
#include "plot_with_mgl.h"

#include <eigen3/Eigen/Dense>
#include <dlib/optimization.h>


using namespace std;
using namespace Eigen;

typedef dlib::matrix<double,2,1> input_vector;
typedef dlib::matrix<double,4,1> parameter_vector;

class lmfitter
{

public:
	lmfitter(MatrixXd dens, MatrixData::MetaData& m) : density(dens) , meta(m) {};

	vector<double> optimize();

private:
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