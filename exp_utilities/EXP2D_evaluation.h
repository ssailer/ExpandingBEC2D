#ifndef EXP2D_EVALUATION_H__
#define EXP2D_EVALUATION_H__

#include <iostream>
#include <complex>
#include <cmath>
#include <numeric>
#include <complexgrid.h>
#include <realgrid.h>
#include <bh3binaryfile.h>
#include <vector>
#include <omp.h>
#include <string>
#include <plot_with_mgl.h>
#include <EXP2D_tools.h>
#include <EXP2D_observables.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

class Eval{
public:
	Eval();
	~Eval();

	// wrapperfunctions 
	void saveData(vector<MatrixXcd> &wavefctVec,Options &externalopt,int &external_snapshot_time,string runname_external); // If data comes as a vector of matrices (from statistics RTE)
	void saveData(MatrixXcd &wavefct,Options &externalopt,int &external_snapshot_time,string runname_external); // If data comes only as a Matrix (from ITP)
	void evaluateData(); // calculate the observables
	void plotData(); // plot Results

	// public total Result of Evaluation
	Observables totalResult;


private:

	// data savefiles
	string runname;
	vector<ComplexGrid> PsiVec;
	Options opt;
	int snapshot_time;
	vector<RealGrid> vortexLocationMap;
	vector<RealGrid> densityLocationMap;
	vector<vector<Coordinate<int32_t>>> vortexCoordinates;
	vector<vector<Coordinate<int32_t>>> densityCoordinates;
	vector<double> x_dist;
	vector<double> y_dist;

	// doing functinos
	Observables calculator(ComplexGrid data,int sampleindex);
	void findVortices(ComplexGrid data, RealGrid &vortexLocationMap_local, vector<Coordinate<int32_t>> &vortexCoordinates);
	void findDensity(ComplexGrid data, RealGrid &densityLocationMap_local, vector<Coordinate<int32_t>> &densityCoordinates);




};


#endif // EXP2D_EVALUATION_H__