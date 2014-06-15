#ifndef EXP2D_EVALUATION_H__
#define EXP2D_EVALUATION_H__

#include <iostream>
#include <sys/stat.h>
#include <complex>
#include <cmath>
#include <numeric>
#include <algorithm>
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



	PathResults pres;

	


private:

	// data savefiles

	RealGrid *phase, *zeros;
	string runname;
	vector<ComplexGrid> PsiVec;
	Options opt;
	int snapshot_time;
	// vector<RealGrid> vortexLocationMap;
	vector<RealGrid> densityLocationMap;
	// vector<vector<Coordinate<int32_t>>> vortexCoordinates;
	vector<vector<Coordinate<int32_t>>> densityCoordinates;
	vector<double> x_dist,y_dist,x_dist_grad,y_dist_grad;
	int densityCounter;

	// doing functinos
	Observables calculator(ComplexGrid data,int sampleindex);
	void getVortices(ComplexGrid data, vector<Coordinate<int32_t>> &densityCoordinates);
	void getDensity(ComplexGrid data, RealGrid &densityLocationMap_local, vector<Coordinate<int32_t>> &densityCoordinates);

	int get_phase_jump(const Coordinate<int32_t> &c, const Vector<int32_t> &v, const RealGrid *phase);
	void find_vortices(const RealGrid *phase, vector<Coordinate<int32_t>> &densityCoordinates, list<VortexData> &vlist);
	void calc_fields(const ComplexGrid &data);




};


#endif // EXP2D_EVALUATION_H__