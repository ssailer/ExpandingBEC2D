#ifndef EXP2D_OBSERVABLES_H__
#define EXP2D_OBSERVABLES_H__

#include <iostream>
#include <complex>
#include <math.h>
#include <complexgrid.h>
#include <realgrid.h>
#include <bh3binaryfile.h>
#include <vector>
#include <omp.h>
#include <string>
#include <plot_with_mgl.h>
#include <EXP2D_tools.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

class Eval{
public:
	Eval();
	~Eval();

	// wrapperfunctions 
	void saveData(vector<MatrixXcd> &wavefctVec,Options &externalopt,int &external_snapshot_time); // If data comes as a vector(from statistics RTE)
	void saveData(MatrixXcd &wavefct,Options &externalopt,int &external_snapshot_time); // If data comes only as on Matrix(from ITP)
	void evaluateData(); // calculate the observables
	void plotData(); // plot totalResult 
	// public total Result of Evaluation
	Observables totalResult;

		// doing functinos
	Observables evaluate(ComplexGrid data);
	RealGrid findVortices(ComplexGrid data);


private:

	// data savefiles
	vector<ComplexGrid> PsiVec;
	Options opt;
	int snapshot_time;
	vector<RealGrid> vortexLocationMap;


};


#endif // EXP2D_OBSERVABLES_H__