#ifndef EXP2D_OBSERVABLES_H__
#define EXP2D_OBSERVABLES_H__

#include <iostream>
#include <complex>
#include <math.h>
#include <complexgrid.h>
#include <bh3binaryfile.h>
#include <vector>
#include <omp.h>
#include <string>
#include <plot_with_mgl.h>
#include <EXP2D_tools.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

class Averages{
public:
	Averages();
	~Averages();

	// wrapperfunctions 
	void saveData(vector<MatrixXcd> &wavefctVec,Options &externalopt,int &external_snapshot_time); // If data comes as a vector(from statistics RTE)
	void saveData(MatrixXcd &wavefct,Options &externalopt,int &external_snapshot_time); // If data comes only as on Matrix(from ITP)
	void evaluateData(); // calculate the observables
	void plotTotalResult(); // plot totalResult 
	// public total Result of Evaluation
	Evaluation totalResult;

		// doing functinos
	Evaluation evaluate(ComplexGrid &data);
	void plot(const int &snapshot_time,Evaluation &eval);


private:

	// data savefiles
	vector<ComplexGrid> PsiVec;
	Options opt;
	int snapshot_time;


};


#endif // EXP2D_OBSERVABLES_H__