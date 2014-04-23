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
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

class Averages {
public:
	Averages();
	~Averages();
	void saveData(vector<MatrixXcd> &wavefctVec,Options &externalopt);
private:
	vector<ComplexGrid> PsiVec;
	Options opt;

};

class Evaluation {
public:
	
	double Ekin, particle_count;
	ArrayXd number;
	ArrayXd k;
	
	Evaluation() {};
	Evaluation(int avgrid);
	~Evaluation() {};

};

#endif // EXP2D_OBSERVABLES_H__