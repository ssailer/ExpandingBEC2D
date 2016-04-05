#ifndef EXP2D_PLOTTER_H__
#define EXP2D_PLOTTER_H__

#include <iostream>
#include <complex>
#include <eigen3/Eigen/Dense>
#include <mgl2/mgl.h>
#include <math.h>
// #include <complexgrid.h>
// #include <realgrid.h>
// #include <bh3binaryfile.h>

#include <vector>
#include <unordered_set>
#include <omp.h>
#include <cstring>


#include "coordinate.h"
#include "tools.h"
#include "observables.h"
#include "evaluation.h"

using namespace std;

class Plotter
{
public:
	Plotter(shared_ptr<Eval> e, Options &o);
	// ~Plotter();

	void plotEval();
	void control();
	void spectrum();
	void alphas();
	void contour();
	void densityMap();
	void vortices();
	void prepareData();
	void combinedControl();

	void writeTexData(string filename,vector<double> x, vector<double> y);
private:


	shared_ptr<Eval> eval;
	Options opt;

	double xrange,yrange;
	string stepsString, dirname;
	string title;

	mglData density, phase, densitymap, contour_x, contour_y, vortex_x, vortex_y, cover_x, cover_y, k, number, ableitung;
	float kmin, kmax, numbermin, numbermax;
	mglPoint major_1, minor_1,major_2, minor_2, origin;
	mglPoint reg_1, reg_2;


	
};


#endif // EXP2D_PLOTTER_H__