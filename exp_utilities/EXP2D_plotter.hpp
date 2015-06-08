#ifndef EXP2D_PLOTTER_H__
#define EXP2D_PLOTTER_H__

#include <iostream>
#include <complex>
#include <eigen3/Eigen/Dense>
#include <mgl2/mgl.h>
#include <math.h>
#include <complexgrid.h>
#include <realgrid.h>
#include <bh3binaryfile.h>
#include <coordinate.h>
#include <vector>
#include <unordered_set>
#include <omp.h>
#include <cstring>
#include <EXP2D_tools.h>
#include <EXP2D_observables.h>
#include <EXP2D_evaluation.h>

using namespace std;

class Plotter
{
public:
	Plotter(Eval &e, Options &o);
	// ~Plotter();

	void plotEval();
private:
	void control();
	void spectrum();
	void contour();
	void densityMap();
	void vortices();
	void prepareData();
	void combinedControl();

	Eval eval;
	Options opt;

	double xrange,yrange;
	string stepsString, dirname;
	string title;

	mglData density, phase, densitymap, contour_x, contour_y, vortex_x, vortex_y, k, number;
	float kmin, kmax, numbermin, numbermax;
	mglPoint major_1, minor_1,major_2, minor_2, origin;


	
};


#endif // EXP2D_PLOTTER_H__