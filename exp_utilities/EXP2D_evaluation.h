#ifndef EXP2D_EVALUATION_H__
#define EXP2D_EVALUATION_H__

#define EIGEN_FFTW_DEFAULT

#include <iostream>
#include <sys/stat.h>
#include <complex>
#include <cmath>
#include <numeric>
#include <stack>
#include <algorithm>
#include <complexgrid.h>
#include <realgrid.h>
#include <bh3binaryfile.h>
#include <coordinate.h>
#include <vector>
#include <unordered_set>
#include <omp.h>
#include <string>
#include <sstream>
#include <EXP2D_Contour.h>
#include <plot_with_mgl.h>
#include <EXP2D_tools.h>
#include <EXP2D_MatrixData.h>
#include <EXP2D_observables.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;



class Eval{
public:
	Eval(MatrixData d,Options o);
	// Eval();
	// ~Eval();

	// wrapperfunctions 
	// void saveData(vector<MatrixXcd> &wavefctVec,Options &external_opt,int external_snapshot_time,string external_runname); // If data comes as a vector of matrices (from statistics RTE)
	// void saveData(MatrixXcd &wavefct,Options &external_opt,int external_snapshot_time,string external_runname); // If data comes only as a Matrix (from ITP)
	// void saveData2DSlice(vector<ComplexGrid> &wavefctVec, Options & external_opt, int external_snapshot_time, string external_runname, int sliceNumber); // if data comes as a vector of ComplexGrids, just eval a sclice of the 3D data.
	// void saveDataFromEval(Options &external_opt,int &external_snapshot_time,string &external_runname,vector<Eval> &extEval);
	void process(); // calculate the observables
	void save();
	// void evaluateDataITP();
	// void plot(); // plot Results
	bool checkResizeCondition();
	int getVortexNumber();
	bool checkResizeCondition(vector<int> &edges);


	// Observables.h
	Observables totalResult;
	// vector<PathResults> pres;
	vector<list<VortexData>> vlist;
	vector<c_set> contour;
	vector<MatrixXi> densityLocationMap;

	Ellipse ellipse;

	MatrixData data;
	Options opt;

private:

	typedef struct{
		int32_t length;
		Coordinate<int32_t> start;
		Coordinate<int32_t> stop;
	} lineData;

	typedef struct{
		Coordinate<int32_t> c;
		double phi;
		double r;
	} contourData;

	
	// data savefiles

	// RealGrid phase, zeros;
	MatrixXd phase;
	MatrixXd density;

	// vector<ComplexGrid> PsiVec;
	
	int snapshot_time;

	// vector<RealGrid> densityLocationMap;
	

	vector<vector<Coordinate<int32_t>>> densityCoordinates;
	vector<double> x_dist,y_dist,x_dist_grad,y_dist_grad;
	vector<int> densityCounter;

	void CombinedEval();
	void CombinedSpectrum();

	// doing functinos
	Observables calculator(MatrixXcd DATA,int sampleindex);
	// Observables calculatorITP(ComplexGrid data,int sampleindex);
	void aspectRatio(Observables &obs, int &sampleindex);
	void getVortices(MatrixXcd &DATA, vector<Coordinate<int32_t>> &densityCoordinates,list<VortexData> &vlist);
	// void getDensity(ComplexGrid &data, RealGrid &densityLocationMap, vector<Coordinate<int32_t>> &densityCoordinates,int &densityCounter);
	void getDensity();

	int checkSum(MatrixXi &d,int &i, int &j);
	void erosion(MatrixXi &d);
	void dilation(MatrixXi &d);

	void floodFill(MatrixXi &dens);
	void fillHoles(MatrixXi &dens);

	vector<double> polarDensity();
	vector<int> findMajorMinor();
	Ellipse fitEllipse(c_set &Data);
	c_set generateContour(Ellipse &ellipse);
	

	int get_phase_jump(const Coordinate<int32_t> &c, const Vector<int32_t> &v);
	void findVortices(vector<Coordinate<int32_t>> &densityCoordinates, list<VortexData> &vlist);

	inline double norm(Coordinate<double> &a, Coordinate<double> &b, double &h_x, double &h_y);
	// inline void pairDistanceHistogram(list<Vo &pres, double &distance, double &coordDistance);
	// void getVortexDistance(PathResults &pres);
	void calc_fields(MatrixXcd &DATA, Options &opt);
	// void checkEdges();

	// Contour Tracking Algorithm

	// c_set trackContour(const RealGrid &data);
	
	// inline Coordinate<int32_t> nextClockwise(Coordinate<int32_t> &s, int32_t &direction);
	// inline void setDirection(int32_t &direction);
	// void findInitialP(RealGrid &data,Coordinate<int32_t> &p,Coordinate<int32_t> &s, Coordinate<int32_t> *initial);
	// void findMostRightP(c_set &contour, Coordinate<int32_t> &p);


	// Helpermethods to access radial vectors:
	void checkNextAngles(vector<double> &r, int &i);
	void cyclicAssignment(vector<double> &r, int i, double rvalue);
	double cyclicReadout(vector<double> &r, int i);

};





#endif // EXP2D_EVALUATION_H__