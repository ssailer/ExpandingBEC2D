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
#include <unordered_set>
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
	void getVortices(const ComplexGrid &data, vector<Coordinate<int32_t>> &densityCoordinates);
	void getDensity(const ComplexGrid &data, RealGrid &densityLocationMap_local, vector<Coordinate<int32_t>> &densityCoordinates);
	

	int get_phase_jump(const Coordinate<int32_t> &c, const Vector<int32_t> &v, const RealGrid *phase);
	void find_vortices(const RealGrid *phase, const RealGrid *zeros, vector<Coordinate<int32_t>> &densityCoordinates, list<VortexData> &vlist);
	void calc_fields(const ComplexGrid &data, Options &opt);

	// Contour Tracking Algorithm
	std::unordered_set<Coordinate<int32_t>> trackContour(const ComplexGrid &data, const Options &opt);
	Vector<int32_t> v_left,v_right, v_up, v_down;
	inline Coordinate<int32_t> nextClockwise(Coordinate<int32_t> &s, int32_t &direction);
	inline void setDirection(int32_t &direction);
	void findInitialP(Coordinate<int32_t> &p,Coordinate<int32_t> &s, Coordinate<int32_t> *initial,const Options &opt);


};



#endif // EXP2D_EVALUATION_H__