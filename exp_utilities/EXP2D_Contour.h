#ifndef EXP2D_CONTOUR_H__
#define EXP2D_CONTOUR_H__

#include <iostream>
#include <vector>
#include <unordered_set>

// #include <realgrid.h>
#include <coordinate.h>
#include <EXP2D_tools.h>
#include <plot_with_mgl.h>
#include <EXP2D_MatrixData.h>
#include <eigen3/Eigen/Dense>


using namespace std;

typedef std::unordered_set<Coordinate<int32_t>,Hash> c_set;

class Contour{
public:
	Contour(Options &opt);
	c_set trackContour(MatrixXi &data);
private:
	Options opt;
	Vector<int32_t> v_left,v_right, v_up, v_down;
	inline void findMostRightP(c_set &contour, Coordinate<int32_t> &p);
	inline void findInitialP(MatrixXi &data,Coordinate<int32_t> &p,Coordinate<int32_t> &s);
	inline void findSecondP(MatrixXi &data,Coordinate<int32_t> &p,Coordinate<int32_t> &s);
	inline Coordinate<int32_t> nextClockwise(Coordinate<int32_t> &s, int32_t &direction);
	inline void setDirection(int32_t &direction);
	
};

#endif // EXP2D_CONTOUR_H__
