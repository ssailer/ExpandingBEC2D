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

class Averages{
public:
	Averages();
	~Averages();
	class Evaluation {
	public:
	
	double Ekin, particle_count;
	ArrayXd number;
	ArrayXd k;
	
	Evaluation() {};
	Evaluation(int avgrid);
	~Evaluation() {};

	Evaluation operator+ (const Evaluation &a) const;
	Evaluation operator- (const Evaluation &a) const;
	Evaluation operator* (const Evaluation &a) const;

	Evaluation operator* (double d) const;
	Evaluation operator/ (double d) const;

	Evaluation &operator+= (const Evaluation &a);
	Evaluation &operator-= (const Evaluation &a);
	Evaluation &operator*= (const Evaluation &a);
	
	Evaluation &operator/= (double d);

	Evaluation &operator*= (double d);

};
	void saveData(vector<MatrixXcd> &wavefctVec,Options &externalopt,int &external_snapshot_time);
	void evaluateData();
	Evaluation evaluate(ComplexGrid &data);
	void plot(const int &snapshot_time,Evaluation &eval);

private:
	vector<ComplexGrid> PsiVec;
	Options opt;
	int snapshot_time;

};

inline Averages::Evaluation Averages::Evaluation::operator+ (const Evaluation &a) const
{
	Evaluation ret(number.size());	

	ret.particle_count = particle_count + a.particle_count;	
	ret.Ekin = Ekin + a.Ekin;
	ret.number = number + a.number;	
	ret.k = k + a.k;
	
	return ret;
}

inline Averages::Evaluation Averages::Evaluation::operator- (const Evaluation &a) const
{
	Evaluation ret(number.size());	

	ret.particle_count = particle_count - a.particle_count;	
	ret.Ekin = Ekin - a.Ekin;
	ret.number = number - a.number;	
	ret.k = k - a.k;
	
	return ret;
}

inline Averages::Evaluation Averages::Evaluation::operator* (const Evaluation &a) const
{
	Evaluation ret(number.size());	

	ret.particle_count = particle_count * a.particle_count;	
	ret.Ekin = Ekin * a.Ekin;
	ret.number = number * a.number;	
	ret.k = k * a.k;
	
	return ret;
}

inline Averages::Evaluation Averages::Evaluation::operator* (double d) const
{  
	Evaluation ret(number.size());

	ret.particle_count = particle_count * d;	
	ret.Ekin = Ekin * d;
	ret.number = number * d;	
	ret.k = k * d;
}

inline Averages::Evaluation Averages::Evaluation::operator/ (double d) const
{  
	Evaluation ret(number.size());

	ret.particle_count = particle_count / d;	
	ret.Ekin = Ekin / d;
	ret.number = number / d;	
	ret.k = k / d;
}



inline Averages::Evaluation & Averages::Evaluation::operator+= (const Evaluation &a)
{
	*this = *this + a;
	return *this;
}

inline Averages::Evaluation & Averages::Evaluation::operator-= (const Evaluation &a)
{
	*this = *this - a;
	return *this;
}

inline Averages::Evaluation & Averages::Evaluation::operator*= (const Evaluation &a)
{
	*this = *this * a;
	return *this;
}

inline Averages::Evaluation & Averages::Evaluation::operator/= (double d)
{
	*this = *this / d;
	return *this;
}

inline Averages::Evaluation & Averages::Evaluation::operator*= (double d)
{
	*this = *this * d;
	return *this;
}

#endif // EXP2D_OBSERVABLES_H__