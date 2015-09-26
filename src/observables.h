#ifndef EXP2D_OBSERVABLE_H__
#define EXP2D_OBSERVABLE_H__


#include <iostream>
#include <complex>
#include <math.h>
#include <vector>
#include <omp.h>
#include <string>
#include <tools.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

class Observables {
        public:
        
        double Ekin, particle_count, healing_length, volume, density, aspectRatio, alpha, r_max, r_min, r_max_phi, r_min_phi, Rx, Ry, n0;
        ArrayXd number;
        ArrayXd k,r;
        ArrayXd angularDensity,radialDensity;
        ArrayXd fixedAspectRatio;
        
        Observables();
        Observables(int avgrid);
    
        Observables operator+ (const Observables &a) const;
        Observables operator- (const Observables &a) const;
        Observables operator* (const Observables &a) const;
    
        Observables operator* (double d) const;
        Observables operator/ (double d) const;
    
        Observables &operator+= (const Observables &a);
        Observables &operator-= (const Observables &a);
        Observables &operator*= (const Observables &a);
        
        Observables &operator/= (double d);    
        Observables &operator*= (double d);
};

typedef struct {
    int32_t n;
    Coordinate<double> c;
    // vector<Vector<double>> velocity;
    // list<Coordinate<int32_t>> points;
    // int32_t num_points;
    // double pair_distance;
    double surroundDens;
    double zeroDensity;
} VortexData;

struct Ellipse {
	Matrix<double, 6,1> coef;
	vector<double> center;
	double major;
	double minor;
	double angle;
};

// struct PathResults {
//     list<VortexData> vlist;
//     vector<double> histogram;
//     vector<double> distance;
// };

inline Observables::Observables() :
        number(),
        k(),
        r(),
        radialDensity(),
        angularDensity(360),
        fixedAspectRatio(360)
{
    Ekin = particle_count = healing_length = volume = density = aspectRatio = alpha = r_max = r_min = r_max_phi = r_min_phi = Rx = Ry = 0.0;
    number.setZero();
    k.setZero();
    r.setZero();
    radialDensity.setZero();
    angularDensity.setZero();
    fixedAspectRatio.setZero();
}



inline Observables::Observables(int avgrid) :
        number(avgrid),
        k(avgrid),
        r(avgrid),
        radialDensity(avgrid),
        angularDensity(360),
        fixedAspectRatio(360)
{
    Ekin = particle_count = healing_length = volume = density = aspectRatio = alpha = r_max = r_min = r_max_phi = r_min_phi = Rx = Ry = n0 = 0.0;
    number.setZero();
    k.setZero();
     r.setZero();
    radialDensity.setZero();
    angularDensity.setZero();
    fixedAspectRatio.setZero();
}

inline Observables Observables::operator+ (const Observables &a) const
{   
    Observables ret(number.size());  

    ret.particle_count = particle_count + a.particle_count;
    ret.healing_length = healing_length + a.healing_length; 
    ret.Ekin = Ekin + a.Ekin;
    ret.aspectRatio = aspectRatio + a.aspectRatio;
    ret.fixedAspectRatio = fixedAspectRatio + a.fixedAspectRatio;
    ret.alpha = alpha + a.alpha;
    ret.r_max = r_max + a.r_max;
    ret.r_min = r_min + a.r_min;
    ret.r_max_phi = r_max_phi + a.r_max_phi;    
    ret.r_min_phi = r_min_phi + a.r_min_phi;
    ret.Rx = Rx + a.Rx;
    ret.Ry = Ry + a.Ry;
    ret.density = density + a.density;
    ret.volume = volume + a.volume;
    ret.angularDensity = angularDensity + a.angularDensity;
    ret.radialDensity = radialDensity + a.radialDensity;
    ret.n0 = n0 + a.n0;


    ret.number = number + a.number; 
    ret.k = k + a.k;
    ret.r = r + a.r;
    
    return ret;
}

inline Observables Observables::operator- (const Observables &a) const
{
    Observables ret(number.size());  

    ret.particle_count = particle_count - a.particle_count;
    ret.healing_length = healing_length - a.healing_length;     
    ret.Ekin = Ekin - a.Ekin;
    ret.aspectRatio = aspectRatio - a.aspectRatio;
    ret.fixedAspectRatio = fixedAspectRatio - a.fixedAspectRatio;
    ret.alpha = alpha - a.alpha;
    ret.r_max = r_max - a.r_max;
    ret.r_min = r_min - a.r_min;
    ret.r_max_phi = r_max_phi - a.r_max_phi; 
    ret.r_min_phi = r_min_phi - a.r_min_phi;
    ret.Rx = Rx - a.Rx;
    ret.Ry = Ry - a.Ry;
    ret.density = density - a.density;
    ret.volume = volume - a.volume;
    ret.angularDensity = angularDensity - a.angularDensity;
    ret.radialDensity = radialDensity - a.radialDensity;
    ret.n0 = n0 - a.n0;

    ret.number = number - a.number; 
    ret.k = k - a.k;
    ret.r = r - a.r;
    
    return ret;
}

inline Observables Observables::operator* (const Observables &a) const
{
    Observables ret(number.size());  

    ret.particle_count = particle_count * a.particle_count;
    ret.healing_length = healing_length * a.healing_length;     
    ret.Ekin = Ekin * a.Ekin;
    ret.aspectRatio = aspectRatio * a.aspectRatio;
    ret.fixedAspectRatio = fixedAspectRatio * a.fixedAspectRatio;
    ret.alpha = alpha * a.alpha;
    ret.r_max = r_max * a.r_max;
    ret.r_min = r_min * a.r_min;
    ret.r_max_phi = r_max_phi * a.r_max_phi; 
    ret.r_min_phi = r_min_phi * a.r_min_phi;
    ret.Rx = Rx * a.Rx;
    ret.Ry = Ry * a.Ry;
    ret.density = density * a.density;
    ret.volume = volume * a.volume;
    ret.angularDensity = angularDensity * a.angularDensity;
    ret.radialDensity = radialDensity * a.radialDensity;
    ret.n0 = n0 * a.n0;

    ret.number = number * a.number; 
    ret.k = k * a.k;
    ret.r = r * a.r;
    
    return ret;
}

inline Observables Observables::operator* (double d) const
{  
    Observables ret(number.size());

    ret.particle_count = particle_count * d;
    ret.healing_length = healing_length * d;    
    ret.Ekin = Ekin * d;
    ret.aspectRatio = aspectRatio * d;
    ret.fixedAspectRatio = fixedAspectRatio * d;
    ret.alpha = alpha * d;
    ret.r_max = r_max * d;
    ret.r_min = r_min * d;
    ret.r_max_phi = r_max_phi * d; 
    ret.r_min_phi = r_min_phi * d;
    ret.Rx = Rx * d;
    ret.Ry = Ry * d;
    ret.density = density * d;
    ret.volume = volume * d;
    ret.angularDensity = angularDensity * d;
    ret.radialDensity = radialDensity * d;
    ret.n0 = n0 * d;

    ret.number = number * d;    
    ret.k = k * d;
    ret.r = r * d;

    return ret;
}

inline Observables Observables::operator/ (double d) const
{  
    Observables ret(number.size());

    ret.particle_count = particle_count / d;
    ret.healing_length = healing_length / d;    
    ret.Ekin = Ekin / d;
    ret.aspectRatio = aspectRatio / d;
    ret.fixedAspectRatio = fixedAspectRatio / d;
    ret.alpha = alpha / d;
    ret.r_max = r_max / d;
    ret.r_min = r_min / d;
    ret.r_max_phi = r_max_phi / d; 
    ret.r_min_phi = r_min_phi / d;
    ret.Rx = Rx / d;
    ret.Ry = Ry / d;
    ret.density = density / d;
    ret.volume = volume / d;
    ret.angularDensity = angularDensity / d;
    ret.radialDensity = radialDensity / d;
    ret.n0 = n0 / d;

    ret.number = number / d;    
    ret.k = k / d;
    ret.r = r / d;
    
    return ret;
}

inline Observables & Observables::operator+= (const Observables &a)
{
    *this = *this + a;
    return *this;
}

inline Observables & Observables::operator-= (const Observables &a)
{
    *this = *this - a;
    return *this;
}

inline Observables & Observables::operator*= (const Observables &a)
{
    *this = *this * a;
    return *this;
}

inline Observables & Observables::operator/= (double d)
{
    *this = *this / d;
    return *this;
}

inline Observables & Observables::operator*= (double d)
{
    *this = *this * d;
    return *this;
}

#endif // EXP2D_OBSERVABLE_H__