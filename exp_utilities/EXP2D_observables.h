#ifndef EXP2D_OBSERVABLE_H__
#define EXP2D_OBSERVABLE_H__


#include <iostream>
#include <complex>
#include <math.h>
#include <vector>
#include <omp.h>
#include <string>
#include <EXP2D_tools.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

class Observables {
        public:
        
        double Ekin, particle_count, healing_length, volume, density, aspectRatio, aspectRatioAngle, r_max, r_min, r_max_phi, r_min_phi;
        ArrayXd number;
        ArrayXd k;
        ArrayXd angularDensity;
        
        Observables() {};
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
    Coordinate<double> x;
    vector<Vector<double>> velocity;
    list<Coordinate<int32_t>> points;
    int32_t num_points;
    double pair_distance;
    double surroundDens;
    double zeroDensity;
} VortexData;

struct PathResults {
    list<VortexData> vlist;
};





inline Observables::Observables(int avgrid) :
        number(avgrid),
        k(avgrid),
        angularDensity(360)
{
    Ekin = particle_count = healing_length = volume = density = aspectRatio = r_max = r_min = r_max_phi = r_min_phi = 0.0;
    number.setZero();
    k.setZero();
    angularDensity.setZero();
}

inline Observables Observables::operator+ (const Observables &a) const
{
    Observables ret(number.size());  

    ret.particle_count = particle_count + a.particle_count;
    ret.healing_length = healing_length + a.healing_length; 
    ret.Ekin = Ekin + a.Ekin;
    ret.aspectRatio = aspectRatio + a.aspectRatio;
    ret.aspectRatioAngle = aspectRatioAngle + a.aspectRatioAngle;
    ret.r_max = r_max + a.r_max;
    ret.r_min = r_min + a.r_min;
    ret.r_max_phi = r_max_phi + a.r_max_phi;    
    ret.r_min_phi = r_min_phi + a.r_min_phi;
    ret.density = density + a.density;
    ret.number = number + a.number; 
    ret.k = k + a.k;
    ret.volume = volume + a.volume;
    ret.angularDensity = angularDensity + a.angularDensity;
    
    return ret;
}

inline Observables Observables::operator- (const Observables &a) const
{
    Observables ret(number.size());  

    ret.particle_count = particle_count - a.particle_count;
    ret.healing_length = healing_length - a.healing_length;     
    ret.Ekin = Ekin - a.Ekin;
    ret.aspectRatio = aspectRatio - a.aspectRatio;
    ret.aspectRatioAngle = aspectRatioAngle - a.aspectRatioAngle;
    ret.r_max = r_max - a.r_max;
    ret.r_min = r_min - a.r_min;
    ret.r_max_phi = r_max_phi - a.r_max_phi; 
    ret.r_min_phi = r_min_phi - a.r_min_phi;
    ret.density = density - a.density;
    ret.number = number - a.number; 
    ret.k = k - a.k;
    ret.volume = volume - a.volume;
    ret.angularDensity = angularDensity - a.angularDensity;
    
    return ret;
}

inline Observables Observables::operator* (const Observables &a) const
{
    Observables ret(number.size());  

    ret.particle_count = particle_count * a.particle_count;
    ret.healing_length = healing_length * a.healing_length;     
    ret.Ekin = Ekin * a.Ekin;
    ret.aspectRatio = aspectRatio * a.aspectRatio;
    ret.aspectRatioAngle = aspectRatioAngle * a.aspectRatioAngle;
    ret.r_max = r_max * a.r_max;
    ret.r_min = r_min * a.r_min;
    ret.r_max_phi = r_max_phi * a.r_max_phi; 
    ret.r_min_phi = r_min_phi * a.r_min_phi;
    ret.density = density * a.density;
    ret.number = number * a.number; 
    ret.k = k * a.k;
    ret.volume = volume * a.volume;
    ret.angularDensity = angularDensity * a.angularDensity;
    
    return ret;
}

inline Observables Observables::operator* (double d) const
{  
    Observables ret(number.size());

    ret.particle_count = particle_count * d;
    ret.healing_length = healing_length * d;    
    ret.Ekin = Ekin * d;
    ret.aspectRatio = aspectRatio * d;
    ret.aspectRatioAngle = aspectRatioAngle * d;
    ret.r_max = r_max * d;
    ret.r_min = r_min * d;
    ret.r_max_phi = r_max_phi * d; 
    ret.r_min_phi = r_min_phi * d;
    ret.density = density * d;
    ret.number = number * d;    
    ret.k = k * d;
    ret.volume = volume * d;
    ret.angularDensity = angularDensity * d;

    return ret;
}

inline Observables Observables::operator/ (double d) const
{  
    Observables ret(number.size());

    ret.particle_count = particle_count / d;
    ret.healing_length = healing_length / d;    
    ret.Ekin = Ekin / d;
    ret.aspectRatio = aspectRatio / d;
    ret.aspectRatioAngle = aspectRatioAngle / d;
    ret.r_max = r_max / d;
    ret.r_min = r_min / d;
    ret.r_max_phi = r_max_phi / d; 
    ret.r_min_phi = r_min_phi / d;
    ret.density = density / d;
    ret.number = number / d;    
    ret.k = k / d;
    ret.volume = volume / d;
    ret.angularDensity = angularDensity / d;
    
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