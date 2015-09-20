#ifndef PLOT_WITH_MGL_H__
#define PLOT_WITH_MGL_H__

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

#include <dlib/optimization.h>




#define dual std::complex<double>


using namespace std;

class mglComplex : public mglDataA
{
public:
  long nx;      ///< number of points in 1st dimensions ('x' dimension)
  long ny;      ///< number of points in 2nd dimensions ('y' dimension)
  long nz;      ///< number of points in 3d dimensions ('z' dimension)
  dual *a;      ///< data array
  bool use_abs; ///< flag to use abs() or arg()

  inline mglComplex(long xx=1,long yy=1,long zz=1)
  { a=0;  use_abs=true; Create(xx,yy,zz); }
  virtual ~mglComplex()  { if(a)  delete []a; }

  /// Get sizes
  inline long GetNx() const { return nx;  }
  inline long GetNy() const { return ny;  }
  inline long GetNz() const { return nz;  }
  /// Create or recreate the array with specified size and fill it by zero
  inline void Create(long mx,long my=1,long mz=1)
  { nx=mx;  ny=my;  nz=mz;  if(a) delete []a;
  a = new dual[nx*ny*nz]; }
  /// Get maximal value of the data
  inline mreal Maximal() const  { return mgl_data_max(this);  }
  /// Get minimal value of the data
  inline mreal Minimal() const  { return mgl_data_min(this);  }

protected:
  inline mreal v(long i,long j=0,long k=0) const
  { return use_abs ? abs(a[i+nx*(j+ny*k)]) : arg(a[i+nx*(j+ny*k)]);  }
  inline mreal vthr(long i) const
  { return use_abs ? abs(a[i]) : arg(a[i]);  }
  inline mreal dvx(long i,long j=0,long k=0) const
  { register long i0=i+nx*(j+ny*k);
    std::complex<double> res=i>0? (i<nx-1? (a[i0+1]-a[i0-1])/2.:a[i0]-a[i0-1]) : a[i0+1]-a[i0];
    return use_abs? abs(res) : arg(res);  }
  inline mreal dvy(long i,long j=0,long k=0) const
  { register long i0=i+nx*(j+ny*k);
    std::complex<double> res=j>0? (j<ny-1? (a[i0+nx]-a[i0-nx])/2.:a[i0]-a[i0-nx]) : a[i0+nx]-a[i0];
    return use_abs? abs(res) : arg(res);  }
  inline mreal dvz(long i,long j=0,long k=0) const
  { register long i0=i+nx*(j+ny*k), n=nx*ny;
    std::complex<double> res=k>0? (k<nz-1? (a[i0+n]-a[i0-n])/2.:a[i0]-a[i0-n]) : a[i0+n]-a[i0];
    return use_abs? abs(res) : arg(res);  }
};

void plotSpectrum(string name,string title,Observables& eval);
void plotRadialDensity(string name,string title, Observables &ares);
void plotPairDistance(string name,string title,list<VortexData> vlist);
void plotVortexList(string name,string title,const RealGrid &phase,list<VortexData> &vlist,Options &opt);
void plotContour(string filename, MatrixXcd &data, std::unordered_set<Coordinate<int32_t>> &contour, Options &opt);

void plotDataToPng(string filename,string title,ComplexGrid* &g,Options &opt);
void plotDataToPng(string filename,string title,ComplexGrid &g,Options &opt);
void plotDataToPng(string filename,string title,RealGrid g,Options &opt);
void plotDataToPng(string filename,MatrixXi &g,Options &opt);
void plotDataToPngExpanding(string filename,string title,ComplexGrid &g,Options &opt);
void plotWithExpandingFrame(string filename,string title, ComplexGrid &Psi,vector<double> &ranges,vector<double> &Xexpanding,vector<double> &Yexpanding,Options &opt);
void plotVector(string filename,string title,vector<double> v,Options &opt);
void plotVector(string filename,vector<double> v,vector<double> w,Options &opt);
void plotVector(string filename,string title,ArrayXd v,ArrayXi w);


void plotContourSurround(string name, RealGrid &Psi, std::unordered_set<Coordinate<int32_t>> &contour, Options &opt);
void plotAngularDensity(string filename,vector<double> phi,vector<double> density,Options &opt);
void plotDataToPngEigen(string filename, Eigen::MatrixXcd& mPsi,Options opt);
void plotDataToPngEigenExpanding(string filename, Eigen::MatrixXcd& mPsi,vector<double> &ranges,Eigen::VectorXd &Xexpanding,Eigen::VectorXd &Yexpanding,Options opt);


void plotPair(std::vector<std::pair<dlib::matrix<double,2,1>, double> > data_samples);
void plotGauss(dlib::matrix<double,3,1> params,int n,double coordinate_axis);
void plotPairAndGauss(std::vector<std::pair<dlib::matrix<double,2,1>, double> > data_samples, dlib::matrix<double,3,1> params,int n, double coordinate_axis);


#endif // PLOT_WITH_MGL_H__