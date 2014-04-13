#ifndef PLOT_WITH_MGL_H__
#define PLOT_WITH_MGL_H__

#include <iostream>
#include <mgl2/mgl.h>
#include <complex>
#include <math.h>
#include <complexgrid.h>
#include <bh3binaryfile.h>
#include <vector>
#include <omp.h>
#include <cstring>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

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

inline void plotdatatopng(ComplexGrid* &g,Options &opt)
{
	

	int n = opt.grid[1];
	int m = opt.grid[2];

	mglComplex data(n,m);

	int i,j,k;

	// data.Create(n,m);

	// complex<double> data1;

	for(i=0;i<n;i++) for(j=0;j<m;j++)
	{	
		k = i+n*j;
		// data1 = g->at(0,i,j,0);
		data.a[k] = abs2(g->at(0,i,j,0));
	}

	mglGraph gr;

		
		// gr.Light(0,true);
		// gr.Alpha(true);

	string filename = opt.name + ".png";

	gr.SetSize(1800,1800);
	gr.SetQuality(3);
	gr.Title(opt.name.c_str());
	// gr.Alpha(true);



	data.use_abs=false;
	gr.SetRange('x',-opt.min_x,opt.min_x);
	gr.SetRange('y',-opt.min_y,opt.min_y);
	gr.SetRange('z',data);
	gr.SetRange('c',data);

	gr.SubPlot(2,2,0);

	gr.Rotate(40,40);
	gr.Box();
	gr.Axis();
	gr.Surf(data);


	gr.SubPlot(2,2,2);
	gr.Axis();
	gr.Colorbar("_");
	gr.Dens(data);


	data.use_abs=true;
	gr.SetRange('x',-opt.min_x,opt.min_x);
	gr.SetRange('y',-opt.min_y,opt.min_y);
	gr.SetRange('z',data);
	gr.SetRange('c',data);

	gr.SubPlot(2,2,1);

	// gr.Light(true);
	gr.Rotate(40,40);
	gr.Box();
	gr.Axis();

	gr.Surf(data);

	gr.SubPlot(2,2,3);
	gr.Axis();
	gr.Colorbar("_");
	gr.Dens(data);

	gr.WritePNG(filename.c_str(),"ExpandingVortexGas2D",false);

}

inline void plotdatatopngEigen(Eigen::MatrixXcd& mPsi,Options &opt)
{
	

	int n = opt.grid[1];
	int m = opt.grid[2];

	mglComplex data(n,m);

	int i,j,k;

	// data.Create(n,m);

	// complex<double> data1;

	for(i=0;i<n;i++) for(j=0;j<m;j++)
	{	
		k = i+n*j;
		// data1 = g->at(0,i,j,0);
		data.a[k] = abs2(mPsi(i,j));
	}

	mglGraph gr;

		
		// gr.Light(0,true);
		// gr.Alpha(true);

	
	gr.SetRange('x',-opt.min_x,opt.min_x);
	gr.SetRange('y',-opt.min_y,opt.min_y);
	gr.SetSize(1800,1800);
	gr.SetQuality(3);
	gr.Title(opt.name.c_str());
	// gr.Alpha(true);



	// data.use_abs=false;
	// string filename = "PHASE-" + opt.name + ".png";

	// gr.SetRange('z',data);
	// gr.SetRange('c',data);

	// // gr.SubPlot(1,2,0);

	// // gr.Rotate(40,40);
	// // gr.Box();
	// // gr.Axis();
	// // gr.Surf(data);


	// // gr.SubPlot(2,2,2);
	// gr.Axis();
	// gr.Colorbar("_");
	// gr.Dens(data);
	// gr.WritePNG(filename.c_str(),"ExpandingVortexGas2D",false);


	data.use_abs=true;
	string filename = "DENS-" + opt.name + ".png";
	gr.SetRange('z',data);
	gr.SetRange('c',data);

	// gr.SubPlot(1,2,1);

	// // gr.Light(true);
	// gr.Rotate(40,40);
	// gr.Box();
	// gr.Axis();

	// gr.Surf(data);

	// gr.SubPlot(2,2,3);
	gr.Axis();
	gr.Colorbar("_");
	gr.Dens(data);

	gr.WritePNG(filename.c_str(),"ExpandingVortexGas2D",false);

}

inline void plotdatatopngEigenExpanding(Eigen::MatrixXcd& mPsi,vector<double> &ranges,VectorXd &Xexpanding,VectorXd &Yexpanding,Options &opt)
{
	

	int n = opt.grid[1];
	int m = opt.grid[2];



	mglData data(n,m);
	mglData xaxis(n);
	mglData yaxis(m);

	int i,j,k;

	// data.Create(n,m);

	for(i=0;i<n;i++) for(j=0;j<m;j++)
	{	
		k = i+n*j;

		data.a[k] = abs2(mPsi(i,j));		
	}

	for( i = 0; i < n; i++){ xaxis.a[i] = Xexpanding(i); }
	for( j = 0; j < m; j++){ yaxis.a[j] = Yexpanding(j); }



	mglGraph gr;

		
		// gr.Light(0,true);
		// gr.Alpha(true);

	string filename = opt.name + "-DENS.png";

	gr.SetSize(1800,1800);
	gr.SetQuality(3);
	gr.Title(opt.name.c_str());
	gr.SetRange('x',-ranges[0],ranges[0]);
	gr.SetRange('y',-ranges[1],ranges[1]);
	// gr.Alpha(true);



	// data.use_abs=false;

	// gr.SetRange('z',data);
	// gr.SetRange('c',data);

	// gr.SubPlot(2,2,0);

	// gr.Rotate(40,40);
	// gr.Box();
	// gr.Axis();
	// gr.Surf(xaxis,yaxis,data);


	// gr.SubPlot(2,2,2);
	// gr.Axis();
	// gr.Colorbar("_");
	// gr.Dens(xaxis,yaxis,data);


	// data.use_abs=true;
	gr.SetRange('z',data);
	gr.SetRange('c',data);

	// gr.SubPlot(2,2,1);

	// // gr.Light(true);
	// gr.Rotate(40,40);
	// gr.Box();
	// gr.Axis();

	// gr.Surf(xaxis,yaxis,data);

	// gr.SubPlot(2,2,3);
	gr.Axis();
	gr.Colorbar("_");
	gr.Dens(xaxis,yaxis,data);

	gr.WritePNG(filename.c_str(),"ExpandingVortexGas2D",false);

}

#endif // PLOT_WITH_MGL_H__