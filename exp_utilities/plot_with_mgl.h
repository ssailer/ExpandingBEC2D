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
		data.a[k] = g->at(0,i,j,0);
	}

	mglGraph gr;

		
		// gr.Light(0,true);
		// gr.Alpha(true);

	string filename = "OBDM_" + opt.name + ".png";

	gr.SetSize(1600,800);
	gr.SetQuality(3);
	gr.Title(opt.name.c_str());
	// gr.Alpha(true);



	data.use_abs=false;
	gr.SubPlot(2,1,0);
	gr.SetRange('x',-opt.min_x,opt.min_x);
	gr.SetRange('y',-opt.min_y,opt.min_y);
	gr.SetRange('c',data);	
	
	gr.Axis();
	gr.Colorbar("<");
	gr.Dens(data);
	

	data.use_abs=true;
	gr.SubPlot(2,1,1);
	gr.SetRange('x',-opt.min_x,opt.min_x);
	gr.SetRange('y',-opt.min_y,opt.min_y);
	gr.SetRange('c',data);	
	
	gr.Light(true);
	gr.Rotate(60,40);
	gr.Box();
	gr.SetRange('z',data);
	gr.Surf(data);
	gr.Axis();
	gr.Colorbar("<");
	// gr.Dens(data);
	

	


	gr.WritePNG(filename.c_str());

}

#endif // PLOT_WITH_MGL_H__