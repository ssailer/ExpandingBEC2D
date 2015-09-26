// File Name:     ComplexComplexGridTimeVariable2.h
// Author:        Matthias Kronenwett and Mr. B
// Email Address: m.kronenwett@thphys.uni-heidelberg.de
//						b.nowak@thphys.uni-heidelberg.de
// Description:   This is the header file for the class ComplexGridTimeVariable whose
//                values are gsl_complexs. An object is declared as follows:
//                ComplexGridTimeVariable the_object(Gridzize
//                                         Gridsize,
//                                         Gridsize,
//                                         2);
//                where Gridsize is the number of grid points.
//                The implementation for the class ComplexGridTimeVariable is in the
//                implementation file ComplexGridTimeVariable.cpp.
// Last Change:   Feb 16, 2010

#ifndef REALGRID_H__
#define REALGRID_H__
#include <iostream>
#include <stdint.h>
#include <cmath>
#include <complex>
#include <fftw3.h>


using namespace std;

#include "coordinate.h"

class CudaRealGrid;

class RealGrid
{
	public:
		RealGrid(uint32_t int_dim, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid);
		RealGrid();
		RealGrid(const RealGrid& c);
		
		~RealGrid();
		
		RealGrid& operator =(const RealGrid& c);
		friend CudaRealGrid &copyHostToDevice_as_complex (CudaRealGrid &cgrid, const RealGrid &grid);
		friend CudaRealGrid &copyHostToDevice3D (CudaRealGrid &cgrid, const RealGrid &grid);
		friend CudaRealGrid &copyHostToDevice2D (CudaRealGrid &cgrid, const RealGrid &grid);
		friend CudaRealGrid &copyHostToDevice1D (CudaRealGrid &cgrid, const RealGrid &grid);

		friend RealGrid &copyDeviceToHost_as_complex (RealGrid &grid, const CudaRealGrid &cgrid);
		friend RealGrid &copyDeviceToHost3D (RealGrid &grid, const CudaRealGrid &cgrid);
		friend RealGrid &copyDeviceToHost2D (RealGrid &grid, const CudaRealGrid &cgrid);
		friend RealGrid &copyDeviceToHost1D (RealGrid &grid, const CudaRealGrid &cgrid);		

		inline RealGrid& operator =(const CudaRealGrid& right_side)
		{
			return copyDeviceToHost_as_complex(*this, right_side);
		}
		
		RealGrid& set_component(const RealGrid& c, uint32_t mu);
		RealGrid& get_component(RealGrid& c, uint32_t mu);
		RealGrid& multiply_as_real_with_real(const RealGrid& grid);
		RealGrid& multiply_as_complex_with_complex(const RealGrid& grid);
		
		inline Coordinate<int32_t> make_coord(int x, int y, int z) const {return Coordinate<int32_t>(x,y,z,width(),height(),depth());}
		inline Vector<int32_t> make_vector(int x, int y, int z) const {return Vector<int32_t>(x,y,z,width(),height(),depth());}
		
		double at(uint32_t mu, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid)const;
		double & at(uint32_t mu, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid);
		double at(uint32_t mu, const Coordinate<int32_t> &c)const;
		double & at(uint32_t mu, const Coordinate<int32_t> &c);
		complex<double> fat(uint32_t mu, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid)const;
		complex<double> & fat(uint32_t mu, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid);
		complex<double> fat(uint32_t mu, const Coordinate<int32_t> &c)const;
		complex<double> & fat(uint32_t mu, const Coordinate<int32_t> &c);
		
		double operator() (uint32_t mu, int32_t xGrid, int32_t yGrid, int32_t zGrid)const;
		double & operator() (uint32_t mu, int32_t xGrid, int32_t yGrid, int32_t zGrid);
		double operator() (uint32_t mu, const Coordinate<int32_t> &c)const;
		double & operator() (uint32_t mu, const Coordinate<int32_t> &c);
		
		static void fft(RealGrid &in, RealGrid &out, bool forward = true);
		static void fft_unnormalized(RealGrid &in, RealGrid &out, bool forward = true);
		static void set_fft_threads(int num_threads);
		static void set_fft_planning_rigorosity(unsigned flags);
		
		double *get_address() const {return store;}
		
		void resize(uint32_t int_dim, uint32_t w, uint32_t h, uint32_t d);
		inline uint32_t width() const {return dim[0];}
		inline uint32_t height() const {return dim[1];}
		inline uint32_t depth() const {return dim[2];}
		inline uint32_t int_dim() const {return internal;}
		inline uint32_t fft_width() const {return fft_dim[0];}
		inline uint32_t fft_height() const {return fft_dim[1];}
		inline uint32_t fft_depth() const {return fft_dim[2];}
		inline uint32_t cfft_width() const {return cfft_dim[0];}
		inline uint32_t cfft_height() const {return cfft_dim[1];}
		inline uint32_t cfft_depth() const {return cfft_dim[2];}
		
		friend ostream& operator<< (ostream& out, const RealGrid &grid);
		friend istream& operator>> (istream& in, RealGrid &grid);
	private:		
		static fftw_plan plan_r2c_out_of_place;
		static fftw_plan plan_c2r_out_of_place;
		static fftw_plan plan_r2c_in_place;
		static fftw_plan plan_c2r_in_place;
		static int plan_dim[4];
		static bool threads_initialized;
		static unsigned planning_rigor;
		
		uint32_t dim[3];   // dimensions of grid points in real space
		uint32_t fft_dim[3]; // dimesnions of grid points in Fourier space
		uint32_t cfft_dim[3];
		uint32_t internal; // internal dimension of field
		
		double *store;   // Pointer to dynamic array that stores gsl_complex.
		void allocate(uint32_t internal, uint32_t fft_xGrid, uint32_t fft_yGrid, uint32_t fft_zGrid);
		
		void deallocate();
		
		static void check_plan(uint32_t internal, uint32_t x, uint32_t y, uint32_t z);
};

typedef RealGrid* RealGridPtr;

// Inlined function included here for optimization:
inline double RealGrid::at(uint32_t mu, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid)const
{
	return store[mu*fft_dim[0]*fft_dim[1]*fft_dim[2] + xGrid*fft_dim[1]*fft_dim[2] + yGrid*fft_dim[2] + zGrid];
}

inline double& RealGrid::at(uint32_t mu, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid)
{
	return store[mu*fft_dim[0]*fft_dim[1]*fft_dim[2] + xGrid*fft_dim[1]*fft_dim[2] + yGrid*fft_dim[2] + zGrid];
}

inline complex<double> RealGrid::fat(uint32_t mu, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid)const
{
	return (reinterpret_cast<complex<double> *>(store))[mu*cfft_dim[0]*cfft_dim[1]*cfft_dim[2] + xGrid*cfft_dim[1]*cfft_dim[2] + yGrid*cfft_dim[2] + zGrid];
}

inline complex<double>& RealGrid::fat(uint32_t mu, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid)
{
        return (reinterpret_cast<complex<double> *>(store))[mu*cfft_dim[0]*cfft_dim[1]*cfft_dim[2] + xGrid*cfft_dim[1]*cfft_dim[2] + yGrid*cfft_dim[2] + zGrid];
}


inline double RealGrid::operator() (uint32_t mu, int32_t xGrid, int32_t yGrid, int32_t zGrid)const
{
	return at(mu, xGrid, yGrid, zGrid);
}

inline double& RealGrid::operator() (uint32_t mu, int32_t xGrid, int32_t yGrid, int32_t zGrid)
{
	return at(mu, xGrid, yGrid, zGrid);
}

inline double RealGrid::at(uint32_t mu, const Coordinate<int32_t> &c)const
{
	return at(mu, c.x(), c.y(), c.z());
}

inline double& RealGrid::at(uint32_t mu, const Coordinate<int32_t> &c)
{
	return at(mu, c.x(), c.y(), c.z());
}

inline complex<double> RealGrid::fat(uint32_t mu, const Coordinate<int32_t> &c)const
{
        return fat(mu, c.x(), c.y(), c.z());
}

inline complex<double>& RealGrid::fat(uint32_t mu, const Coordinate<int32_t> &c)
{
        return fat(mu, c.x(), c.y(), c.z());
}


inline double RealGrid::operator() (uint32_t mu, const Coordinate<int32_t> &c)const
{
	return at(mu, c.x(), c.y(), c.z());
}

inline double& RealGrid::operator() (uint32_t mu, const Coordinate<int32_t> &c)
{
	return at(mu, c.x(), c.y(), c.z());
}

#endif //GRID_H__.
