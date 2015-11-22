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

#ifndef COMPLEX_GRID_H__
#define COMPLEX_GRID_H__
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <iostream>
#include <stdint.h>
using namespace std;

// Project includes
#include <coordinate.h>

class CudaComplexGrid;

class ComplexGrid
{
	public:
		ComplexGrid(uint32_t int_dim, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid);
		ComplexGrid();
		ComplexGrid(const ComplexGrid& c);
  
		~ComplexGrid();
  
		ComplexGrid& operator =(const ComplexGrid& right_side);
		
		ComplexGrid& set_component(const ComplexGrid& c, uint32_t mu);
		ComplexGrid& get_component(ComplexGrid& c, uint32_t mu);

		friend CudaComplexGrid &copyHostToDevice (CudaComplexGrid &cgrid, const ComplexGrid &grid);
		friend ComplexGrid &copyDeviceToHost (ComplexGrid &grid, const CudaComplexGrid &cgrid);

		inline ComplexGrid& operator =(const CudaComplexGrid& right_side)
		{
			return copyDeviceToHost(*this, right_side);
		}
		
		inline Coordinate<int32_t> make_coord(int x, int y, int z) const {return Coordinate<int32_t>(x,y,z,width(),height(),depth());}
		inline Vector<int32_t> make_vector(int x, int y, int z) const {return Vector<int32_t>(x,y,z,width(),height(),depth());}
		
		const complex<double> & at(uint32_t mu, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid)const;
		complex<double> & at(uint32_t mu, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid);
		const complex<double> & at(uint32_t mu, const Coordinate<int32_t> &c)const;
		complex<double> & at(uint32_t mu, const Coordinate<int32_t> &c);
		
		const complex<double> & operator() (uint32_t mu, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid)const;
		complex<double> & operator() (uint32_t mu, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid);
		const complex<double> & operator() (uint32_t mu, const Coordinate<int32_t> &c)const;
		complex<double> & operator() (uint32_t mu, const Coordinate<int32_t> &c);
		
		static void fft(ComplexGrid &in, ComplexGrid &out, bool forward = true);
		static void fft_unnormalized(ComplexGrid &in, ComplexGrid &out, bool forward = true);
		static void set_fft_threads(int num_threads);
		static void set_fft_planning_rigorosity(unsigned flags);
		
		complex<double> *get_address() const {return store;}
		
		void resize(uint32_t int_dim, uint32_t w, uint32_t h, uint32_t d);
		inline uint32_t width() const {return dim[0];}
		inline uint32_t height() const {return dim[1];}
		inline uint32_t depth() const {return dim[2];}
		inline uint32_t int_dim() const {return internal;}
		
		friend void write(ostream &o, const ComplexGrid &grid);
		friend void read(istream &i, ComplexGrid &grid);
	private:
		static fftw_plan plan_forward_out_of_place;
		static fftw_plan plan_backward_out_of_place;
		static fftw_plan plan_forward_in_place;
		static fftw_plan plan_backward_in_place;
		static int plan_dim[4];
		static bool threads_initialized;
		static unsigned planning_rigor;

		uint32_t dim[3]; // dimensions of grid points.
		uint32_t internal;
  
		complex<double> *store;   // Pointer to dynamic array that stores gsl_complex.
		complex<double>* allocate(uint32_t int_dim, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid);

		void deallocate();

		static void check_plan(uint32_t int_dim, uint32_t x, uint32_t y, uint32_t z);
};

typedef ComplexGrid* ComplexGridPtr;


// for shortening
template <class T> inline T abs2(const complex<T> &c)
{
	return c.real()*c.real() + c.imag()*c.imag();
}

// Inlined function included here for optimization:
inline const complex<double>& ComplexGrid::at(uint32_t mu, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid)const
{
	return store[mu*dim[0]*dim[1]*dim[2] + xGrid*dim[1]*dim[2] + yGrid*dim[2] + zGrid];
}

inline complex<double>& ComplexGrid::at(uint32_t mu, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid)
{
	return store[mu*dim[0]*dim[1]*dim[2] + xGrid*dim[1]*dim[2] + yGrid*dim[2] + zGrid];
}

inline const complex<double>& ComplexGrid::operator() (uint32_t mu, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid)const
{
	return at(mu, xGrid, yGrid, zGrid);
}

inline complex<double>& ComplexGrid::operator() (uint32_t mu, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid)
{
	return at(mu, xGrid, yGrid, zGrid);
}

inline const complex<double>& ComplexGrid::at(uint32_t mu, const Coordinate<int32_t> &c)const
{
	return at(mu, c[0], c[1], c[2]);
}

inline complex<double>& ComplexGrid::at(uint32_t mu, const Coordinate<int32_t> &c)
{
	return at(mu, c[0], c[1], c[2]);
}

inline const complex<double>& ComplexGrid::operator() (uint32_t mu, const Coordinate<int32_t> &c)const
{
	return at(mu, c[0], c[1], c[2]);
}

inline complex<double>& ComplexGrid::operator() (uint32_t mu, const Coordinate<int32_t> &c)
{
	return at(mu, c[0], c[1], c[2]);
}

inline ostream& operator<< (ostream& out, const ComplexGrid &grid)
{
	write(out, grid);
	return out;
}

inline istream& operator>> (istream& in, ComplexGrid &grid)
{
	read(in, grid);
	return in;
}

#endif //COMPLEX_GRID_H__.

