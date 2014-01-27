// File Name:     ComplexGrid.cpp
// Author:        Matthias Kronenwett, Jan Schole
// Email Address: m.kronenwett@thphys.uni-heidelberg.de
//
// Description:   This is the implementation file for the class ComplexGrid.
//                The interface for the class ComplexGrid is in the header
//                file complexgrid.h.
// Last Change:   February 06, 2011

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <complexgrid.h>
#include <bh3binaryfile.h>

using namespace std;


int ComplexGrid::plan_dim[4] = {0,0,0,0};
fftw_plan ComplexGrid::plan_forward_in_place;
fftw_plan ComplexGrid::plan_backward_in_place;
fftw_plan ComplexGrid::plan_forward_out_of_place;
fftw_plan ComplexGrid::plan_backward_out_of_place;
bool ComplexGrid::threads_initialized = false;
unsigned ComplexGrid::planning_rigor = FFTW_MEASURE;

// Allocation and initializing, Constructor
ComplexGrid::ComplexGrid(uint32_t int_dim, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid)
{
	dim[0]=xGrid;
	dim[1]=yGrid;
	dim[2]=zGrid;
	internal = int_dim;

	store = allocate(int_dim, xGrid, yGrid, zGrid);
}

//Default constructor
ComplexGrid::ComplexGrid()
{
	dim[0]=0;
	dim[1]=0;
	dim[2]=0;
	internal = 0;
	
	store = NULL;
}


// Copy constructor
ComplexGrid::ComplexGrid(const ComplexGrid& c)
{
	dim[0]=c.dim[0];
	dim[1]=c.dim[1];
	dim[2]=c.dim[2];
	internal = c.internal;
	
	store = allocate(internal, dim[0], dim[1], dim[2]);

	// Copy the data:
	memcpy(store, c.store, internal*dim[0]*dim[1]*dim[2]*sizeof(complex<double>));
}


// Constructing '=' operation
// Uses iostream:
ComplexGrid& ComplexGrid::operator =(const ComplexGrid& c)
{
	if (this == &c)
	{
		return *this;
	}
	else if ((dim[0] != c.dim[0]) || (dim[1] != c.dim[1]) || (dim[2] != c.dim[2]))
	{
		// Delete the old store:
		deallocate();

		// Adapt the private data:
		dim[0] = c.dim[0];
		dim[1] = c.dim[1];
		dim[2] = c.dim[2];
		internal = c.internal;

		// Allocate the new store:
		allocate(internal, dim[0], dim[1], dim[2]);
	}
	
	// Assign data:    
	memcpy(store, c.store, internal*dim[0]*dim[1]*dim[2]*sizeof(complex<double>));
	
	return *this;
}


// Allocation function
complex<double>* ComplexGrid::allocate(uint32_t int_dim, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid)

{
	store = (complex<double> *) fftw_malloc(sizeof(complex<double>) * int_dim*xGrid*yGrid*zGrid);

	if (!store)
	{
		cout << "Allocation failure 1 in ComplexGrid::allocate'"
				<< endl;
		exit(1);
	}

	memset(store, 0, int_dim*xGrid*yGrid*zGrid*sizeof(complex<double>));

	return store;
}


void ComplexGrid::deallocate()
{
	if(store != NULL)
		fftw_free(store);
	
	dim[0] = 0;
	dim[1] = 0;
	dim[2] = 0;
	internal = 0;
}



ComplexGrid::~ComplexGrid()
{
	deallocate();
}

void ComplexGrid::resize(uint32_t int_dim, uint32_t w, uint32_t h, uint32_t d)
{
	deallocate();
	dim[0] = w;
	dim[1] = h;
	dim[2] = d;
	internal = int_dim;
	store = allocate(int_dim, w, h, d);
}

void ComplexGrid::set_fft_threads(int num_threads)
{
	if(!threads_initialized)
	{
		fftw_init_threads();
		threads_initialized = true;
	}
	fftw_plan_with_nthreads(num_threads);
}

void ComplexGrid::check_plan(uint32_t int_dim, uint32_t x, uint32_t y, uint32_t z)
{
	#pragma omp critical
	{
		if((plan_dim[0] != int_dim) || (plan_dim[1] != x) || (plan_dim[2] != y) || (plan_dim[3] != z))
		{
			ComplexGrid in(int_dim,x,y,z);
			ComplexGrid out(int_dim,x,y,z);
			int rank = z > 1 ? 3 : ((z==1 && y > 1) ? 2 : ((z==1 && y==1 && x>1) ? 1 : 0));
			int n[rank];

			if (rank == 1)
			{
				n[0] = x;
			}
			else if (rank == 2)
			{
				n[0] = x;
				n[1] = y;
			}
			else if (rank == 3)
			{
				n[0] = x;
				n[1] = y;
				n[2] = z;
			}

			plan_forward_out_of_place = fftw_plan_many_dft(rank, n, int_dim, reinterpret_cast<fftw_complex*>(in.store), NULL, 1, x*y*z, reinterpret_cast<fftw_complex*>(out.store), NULL, 1, x*y*z, FFTW_FORWARD, planning_rigor);
			plan_backward_out_of_place = fftw_plan_many_dft(rank, n, int_dim, reinterpret_cast<fftw_complex*>(in.store), NULL, 1, x*y*z, reinterpret_cast<fftw_complex*>(out.store), NULL, 1, x*y*z, FFTW_BACKWARD, planning_rigor);
			plan_forward_in_place = fftw_plan_many_dft(rank, n, int_dim, reinterpret_cast<fftw_complex*>(in.store), NULL, 1, x*y*z, reinterpret_cast<fftw_complex*>(in.store), NULL, 1, x*y*z, FFTW_FORWARD, planning_rigor);
			plan_backward_in_place = fftw_plan_many_dft(rank, n, int_dim, reinterpret_cast<fftw_complex*>(in.store), NULL, 1, x*y*z, reinterpret_cast<fftw_complex*>(in.store), NULL, 1, x*y*z, FFTW_BACKWARD, planning_rigor);

			
			plan_dim[0] = int_dim;
			plan_dim[1] = x;
			plan_dim[2] = y;
			plan_dim[3] = z;
		}
	}
}

void ComplexGrid::set_fft_planning_rigorosity(unsigned flags)
{
	planning_rigor = flags;
	if((plan_dim[0] != 0) || (plan_dim[1] != 0) || (plan_dim[2] != 0) || (plan_dim[3] != 0))
	{
		int x = plan_dim[0];
		plan_dim[0]++;
		check_plan(x, plan_dim[1], plan_dim[2], plan_dim[3]);
	}
}

void ComplexGrid::fft_unnormalized(ComplexGrid &in, ComplexGrid &out, bool forward)
{
	check_plan(in.int_dim(), in.width(),in.height(), in.depth());
	
	if(&in == &out)
	{
		if(forward)
			fftw_execute_dft(plan_forward_in_place, reinterpret_cast<fftw_complex*>(in.store), reinterpret_cast<fftw_complex*>(in.store));
		else
			fftw_execute_dft(plan_backward_in_place, reinterpret_cast<fftw_complex*>(in.store), reinterpret_cast<fftw_complex*>(in.store));
	}
	else
	{
		if((in.width() != out.width()) || (in.height() != out.height()) || (in.depth() != out.depth()))
			out.resize(in.int_dim(), in.width(), in.height(), in.depth());

		if(forward)
			fftw_execute_dft(plan_forward_out_of_place, reinterpret_cast<fftw_complex*>(in.store), reinterpret_cast<fftw_complex*>(out.store));
		else
			fftw_execute_dft(plan_backward_out_of_place, reinterpret_cast<fftw_complex*>(in.store), reinterpret_cast<fftw_complex*>(out.store));
	}
}

void ComplexGrid::fft(ComplexGrid &in, ComplexGrid &out, bool forward)
{
	ComplexGrid::fft_unnormalized(in, out, forward);
	double fft_factor = 1.0 / sqrt(in.width()*in.height()*in.depth());
	
	for(int i = 0; i < out.int_dim()*out.width()*out.height()*out.depth(); i++)
	{
		out.store[i] *= fft_factor;
	}
}

void write(ostream& out, const ComplexGrid &grid)
{
	write(out, grid.int_dim());  
	write(out, grid.width());
	write(out, grid.height());
	write(out, grid.depth());

	out.write((const char *)grid.store, grid.int_dim()*grid.width()*grid.height()*grid.depth()*sizeof(complex<double>));
}

void read(istream& in, ComplexGrid &grid)
{
	uint32_t newdim[4];
	read(in, newdim[0]);
	read(in, newdim[1]);
	read(in, newdim[2]);
	read(in, newdim[3]);
	
	if((newdim[0] != grid.internal) || (newdim[1] != grid.dim[1]) || (newdim[2] != grid.dim[2]) || (newdim[3] != grid.dim[3]))
	{
		grid.resize(newdim[0], newdim[1], newdim[2], newdim[3]);
	}
	in.read((char *)grid.store, grid.int_dim()*grid.width()*grid.height()*grid.depth()*sizeof(complex<double>));
}

ComplexGrid& ComplexGrid::set_component(const ComplexGrid& c, uint32_t mu)
{
	if (c.internal != 1)
	{
		cout << "==> Warning in 'RealGrid::set_component':" << endl
		<< "    r.h.s field must have internal dimension of 1!" << endl;
		return *this;
	}
	
	else if (mu >= internal)
	{
		cout << "==> Warning in 'RealGrid::operator=':" << endl
		<< "    l.h.s does not have a field component "<< mu << endl;
		return *this;
	}
	
	else if ((dim[0] != c.dim[0]) || (dim[1] != c.dim[1]) || (dim[2] != c.dim[2]))
	{
		cout << "==> Warning in 'RealGrid::set_component':" << endl
		<< "    The fields are of different dimensions," << endl
		<< "    the field on the l.h.s. will be resized." << endl;
		
		// Delete the old store:
		deallocate();
		
		// Adapt the private data:
		dim[0] = c.dim[0];
		dim[1] = c.dim[1];
		dim[2] = c.dim[2];	

		
		// Allocate the new store:
		allocate(internal, dim[0], dim[1], dim[2]);
	}
	else if (this == &c)
	{
		return *this;
	}
	
	// Assign data:
	memcpy(store + mu*dim[0]*dim[1]*dim[2], c.store, dim[0]*dim[1]*dim[2]*sizeof(complex<double>));
	
	return *this;
}

ComplexGrid& ComplexGrid::get_component(ComplexGrid& c, uint32_t mu)
{
	
	if (mu >= internal)
	{
		cout << "==> Warning in 'Grid::get_component':" << endl
		<< "    field does not have a component "<< mu << endl;
		return *this;
	}
	else if ((dim[0] != c.dim[0]) || (dim[1] != c.dim[1]) || (dim[2] != c.dim[2]))
	{

	
		// resize target field to fit spatial dimensions
		c.resize(1,dim[0],dim[1],dim[2]);

	}
	
	// Assign data:
	memcpy(c.store, store + mu*dim[0]*dim[1]*dim[2], dim[0]*dim[1]*dim[2]*sizeof(complex<double>));
	
	return *this;
}
