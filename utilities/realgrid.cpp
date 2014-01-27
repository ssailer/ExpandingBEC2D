
#include <iostream>
#include <cstring>
#include <cstdlib>
#include "bh3binaryfile.h"
#include "realgrid.h"
#include <omp.h>

using namespace std;

int RealGrid::plan_dim[4] = {0,0,0,0};
fftw_plan RealGrid::plan_r2c_in_place;
fftw_plan RealGrid::plan_c2r_in_place;
fftw_plan RealGrid::plan_r2c_out_of_place;
fftw_plan RealGrid::plan_c2r_out_of_place;
bool RealGrid::threads_initialized = false;
unsigned RealGrid::planning_rigor = FFTW_MEASURE;

// Allocation and initializing, Constructor
RealGrid::RealGrid(uint32_t int_dim, uint32_t xGrid, uint32_t yGrid, uint32_t zGrid)
{
	dim[0] = xGrid;
	dim[1] = yGrid;
	dim[2] = zGrid;
	internal = int_dim;
	if (zGrid > 1)
	{
		fft_dim[0] = dim[0];
		fft_dim[1] = dim[1];
		fft_dim[2] = 2*(int(dim[2]/2) +1);
		cfft_dim[0] = dim[0];
		cfft_dim[1] = dim[1];
		cfft_dim[2] = fft_dim[2]/2;
	}
	else if (zGrid == 1 && yGrid > 1)
	{
		fft_dim[0] = dim[0];
		fft_dim[1] = 2*(int(dim[1]/2) +1);	
		fft_dim[2] = dim[2];
		cfft_dim[0] = dim[0];
		cfft_dim[1] = fft_dim[1]/2;
		cfft_dim[2] = dim[2];
	}
	else if (zGrid == 1 && yGrid == 1 && xGrid > 1)
	{
		fft_dim[0] = 2*(int(dim[0]/2) +1);
		fft_dim[1] = dim[1];	
		fft_dim[2] = dim[2];
		cfft_dim[0] = fft_dim[0]/2;
		cfft_dim[1] = dim[1];
		cfft_dim[2] = dim[2];
	}
	else
	{
		fft_dim[0] = dim[0];
		fft_dim[1] = dim[1];	
		fft_dim[2] = dim[2];
		cfft_dim[0] = dim[0];
		cfft_dim[1] = dim[1];
		cfft_dim[2] = dim[2];
	}
	
	allocate(int_dim, fft_dim[0], fft_dim[1], fft_dim[2]);
	

		
}

//Default constructor
RealGrid::RealGrid()
{
	dim[0] = 0;
	dim[1] = 0;
	dim[2] = 0;
	internal = 0;
	fft_dim[0] = 0;
	fft_dim[1] = 0;
	fft_dim[2] = 0;
	cfft_dim[0] = 0;
	cfft_dim[1] = 0;
	cfft_dim[2] = 0;
	
	store = NULL;
}


// Copy constructor
RealGrid::RealGrid(const RealGrid& c)
{
	dim[0] = c.dim[0];
	dim[1] = c.dim[1];
	dim[2] = c.dim[2];
	internal = c.internal;
	fft_dim[0] = c.fft_dim[0];
	fft_dim[1] = c.fft_dim[1];
	fft_dim[2] = c.fft_dim[2];
	cfft_dim[0] = c.cfft_dim[0];
	cfft_dim[1] = c.cfft_dim[1];
	cfft_dim[2] = c.cfft_dim[2];
	
	allocate(internal, fft_dim[0], fft_dim[1], fft_dim[2]);
	
	// Copy the data:
	memcpy(store, c.store, internal*fft_dim[0]*fft_dim[1]*fft_dim[2]*sizeof(double));
}


// Constructing '=' operation
// Uses iostream:
RealGrid& RealGrid::operator =(const RealGrid& c)
{
	if ((dim[0] != c.dim[0]) || (dim[1] != c.dim[1]) || (dim[2] != c.dim[2]) || (internal != c.internal))
	{
		// cout << "==> Warning in 'RealGrid::operator=':" << endl
		// << "    The objects are of different sizes," << endl
		// << "    the object on the l.h.s. will be resized." << endl;
		
		// Delete the old store:
		deallocate();
		
		// Adapt the private data:
		dim[0] = c.dim[0];
		dim[1] = c.dim[1];
		dim[2] = c.dim[2];
		internal = c.internal;
		fft_dim[0] = c.fft_dim[0];
		fft_dim[1] = c.fft_dim[1];
		fft_dim[2] = c.fft_dim[2];
		cfft_dim[0] = c.cfft_dim[0];
		cfft_dim[1] = c.cfft_dim[1];
		cfft_dim[2] = c.cfft_dim[2];		
		
		// Allocate the new store:
		allocate(internal, fft_dim[0], fft_dim[1], fft_dim[2]);
	}
	else if (this == &c)
	{
		return *this;
	}
	
	// Assign data:    
	memcpy(store, c.store, internal*fft_dim[0]*fft_dim[1]*fft_dim[2]*sizeof(double));
	
	return *this;
}

RealGrid& RealGrid::set_component(const RealGrid& c, uint32_t mu)
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
		// cout << "==> Warning in 'RealGrid::set_component':" << endl
		// << "    The fields are of different dimensions," << endl
		// << "    the field on the l.h.s. will be resized." << endl;
		
		// Delete the old store:
		deallocate();
		
		// Adapt the private data:
		dim[0] = c.dim[0];
		dim[1] = c.dim[1];
		dim[2] = c.dim[2];
		fft_dim[0] = c.fft_dim[0];
		fft_dim[1] = c.fft_dim[1];
		fft_dim[2] = c.fft_dim[2];		
		cfft_dim[0] = c.cfft_dim[0];
		cfft_dim[1] = c.cfft_dim[1];
		cfft_dim[2] = c.cfft_dim[2];
		
		// Allocate the new store:
		allocate(internal, fft_dim[0], fft_dim[1], fft_dim[2]);
	}
	else if (this == &c)
	{
		return *this;
	}
	
	// Assign data:
	memcpy(store + mu*fft_dim[0]*fft_dim[1]*fft_dim[2], c.store, fft_dim[0]*fft_dim[1]*fft_dim[2]*sizeof(double));
	
	return *this;
}

RealGrid& RealGrid::get_component(RealGrid& c, uint32_t mu)
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
	memcpy(c.store, store + mu*fft_dim[0]*fft_dim[1]*fft_dim[2], fft_dim[0]*fft_dim[1]*fft_dim[2]*sizeof(double));
	
	return *this;
}


// Allocation function
void RealGrid::allocate(uint32_t internal, uint32_t fft_xGrid, uint32_t fft_yGrid, uint32_t fft_zGrid)
{
	
	store = new double[internal*fft_xGrid*fft_yGrid*fft_zGrid];  //allocate memory for Grid
	
	if (!store)
	{
		cout << "Allocation failure 1 in Grid::allocate'"
		<< endl;
		exit(1);
	}
	
	memset(store, 0, internal*fft_xGrid*fft_yGrid*fft_zGrid*sizeof(double)); //initialise Grid memory with zeros
}


void RealGrid::deallocate()
{
	if(store)
		delete store;
}



RealGrid::~RealGrid()
{
	deallocate();
}

void RealGrid::resize(uint32_t int_dim, uint32_t w, uint32_t h, uint32_t d)
{
	if(store)
		deallocate();
	
	dim[0] = w;
	dim[1] = h;
	dim[2] = d;
	internal = int_dim;
	if (d > 1)
	{
		fft_dim[0] = dim[0];
		fft_dim[1] = dim[1];
		fft_dim[2] = 2*(int(dim[2]/2) +1);
		cfft_dim[0] = dim[0];
		cfft_dim[1] = dim[1];
		cfft_dim[2] = fft_dim[2]/2;		
	}
	else if (d == 1 && h > 1)
	{
		fft_dim[0] = dim[0];
		fft_dim[1] = 2*(int(dim[1]/2) +1);	
		fft_dim[2] = dim[2];
		cfft_dim[0] = dim[0];
		cfft_dim[1] = fft_dim[1]/2;
		cfft_dim[2] = dim[2];		
	}
	else if (d == 1 && h == 1 && w > 1)
	{
		fft_dim[0] = 2*(int(dim[0]/2) +1);
		fft_dim[1] = dim[1];	
		fft_dim[2] = dim[2];
		cfft_dim[0] = fft_dim[0]/2;
		cfft_dim[1] = dim[1];
		cfft_dim[2] = dim[2];		
	}
	else
	{
		fft_dim[0] = dim[0];
		fft_dim[1] = dim[1];	
		fft_dim[2] = dim[2];
		cfft_dim[0] = dim[0];
		cfft_dim[1] = dim[1];
		cfft_dim[2] = dim[2];		
	}
	
	allocate(int_dim, fft_dim[0], fft_dim[1], fft_dim[2]);
}

void RealGrid::set_fft_threads(int num_threads)
{
	if(!threads_initialized)
	{
		fftw_init_threads();
		threads_initialized = true;
	}
	fftw_plan_with_nthreads(num_threads);
}

void RealGrid::check_plan(uint32_t internal, uint32_t x, uint32_t y, uint32_t z)
{
	#pragma omp critical
	{
		
		if((plan_dim[0] != internal) || (plan_dim[1] != x) || (plan_dim[2] != y) || (plan_dim[3] != z))
		{
			RealGrid in(internal, x,y,z);
			RealGrid out(internal, x,y,z);
			uint32_t rank;
			uint32_t real_dist = 0;
			
			if (z > 1)
			{
				rank = 3;
				real_dist = x*y*2*(int(z/2) +1);
			}
			else if (z == 1 && y > 1)
			{
				rank = 2;
				real_dist = x*2*(int(y/2) +1);
			}
			else if (z == 1 && y == 1 && x > 1)
			{
				rank = 1;
				real_dist = 2*(int(x/2) +1);
			}
			else
				rank = 0;
			
			int n[rank];
			
			if (rank==1)
				n[0] = x;
			else if (rank ==2)
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

			plan_r2c_out_of_place = fftw_plan_many_dft_r2c(rank, n, internal, in.store, NULL, 1, x*y*z, reinterpret_cast<fftw_complex*>(out.store), NULL,1, real_dist/2, planning_rigor);
			plan_c2r_out_of_place = fftw_plan_many_dft_c2r(rank, n, internal, reinterpret_cast<fftw_complex*>(in.store), NULL,1,real_dist/2, out.store, NULL, 1, x*y*z, planning_rigor);
			plan_r2c_in_place = fftw_plan_many_dft_r2c(rank, n, internal, in.store, NULL ,1,real_dist, reinterpret_cast<fftw_complex*>(in.store), NULL, 1, real_dist/2, planning_rigor);
			plan_c2r_in_place = fftw_plan_many_dft_c2r(rank, n, internal, reinterpret_cast<fftw_complex*>(in.store), NULL,1,real_dist/2, in.store, NULL, 1, real_dist, planning_rigor);
	

			plan_dim[0] = internal;
			plan_dim[1] = x;
			plan_dim[2] = y;
			plan_dim[3] = z;
			
		}
	}
}

void RealGrid::set_fft_planning_rigorosity(unsigned flags)
{
	planning_rigor = flags;
	if((plan_dim[0] != 0) || (plan_dim[1] != 0) || (plan_dim[2] != 0) || (plan_dim[3] != 0))
	{
		int x = plan_dim[0];
		plan_dim[0]++;
		check_plan(x, plan_dim[1], plan_dim[2], plan_dim[3]);
	}
}

void RealGrid::fft_unnormalized(RealGrid &in, RealGrid &out, bool forward)
{

	check_plan(in.int_dim(), in.width(),in.height(), in.depth());

	if(&in == &out)
	{
		if(forward)
			fftw_execute_dft_r2c(plan_r2c_in_place, in.store, reinterpret_cast<fftw_complex*>(in.store));
		else
			fftw_execute_dft_c2r(plan_c2r_in_place, reinterpret_cast<fftw_complex*>(in.store), in.store);
	}
	else
	{
		cout << "out of place r2c or c2r FFTs are not supported at this stage" << endl;
		/*if((in.int_dim() != out.int_dim()) || (in.width() != out.width()) || (in.height() != out.height()) || (in.depth() != out.depth()))
			out.resize(in.int_dim(), in.width(), in.height(), in.depth());

		if(forward)
			fftw_execute_dft_r2c(plan_r2c_out_of_place, in.store, reinterpret_cast<fftw_complex*>(out.store));
		else
			fftw_execute_dft_c2r(plan_c2r_out_of_place, reinterpret_cast<fftw_complex*>(in.store), out.store);*/
	}
}

void RealGrid::fft(RealGrid &in, RealGrid &out, bool forward)
{
	RealGrid::fft_unnormalized(in, out, forward);
	double fft_factor = 1.0 / sqrt(in.width()*in.height()*in.depth());
	
	for(int i = 0; i < out.int_dim()*out.fft_width()*out.fft_height()*out.fft_depth(); i++)
	{
		out.store[i] *= fft_factor;
	}
}

RealGrid& RealGrid::multiply_as_real_with_real(const RealGrid& grid)
{
	if ((grid.int_dim()!= internal) || (grid.width()!= dim[0]) || (grid.height()!= dim[1]) || (grid.depth()!= dim[2]))
	{
		cout << "Dimensions do not fit. No multiplication is carried out" << endl;
		return *this;
	}
	
	
	#pragma omp for schedule(guided, 1)
	for (int mu = 0; mu < internal; mu++)
	{
		for (int x = 0; x < dim[0]; x++)
		{
			for (int y = 0; y < dim[1]; y++)
			{
				for (int z = 0; z < dim[2]; z++)
				{
					this->at(mu,x,y,z) *= grid(mu,x,y,z);
				}
			}
		}	  
	}
	
	return *this;
	  
}

RealGrid& RealGrid::multiply_as_complex_with_complex(const RealGrid& grid)
{
	if ((grid.int_dim()!= internal) || (grid.width()!= dim[0]) || (grid.height()!= dim[1]) || (grid.depth()!= dim[2]))
	{
		cout << "Dimensions do not fit. No multiplication is carried out" << endl;
		return *this;
	}
	
	int last_dim =  dim[2]!=fft_dim[2] ? 2 : (dim[1]!=fft_dim[1] ? 1 : (dim[0]!=fft_dim[0] ? 0 : -1));
	
	if (last_dim == 2)
	{
		#pragma omp for schedule(guided, 1)
		for (int mu = 0; mu < internal; mu++)
		{
			for (int x = 0; x < dim[0]; x++)
			{
				for (int y = 0; y < dim[1]; y++)
				{
					for (int z = 0; z < dim[2]; z+=2)
					{
						complex<double> value1 = complex<double>(grid(mu,x,y,z), grid(mu,x,y,z+1));
						complex<double> value2 = complex<double>(this->at(mu,x,y,z), this->at(mu,x,y,z+1));
						this->at(mu,x,y,z) = real(value1*value2);
						this->at(mu,x,y,z+1) = imag(value1*value2);
					}
				}
			}	  
		}
	}
	
	else if (last_dim == 1)
	{
		#pragma omp for schedule(guided, 1)
		for (int mu = 0; mu < internal; mu++)
		{
			for (int x = 0; x < dim[0]; x++)
			{
				for (int y = 0; y < dim[1]; y +=2)
				{
					complex<double> value1 = complex<double>(grid(mu,x,y,1), grid(mu,x,y+1,1));
					complex<double> value2 = complex<double>(this->at(mu,x,y,1), this->at(mu,x,y+1,1));
					this->at(mu,x,y,1) = real(value1*value2);
					this->at(mu,x,y+1,1) = imag(value1*value2);
				}
			}	  
		}
	}
	
	else if (last_dim == 0)
	{
		#pragma omp for schedule(guided, 1)
		for (int mu = 0; mu < internal; mu++)
		{
			for (int x = 0; x < dim[0]; x+=2)
			{
				complex<double> value1 = complex<double>(grid(mu,x,1,1), grid(mu,x+1,1,1));
				complex<double> value2 = complex<double>(this->at(mu,x,1,1), this->at(mu,x+1,1,1));
				this->at(mu,x,1,1) = real(value1*value2);
				this->at(mu,x+1,1,1) = imag(value1*value2);
			}	  
		}
	}
	
	return *this;
	  
}


ostream& operator<< (ostream& out, const RealGrid &grid)
{
	write(out, grid.fft_width());
	write(out, grid.fft_height());
	write(out, grid.fft_depth());
	write(out, grid.int_dim());
	out.write((const char *)grid.store, grid.int_dim()*grid.fft_width()*grid.fft_height()*grid.fft_depth()*sizeof(double));
	return out;
}

istream& operator>> (istream& in, RealGrid &grid)
{
	uint32_t newdim[3];
	uint32_t newintdim;
	read(in, newdim[0]);
	read(in, newdim[1]);
	read(in, newdim[2]);
	read(in, newintdim);
	if((newdim[0] != grid.dim[0]) || (newdim[1] != grid.dim[1]) || (newdim[2] != grid.dim[2]) || (newintdim != grid.internal))
	{
		grid.resize(newdim[0], newdim[1], newdim[2], newintdim);
	}
	in.read((char *)grid.store, grid.int_dim()*grid.fft_width()*grid.fft_height()*grid.fft_depth()*sizeof(double));
	return in;
}

