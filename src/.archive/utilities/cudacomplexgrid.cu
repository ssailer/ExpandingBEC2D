#ifndef CUDACOMPLEXGRID_CU__
#define CUDACOMPLEXGRID_CU__

#include <iostream>
#include <omp.h>

#include "cudacomplexgrid.h"
#include "complexgrid.h"
#include "wrapped_cuda_functions.h"

CudaComplexGrid::CudaComplexGrid()
{
	dim[0] = dim[1] = dim[2] = internal = 0;
	store = NULL;
}

CudaComplexGrid::CudaComplexGrid(uint32_t int_dim, uint32_t w, uint32_t h, uint32_t d)
{
	initialize(int_dim, w, h, d);
}

CudaComplexGrid::CudaComplexGrid(const CudaComplexGrid &grid)
{
	initialize(grid.internal, grid.dim[0], grid.dim[1], grid.dim[2]);
	cudaMemcpy2D(store, pitch, grid.store, grid.pitch, dim[0]*dim[1]*dim[2]*sizeof(cuDoubleComplex), internal, cudaMemcpyDeviceToDevice);
}

CudaComplexGrid::~CudaComplexGrid()
{
	deallocate();
}

void CudaComplexGrid::createPlan(const size_t ipitch, const size_t opitch)
{
	cufftResult res;
			
	int rank = dim[2] > 1 ? 3 : ((dim[2]==1 && dim[1] > 1) ? 2 : ((dim[2]==1 && dim[1]==1 && dim[0]>1) ? 1 : 0));
	int n[rank];

	if (rank == 1)
	{
		n[0] = dim[0];
	}
	else if (rank == 2)
	{
		n[0] = dim[0];
		n[1] = dim[1];
	}
	else if (rank == 3)
	{
		n[0] = dim[0];
		n[1] = dim[1];
		n[2] = dim[2];
	}

	res = cufftPlanMany(&plan, rank, n, NULL, 1, (size_t) ipitch/sizeof(cuDoubleComplex), NULL, 1, (size_t) opitch/sizeof(cuDoubleComplex), CUFFT_Z2Z, internal);
	

	if(res != CUFFT_SUCCESS)
	{
		cout << "Creating CUFFT-Plan failed in thread " << omp_get_thread_num() << "!" << endl << "Reason:  ";
		switch(res) {
			case CUFFT_SETUP_FAILED:
				cout << "Setup failed!" << endl;
			case CUFFT_INVALID_SIZE:
				cout << "Invalid size!" << endl;
			case CUFFT_INVALID_TYPE:
				cout << "Invalid type!" << endl;
			case CUFFT_ALLOC_FAILED:
				cout << "Allocation failed!" << endl;
			default:
				cout << "Unknown error!" << endl;
				break;
		}
	}
}

void CudaComplexGrid::initialize(uint32_t int_dim, uint32_t w, uint32_t h, uint32_t d)
{
	dim[0] = w;
	dim[1] = h;
	dim[2] = d;
	internal = int_dim;

	if(cudaMallocPitch(&store, &pitch,  sizeof(cuDoubleComplex)*dim[0]*dim[1]*dim[2], internal) != cudaSuccess)
	{
		cout << "Error allocating Cuda-Memory in thread " << omp_get_thread_num() << "!" << endl;
	}
	if (pitch/(sizeof(cuDoubleComplex)) != floor(pitch/(sizeof(cuDoubleComplex))))
	{
		cout << "development error: pitch for cudacomplexgrid is not an integer multiple of cuDoubleComplex. CUFFT will most likely fail" << endl;
	}
	createPlan(pitch, pitch);
}

void CudaComplexGrid::deallocate()
{
	cufftDestroy(plan);
	cudaFree(store);
}

inline const char *cufft_error_string(int error_id)
{
	switch(error_id) {
		case CUFFT_SUCCESS:					// Should never occur
			return "Success!";
		case CUFFT_SETUP_FAILED:
			return "Setup failed!";
		case CUFFT_INVALID_PLAN:
			return "Invalid plan!";
		case CUFFT_INVALID_VALUE:
			return "Invalid value!";
		case CUFFT_EXEC_FAILED:
			return "Execution on GPU failed!";
		default:
			return "Unknown error!";
	}
}

bool CudaComplexGrid::fft(CudaComplexGrid &i, CudaComplexGrid &o, int direction)
{
	int res[4] = {CUFFT_SUCCESS, CUFFT_SUCCESS, CUFFT_SUCCESS, CUFFT_SUCCESS};
	int j = -1;

	if(i.pitch != o.pitch)
	{
		cout << "development error: input pitch does not match output pitch. out of place CUFFT will most likely fail" << endl;
	}

	do
	{
		j++;
		res[j] = cufftExecZ2Z(i.plan, i.store, o.store, direction);
	} while((res[j] != CUFFT_SUCCESS) && (j < 3));
	if(res[0] != CUFFT_SUCCESS)
	{
		string dir = (direction == CUFFT_INVERSE) ? "inverse" : ((direction == CUFFT_FORWARD) ? "forward" : "Invalid direction");
		cout << "Thread " << omp_get_thread_num() << ", Direction " << dir << ":" << endl;
		for(int k = 0; k < 4; k++)
		{
			if(res[k] != CUFFT_SUCCESS)
				cout << "Warning: Try " << k << " of fourier transform failed: " << cufft_error_string(res[k]) << endl;
		}
		if((res[1] != CUFFT_SUCCESS) && (res[2] != CUFFT_SUCCESS) && (res[3] != CUFFT_SUCCESS))
		{
			cout << "Error: All trys of fourier transforming failed in thread " << omp_get_thread_num() << "!" << endl;
			return false;
		}
	}
	return true;
}

void CudaComplexGrid::resize (uint32_t int_dim, uint32_t w, uint32_t h, uint32_t d)
{
	deallocate();
	initialize(int_dim, w, h, d);
}

CudaComplexGrid &CudaComplexGrid::operator= (const CudaComplexGrid &grid)
{
	if(&grid != this)
	{
		if((internal != grid.internal) || (dim[0] != grid.dim[0]) || (dim[1] != grid.dim[1]) || (dim[2] != grid.dim[2]))
			resize(grid.internal, grid.dim[0], grid.dim[1], grid.dim[2]);
		cudaMemcpy2D(store, pitch, grid.store, grid.get_pitch(), dim[0]*dim[1]*dim[2]*sizeof(cuDoubleComplex), internal, cudaMemcpyDeviceToDevice);
	}
	return *this;
}

CudaComplexGrid &copyHostToDevice (CudaComplexGrid &cgrid, const ComplexGrid &grid)
{
	if ((cgrid.int_dim() != grid.int_dim()) || (cgrid.width() != grid.width()) || (cgrid.height() != grid.height()) || (cgrid.depth() != grid.depth()))
	{
		cout << "Warning resizing CudaComplexGrid for copying!" << endl;
		cgrid.resize(grid.int_dim(), grid.width(), grid.height(), grid.depth());
	}
	memcpypitch_host_to_device(cgrid.store, cgrid.get_pitch(), grid.store, sizeof(cuDoubleComplex)*grid.width()*grid.height()*grid.depth(), sizeof(cuDoubleComplex)*grid.width()*grid.height()*grid.depth(), grid.int_dim());
	return cgrid;
}

ComplexGrid &copyDeviceToHost (ComplexGrid &grid, const CudaComplexGrid &cgrid)
{
	if ((cgrid.int_dim() != grid.int_dim()) || (cgrid.width() != grid.width()) || (cgrid.height() != grid.height()) || (cgrid.depth() != grid.depth()))
	{
		cout << "Warning resizing ComplexGrid for copying!" << endl;
		grid.resize(cgrid.int_dim(), cgrid.width(), cgrid.height(), cgrid.depth());
	}
	memcpypitch_device_to_host(grid.store, sizeof(cuDoubleComplex)*grid.width()*grid.height()*grid.depth(), cgrid.store, cgrid.get_pitch(), cgrid.width()*cgrid.height()*cgrid.depth()*sizeof(cuDoubleComplex), cgrid.int_dim());
	return grid;
}

#endif
