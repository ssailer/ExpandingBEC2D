#ifndef CUDAREALGRID_CU__
#define CUDAREALGRID_CU__

#include <iostream>
#include <omp.h>

#include "cudarealgrid.h"
#include "realgrid.h"
#include "wrapped_cuda_functions.h"

CudaRealGrid::CudaRealGrid()
{
	dim[0] = dim[1] = dim[2] = internal = fft_dim[0] = fft_dim[1] = fft_dim[2] = 0;
	store = NULL;
}

CudaRealGrid::CudaRealGrid(uint32_t internal, uint32_t w, uint32_t h, uint32_t d)
{
	initialize(internal, w, h, d);
}

CudaRealGrid::CudaRealGrid(const CudaRealGrid &grid)
{
	initialize(internal, grid.dim[0], grid.dim[1], grid.dim[2]);
	cudaMemcpy(store, grid.store, internal*fft_dim[0]*fft_dim[1]*fft_dim[2]/2*sizeof(cufftDoubleComplex), cudaMemcpyDeviceToDevice);
}

CudaRealGrid::~CudaRealGrid()
{
	deallocate();
}

void CudaRealGrid::createPlan()
{
	cufftResult res_forward;
	cufftResult res_backward;

	uint32_t rank = dim[2]>1 ? 3 : ((dim[2]==1 && dim[1]>1) ? 2 : ((dim[2] == 1 && dim[1]==1 && dim[0] > 1) ? 1 : 3));
			
	int n[rank];
			
	if (rank==1)
		n[0] = dim[0];
	else if (rank ==2)
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

	res_forward = cufftPlanMany(&plan_forward, rank, n, NULL, 1, dim[0]*dim[1]*dim[2], NULL, 1, fft_dim[0]*fft_dim[1]*fft_dim[2]/2, CUFFT_D2Z, internal);
	res_backward = cufftPlanMany(&plan_backward, rank, n, NULL, 1, fft_dim[0]*fft_dim[1]*fft_dim[2]/2, NULL, 1, dim[0]*dim[1]*dim[2], CUFFT_Z2D, internal);

	if(res_forward != CUFFT_SUCCESS)
	{
		cout << "Creating CUFFT-Plan D2Z failed in thread " << omp_get_thread_num() << "!" << endl << "Reason:  ";
		switch(res_forward) {
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
	
	if(res_backward != CUFFT_SUCCESS)
	{
		cout << "Creating CUFFT-Plan Z2D failed in thread " << omp_get_thread_num() << "!" << endl << "Reason:  ";
		switch(res_backward) {
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
	
	if (cufftSetCompatibilityMode(plan_forward, CUFFT_COMPATIBILITY_NATIVE) != CUFFT_SUCCESS)
	{
		cout << "error occured while setting FFT compatibility to NATIVE for forward plan" << endl;
	}
	
	if (cufftSetCompatibilityMode(plan_backward, CUFFT_COMPATIBILITY_NATIVE) != CUFFT_SUCCESS)
	{
		cout << "error occured while setting FFT compatibility to NATIVE for backward plan" << endl;
	}
}

void CudaRealGrid::initialize(uint32_t int_dim, uint32_t w, uint32_t h, uint32_t d)
{
	dim[0] = w;
	dim[1] = h;
	dim[2] = d;
	internal = int_dim;
	if (d > 1)
	{
		fft_dim[0] = w;
		fft_dim[1] = h;
		fft_dim[2] = 2*(int(d/2) + 1);
 	}
	else if (d == 1 && h > 1)
	{
		fft_dim[0] = w;
		fft_dim[1] = 2*(int(h/2) + 1);
		fft_dim[2] = 1;
	}
	else if (d == 1 && h == 1 && w > 1)
	{
		fft_dim[0] = 2*(int(w/2) + 1);
		fft_dim[1] = 1;
		fft_dim[2] = 1;
	}
	else
	{
		fft_dim[0] = 1;
		fft_dim[1] = 1;
		fft_dim[2] = 1;
	}
	
	if(cudaMalloc(&store, sizeof(cufftDoubleComplex)*fft_dim[0]*fft_dim[1]*fft_dim[2]*internal/2) != cudaSuccess)
	{
		cout << "Error allocating Cuda-Memory in thread " << omp_get_thread_num() << "!" << endl;
	}
	createPlan();
}

void CudaRealGrid::deallocate()
{
	cufftDestroy(plan_forward);
	cufftDestroy(plan_backward);
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

bool CudaRealGrid::fft(CudaRealGrid &i, CudaRealGrid &o, int direction)
{
	int res_forward[4] = {CUFFT_SUCCESS, CUFFT_SUCCESS, CUFFT_SUCCESS, CUFFT_SUCCESS};
	int res_backward[4] = {CUFFT_SUCCESS, CUFFT_SUCCESS, CUFFT_SUCCESS, CUFFT_SUCCESS};
	
	if (direction == CUFFT_FORWARD)
	{
		int j = -1;
		do
		{
			j++;
			res_forward[j] = cufftExecD2Z(i.plan_forward, (double*)i.store, o.store);
		} while((res_forward[j] != CUFFT_SUCCESS) && (j < 3));
	}
	
	else if (direction == CUFFT_INVERSE)
	{
		int j = -1;
		do
		{
			j++;
			res_backward[j] = cufftExecZ2D(i.plan_backward, i.store, (double*)o.store);
		} while((res_backward[j] != CUFFT_SUCCESS) && (j < 3));
	}
	
	else
	{
		cout << "no valid direction chosen" << endl;
		return false;
	}
	
	if(res_forward[0] != CUFFT_SUCCESS)
	{
		cout << "Thread " << omp_get_thread_num() << ", Direction forward" << ":" << endl;
		for(int k = 0; k < 4; k++)
		{
			if(res_forward[k] != CUFFT_SUCCESS)
				cout << "Warning: Try " << k << " of fourier transform failed: " << cufft_error_string(res_forward[k]) << endl;
		}
		if((res_forward[1] != CUFFT_SUCCESS) && (res_forward[2] != CUFFT_SUCCESS) && (res_forward[3] != CUFFT_SUCCESS))
		{
			cout << "Error: All trys of fourier transforming failed in thread " << omp_get_thread_num() << "!" << endl;
			return false;
		}
	}
	
	if(res_backward[0] != CUFFT_SUCCESS)
	{
		cout << "Thread " << omp_get_thread_num() << ", Direction backward" << ":" << endl;
		for(int k = 0; k < 4; k++)
		{
			if(res_backward[k] != CUFFT_SUCCESS)
				cout << "Warning: Try " << k << " of fourier transform failed: " << cufft_error_string(res_backward[k]) << endl;
		}
		if((res_backward[1] != CUFFT_SUCCESS) && (res_backward[2] != CUFFT_SUCCESS) && (res_backward[3] != CUFFT_SUCCESS))
		{
			cout << "Error: All trys of fourier transforming failed in thread " << omp_get_thread_num() << "!" << endl;
			return false;
		}
	}
	
	return true;
}

void CudaRealGrid::resize (uint32_t internal, uint32_t w, uint32_t h, uint32_t d)
{
	deallocate();
	initialize(internal, w, h, d);
}

CudaRealGrid &CudaRealGrid::operator= (const CudaRealGrid &grid)
{
	if(&grid != this)
	{
		if((internal != grid.internal) || (dim[0] != grid.dim[0]) || (dim[1] != grid.dim[1]) || (dim[2] != grid.dim[2]))
			resize(grid.internal, grid.dim[0], grid.dim[1], grid.dim[2]);

		cudaMemcpy(store, grid.store, fft_dim[0]*fft_dim[1]*fft_dim[2]/2*sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice);
	}
	return *this;
}

CudaRealGrid &copyHostToDevice_as_complex (CudaRealGrid &cgrid, const RealGrid &grid)
{
	if ((cgrid.int_dim() != grid.int_dim()) || (cgrid.width() != grid.width()) || (cgrid.height() != grid.height()) || (cgrid.depth() != grid.depth()))
	{
		cout << "Warning resizing CudaRealGrid for copying!" << endl;
		cgrid.resize(grid.int_dim(), grid.width(), grid.height(), grid.depth());
	}
	double *copygrid = new double[grid.int_dim()*grid.width()*grid.height()*grid.depth()];


	memcpy_host_to_device(cgrid.store, grid.store, sizeof(complex<double>)*grid.int_dim()*grid.fft_width()*grid.fft_height()*grid.fft_depth()/2);

	
	return cgrid;
}


CudaRealGrid &copyHostToDevice3D (CudaRealGrid &cgrid, const RealGrid &grid)
{
	if ((cgrid.int_dim() != grid.int_dim()) || (cgrid.width() != grid.width()) || (cgrid.height() != grid.height()) || (cgrid.depth() != grid.depth()))
	{
		cout << "Warning resizing CudaRealGrid for copying!" << endl;
		cgrid.resize(grid.int_dim(), grid.width(), grid.height(), grid.depth());
	}
	double *copygrid = new double[grid.int_dim()*grid.width()*grid.height()*grid.depth()];
	//remove padding in grid by copying row by row to copygrid
	for(int mu = 0; mu < grid.int_dim(); mu++)
	{
		for(int x = 0; x < grid.width(); x++)
		{
			for(int y = 0; y < grid.height(); y++)
			{
				memcpy(copygrid + mu*cgrid.width()*cgrid.height()*cgrid.depth() + x*cgrid.height()*cgrid.depth() + y*cgrid.depth(), grid.store + mu*grid.fft_width()*grid.fft_height()*grid.fft_depth() + x*grid.fft_height()*grid.fft_depth() + y*grid.fft_depth(), grid.depth()*sizeof(double));
			}
		}
	}
	//finally copy unpadded grid to CUDA device
	memcpy_host_to_device(cgrid.store, copygrid, sizeof(double)*grid.int_dim()*grid.width()*grid.height()*grid.depth());
	//clean up
	delete copygrid;
	
	return cgrid;
}

CudaRealGrid &copyHostToDevice2D (CudaRealGrid &cgrid, const RealGrid &grid)
{
	if ((cgrid.int_dim() != grid.int_dim()) || (cgrid.width() != grid.width()) || (cgrid.height() != grid.height()) || (cgrid.depth() != grid.depth()))
	{
		cout << "Warning resizing CudaRealGrid for copying!" << endl;
		cgrid.resize(grid.int_dim(), grid.width(), grid.height(), grid.depth());
	}
	double *copygrid = new double[grid.int_dim()*grid.width()*grid.height()*grid.depth()];
	//remove padding in grid by copying row by row to copygrid
	for(int mu = 0; mu < grid.int_dim(); mu++)
	{
		for(int x = 0; x < grid.width(); x++)
		{
			memcpy(copygrid + mu*cgrid.width()*cgrid.height()*cgrid.depth() + x*cgrid.height()*cgrid.depth(), grid.store + mu*grid.fft_width()*grid.fft_height()*grid.fft_depth() + x*grid.fft_height()*grid.fft_depth(), grid.height()*sizeof(double));			
		}
	}
	//finally copy unpadded grid to CUDA device
	memcpy_host_to_device(cgrid.store, copygrid, sizeof(double)*grid.int_dim()*grid.width()*grid.height()*grid.depth());
	//clean up
	delete copygrid;
	
	return cgrid;
}

CudaRealGrid &copyHostToDevice1D (CudaRealGrid &cgrid, const RealGrid &grid)
{
	if ((cgrid.int_dim() != grid.int_dim()) || (cgrid.width() != grid.width()) || (cgrid.height() != grid.height()) || (cgrid.depth() != grid.depth()))
	{
		cout << "Warning resizing CudaRealGrid for copying!" << endl;
		cgrid.resize(grid.int_dim(), grid.width(), grid.height(), grid.depth());
	}
	double *copygrid = new double[grid.int_dim()*grid.width()*grid.height()*grid.depth()];
	//remove padding in grid by copying row by row to copygrid
	for(int mu = 0; mu < grid.int_dim(); mu++)
	{
		memcpy(copygrid + mu*cgrid.width()*cgrid.height()*cgrid.depth(), grid.store + mu*grid.fft_width()*grid.fft_height()*grid.fft_depth(), grid.width()*sizeof(double));
	}
	//finally copy unpadded grid to CUDA device
	memcpy_host_to_device(cgrid.store, copygrid, sizeof(double)*grid.int_dim()*grid.width()*grid.height()*grid.depth());
	//clean up
	delete copygrid;
	
	return cgrid;
}

RealGrid &copyDeviceToHost_as_complex (RealGrid &grid, const CudaRealGrid &cgrid)
{
	if ((cgrid.int_dim() != grid.int_dim()) || (cgrid.width() != grid.width()) || (cgrid.height() != grid.height()) || (cgrid.depth() != grid.depth()))
	{
		cout << "Warning resizing RealGrid for copying!" << endl;
		grid.resize(cgrid.int_dim(), cgrid.width(), cgrid.height(), cgrid.depth());
	}
	
	memcpy_device_to_host(grid.store, cgrid.store, sizeof(complex<double>)*grid.int_dim()*grid.fft_width()*grid.fft_height()*grid.fft_depth()/2);
	
	return grid;
}



RealGrid &copyDeviceToHost3D (RealGrid &grid, const CudaRealGrid &cgrid)
{
	if ((cgrid.int_dim() != grid.int_dim()) || (cgrid.width() != grid.width()) || (cgrid.height() != grid.height()) || (cgrid.depth() != grid.depth()))
	{
		cout << "Warning resizing RealGrid for copying!" << endl;
		grid.resize(cgrid.int_dim(), cgrid.width(), cgrid.height(), cgrid.depth());
	}
	
	double *copygrid = new double[cgrid.int_dim()*cgrid.width()*cgrid.height()*cgrid.depth()];
	//copy unpadded grid from CUDA device to copygrid on host
	memcpy_device_to_host(copygrid, cgrid.store, cgrid.int_dim()*cgrid.width()*cgrid.height()*cgrid.depth()*sizeof(double));
	//add padding by copying row by row to grid.store
	for(int mu = 0; mu < grid.int_dim(); mu++)
	{
		for(int x = 0; x < grid.width(); x++)
		{
			for(int y = 0; y < grid.height(); y++)
			{
				memcpy(grid.store + mu*grid.fft_width()*grid.fft_height()*grid.fft_depth() + x*grid.fft_height()*grid.fft_depth() + y*grid.fft_depth(), copygrid + mu*cgrid.width()*cgrid.height()*cgrid.depth() + x*cgrid.height()*cgrid.depth() + y*cgrid.depth(), cgrid.depth()*sizeof(double));
			}
		}
	}	
	//clean up
	delete copygrid;
	
	return grid;
}

RealGrid &copyDeviceToHost2D (RealGrid &grid, const CudaRealGrid &cgrid)
{
	if ((cgrid.int_dim() != grid.int_dim()) || (cgrid.width() != grid.width()) || (cgrid.height() != grid.height()) || (cgrid.depth() != grid.depth()))
	{
		cout << "Warning resizing RealGrid for copying!" << endl;
		grid.resize(cgrid.int_dim(), cgrid.width(), cgrid.height(), cgrid.depth());
	}
	
	double *copygrid = new double[cgrid.int_dim()*cgrid.width()*cgrid.height()*cgrid.depth()];
	//copy unpadded grid from CUDA device to copygrid on host
	memcpy_device_to_host(copygrid, cgrid.store, cgrid.int_dim()*cgrid.width()*cgrid.height()*cgrid.depth()*sizeof(double));
	//add padding by copying row by row to grid.store
	for(int mu = 0; mu < grid.int_dim(); mu++)
	{
		for(int x = 0; x < grid.width(); x++)
		{
			memcpy(grid.store + mu*grid.fft_width()*grid.fft_height()*grid.fft_depth() + x*grid.fft_height()*grid.fft_depth(), copygrid + mu*cgrid.width()*cgrid.height()*cgrid.depth() + x*cgrid.height()*cgrid.depth(), cgrid.height()*sizeof(double));
		}
	}	
	//clean up
	delete copygrid;
	
	return grid;
}

RealGrid &copyDeviceToHost1D (RealGrid &grid, const CudaRealGrid &cgrid)
{
	if ((cgrid.int_dim() != grid.int_dim()) || (cgrid.width() != grid.width()) || (cgrid.height() != grid.height()) || (cgrid.depth() != grid.depth()))
	{
		cout << "Warning resizing RealGrid for copying!" << endl;
		grid.resize(cgrid.int_dim(), cgrid.width(), cgrid.height(), cgrid.depth());
	}
	
	double *copygrid = new double[cgrid.int_dim()*cgrid.width()*cgrid.height()*cgrid.depth()];
	//copy unpadded grid from CUDA device to copygrid on host
	memcpy_device_to_host(copygrid, cgrid.store, cgrid.int_dim()*cgrid.width()*cgrid.height()*cgrid.depth()*sizeof(double));
	//add padding by copying row by row to grid.store
	for(int mu = 0; mu < grid.int_dim(); mu++)
	{
		memcpy(grid.store + mu*grid.fft_width()*grid.fft_height()*grid.fft_depth(), copygrid + mu*cgrid.width()*cgrid.height()*cgrid.depth(), cgrid.width()*sizeof(double));
	}	
	//clean up
	delete copygrid;
	
	return grid;
}

#endif
