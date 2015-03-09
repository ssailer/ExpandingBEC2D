/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and 
 * proprietary rights in and to this software and related documentation. 
 * Any use, reproduction, disclosure, or distribution of this software 
 * and related documentation without an express license agreement from
 * NVIDIA Corporation is strictly prohibited.
 *
 * Please refer to the applicable NVIDIA end user license agreement (EULA) 
 * associated with this source code for terms and conditions that govern 
 * your use of this NVIDIA software.
 * 
 */

/* Template project which demonstrates the basics on how to setup a project 
* example application.
* Host code.
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
//#include <cuda_runtime.h>

// includes, kernels
#include "bh3fastcudapropagator.h"
#include <complexgrid.h>
#include <wrapped_cuda_functions.h>

#define BLOCK_1D_LENGTH 512
#define KP_LENGTH 128
#define fix_N 15

// Kernels for r- and k-propagation
#include <cuComplex.h>

using namespace std;

typedef struct {
	cuDoubleComplex kx[KP_LENGTH];
	cuDoubleComplex ky[KP_LENGTH];
	cuDoubleComplex kz[KP_LENGTH];
	double timestepsize;
	double Ut;
	int N;
} KParam;

__constant__ __device__ KParam kp;

__device__ static __inline__ double cuCabs2 (cuDoubleComplex x)
{
	return cuCreal(x)*cuCreal(x) + cuCimag(x)*cuCimag(x);
}

extern __shared__ double base[];

#define result_real(index) base[2 * blockDim.x + index]
#define result_imag(index) base[3 * blockDim.x + index]
#define temp_real(index) base[index]
#define temp_imag(index) base[blockDim.x + index]

__inline__ __device__ void fastprop()
{
	result_real(threadIdx.x) = temp_real(threadIdx.x);
	result_imag(threadIdx.x) = temp_imag(threadIdx.x);
	__syncthreads();
	
	for(int i = 1; i <= kp.N; i++)
	{
		cuDoubleComplex t;
		int ind;
		ind = threadIdx.x - 1;
		if(threadIdx.x == 0)
			ind = blockDim.x - 1;
		t.x = temp_real(ind);
		t.y = temp_imag(ind);
		t.x += -2.0 * temp_real(threadIdx.x);
		t.y += -2.0 * temp_imag(threadIdx.x);
		ind = threadIdx.x + 1;
		if(threadIdx.x == blockDim.x - 1)
			ind = 0;
		t.x += temp_real(ind);
		t.y += temp_imag(ind);
		t.x *= kp.timestepsize / i;
		t.y *= kp.timestepsize / i;
		__syncthreads();
		temp_real(threadIdx.x) = - t.y;
		temp_imag(threadIdx.x) = t.x;
		result_real(threadIdx.x) -= t.y;
		result_imag(threadIdx.x) += t.x;
		__syncthreads();
	}
}

////////////////////////////////////////////////////////////////////////////////
//! propagator in r-space
//! @param g_idata  input data in global memory
//! @param g_odata  output data in global memory
////////////////////////////////////////////////////////////////////////////////
__global__ void short_fast_zpropagate(cuDoubleComplex *grid)
{
	// grid.z * grid.y * x + grid.z * y + z		: z = threadIdx.x; y = blockIdx.y; x = blockIdx.x;
	// => blockDim.y = 1; blockDim.x = grid.z; gridDim.y = grid.y; gridDim.x = grid.x;
	const unsigned int index = blockDim.x * gridDim.y * blockIdx.x + blockDim.x * blockIdx.y + threadIdx.x;
	
	temp_real(threadIdx.x) = grid[index].x;
	temp_imag(threadIdx.x) = grid[index].y;
	fastprop();
	grid[index].x = result_real(threadIdx.x);
	grid[index].y = result_imag(threadIdx.x);
}

__global__ void short_fast_ypropagate(cuDoubleComplex *grid)
{
	// grid.z * grid.y * x + grid.z * y + z		: z = blockIdx.y; y = threadIdx.x; x = blockIdx.x;
	// => blockDim.y = 1; blockDim.x = grid.y; gridDim.y = grid.z; gridDim.x = grid.x;
	const unsigned int index = gridDim.y * blockDim.x * blockIdx.x + gridDim.y * threadIdx.x + blockIdx.y;
	
	temp_real(threadIdx.x) = grid[index].x;
	temp_imag(threadIdx.x) = grid[index].y;
	fastprop();
	grid[index].x = result_real(threadIdx.x);
	grid[index].y = result_imag(threadIdx.x);
}

__global__ void short_fast_xpropagate(cuDoubleComplex *grid)
{
	// grid.z * grid.y * x + grid.z * y + z		: z = blockIdx.y; y = blockIdx.x; x = threadIdx.x;
	// => blockDim.y = 1; blockDim.x = grid.x; gridDim.y = grid.z; gridDim.x = grid.y;
	const unsigned int index = gridDim.y * gridDim.x * threadIdx.x + gridDim.y * blockIdx.x + blockIdx.y;
	
	temp_real(threadIdx.x) = grid[index].x;
	temp_imag(threadIdx.x) = grid[index].y;
	fastprop();
	// calculate the propagator and propagate
	cuDoubleComplex t;
	cuDoubleComplex t2;
	t.x = result_real(threadIdx.x);
	t.y = result_imag(threadIdx.x);
	t2 = t;
	t.x *= t.x;
	t.y *= t.y;
	t.x += t.y;
	t.x *= kp.Ut;
	sincos(t.x, &t.y, &t.x);
	grid[index] = cuCmul(t, t2);
}

// Bh3CudaPropagator - Class implementation

Bh3FastCudaPropagator::Bh3FastCudaPropagator(const PathOptions &opt, const ComplexGrid &start) :
			Bh3Propagator(opt, start)
{
	if((opt.grid[0] > BLOCK_1D_LENGTH) ||
		(opt.grid[1] > BLOCK_1D_LENGTH) ||
		(opt.grid[2] > BLOCK_1D_LENGTH))
	{
		cout << "Warning: invalid grid-sizes: Must be <= BLOCK_1D_LENGTH=" << BLOCK_1D_LENGTH  << " !" << endl;
	}
	
	dev_rgrid = new CudaComplexGrid(opt.grid[0], opt.grid[1], opt.grid[2]);
	
	KParam p;
	
	p.timestepsize = options.timestepsize;
	p.Ut = - options.U * options.timestepsize;
	p.N = fix_N;
	memcpy_host_to_symbol("kp", &p, sizeof(KParam), 0);
	
	*dev_rgrid = rgrid[0];
}

Bh3FastCudaPropagator::~Bh3FastCudaPropagator()
{
	delete dev_rgrid;
}

bool Bh3FastCudaPropagator::propagate1()
{
	if(options.grid[2] > 1)
	{
		dim3 dimBlock(options.grid[2], 1);
		dim3 dimGrid(options.grid[0],options.grid[1]);
		thread_synchronize("before z");
		short_fast_zpropagate<<<dimGrid, dimBlock, options.grid[2] * 2 * sizeof(cuDoubleComplex)>>>(dev_rgrid->getDevicePointer());
		thread_synchronize("after z");
	}
	if(options.grid[1] > 1)
	{
		dim3 dimBlock(options.grid[1], 1);
		dim3 dimGrid(options.grid[0], options.grid[2]);
		thread_synchronize("before y");
		short_fast_ypropagate<<<dimGrid, dimBlock, options.grid[1] * 2 * sizeof(cuDoubleComplex)>>>(dev_rgrid->getDevicePointer());
		thread_synchronize("after y");
	}
	if(options.grid[0] > 1)
	{
		dim3 dimBlock(options.grid[0], 1);
		dim3 dimGrid(options.grid[1], options.grid[2]);
		thread_synchronize("before x");
		short_fast_xpropagate<<<dimGrid, dimBlock, options.grid[0] * 2 * sizeof(cuDoubleComplex)>>>(dev_rgrid->getDevicePointer());
		thread_synchronize("after x");
	}
	return true;
}

bool Bh3FastCudaPropagator::propagateN(int N)
{
	cout << "Starting propagation to " << N << endl;
	int steps = N;
	for(int n = delta_N.size() - 1; n >= 0; n--)
	{
		steps -= delta_N[n];
		if(steps < 0)
			steps = 0;
		for(int i = 0; i < steps; i++)
		{
			if(!propagate1())
				return false;
		}
		
		rgrid[n] = *dev_rgrid;
		N = steps = N - steps;
	}
	cout << "Finishing propagation to " << N << endl;
	return true;
}
