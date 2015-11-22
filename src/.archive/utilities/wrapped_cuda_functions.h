#ifndef WRAPPED_CUDA_FUNCTIONS_H__
#define WRAPPED_CUDA_FUNCTIONS_H__

#include <cuda_runtime_api.h>
#include <omp.h>


inline bool init_cuda_device(int i = 0)
{
	cudaError_t result;
	switch(result = cudaSetDevice(i))
	{
		case cudaSuccess:
			cout << "Cuda Device " << i << " chosen for thread " << omp_get_thread_num() << "." << endl;
			if (cudaSetDeviceFlags(cudaDeviceMapHost) != cudaSuccess)
				cout << "Activation of Host to Device memory mapping failed. Propagation algorithm will most likely fail." << endl;
			return true;
		case cudaErrorInvalidDevice:
			cout << "Error: Could not set device " << i << " for thread " << omp_get_thread_num() << "! Invalid device!" << endl;
			return false;
		case cudaErrorSetOnActiveProcess:
			cout << "Error: Could not set device " << i << " for thread " << omp_get_thread_num() << "! Cuda already initialized!" << endl;
			return false;
		default:
			cout << "Internal Bug when setting CUDA-Device!" << endl;
			cout << "Error: " << cudaGetErrorString(result) << endl;
			cout << "Last Error: " << cudaGetErrorString(cudaGetLastError()) << endl;
			return false;
	}
}

inline void thread_synchronize(const char *error_string)
{
	cudaError_t res;
	if(cudaThreadSynchronize() != cudaSuccess)
	{
		cout << "Error synchronizing thread " << omp_get_thread_num() << ": " << error_string << "!" << endl;
	}
	if((res = cudaGetLastError()) != cudaSuccess)
	{
		cout << "Error catched when synchronizing thread " << omp_get_thread_num() << ": " << error_string << ": " << cudaGetErrorString(res) << endl;
	}
}

inline bool memcpy_host_to_device(void *dst, const void *src, size_t n)
{
	switch(cudaMemcpy(dst, src, n, cudaMemcpyHostToDevice))
	{
		case cudaSuccess:
			break;
		case cudaErrorInvalidValue:
			cout << "Error copying from host to device in thread " << omp_get_thread_num() << ": Invalid value!" << endl;
			return false;
		case cudaErrorInvalidDevicePointer:
			cout << "Error copying from host to device in thread " << omp_get_thread_num() << ": Invalid device pointer!" << endl;
			return false;
		case cudaErrorInvalidMemcpyDirection:
			cout << "Error copying from host to device in thread " << omp_get_thread_num() << ": Invalid memcpy direction!" << endl;
			return false;
		default:
			cout << "Internal bug in thread " << omp_get_thread_num() << " while copying from host to device!" << endl;
			return false;
	}
	return true;
}

inline bool memcpypitch_host_to_device(void *dst, size_t dstpitch, const void *src, size_t srcpitch, size_t width, size_t height)
{
	switch(cudaMemcpy2D(dst, dstpitch, src, srcpitch, width, height, cudaMemcpyHostToDevice))
	{
		case cudaSuccess:
			break;
		case cudaErrorInvalidValue:
			cout << "Error copying from host to device in thread " << omp_get_thread_num() << ": Invalid value!" << endl;
			return false;
		case cudaErrorInvalidDevicePointer:
			cout << "Error copying from host to device in thread " << omp_get_thread_num() << ": Invalid device pointer!" << endl;
			return false;
		case cudaErrorInvalidMemcpyDirection:
			cout << "Error copying from host to device in thread " << omp_get_thread_num() << ": Invalid memcpy direction!" << endl;
			return false;
		default:
			cout << "Internal bug in thread " << omp_get_thread_num() << " while copying from host to device!" << endl;
			return false;
	}
	return true;
}

inline bool memcpy_device_to_host(void *dst, const void *src, size_t n)
{
	switch(cudaMemcpy(dst, src, n, cudaMemcpyDeviceToHost))
	{
		case cudaSuccess:
			break;
		case cudaErrorInvalidValue:
			cout << "Error copying from device to host in thread " << omp_get_thread_num() << ": Invalid value!" << endl;
			return false;
		case cudaErrorInvalidDevicePointer:
			cout << "Error copying from device to host in thread " << omp_get_thread_num() << ": Invalid device pointer!" << endl;
			return false;
		case cudaErrorInvalidMemcpyDirection:
			cout << "Error copying from device to host in thread " << omp_get_thread_num() << ": Invalid memcpy direction!" << endl;
			return false;
		default:
			cout << "Internal bug in thread " << omp_get_thread_num() << " while copying from device to host!" << endl;
			return false;
	}
	return true;
}

inline bool memcpypitch_device_to_host(void *dst, size_t dstpitch, const void *src, size_t srcpitch, size_t width, size_t height)
{
	switch(cudaMemcpy2D(dst, dstpitch, src, srcpitch, width, height, cudaMemcpyDeviceToHost))
	{
		case cudaSuccess:
			break;
		case cudaErrorInvalidValue:
			cout << "Error copying from device to host in thread " << omp_get_thread_num() << ": Invalid value!" << endl;
			return false;
		case cudaErrorInvalidDevicePointer:
			cout << "Error copying from device to host in thread " << omp_get_thread_num() << ": Invalid device pointer!" << endl;
			return false;
		case cudaErrorInvalidMemcpyDirection:
			cout << "Error copying from device to host in thread " << omp_get_thread_num() << ": Invalid memcpy direction!" << endl;
			return false;
		default:
			cout << "Internal bug in thread " << omp_get_thread_num() << " while copying from device to host!" << endl;
			return false;
	}
	return true;
}

template <class T>
inline bool memcpy_host_to_symbol(const T & dst, const void *src, size_t n, size_t offset)
{
	switch(cudaMemcpyToSymbol(dst, src, n, offset, cudaMemcpyHostToDevice))
	{
		case cudaSuccess:
			break;
		case cudaErrorInvalidValue:
			cout << "Error copying from host to symbol in thread " << omp_get_thread_num() << ": Invalid value!" << endl;
			return false;
		case cudaErrorInvalidSymbol:
			cout << "Error copying from host to symbol in thread " << omp_get_thread_num() << ": Invalid symbol!" << endl;
			return false;
		case cudaErrorInvalidDevicePointer:
			cout << "Error copying from host to symbol in thread " << omp_get_thread_num() << ": Invalid device pointer!" << endl;
			return false;
		case cudaErrorInvalidMemcpyDirection:
			cout << "Error copying from host to symbol in thread " << omp_get_thread_num() << ": Invalid memcpy direction!" << endl;
			return false;
		default:
			cout << "Internal bug in thread " << omp_get_thread_num() << " while copying from host to symbol!" << endl;
			return false;
	}
	return true;
}

#endif
