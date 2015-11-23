#ifndef CUDACOMPLEXGRID_H__
#define CUDACOMPLEXGRID_H__


#include <cufft.h>
#include <stdint.h>

class ComplexGrid;

class CudaComplexGrid {
	private:
		uint32_t dim[3];
		uint32_t internal;
		cuDoubleComplex *store;
		size_t pitch;
		cufftHandle plan;
		
		void createPlan(const size_t ipitch, const size_t opitch);
		void initialize(uint32_t int_dim, uint32_t w, uint32_t h, uint32_t d);
		void deallocate();
	public:
		CudaComplexGrid();
		CudaComplexGrid(const CudaComplexGrid &grid);
		CudaComplexGrid(uint32_t int_dim, uint32_t w, uint32_t h, uint32_t d);
		~CudaComplexGrid();
		
		static bool fft(CudaComplexGrid &i, CudaComplexGrid &o, int direction);
		
		inline cuDoubleComplex* getDevicePointer() const
		{
			return store;
		}
		
		inline size_t get_pitch() const
		{
			return pitch;
		}
		
		inline uint32_t width() const {return dim[0];}
		inline uint32_t height() const {return dim[1];}
		inline uint32_t depth() const {return dim[2];}
		inline uint32_t int_dim() const {return internal;}
		
		friend CudaComplexGrid &copyHostToDevice (CudaComplexGrid &cgrid, const ComplexGrid &grid);
		friend ComplexGrid &copyDeviceToHost (ComplexGrid &grid, const CudaComplexGrid &cgrid);
		
		inline CudaComplexGrid &operator= (const ComplexGrid &grid)
		{
			return copyHostToDevice(*this, grid);
		}
		
		CudaComplexGrid &operator= (const CudaComplexGrid &grid);
		void resize(uint32_t int_dim, uint32_t w, uint32_t h, uint32_t d);
};

#endif
