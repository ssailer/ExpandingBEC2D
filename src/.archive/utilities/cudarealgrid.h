#ifndef CUDAREALGRID_H__
#define CUDAREALGRID_H__


#include <cufft.h>
#include <stdint.h>

class RealGrid;

class CudaRealGrid {
	private:
		uint32_t dim[3];
		uint32_t fft_dim[3];
		uint32_t internal;
		cufftDoubleComplex *store;
		cufftHandle plan_forward;
		cufftHandle plan_backward;
		
		void createPlan();
		void initialize(uint32_t int_dim, uint32_t w, uint32_t h, uint32_t d);
		void deallocate();
	public:
		CudaRealGrid();
		CudaRealGrid(const CudaRealGrid &grid);
		CudaRealGrid(uint32_t internal, uint32_t w, uint32_t h, uint32_t d);
		~CudaRealGrid();
		
		static bool fft(CudaRealGrid &i, CudaRealGrid &o, int direction);
		
		inline cufftDoubleComplex* getDevicePointer()
		{
			return store;
		}
		
		inline uint32_t width() const {return dim[0];}
		inline uint32_t height() const {return dim[1];}
		inline uint32_t depth() const {return dim[2];}
		inline uint32_t int_dim() const {return internal;}
		
		friend CudaRealGrid &copyHostToDevice_as_complex (CudaRealGrid &cgrid, const RealGrid &grid);
		friend CudaRealGrid &copyHostToDevice3D (CudaRealGrid &cgrid, const RealGrid &grid);
		friend CudaRealGrid &copyHostToDevice2D (CudaRealGrid &cgrid, const RealGrid &grid);
		friend CudaRealGrid &copyHostToDevice1D (CudaRealGrid &cgrid, const RealGrid &grid);

		friend RealGrid &copyDeviceToHost_as_complex (RealGrid &grid, const CudaRealGrid &cgrid);
		friend RealGrid &copyDeviceToHost3D (RealGrid &grid, const CudaRealGrid &cgrid);
		friend RealGrid &copyDeviceToHost2D (RealGrid &grid, const CudaRealGrid &cgrid);
		friend RealGrid &copyDeviceToHost1D (RealGrid &grid, const CudaRealGrid &cgrid);
		
		inline CudaRealGrid &operator= (const RealGrid &grid)
		{
			return copyHostToDevice_as_complex(*this, grid);
		}
		
		CudaRealGrid &operator= (const CudaRealGrid &grid);
		void resize(uint32_t internal, uint32_t w, uint32_t h, uint32_t d);
};

#endif  //cudarealgrid
