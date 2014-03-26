#ifndef BH2CUDAPROPAGATOR_H__
#define BH2CUDAPROPAGATOR_H__

#include "cudacomplexgrid.h"
#include "bh3propagator.h"

class Bh3FastCudaPropagator : public Bh3Propagator {
	protected:
		CudaComplexGrid *dev_rgrid;
	public:
		Bh3FastCudaPropagator(const PathOptions &opt, const ComplexGrid &start);
		virtual ~Bh3FastCudaPropagator();
		
	protected:
		bool propagateN(int N);
		bool propagate1();
};

#endif
 
