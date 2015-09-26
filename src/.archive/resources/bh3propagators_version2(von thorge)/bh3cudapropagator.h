#ifndef BH2CUDAPROPAGATOR_H__
#define BH2CUDAPROPAGATOR_H__

#include "cudacomplexgrid.h"
#include "bh3propagator.h"

class Bh3CPUPropagator : public Bh3Propagator {
	protected:
		CudaComplexGrid *dev_rgrid1;
		CudaComplexGrid *dev_rgrid2;
		CudaComplexGrid *dev_kgrid1;
		CudaComplexGrid *dev_kgrid2;
		CudaComplexGrid *dev_kprop;
		CudaComplexGrid *dev_kprop_imag;
		double Omega_cr;
		double tau_q;

	public:
		Bh3CPUPropagator(const PathOptions &opt, const ComplexGrid &start1, const ComplexGrid &start2, const double t_q);
		virtual ~Bh3CPUPropagator();
              
		
	protected:
		bool propagateN(int N);
		bool propagate1();
                bool imag_Time_Prop();
                bool Normfaktor();
		bool noisefunction();
};



#endif
 
