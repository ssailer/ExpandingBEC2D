#ifndef BH2CPUOMPPROPAGATOR_H__
#define BH2CPUOMPPROPAGATOR_H__

#include <math.h>
#include <complex>
#include "bh3propagator.h"

class Bh3CPUOMPPropagator : public Bh3Propagator {
	public:
		Bh3CPUOMPPropagator(const PathOptions &opt, const ComplexGrid &start);
		virtual ~Bh3CPUOMPPropagator();
	protected:
		bool propagate1();
		ComplexGrid *kprop;
};

#endif

