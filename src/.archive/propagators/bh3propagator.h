#ifndef BH3PROPAGATOR_H__
#define BH3PROPAGATOR_H__

#include <vector>
#include <string>
#include <ostream>
#include <stdint.h>

#include <bh3binaryfile.h>

class ComplexGrid;

class Bh3Propagator {
	protected:
		std::vector<ComplexGrid> rgrid;
		std::vector<ComplexGrid> kgrid;
		
        
        PathOptions options;
		double current_time;
		
		std::vector<int> delta_N;
		
	public:
		Bh3Propagator(const PathOptions &opt, const ComplexGrid &start);
		virtual ~Bh3Propagator();
		
		virtual bool start(const std::vector<double> &snapshots, const std::string &filename); //==============was bringt virtual hier??? Polymorphie !!!============================//
		virtual bool propagateToTime(double time);
		virtual const vector<ComplexGrid> &getRGrid() const;
				
  protected:
		virtual bool propagateN(int N);
		virtual bool propagate1() = 0;
};

inline const vector<ComplexGrid> &Bh3Propagator::getRGrid() const
{
    return rgrid;
}

#endif

