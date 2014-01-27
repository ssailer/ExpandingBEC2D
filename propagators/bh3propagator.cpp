// includes, system
#include <fstream>
#include <omp.h>

// includes, kernels
#include "bh3propagator.h"
#include "complexgrid.h"

using namespace std;

Bh3Propagator::Bh3Propagator(const PathOptions &opt, const ComplexGrid &start)
{
	options = opt;
	rgrid.resize(options.delta_t.size() + 1);

	// start-grid kopieren
	rgrid[0] = start;
	
	current_time = 0.0;
	
	delta_N.resize(options.delta_t.size() + 1, 0);
	for(int i = 1; i < delta_N.size(); i++)
	{
		delta_N[i] = options.delta_t[i-1] / options.timestepsize;
	
	}
}

Bh3Propagator::~Bh3Propagator()
{
}

bool Bh3Propagator::start(const std::vector<double> &snapshots, const string &filename)
{
	Bh3BinaryFile file(filename, options, Bh3BinaryFile::out);
        	
	for(std::vector<double>::const_iterator it = snapshots.begin(); it != snapshots.end(); ++it)
	{
		//cout << "Thread " << omp_get_thread_num() << ": Propagation to " << *it << "." << endl;
		if(!propagateToTime(*it))
			return false;
		
		file.append_snapshot(*it, rgrid);
    }
	return true;
}

bool Bh3Propagator::propagateToTime(double time)
{
  	if(!propagateN(((int) (time / options.timestepsize)) - ((int) (current_time / options.timestepsize))))
		return false;

	current_time = time;

	return true;
}

bool Bh3Propagator::propagateN(int N)
{ 
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
		
		rgrid[n] = rgrid[0];
	
		N = steps = N - steps;
	}
	return true;
}

