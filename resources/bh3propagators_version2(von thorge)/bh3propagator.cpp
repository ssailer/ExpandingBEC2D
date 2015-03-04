// includes, system
#include <fstream>
#include <omp.h>

// includes, kernels
#include "bh3propagator.h"
#include "complexgrid.h"

using namespace std;

Bh3Propagator::Bh3Propagator(const PathOptions &opt, const ComplexGrid &start1, const ComplexGrid &start2)
{
	options = opt;
	rgrid1.resize(options.delta_t.size() + 1);
	kgrid1.resize(options.delta_t.size() + 1);
	rgrid2.resize(options.delta_t.size() + 1);
	kgrid2.resize(options.delta_t.size() + 1);
	
	annealing = true;
	
	// start-grid kopieren
	rgrid1[0] = start1;
	rgrid2[0] = start2;

	current_time = 0.0;
	Omega_t = 0.0;
	
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
	double offset_time = 0.0;
	
	if(!propagateToTime(offset_time)) //evolution to properly mixed system;
		return false;
	
	
	for(std::vector<double>::const_iterator it = snapshots.begin(); it != snapshots.end(); ++it)
	{
		//cout << "Thread " << omp_get_thread_num() << ": Propagation to " << *it << "." << endl;
		if(!propagateToTime(*it + offset_time))
			return false;
		
		file.append_snapshot(*it, rgrid1);
		file.append_snapshot(Omega_t, rgrid2);
	}
	return true;
}

bool Bh3Propagator::propagateToTime(double time)
{
	if(!propagateN(((int) (time / options.timestepsize)) - ((int) (current_time / options.timestepsize))))
		return false;
	current_time = time;
	//cout<<current_time<<endl;
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
		
		rgrid1[n] = rgrid1[0];
		kgrid1[n] = kgrid1[0];
		rgrid2[n] = rgrid2[0];
		kgrid2[n] = kgrid2[0];
		N = steps = N - steps;
	}
	return true;
}

