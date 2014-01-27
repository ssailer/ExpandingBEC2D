#ifndef BH2PROPAGATOR_H__
#define BH2PROPAGATOR_H__

#include <vector>
#include <string>
#include <ostream>
#include <stdint.h>

#include <bh3binaryfile.h>

class ComplexGrid;

class Bh3Propagator {
	protected:
		std::vector<ComplexGrid> rgrid1;
		std::vector<ComplexGrid> kgrid1;
		std::vector<ComplexGrid> rgrid2;
		std::vector<ComplexGrid> kgrid2;
		PathOptions options;
		double current_time;
		double Omega_t;
		std::vector<int> delta_N;
		bool annealing;
	public:
		Bh3Propagator(const PathOptions &opt, const ComplexGrid &start1, const ComplexGrid &start2);
		virtual ~Bh3Propagator();
		
		virtual bool start(const std::vector<double> &snapshots, const std::string &filename); //==============was bringt virtual hier???============================//
		virtual bool propagateToTime(double time);
		virtual const vector<ComplexGrid> &getRGrids(int d) const;
		virtual const vector<ComplexGrid> &getKGrids(int d) const;
		
		inline void switch_annealing(){annealing = !annealing;};
		inline void set_Omega(double Omega_start) {Omega_t = Omega_start;
									annealing = false;};
		inline double get_Omega(){return Omega_t;};
	protected:
		virtual bool propagateN(int N);
		virtual bool propagate1() = 0;
};

inline const vector<ComplexGrid> &Bh3Propagator::getRGrids(int d) const
{
	if (d == 1)
		return rgrid1;
	else 
		return rgrid2;
}

inline const vector<ComplexGrid> &Bh3Propagator::getKGrids(int d) const
{
	if (d == 1)
		return kgrid1;
	else 
		return kgrid2;
}

#endif

