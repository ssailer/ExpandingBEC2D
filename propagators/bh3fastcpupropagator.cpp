// includes, system
#include <string.h>
#include <math.h>
#include <iostream>

// includes, kernels
#include "bh3fastcpupropagator.h"
#include "complexgrid.h"
#include <bench_time.h>

using namespace std;

Bh3FastCPUPropagator::Bh3FastCPUPropagator(const PathOptions &opt, const ComplexGrid &start) :
							Bh3Propagator(opt,start)
{
	bench = 0.0;
}

Bh3FastCPUPropagator::~Bh3FastCPUPropagator()
{
	cout << "Benchmark results: " << bench << endl;
}

void Bh3FastCPUPropagator::kprop(ComplexGrid *g1, ComplexGrid *g2)
{
	const int depth = 15;
	*g2 = *g1;
	if(g1->depth() > 1)
	{
		for(int i = 1; i <= depth; i++)
		{
			double factor = options.timestepsize / (double) i * options.klength[2] / 2.0;
			for(int x = 0; x < g1->width(); x++)
			{
				for(int y = 0; y < g1->height(); y++)
				{
					complex<double> first = g1->at(x,y,0);
					complex<double> previous = g1->at(x,y,g1->depth()-1);
					complex<double> now = first;
					complex<double> temp;
					for(int z = 0; z < g1->depth()-1; z++)
					{
						complex<double> next = g1->at(x,y,z+1);
						temp = (previous + next - 2.0*now) * factor;
						g1->at(x,y,z).real() = - temp.imag();
						g1->at(x,y,z).imag() = temp.real();
						g2->at(x,y,z) += g1->at(x,y,z);
						previous = now;
						now = next;
					}
					temp = (previous + first - 2.0*now) * factor;
					g1->at(x,y,g1->depth()-1).real() = - temp.imag();
					g1->at(x,y,g1->depth()-1).imag() = temp.real();
					g2->at(x,y,g1->depth()-1) += g1->at(x,y,g1->depth()-1);
				}
			}
		}
		*g1 = *g2;
	}
	
	if(g1->height() > 1)
	{
		for(int i = 1; i <= depth; i++)
		{
			double factor = options.timestepsize / (double) i * options.klength[1] / 2.0;
			for(int z = 0; z < g1->depth(); z++)
			{
				for(int x = 0; x < g1->width(); x++)
				{
					complex<double> first = g1->at(x,0,z);
					complex<double> previous = g1->at(x,g1->height()-1,z);
					complex<double> now = first;
					complex<double> temp;
					for(int y = 0; y < g1->height()-1; y++)
					{
						complex<double> next = g1->at(x,y+1,z);
						temp = (previous + next - 2.0*now) * factor;
						g1->at(x,y,z).real() = - temp.imag();
						g1->at(x,y,z).imag() = temp.real();
						g2->at(x,y,z) += g1->at(x,y,z);
						previous = now;
						now = next;
					}
					temp = (previous + first - 2.0*now) * factor;
					g1->at(x,g1->height()-1,z).real() = - temp.imag();
					g1->at(x,g1->height()-1,z).imag() = temp.real();
					g2->at(x,g2->height()-1,z) += g1->at(x,g1->height()-1,z);
				}
			}
		}
		*g1 = *g2;
	}
	
	if(g1->width() > 1)
	{
		for(int i = 1; i <= depth; i++)
		{
			double factor = options.timestepsize / (double) i * options.klength[0] / 2.0;
			for(int z = 0; z < g1->depth(); z++)
			{
				for(int y = 0; y < g1->height(); y++)
				{
					complex<double> first = g1->at(0,y,z);
					complex<double> previous = g1->at(g1->width()-1,y,z);
					complex<double> now = first;
					complex<double> temp;
					for(int x = 0; x < g1->width()-1; x++)
					{
						complex<double> next = g1->at(x+1,y,z);
						temp = (previous + next - 2.0*now) * factor;
						g1->at(x,y,z).real() = - temp.imag();
						g1->at(x,y,z).imag() = temp.real();
						g2->at(x,y,z) += g1->at(x,y,z);
						previous = now;
						now = next;
					}
					temp = (previous + first - 2.0*now) * factor;
					g1->at(g1->width()-1,y,z).real() = - temp.imag();
					g1->at(g1->width()-1,y,z).imag() = temp.real();
					g2->at(g2->width()-1,y,z) += g1->at(g1->width()-1,y,z);
				}
			}
		}
		*g1 = *g2;
	}
}

bool Bh3FastCPUPropagator::propagate1()
{
	ComplexGrid *rg, *kg;
	rg = &rgrid[0];
	kg = &kgrid[0];
	double utime1, stime1, utime2, stime2;
	bench_time(utime1, stime1);
	kprop(rg, kg);
	double factor = -options.U * options.timestepsize;
	for(int z = 0; z < rg->depth(); z++)
	{
		for(int x = 0; x < rg->width(); x++)
		{
			for(int y = 0; y < rg->height(); y++)
			{
				
				
				complex<double> value = rg->at(x,y,z);
				double V = abs2(value) * factor;
				rg->at(x,y,z) = complex<double>(cos(V), sin(V)) * value;
			}
		}
	}
	bench_time(utime2, stime2);
	bench += utime2 - utime1;
	
	return true;
}


