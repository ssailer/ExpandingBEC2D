// includes, system
#include <string.h>
#include <math.h>
#include <iostream>
#include <omp.h>

// includes, kernels
#include "bh3cpuomppropagator.h"
#include "complexgrid.h"

using namespace std;

Bh3CPUOMPPropagator::Bh3CPUOMPPropagator(const PathOptions &opt, const ComplexGrid &start) :
Bh3Propagator(opt,start)
{
	kprop = new ComplexGrid(opt.grid[0], opt.grid[1], opt.grid[2]);
	
	for(int x = 0; x < kprop->width(); x++)
	{
		for(int y = 0; y < kprop->height(); y++)
		{
			for(int z = 0; z < kprop->depth(); z++)
			{
				double k[3];
				k[0] = opt.klength[0] * sin(M_PI * x / (double) opt.grid[0]);
				k[1] = opt.klength[1] * sin(M_PI * y / (double) opt.grid[1]);
				k[2] = opt.klength[2] * sin(M_PI * z / (double) opt.grid[2]);
				
				double T = - (k[0] * k[0] + k[1] * k[1] + k[2] * k[2]) * opt.timestepsize;
				kprop->at(x,y,z) = exp(complex<double>(0,T)) / (double) (options.grid[0]*options.grid[1]*options.grid[2]);
				
			}
		}
	}
}

Bh3CPUOMPPropagator::~Bh3CPUOMPPropagator()
{
	if(kprop)
		delete kprop;
}

bool Bh3CPUOMPPropagator::propagate1()
{
	ComplexGrid::fft_unnormalized(rgrid[0], kgrid[0], true);
	for(int z = 0; z < kgrid[0].depth(); z++)
	{
		#pragma omp parallel for schedule(guided,4)
		for(int x = 0; x < kgrid[0].width(); x++)
		{
			for(int y = 0; y < kgrid[0].height(); y++)
			{
				kgrid[0](x,y,z) = kprop->at(x,y,z) * kgrid[0](x,y,z);

			}
		}
	}
	
	ComplexGrid::fft_unnormalized(kgrid[0], rgrid[0], false);
	double factor = -options.U * options.timestepsize;
	for(int z = 0; z < rgrid[0].depth(); z++)
	{
		#pragma omp parallel for schedule(guided,4)
		for(int x = 0; x < rgrid[0].width(); x++)
		{
			for(int y = 0; y < rgrid[0].height(); y++)
			{
				complex<double> value = rgrid[0](x,y,z);
				double V = abs2(value) * factor;
				rgrid[0](x,y,z) = complex<double>(cos(V), sin(V)) * value;
			}
		}
	}
	
	return true;
}


