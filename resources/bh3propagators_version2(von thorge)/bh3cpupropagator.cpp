// includes, system
#include <string.h>
#include <math.h>
#include <iostream>

// includes, kernels
#include "bh3cpupropagator.h"
#include "complexgrid.h"
#include <bench_time.h>

using namespace std;

Bh3CPUPropagator::Bh3CPUPropagator(const PathOptions &opt, const ComplexGrid &start) :
							Bh3Propagator(opt,start)
{
	bench_fft = bench_rprop = bench_kprop = 0.0;
	bench_fft_sys = bench_rprop_sys = bench_kprop_sys = 0.0;
	
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

Bh3CPUPropagator::~Bh3CPUPropagator()
{
	cout << "Benchmark results:" << endl
		<< "FFT: " << bench_fft << "s" << " (" << bench_fft_sys << "s in sys)" << endl
		<< "R-Propagation" << bench_rprop << "s" << " (" << bench_rprop_sys << "s in sys)" << endl
		<< "K-Propagation" << bench_kprop << "s" << " (" << bench_kprop_sys << "s in sys)" << endl;
	if(kprop)
		delete kprop;
}

bool Bh3CPUPropagator::propagate1()
{
	double utime1, stime1, utime2, stime2;
	bench_time(utime1, stime1);
	ComplexGrid::fft_unnormalized(rgrid[0], kgrid[0], true);
	bench_time(utime2, stime2);
	bench_fft += utime2 - utime1;
	bench_fft_sys += stime2 - stime1;
	
	bench_time(utime1, stime1);
	for(int z = 0; z < kgrid[0].depth(); z++)
	{
		for(int x = 0; x < kgrid[0].width(); x++)
		{
			for(int y = 0; y < kgrid[0].height(); y++)
			{
				kgrid[0](x,y,z) = kprop->at(x,y,z) * kgrid[0](x,y,z);
			}
		}
	}
	bench_time(utime2, stime2);
	bench_kprop += utime2 - utime1;
	bench_kprop_sys += stime2 - stime1;
	
	bench_time(utime1, stime1);
	ComplexGrid::fft_unnormalized(kgrid[0], rgrid[0], false);
	bench_time(utime2, stime2);
	bench_fft += utime2 - utime1;
	bench_fft_sys += stime2 - stime1;
	
	bench_time(utime1, stime1);
	double factor = -options.U * options.timestepsize;
	for(int z = 0; z < rgrid[0].depth(); z++)
	{
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
	bench_time(utime2, stime2);
	bench_rprop += utime2 - utime1;
	bench_rprop_sys += stime2 - stime1;
	
	return true;
}


