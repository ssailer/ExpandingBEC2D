#include "gauss_random.h"
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <iostream>

GaussRandom::GaussRandom()
{
  r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r, time(NULL));
}

GaussRandom::GaussRandom(unsigned long int seed)
{
  r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r, time(NULL) + seed);
}


GaussRandom::~GaussRandom()
{
  gsl_rng_free(r);
}



std::complex<double> GaussRandom::gauss_random()
{
	double x1, x2, w, y1, y2;
	
	do
	{
		x1 = 2.0 * double_random() - 1.0;
		x2 = 2.0 * double_random() - 1.0;
		w = x1*x1 + x2*x2;
	} while (w >= 1.0);
	
	w = sqrt((-2.0 * log(w)) / w);
	y1 = x1 * w;
	y2 = x2 * w;
	
	return std::complex<double>(y1,y2);
}

std::complex<double> GaussRandom::gauss_random(double mu, double sigma)
{
    double xi1 = double_random();
    double xi2 = double_random();
    
    double x1 = sqrt(-2.*log(xi1))*sin(M_PI*2.*xi2);
    double x2 = sqrt(-2.*log(xi1))*cos(M_PI*2.*xi2);

    return std::complex<double>(mu + sigma*x1, mu + sigma*x2);
}


void init_random()
{
	srandom(time(NULL));
}

unsigned long int get_seed()
{
	return (unsigned long int) random();
}
