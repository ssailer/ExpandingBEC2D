#ifndef _GAUSS_RANDOM_H__
#define _GAUSS_RANDOM_H__

#include <complex>
#include <gsl/gsl_rng.h>


class GaussRandom {

  private:
    gsl_rng *r;

  public:
    GaussRandom();
    GaussRandom(unsigned long int seed);
    ~GaussRandom();
    std::complex<double> gauss_random();
    std::complex<double> gauss_random(double mu, double sigma);
    

    double double_random()
    {
        return gsl_rng_uniform(r);
    };
    

    unsigned long int int_random(const unsigned long int bound)
    {
        return gsl_rng_uniform_int(r,bound);
    };
    

};

void init_random();
unsigned long int  get_seed();

#endif
