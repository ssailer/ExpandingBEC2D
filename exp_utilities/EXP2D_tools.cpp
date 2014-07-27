#define EIGEN_VECTORIZE
#define EIGEN_NO_DEBUG

#include <EXP2D_tools.h>


void noiseTheGrid(ComplexGrid &g){
   GaussRandom r (get_seed());
   double rvalue;
   for(int i = 0;i < g.width();i++){
    for(int j = 0; j < g.height();j++){
        rvalue = real(g(0,i,j,0)) * 0.1;
        g(0,i,j,0) += r.gauss_random(0.0,rvalue);
    }
   }
}
