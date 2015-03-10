#include <inttypes.h>
#include <iostream>
#include <stdio.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

 const complex<double> one = complex<double>(1.0,0);

int main() {
  // printf( "    short int: %zd\n" , sizeof(short int) ) ;
  // printf( "          int: %zd\n" , sizeof(int) ) ;
  // printf( "     long int: %zd\n", sizeof(long int) ) ;
  // printf( "long long int: %zd\n", sizeof(long long int) ) ;
  // printf( "       size_t: %zd\n", sizeof(size_t) ) ;
  // printf( "        void*: %zd\n\n", sizeof(void *) ) ;


  // printf( "PRIu32 usage (see source): %" PRIu32 "\n" , (uint32_t) 42 ) ;

  cout << "TEST" << one << endl;




  return 0;
}