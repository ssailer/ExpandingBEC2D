#include <inttypes.h>
#include <iostream>
#include <stdio.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;


int main() {
  // printf( "    short int: %zd\n" , sizeof(short int) ) ;
  // printf( "          int: %zd\n" , sizeof(int) ) ;
  // printf( "     long int: %zd\n", sizeof(long int) ) ;
  // printf( "long long int: %zd\n", sizeof(long long int) ) ;
  // printf( "       size_t: %zd\n", sizeof(size_t) ) ;
  // printf( "        void*: %zd\n\n", sizeof(void *) ) ;


  // printf( "PRIu32 usage (see source): %" PRIu32 "\n" , (uint32_t) 42 ) ;





  MatrixXcd m1(3,3);
  m1 << 1,2,3,
  		4,5,6,
  		7,8,9;

  MatrixXcd m2(3,3);
  m2 << 10,11,12,
        13,14,15,
        16,17,18;

  VectorXcd v1(3);
  v1 << 1,2,3;

  cout << "m1 * m2 = " << endl << m1*m2 << endl;
  cout << "m1.array() * m2.array() = " << endl << m1.array() * m2.array() << endl;
  cout << "m1.cwiseProduct(m2) = " << endl << m1.cwiseProduct(m2) << endl;

  return 0;
}