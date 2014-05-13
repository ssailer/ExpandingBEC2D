#include <iostream>
#include <EXP2D_observables.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

int main( int argc, char** argv) 
{
	Observables x(2);
	Observables y(2);
	Observables z(2);

	x.Ekin = 2;
	y.Ekin = 2;
	x.particle_count = 3;
	y.particle_count = 3;
	x.healing_length = 4;
	y.healing_length = 4;
	x.number(0) = 5;
	y.number(0) = 5;
	x.number(1) = 6;
	y.number(1) = 6;
	x.k(0) = 7;
	y.k(0) = 7;
	x.k(1) = 8;
	y.k(1) = 8;

	// z = x + y;
	z = x * 3;
	// z = x * 2;
	// z = z * 2.0;
	cout << "Test hier: " << z.Ekin << endl << z.particle_count << endl << z.healing_length << endl << z.number << endl << z.k  << endl << "Ende." << endl;

	// ArrayXd number(2);
	// number(0) = 2;
	// number(1) = 2;
	// number = number / 2;
	// cout << "Array " << endl << number << endl;


	return 0;
}