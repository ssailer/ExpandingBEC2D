#include <boost/program_options.hpp>
#include <iostream>
#include <unistd.h>
#include <cstdlib>
#include <string>
#include <cmath>
#include <complex>
#include <omp.h>
#include <sys/stat.h>

#include <EXP2D_binaryfile.h>

#include <main.h>
#include <EXP2D_tools.h>

#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;


int main( int argc, char** argv) 
{	
	Options opt;

	read_cli_options(argc,argv,opt);
	read_config(argc,argv,opt);
	set_workingdirectory(opt);

	vector<MatrixXcd> testvector(2);
	testvector[0] = MatrixXcd(2,2);
	testvector[1] = MatrixXcd::Zero(3,3);

	testvector[0](0,0) = complex<double>(1.0,2.0);
	testvector[0](1,0) = complex<double>(3.0,4.0);
	testvector[0](0,1) = complex<double>(5.0,6.0);
	testvector[0](1,1) = complex<double>(7.0,8.0);

	cout << testvector[0] << endl;

	string filename = "testData.h5";

	binaryFile *dataFile = new binaryFile(filename,binaryFile::out);
	dataFile->appendSnapshot("RTE",10,testvector,opt);
	delete dataFile;

	cout << "Done with writing to h5" << endl;

	vector<MatrixXcd> testvector1(2);
	testvector1[0] = MatrixXcd(2,2);
	testvector1[1] = MatrixXcd::Zero(3,3);

	binaryFile *dataFile1 = new binaryFile(filename,binaryFile::in);
	dataFile1->getSnapshot("RTE",10,testvector1,opt);
	delete dataFile1;

	cout << testvector1[0] << endl;

	cout << "Done!" << endl;

	return 0;
}