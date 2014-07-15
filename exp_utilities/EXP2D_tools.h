#ifndef EXP2D_TOOLS_H__
#define EXP2D_TOOLS_H__

#include <iostream>
#include <complex>
#include <math.h>
#include <complexgrid.h>
#include <bh3binaryfile.h>
#include <gauss_random.h>
#include <vector>
#include <string>
#include <iomanip>
#include <gauss_random.h>
#include <stdlib.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

typedef struct {

    double N; // Number of particles    
    double klength[3];
    vector<double> stateInformation; // passing information about the state at the absolut time to the observable, lambda(time) FIXME : this is bad, but I don't know how to do it better atm
    complex<double> omega_x,omega_y; // Frequency of the harmonic trap
    complex<double> dispersion_x, dispersion_y; // dispersion relation for the expandion frame
    double min_x,min_y; // Coordinate boundaries    
    complex<double> t_abs; //Absolute time // remove from opt! put into the function, don't need it here
    complex<double> exp_factor; //Expansion factor
    double g; // coupling constant
    double ITP_step, RTE_step; // stepsize for the timeiteration

    int32_t grid[4];  // gridsize
    int n_it_RTE; // number of timesteps
    int samplesize;
    int vortexnumber;
    
    string runmode; // Use this to control the program flow: first char determines if the program is loading from a dataset or using ITP to generate the necessary datafile
                     // second char determines if expanding coordinates are used or not
                     // third char determines if potential is switch on for the differential equation
    string name; // naming of the datafile      // think about that naming system remove it from here
    string config; // name of the config file 
    string workingdirectory;   // remove it from here, only needed in the program itself
    string workingfile;
   

    double scale_factor; //Scale factor
    
} Options;

void optToPath(Options &opt,PathOptions &pathopt);
void pathToOpt(PathOptions &pathopt,Options &opt);
void readDataFromHDF5(ComplexGrid* &g,Options &opt);
void saveDataToHDF5(ComplexGrid* &g, Options &opt);
void noiseTheGrid(ComplexGrid &g);

void saveEigenMatrixToHDF5();
void loadEigenMatrixFromHDF5();

class expException {
public:
    inline expException(std::string const& info);    
    inline void setString(std::string const& info);
    inline void addString(std::string const& info);
    inline std::string printString();
    std::string stringException;
private:
    
};

inline expException::expException(std::string const& info){
    stringException = info;
}

inline void expException::setString(std::string const& info){
    stringException = info;
}

inline void expException::addString(std::string const& info){
    stringException += info;
}

inline std::string expException::printString(){
    cout << stringException.c_str() << endl;
}



#endif // EXP2D_TOOLS_H__