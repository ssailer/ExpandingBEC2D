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

#include <stdlib.h>
#include <time.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

typedef struct Options {

    

    double N; // Number of particles    
    double klength[3];
    vector<double> stateInformation; // passing information about the state at the absolut time to the observable, lambda(time) FIXME : this is bad, but I don't know how to do it better atm
    complex<double> omega_x,omega_y,omega_w; // Frequency of the harmonic trap
    complex<double> dispersion_x, dispersion_y; // dispersion relation for the expandion frame
    double min_x,min_y,min_z; // Coordinate boundaries    
    complex<double> t_abs; //Absolute time // remove from opt! put into the function, don't need it here
    complex<double> exp_factor; //Expansion factor
    double g; // coupling constant
    double ITP_step, RTE_step; // stepsize for the timeiteration

    int32_t grid[4];  // gridsize
    int n_it_RTE; // number of Iterations
    int snapshots; // number of Snapshots
    int samplesize;
    int vortexnumber;
    int vortexspacing;
    double potFactor;
    double Ag;
    double OmegaG;
    
    string runmode; // Use this to control the program flow: first char determines if the program is loading from a dataset or using ITP to generate the necessary datafile
                     // second char determines if expanding coordinates are used or not
                     // third char determines if potential is switch on for the differential equation
    string config; // name of the config file 
    string workingdirectory;   // remove it from here, only needed in the program itself
    bool initialRun;

    Options () : N(100000), stateInformation(2), t_abs(0,0), snapshots(100), vortexnumber(20), vortexspacing(50), potFactor(1), initialRun(true), config("sim.cfg") {}
    
} Options;


class expException {
    public:
        expException() : stringException("expException, nothing happened, yet.") {}
        expException(std::string const& info);
        void setString(std::string const& info);
        void addString(std::string const& info);
        void printString();
    private:
        std::string stringException;
};

inline expException::expException(std::string const& info) : stringException(info) {}    
inline void expException::setString(std::string const& info){
    stringException = info;
};
inline void expException::addString(std::string const& info){
    stringException += info;
};
inline void expException::printString(){
    cerr << stringException.c_str() << endl;
};

inline const std::string currentDate() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d", &tstruct);

    return buf;
}

inline const std::string currentTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%X", &tstruct);

    return buf;
}


inline void progressBar(const float progress){
    while (progress < 1.0) {
        int barWidth = 70;
    
        std::cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();
    }
    std::cout << std::endl;
}



#endif // EXP2D_TOOLS_H__