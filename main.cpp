/**************************************************************************
Title: Simulating the Expansion of Turbulent Bose-Einstein Condensates (2D) 
Author: Simon Sailer (This work is based on the work of Bartholomew Andrews who made this as his master thesis.)
Last Update: 22/07/13
**************************************************************************/

#include <boost/program_options.hpp>
#include <iostream>
#include <unistd.h>
#include <cstdlib>
// #include <cstring>
#include <string>
#include <cmath>
#include <complex>
#include <omp.h>
#include <sys/stat.h>

#include <complexgrid.h>
#include <bh3defaultgrid.h>
#include <averageclass.h>
#include <bh3observables.h>

#include <main.h>
#include <EXP2D_tools.h>
#include <EXP2D_itp.hpp>
#include <EXP2D_rte.hpp>
#include <EXP2D_observables.h>
#include <plot_with_mgl.h>

// #include <typeinfo>

#define SUCCESS 0; 
#define ERROR_IN_COMMAND_LINE 1;
#define ERROR_IN_CONFIG_FILE 2;
#define ERROR_UNHANDLED_EXCEPTION 3;

using namespace std;




int main( int argc, char** argv) 
{	
try{ 

Options opt;

read_cli_options(argc,argv,opt);
read_config(argc,argv,opt);
set_workingdirectory(opt);

// Initialize the needed grid object 
ComplexGrid* data = new ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
ITP* itprun = new ITP(data,opt);

printInitVar(opt); 

if(opt.runmode.compare(0,1,"0") == 0){

double sigma_real[2];
sigma_real[0] = opt.min_x/4;
sigma_real[1] = opt.min_y/4;		
data = set_grid_to_gaussian(data,opt,sigma_real[0],sigma_real[1]);

// set the datafile identifier name and save the initial grid

opt.name = "INIT";
plotdatatopng(opt.name,data,opt);

// //====> Imaginary Time Propagation (ITP)
itprun->propagateToGroundState("ITP1");
opt.name = "ITP1";
plotdatatopng(opt.name,data,opt);

//////////// VORTICES ////////////////

	// if the given value is true, add vortices to the startgrid
if(opt.runmode.compare(3,1,"1") == 0)
{
int sigma_grid[2];
sigma_grid[0] = opt.grid[1]/8;
sigma_grid[1] = opt.grid[2]/8;
double r = (sigma_grid[0]+sigma_grid[1])/2.0; 

// data = add_central_vortex(data,opt);	
// data = add_circle_vortex(data,opt,r,3);
data = add_circle_vortex(data,opt,r/4.0,6);
data = add_circle_vortex(data,opt,r*2.0/4.0,12);
data = add_circle_vortex(data,opt,r*3.0/4.0,24);

cout << "Vortices added." << endl;
}

////// END VORTICES //////////

//====> Imaginary Time Propagation (ITP)
itprun->formVortices("ITP2");
opt.name = "ITP2";
plotdatatopng(opt.name,data,opt);
// saveDataToHDF5(data,opt);
}

delete itprun;
RTE* rterun = new RTE(data,opt);

if(opt.runmode.compare(0,1,"1") == 0)
{
	// readDataFromHDF5(data,opt);	
	rterun->setOptions(opt);
	rterun->RunSetup();
	printInitVar(opt);
}

//====> Real Time Expansion (RTE)
vector<int> snapshot_times(10);
for(int i = 0;i < 10;i++){
	snapshot_times[i] = (i+1) *opt.n_it_RTE / 10.0;
}

ofstream runparameters;
runparameters.open(("runparameters.txt")/*.c_str()*/, ios::out | ios::trunc);
runparameters << "Parameters of this run:" << endl
				<< "Gridsize in x-direction: " << opt.grid[1] << "\t" << "omega_x = " << opt.omega_x.real() << " dispersion_x = " << opt.dispersion_x.real() << endl
				<< "Gridsize in y-direction: " << opt.grid[2] << "\t" << "omega_y = " << opt.omega_y.real() << " dispersion_y = " << opt.dispersion_y.real() << endl
				<< "K-Length in x-direction: " << opt.klength[0] << endl
				<< "K-Length in y-direction: " << opt.klength[1] << endl
				<< "Expansion factor: " << opt.exp_factor.real() << "\t" << "Number of particles: " << opt.N << "\t" << "Interaction constant g: " << opt.g << endl
				<< "Runmode: " << opt.runmode << endl
				<< "Runtime of the RTE: " << opt.n_it_RTE << " steps." << endl << endl;
runparameters.close();

// run
string runname = "RTE";
Eval* eval = new Eval;
rterun->rteToTime(runname,snapshot_times,eval);

cout << "Run finished." << endl;

delete eval;
delete rterun;
delete data;

// Everything finished here
 
}  // exceptions catcher
catch(std::exception& e) 
{ 
  std::cerr << "Unhandled Exception reached the top of main: " 
            << e.what() << ", application will now exit" << std::endl; 
  return ERROR_UNHANDLED_EXCEPTION; 
}  
return SUCCESS; 	
}




