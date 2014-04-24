/**************************************************************************
Title: Simulating the Expansion of Turbulent Bose-Einstein Condensates (2D) 
Author: Simon Sailer (This work is based on the work of Bartholomew Andrews who made this as his master thesis.)
Last Update: 22/07/13
**************************************************************************/

#include <boost/program_options.hpp>
#include <iostream>
#include <unistd.h>
#include <cstdlib>
#include <cstring>
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

// Initialize all option variables

read_config(argc,argv,opt);

set_workingdirectory(opt);

// Initialize the needed grid object 
ComplexGrid* data = new ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
// Initialize the Run Object and load the Grid
// EXP2D* run = new EXP2D(data,opt);
ITP* itprun = new ITP(data,opt);
RTE* rterun = new RTE(data,opt);


// if the given value is true, initialize the startgrid with a gaussian distribution, else load from rundata file

if(opt.startgrid[0]==true)
{
	printInitVar(opt); 

	double sigma_real[2];
	sigma_real[0] = opt.min_x/4;
	sigma_real[1] = opt.min_y/4;		
	data = set_grid_to_gaussian(data,opt,sigma_real[0],sigma_real[1]);
}else
{
	readDataFromHDF5(data,opt);	
	rterun->setOptions(opt);
	rterun->RunSetup();
	printInitVar(opt);
}



if(opt.RTE_only == false){

// set the datafile identifier name and save the initial grid

opt.name = "INIT";
plotdatatopng(data,opt);

// //====> Imaginary Time Propagation (ITP)
itprun->itpToTime("ITP1", true);
opt.name = "ITP1";
plotdatatopng(data,opt);

//////////// VORTICES ////////////////

	// if the given value is true, add vortices to the startgrid
if(opt.startgrid[1]==true)
{
int sigma_grid[2];
sigma_grid[0] = opt.grid[1]/8;
sigma_grid[1] = opt.grid[2]/8;
double r = (sigma_grid[0]+sigma_grid[1])/2.0; 

data = add_central_vortex(data,opt);	
// data = add_circle_vortex(data,opt,r,3);
data = add_circle_vortex(data,opt,r/4.0,6);
data = add_circle_vortex(data,opt,r*2.0/4.0,12);
data = add_circle_vortex(data,opt,r*3.0/4.0,24);

cout << "Vortices added." << endl;
}

////// END VORTICES //////////

//====> Imaginary Time Propagation (ITP)
itprun->itpToTime("ITP2",true);
opt.name = "ITP2";
plotdatatopng(data,opt);
saveDataToHDF5(data,opt);
}

//setting expansion without extending coordinates
// opt.omega_x = 0.0;
// opt.omega_y = 0.0;
rterun->setOptions(opt);
rterun->RunSetup();

//====> Real Time Expansion (RTE)
vector<int> snapshot_times(10);
for(int i = 0;i < 10;i++){
	snapshot_times[i] = (i+1) *opt.n_it_RTE / 10.0;
}

PathOptions pathopt;
	pathopt.timestepsize = opt.RTE_step;
	pathopt.delta_t.resize(0);
	pathopt.g.resize(0);
	pathopt.N =opt.N;
    pathopt.grid[0] = opt.grid[0];
    pathopt.grid[1] = opt.grid[1];
	pathopt.grid[2] = opt.grid[2];
	pathopt.grid[3] = opt.grid[3];
	pathopt.U = opt.g;
	pathopt.klength[0] = 2.0;
	pathopt.klength[1] = 2.0;
	pathopt.klength[2] = 2.0;

ofstream runparameters;
runparameters.open(("runparameters.txt")/*.c_str()*/, ios::out | ios::trunc);
runparameters << "Parameters of this run:" << endl
				<< "Gridsize in x-direction: " << opt.grid[1] << "\t" << "omega_x = " << opt.omega_x.real() << endl
				<< "Gridsize in y-direction: " << opt.grid[2] << "\t" << "omega_y = " << opt.omega_y.real() << endl
				<< "K-Length in x-direction: " << pathopt.klength[0] << endl
				<< "K-Length in y-direction: " << pathopt.klength[1] << endl
				<< "Expansion factor: " << opt.exp_factor.real() << "\t" << "Number of particles: " << opt.N << "\t" << "Interaction constant g: " << opt.g << endl
				<< "Initial gausspacket: " << opt.startgrid[0] << "\t" << "Vortices were be added: " << opt.startgrid[1] << endl
				<< "RTE potential on: " << opt.startgrid[2] << endl
				<< "Runtime of the ITP1: " << opt.n_it_ITP1 << " steps." << endl
				<< "Runtime of the ITP2: " << opt.n_it_ITP2 << " steps." << endl
				<< "Runtime of the RTE: " << opt.n_it_RTE << " steps." << endl << endl;
runparameters.close();

// run
string runname = "RTE";
Averages* eval = new Averages;
rterun->rteToTime(runname,snapshot_times,eval);

// opt.workingfile = "rundata_afterRTE";
// saveDataToHDF5(data,opt);

cout << "Run finished." << endl;

// 	// evalution
// start = omp_get_wtime();
// for(int j = 0; j < snapshot_times.size(); j++){
// run->cli_plot("Evaluating",(j+1)*10,100,start,false);
// evaluate(data_storage[j],pathopt,j+1);
// }

// cout << "Evaluating finished." << endl;
delete eval;

delete itprun;
delete rterun;

delete data;

// Everything finished here, cleanup remaining	




 
}  // exceptions catcher
catch(std::exception& e) 
{ 
  std::cerr << "Unhandled Exception reached the top of main: " 
            << e.what() << ", application will now exit" << std::endl; 
  return ERROR_UNHANDLED_EXCEPTION; 
}  
return SUCCESS; 	
}




