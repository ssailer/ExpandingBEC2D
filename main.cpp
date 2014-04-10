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

#include <EXP2D_tools.h>
#include <main.h>
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

// print the initial values of the run to the console

if(opt.RTE_only == false)
{
	printInitVar(opt); 
}

// Initialize the needed grid object 
ComplexGrid* data = new ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);	

// if the given value is true, initialize the startgrid with a gaussian distribution

if(opt.RTE_only == false){

if(opt.startgrid[0]==true)
{
	double sigma_real[2];
	sigma_real[0] = opt.min_x/4;
	sigma_real[1] = opt.min_y/4;		
	data = set_grid_to_gaussian(data,opt,sigma_real[0],sigma_real[1]);
}else
{
	data = create_noise_Start_Grid(data,opt);
}

// Initialize the Run Object and load the Grid
EXP2D* run = new EXP2D(data,opt);

// set the datafile identifier name and save the initial grid

opt.name = "INIT";
plotdatatopng(data,opt);

// //====> Imaginary Time Propagation (ITP)
run->itpToTime("ITP1", opt.n_it_ITP1,false);

//////////// VORTICES ////////////////

	// if the given value is true, add vortices to the startgrid
if(opt.startgrid[1]==true)
{
  int sigma_grid[2];
	sigma_grid[0] = opt.grid[1]/8;
	sigma_grid[1] = opt.grid[2]/8;
	double r = (sigma_grid[0]+sigma_grid[1])/2.0; 

  data = add_central_vortex(data,opt);	
  data = add_circle_vortex(data,opt,r/4.0,3);
  data = add_circle_vortex(data,opt,r/2.0,6);
  data = add_circle_vortex(data,opt,r*3.0/4.0,12);
	// data = add_circle_vortex(data,opt,r,2);

cout << "Vortices added." << endl;
opt.name = "VORT";
plotdatatopng(data,opt);
}

////// END VORTICES //////////

//====> Imaginary Time Propagation (ITP)
run->itpToTime("ITP2",opt.n_it_ITP2,false);
opt.name = "ITP2";
plotdatatopng(data,opt);
saveDataToHDF5(data,opt);

//====> Real Time Expansion (RTE)
run->rteToTime("RTE",opt.n_it_RTE,true);
delete run;

// endif opt.RTE_only
}else
{
EXP2D* run = new EXP2D(data,opt);
readDataFromHDF5(data,opt);
run->setOptions(opt);
run->RunSetup();
printInitVar(opt);
//====> Real Time Expansion (RTE)
run->rteToTime("RTE",opt.n_it_RTE,true);
delete run;

}

delete data;

// Everything finished here, cleanup remaining	

cout << "Run finished." << endl;


 
}  // exceptions catcher
catch(std::exception& e) 
{ 
  std::cerr << "Unhandled Exception reached the top of main: " 
            << e.what() << ", application will now exit" << std::endl; 
  return ERROR_UNHANDLED_EXCEPTION; 
}  
return SUCCESS; 	
}



