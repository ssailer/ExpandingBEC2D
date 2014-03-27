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

#include <exp_RK4_tools.h>
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
		
// Initialize the needed grid object and run object

ComplexGrid* startgrid = new ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
RK4* run = new RK4(startgrid,opt);	

// if the given value is true, initialize the startgrid with a gaussian distribution

if(opt.RTE_only == false){

if(opt.startgrid[0]==true)
{
	double sigma_real[2];
	sigma_real[0] = opt.min_x/4;
	sigma_real[1] = opt.min_y/4;		
	run->pPsi = set_grid_to_gaussian(run->pPsi,opt,run->x_axis,run->y_axis,sigma_real[0],sigma_real[1]);
}else
{
	run->pPsi = create_noise_Start_Grid(run->pPsi,opt);
}

// set the datafile identifier name and save the initial grid

opt.name = "INIT";
plotdatatopng(run->pPsi,opt);

// //====> Imaginary Time Propagation (ITP)
opt.name = "ITP1";
opt.n_it_ITP = opt.n_it_ITP1;
run->itpToTime(opt,false);

//////////// VORTICES ////////////////

	// if the given value is true, add vortices to the startgrid
if(opt.startgrid[1]==true)
{
  int sigma_grid[2];
	sigma_grid[0] = opt.grid[1]/6;
	sigma_grid[1] = opt.grid[2]/6;
	double r = (sigma_grid[0]+sigma_grid[1])/2.0; 

  run->pPsi = add_central_vortex(run->pPsi,opt);	
  run->pPsi = add_circle_vortex(run->pPsi,opt,r/4.0,6);
  run->pPsi = add_circle_vortex(run->pPsi,opt,r/2.0,12);
  run->pPsi = add_circle_vortex(run->pPsi,opt,r/1.5,24);
	// run->pPsi = add_circle_vortex(run->pPsi,opt,r,2);

  cout << "Vortices added." << endl;
  opt.name = "VORT";

	plotdatatopng(run->pPsi,opt);
}

////// END VORTICES //////////

//====> Imaginary Time Propagation (ITP)
opt.name = "ITP2";
opt.n_it_ITP = opt.n_it_ITP2;

run->itpToTime(opt,false);

plotdatatopng(run->pPsi,opt);
saveDataToHDF5(run->pPsi,opt);

// endif opt.RTE_only
}else
{
readDataFromHDF5(run->pPsi,opt);
printInitVar(opt);
}

//====> Real Time Expansion (RTE)
opt.name = "RTE";
	
run->rteToTime(opt,true);

// Everything finished here, cleanup remaining	

cout << "Run finished." << endl;

delete startgrid;
delete run;
 
}  // exceptions catcher
catch(std::exception& e) 
{ 
  std::cerr << "Unhandled Exception reached the top of main: " 
            << e.what() << ", application will now exit" << std::endl; 
  return ERROR_UNHANDLED_EXCEPTION; 
}  
return SUCCESS; 	
}



