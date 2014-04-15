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

#include <EXP2D_tools.h>
#include <main.h>
#include <plot_with_mgl.h>
// #include <typeinfo>

#define SUCCESS 0; 
#define ERROR_IN_COMMAND_LINE 1;
#define ERROR_IN_CONFIG_FILE 2;
#define ERROR_UNHANDLED_EXCEPTION 3;

using namespace std;

void plot(const string &dirname, const PathOptions& opt, vector<double> &snapshot_times, AverageClass<Bh3Evaluation::Averages> *av)
{  
	stringstream d;
	d << dirname << "/" << "spectrum" << "/";
	string dir = d.str();
	system((string("mkdir -p ") + dir).c_str());
	
    ofstream plotfile;
	
	plotfile.open((dir + string("radial_avgs.dat")).c_str(), ios::out | ios::trunc);
	
	for (int i = 0; i < snapshot_times.size(); i++)
	{
        Bh3Evaluation::Averages means = av[i].av(); 
        
        for (int r = 0; r < means.number.size(); r++)             
		{
            plotfile << r <<"\t"<< means.k(r) <<"\t" << means.number(r) <<"\t";
            plotfile << means.ikinetick(r) << "\t" << means.ckinetick(r) << "\t";
            plotfile << means.kinetick(r) << "\t" << means.pressure(r) << "\t";
            plotfile << means.ikinetick_wo_phase(r) <<  "\t" ;
            plotfile << means.ckinetick_wo_phase(r) << "\t";
            plotfile << means.kinetick_wo_phase(r) << "\t" ;
            plotfile << means.pressure_wo_phase(r) <<"\t";
            plotfile << endl;
		}
		plotfile << endl << endl;
	}
	
	plotfile.close();
}


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
// AverageClass<ComplexGrid>* av = new AverageClass<ComplexGrid>();

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
  // data = add_circle_vortex(data,opt,r,3);
  data = add_circle_vortex(data,opt,r/4.0,6);
  data = add_circle_vortex(data,opt,r*2.0/4.0,12);
	data = add_circle_vortex(data,opt,r*3.0/4.0,24);

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

//====> Noising the Grid before RTE
noiseTest(opt,data);
opt.name = "Noise after ITP";
plotdatatopng(data,opt);


//====> Real Time Expansion (RTE)
vector<double> snapshot_times(10);
for(int i = 0;i<10;i++){
	snapshot_times[i] = i *opt.n_it_RTE / 10.0;
}

PathOptions pathopt;
	pathopt.timestepsize = opt.RTE_step;
	pathopt.delta_t.resize(0);             
        //pathopt.delta_t[0]=0.2;
        //pathopt.delta_t[1]=0.4;

	pathopt.N =opt.N;  //normed for 512*512 N=64*50000
	
    pathopt.grid[0] = opt.grid[0];
    pathopt.grid[1] = opt.grid[1];
	pathopt.grid[2] = opt.grid[2];
	pathopt.grid[3] = opt.grid[3];
	pathopt.U = opt.g;
	
    pathopt.g.resize(0);
    // pathopt.g[0] = 1./4.;
    	
	pathopt.klength[0] = 2.0;
	pathopt.klength[1] = 2.0;
	pathopt.klength[2] = 2.0;


double start = omp_get_wtime();
	AverageClass<Bh3Evaluation::Averages> *av =  new AverageClass<Bh3Evaluation::Averages> [snapshot_times.size()];

for(int j = 0; j < snapshot_times.size(); j++){


run->rteToTime("RTE",snapshot_times[j],false);
run->cli_plot("RTE",j*10,snapshot_times.size()*10,start,true);

Bh3Evaluation ev(pathopt);
								
                ev.setTime(snapshot_times[j]);

                vector<ComplexGrid> eval_data(1);
                eval_data[0] = *data;

                ev.setTime(snapshot_times[j]);
                ev.setData(eval_data,Bh3Evaluation::RSpace);
                                                             
                ev.calc_radial_averages(); 
                
                av->average(ev.get_averageable_results());

                	stringstream dstr;
	dstr << "exp2d_" + std::to_string(j);
	string dirname = dstr.str();
    initialize_binary_dir(dirname, pathopt);
	// mkdir((dirname + "/temp").c_str(), 0755);
    plot(dirname, pathopt, snapshot_times, av);



}
delete [] av;

delete run;

}else 
{
// run RTE only loading data from rundata

EXP2D* run = new EXP2D(data,opt);
readDataFromHDF5(data,opt);
run->setOptions(opt);
run->RunSetup();
printInitVar(opt);

//====> Noising the Grid before RTE
noiseTest(opt,data);
opt.name = "Noise after ITP";
plotdatatopng(data,opt);



//====> Real Time Expansion (RTE)
//====> Real Time Expansion (RTE)
vector<double> snapshot_times(10);
for(int i = 0;i<10;i++){
	snapshot_times[i] = (i+1) *opt.n_it_RTE / 10.0;
}

PathOptions pathopt;
	pathopt.timestepsize = opt.RTE_step;
	pathopt.delta_t.resize(0);             
        //pathopt.delta_t[0]=0.2;
        //pathopt.delta_t[1]=0.4;

	pathopt.N =opt.N;  //normed for 512*512 N=64*50000
	
    pathopt.grid[0] = opt.grid[0];
    pathopt.grid[1] = opt.grid[1];
	pathopt.grid[2] = opt.grid[2];
	pathopt.grid[3] = opt.grid[3];
	pathopt.U = opt.g;
	
    pathopt.g.resize(0);
    // pathopt.g[0] = 1./4.;
    	
	pathopt.klength[0] = 2.0;
	pathopt.klength[1] = 2.0;
	pathopt.klength[2] = 2.0;


double start = omp_get_wtime();
	AverageClass<Bh3Evaluation::Averages> *av =  new AverageClass<Bh3Evaluation::Averages> [snapshot_times.size()];;

for(int j = 0; j < snapshot_times.size(); j++){


run->rteToTime("RTE",snapshot_times[j],false);
run->cli_plot("RTE",j+1,100,start,true);

Bh3Evaluation ev(pathopt);
								
                ev.setTime(snapshot_times[j]);

                vector<ComplexGrid> eval_data(1);
                eval_data[0] = *data;

                ev.setTime(snapshot_times[j]);
                ev.setData(eval_data,Bh3Evaluation::RSpace);
                                                             
                ev.calc_radial_averages(); 
                
                av[j].average(ev.get_averageable_results());

                	stringstream dstr;
	dstr << "exp2d_" + std::to_string(j);
	string dirname = dstr.str();
    initialize_binary_dir(dirname, pathopt);
	mkdir((dirname + "/temp").c_str(), 0755);
    plot(dirname, pathopt, snapshot_times, av);
    delete [] av;
delete run;
}

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




