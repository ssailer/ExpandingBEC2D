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
#include <complexgrid.h>
#include <exp_RK4_tools.h>
#include <bh3defaultgrid.h>
#include <omp.h>
#include <main.h>
#include <plot_with_mgl.h>
// #include <typeinfo>
// #include <vortexcoordinates.h>


using namespace std;


namespace // namespace for program options
{ 
  const size_t ERROR_IN_COMMAND_LINE = 1; 
  const size_t SUCCESS = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
 
}

inline void savedatahdf5(double time,Bh3BinaryFile* &bf, ComplexGrid* &g, Options &opt)
	{
	vector<ComplexGrid> vectork(1);
	vectork[0] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

	for(int i = 0; i < opt.grid[1];i++) for(int j = 0; j < opt.grid[2]; j++)
	vectork.at(0).at(0,i,j,0) = g->at(0,i,j,0);	

	bf->append_snapshot(time, vectork);
	}

//>>>>>main program<<<<< 

int main( int argc, char** argv) 
{	
	Options opt;
	vector<double> snapshot_times;




	

	// Beginning of the options block
try 
{ 
    /** Define and parse the program options 
     */ 
    namespace po = boost::program_options; 
    po::options_description desc("Options"); 
    desc.add_options() 
      ("help,h", "Print help messages.") 
      ("config,c",po::value<string>(&opt.config), "Name of the configfile")
      ("directory,d",po::value<string>(&opt.workingdirectory), "Name of the directory this run saves its data");
      // ("xgrid,x",po::value<int>(&opt.grid[1]),"Gridsize in x direction.")
      // ("ygrid,y",po::value<int>(&opt.grid[2]),"Gridsize in y direction.")
      // ("gauss",po::value<bool>(&opt.startgrid[0]),"Initial Grid has gaussian form.")
      // ("vortices",po::value<bool>(&opt.startgrid[1]),"Add Vortices to the grid.")
      // ("itp",po::value<int>(&opt.n_it_ITP),"Total runtime of the ITP-Step.")
      // ("rte",po::value<int>(&opt.n_it_RTE),"Total runtime of the RTE-Step.")
      // ("number,N",po::value<double>(&opt.N),"Number of particles.")
      // ("expansion,e",po::value<complex<double> >(&opt.exp_factor),"Expansion Factor")
      // ("interaction,g",po::value<double> (&opt.g),"Interaction Constant");

	po::positional_options_description positionalOptions; 
	positionalOptions.add("config", 1);
	positionalOptions.add("directory",1);

    po::variables_map vm; 
    try 
    { 

		po::store(po::command_line_parser(argc, argv).options(desc)
					.positional(positionalOptions).run(), 
          			vm); 
 
      /** --help option 
       */ 
      if ( vm.count("help")  ) 
      { 
        std::cout << "This is a program to simulate the expansion of a BEC with Vortices in 2D" << endl
        		  << "after a harmonic trap has been turned off. The expansion is simulated by" << endl
        		  << "an expanding coordinate system. The implemented algorithm to solve the GPE" << endl
        		  << "is a 4-th order Runge-Kutta Integration." << endl << endl
                  << desc << endl; 
        return SUCCESS; 
      } 
 
      po::notify(vm); // throws on error, so do after help in case 
                      // there are any problems 
    } 
    catch(po::error& e) 
    { 
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
      std::cerr << desc << std::endl; 
      return ERROR_IN_COMMAND_LINE; 
    } 

// Beginning of the main program block

    // Initialize all option variables

    if(init_bh3(argc, argv, opt, snapshot_times) == 1){
		cout << endl << "Could not initialize the values, abort." << endl;
		return 0;
	}else{ 
		cout << endl << "Run parameters initialized. Ready to start the computation." << endl << endl;
	}

    // print the initial values of the run to the console

    printInitVar(opt); 
		
	// Initialize the needed grid object and run object

	ComplexGrid* startgrid = new ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	RK4* run = new RK4(startgrid,opt);
	Bh3BinaryFile *bf = new Bh3BinaryFile("myrun_01", opt, Bh3BinaryFile::out);

	// if the given value is true, initialize the startgrid with a gaussian distribution

	for(int i=0;i<opt.grid[1];i++)for(int j=0;j<opt.grid[2];j++)
		run->pPsi->at(0,i,j,0)=complex<double>(1.0,0.0);


	// if(opt.startgrid[0]==true){
	// double sigma_real[2];
	// sigma_real[0] = opt.min_x/4;
	// sigma_real[1] = opt.min_y/4;		
	// run->pPsi = set_grid_to_gaussian(run->pPsi,opt,run->x_axis,run->y_axis,sigma_real[0],sigma_real[1]);
	// }else{
		
	// 	run->pPsi = create_noise_Start_Grid(run->pPsi,opt);
	// }


	// set the datafile identifier name and save the initial grid

    opt.name = "INIT";
	plotdatatopng(run->pPsi,opt);
	savedatahdf5(1.,bf,run->pPsi,opt);

	//====> Imaginary Time Propagation (ITP)
	opt.name = "ITP1";
	opt.n_it_ITP = opt.n_it_ITP1;
	run->itpToTime(opt,false);

//////////// VORTICES ////////////////

	// if the given value is true, add vortices to the startgrid
	if(opt.startgrid[1]==true)
    {
    int sigma_grid[2];
	sigma_grid[0] = opt.grid[1]/4;
	sigma_grid[1] = opt.grid[2]/4;

    run->pPsi = add_vortex_to_grid(run->pPsi,opt,sigma_grid);
   	cout << "Vortices added." << endl;
   	opt.name = "VORT";

	plotdatatopng(run->pPsi,opt);
   	savedatahdf5(2.,bf,run->pPsi,opt);
   	}
////// END VORTICES //////////

   	//====> Imaginary Time Propagation (ITP)
    opt.name = "ITP2";
	opt.n_it_ITP = opt.n_it_ITP2;

	run->itpToTime(opt,true);

	plotdatatopng(run->pPsi,opt);
	savedatahdf5(3.,bf,run->pPsi,opt);


	//====> Real Time Expansion (RTE)
	opt.name = "RTE";
	
	run->rteToTime(opt,true);

	plotdatatopng(run->pPsi,opt);
	savedatahdf5(4.,bf,run->pPsi,opt);

	// Everything finished here, plots and cleanup remaining	

	cout << "Run finished." << endl;

    delete startgrid;
	delete run;
	delete bf;
 
 
  } 
  // options menu exception catcher
  catch(std::exception& e) 
  { 
    std::cerr << "Unhandled Exception reached the top of main: " 
              << e.what() << ", application will now exit" << std::endl; 
    return ERROR_UNHANDLED_EXCEPTION; 
 
  } 
 
  return SUCCESS; 
	
}



