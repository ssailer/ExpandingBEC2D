#ifndef MAIN_H__
#define MAIN_H__

#include <libconfig.h++>
#include <string>
#include <cstring>
#include <unistd.h>
#include <stdio.h>

#define SUCCESS 0; 
#define ERROR_IN_COMMAND_LINE 1;
#define ERROR_IN_CONFIG_FILE 2;
#define ERROR_UNHANDLED_EXCEPTION 3;

using namespace libconfig;

void printInitVar(Options &opt)
{
	std::cout.setf(std::ios::boolalpha);
	std::cout 	<< "Used configfile: \"" << opt.config << "\"" << endl
				<< "Gridsize in x-direction: " << opt.grid[1] << "\t" << "omega_x = " << opt.omega_x.real() << " dispersion_x = " << opt.dispersion_x.real() << endl
				<< "Gridsize in y-direction: " << opt.grid[2] << "\t" << "omega_y = " << opt.omega_y.real() << " dispersion_y = " << opt.dispersion_y.real() << endl
				<< "Expansion factor: " << opt.exp_factor.real() << "\t" << "Number of particles: " << opt.N << "\t" << "Interaction constant g: " << opt.g << endl
				<< "Reading from Datafile: " << opt.runmode[0] << "\t" << "Vortices will be added: " << opt.runmode[3] << endl
				<< "RTE potential on: " << opt.runmode[2] << endl
				<< "Runmode: " << opt.runmode << endl
				<< "Runtime of the RTE: " << opt.n_it_RTE << " steps." << endl << endl;
}

void set_workingdirectory(Options &opt)
{
	// cout << "Workingdirectory: " << "\"" << opt.workingdirectory << "\"" << endl;
	struct stat wd_stat;
	if(stat(opt.workingdirectory.c_str(),&wd_stat) == 0){
		if(chdir(opt.workingdirectory.c_str()) == 0){
			cout << "Using existing directory: " << "\"" << opt.workingdirectory << "\"." << endl;
			cout << "Check \"run.log\" for output of this run." << endl;
		}
	}else
	{
		char command[256];
		sprintf(command,"mkdir %s",opt.workingdirectory.c_str());
		if(system(command) == 0){
			cout << "Creating directory: " << "\"" << opt.workingdirectory << "\"";
		}
		if(chdir(opt.workingdirectory.c_str()) == 0){
			cout << "Switchting to "<< "\"" << opt.workingdirectory << "\"";
		}
		cout << endl;
		cout << "Check \"run.log\" for output of this run." << endl;
	}
}


int read_cli_options(int argc, char** argv, Options &opt)
{
	// Beginning of the options block

    // Define and parse the program options 

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
      // ("expansion,e",po::value<double>(&exp_factor),"Expansion Factor")
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

    return SUCCESS;

}


int read_config(int argc, char** argv, Options &opt)
{
	libconfig::Config cfg;

	string sConfig = opt.config;

	  // Read the file. If there is an error, report it and exit.
	  try
	  {
	    cfg.readFile(sConfig.c_str());
	  }
	  catch(const libconfig::FileIOException &fioex)
	  {
	    std::cerr << "I/O error while reading file." << std::endl;
	  }
	  catch(const libconfig::ParseException &pex)
	  {
	    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
	              << " - " << pex.getError() << std::endl;
	  }


	const libconfig::Setting & root = cfg.getRoot();

	try
	{

	opt.N                    = root["RunOptions"]["N"];
	opt.min_x                = root["RunOptions"]["min_x"]; 					
	opt.min_y                = root["RunOptions"]["min_y"];
	opt.klength[0] 			 = root["RunOptions"]["klength0"];
	opt.klength[1] 			 = root["RunOptions"]["klength1"];
	opt.klength[2] 			 = root["RunOptions"]["klength2"];
	opt.grid[0]              = root["RunOptions"]["grid0"];				
	opt.grid[1]              = root["RunOptions"]["grid1"];				
	opt.grid[2]              = root["RunOptions"]["grid2"];	   			
	opt.grid[3]              = root["RunOptions"]["grid3"];				
	opt.g                    = root["RunOptions"]["g"]; 						
	opt.n_it_RTE             = root["RunOptions"]["n_it_RTE"]; 				
	opt.ITP_step             = root["RunOptions"]["ITP_step"]; 				
	opt.RTE_step             = root["RunOptions"]["RTE_step"];
	opt.samplesize			 = root["RunOptions"]["samplesize"];
	// opt.workingfile			 = root["RunOptions"]["workingfile"]
	cfg.lookupValue("RunOptions.workingfile",opt.workingfile);
	cfg.lookupValue("RunOptions.runmode",opt.runmode);
	// opt.name

	double exp_factor        = root["RunOptions"]["exp_factor"];
	double omega_x_realValue = root["RunOptions"]["omega_x"];  // cfg.lookup("RunOptions.omega_x");
	double omega_y_realValue = root["RunOptions"]["omega_y"];  // cfg.lookup("RunOptions.omega_y");

	double dispersion_x_realValue = root["RunOptions"]["dispersion_x"]; 
	double dispersion_y_realValue = root["RunOptions"]["dispersion_y"]; 

	opt.exp_factor           = complex<double>(exp_factor,0); //Expansion factor
	opt.omega_x              = complex<double>(omega_x_realValue,0);
	opt.omega_y              = complex<double>(omega_y_realValue,0);
	opt.dispersion_x		 = complex<double>(dispersion_x_realValue,0);
	opt.dispersion_y 		 = complex<double>(dispersion_y_realValue,0);



	}
	catch(const SettingNotFoundException &nfex)
	{
	cerr << endl <<  "Something is wrong here with your config." << endl << endl;

	return ERROR_IN_CONFIG_FILE;

	}

	// runspecific Values, just initilized here
	opt.scale_factor = 0; //Scale factor
	opt.t_abs = complex<double>(0,0); //Absolute time 
	opt.name       = "run";
	opt.stateInformation.resize(2);
	opt.stateInformation[0] = 1;
	opt.stateInformation[1] = 1;
	opt.vortexnumber = 0;

    // Set Parameters manually (default values)
	//opt.timestepsize = 0.2;
	//opt.delta_t.resize(0);   
    //opt.delta_t[0]=0.2;
    //opt.delta_t[1]=0.4;
	// 	opt.klength[0] = 2.0;
	// 	opt.klength[1] = 2.0;
	// 	opt.klength[2] = 2.0;

     return SUCCESS;
	
}





#endif // MAIN_H__