#ifndef INIT_BH3_H__
#define INIT_BH3_H__

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

void saveDataToHDF5(ComplexGrid* &g, Options &opt)
{ 

  PathOptions options;

  	// useless to me
  options.timestepsize = opt.RTE_step;
  for(int i = 0; i<3; i++){options.klength[i] = 2.0;}
  options.delta_t.resize(1);
  options.delta_t[0] = 1.0;
  	// still useless to me

  options.U = opt.g;
  options.N = opt.N;
  for(int i = 0; i<4;i++){ options.grid[i] = opt.grid[i]; }
  options.g.resize(3);
  options.g[0] = real(opt.omega_x);
  options.g[1] = real(opt.omega_y);
  options.g[2] = real(opt.t_abs);

  double time = opt.n_it_ITP1 + opt.n_it_ITP2;

  Bh3BinaryFile *bf = new Bh3BinaryFile(opt.workingfile, options, Bh3BinaryFile::out);

	vector<ComplexGrid> vectork(1);
	vectork[0] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

	for(int i = 0; i < opt.grid[1];i++){for(int j = 0; j < opt.grid[2]; j++){ vectork.at(0).at(0,i,j,0) = g->at(0,i,j,0) ;}}

	bf->append_snapshot(time, vectork);

  delete bf;
}

void readDataFromHDF5(ComplexGrid* &g,Options &opt)
{
	try{
	PathOptions options;

  	double time = opt.n_it_ITP1 + opt.n_it_ITP2;

	Bh3BinaryFile *bf = new Bh3BinaryFile(opt.workingfile, options, Bh3BinaryFile::in);

	options = bf->get_options();

	opt.g = options.U;
	opt.N = options.N;
	for(int i = 0; i<4; i++){ options.grid[i] = opt.grid[i]; }
	options.g.resize(3);
	opt.omega_x = complex<double>(options.g[0],0.0);
	opt.omega_y = complex<double>(options.g[1],0.0);
	opt.t_abs   = complex<double>(options.g[2],0.0);

	vector<ComplexGrid> vectork(1);
	vectork[0] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

	bf->get_snapshot(time, vectork,0);

	for(int i = 0; i < opt.grid[1];i++){for(int j = 0; j < opt.grid[2]; j++){ g->at(0,i,j,0) = vectork.at(0).at(0,i,j,0) ;}}
	
	delete bf;
	}

	catch(std::exception& e) 
	{ 
  		std::cerr << "Reading from HDF5 File failed, whaaaat?: " 
        		  << e.what() << ", application will now exit" << std::endl; 
	} 

}

void printInitVar(Options &opt)
{
	std::cout.setf(std::ios::boolalpha);
	std::cout 	<< "Used configfile: \"" << opt.config << "\"" << endl
				<< "Gridsize in x-direction: " << opt.grid[1] << "\t" << "omega_x = " << opt.omega_x.real() << endl
				<< "Gridsize in y-direction: " << opt.grid[2] << "\t" << "omega_y = " << opt.omega_y.real() << endl
				<< "Expansion factor: " << opt.exp_factor.real() << "\t" << "Number of particles: " << opt.N << "\t" << "Interaction constant g: " << opt.g << endl
				<< "Initial gausspacket: " << opt.startgrid[0] << "\t" << "Vortices will be added: " << opt.startgrid[1] << endl
				<< "RTE potential on: " << opt.startgrid[2] << endl
				<< "Runtime of the ITP1: " << opt.n_it_ITP1 << " steps." << endl
				<< "Runtime of the ITP2: " << opt.n_it_ITP2 << " steps." << endl
				<< "Runtime of the RTE: " << opt.n_it_RTE << " steps." << endl << endl;
}

void set_workingdirectory(Options &opt)
{
	// cout << "Workingdirectory: " << "\"" << opt.workingdirectory << "\"" << endl;
	struct stat wd_stat;
	if(stat(opt.workingdirectory.c_str(),&wd_stat) == 0){
		if(chdir(opt.workingdirectory.c_str()) == 0){
			cout << "Changed to existing directory: " << "\"" << opt.workingdirectory << "\"" << endl;
		}
	}else
	{
		char command[256];
		sprintf(command,"mkdir %s",opt.workingdirectory.c_str());
		if(system(command) == 0){
			cout << "Created directory: " << "\"" << opt.workingdirectory << "\"";
		}
		if(chdir(opt.workingdirectory.c_str()) == 0){
			cout << " and changed into it.";
		}
		cout << endl;
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
	opt.grid[0]              = root["RunOptions"]["grid0"];				
	opt.grid[1]              = root["RunOptions"]["grid1"];				
	opt.grid[2]              = root["RunOptions"]["grid2"];	   			
	opt.grid[3]              = root["RunOptions"]["grid3"];				
	opt.g                    = root["RunOptions"]["g"]; 						
	opt.n_it_RTE             = root["RunOptions"]["n_it_RTE"]; 				
	opt.n_save_RTE           = root["RunOptions"]["n_save_RTE"]; 			
	opt.n_it_ITP1            = root["RunOptions"]["n_it_ITP1"];	
	opt.n_it_ITP2            = root["RunOptions"]["n_it_ITP2"];				
	opt.n_save_ITP           = root["RunOptions"]["n_save_ITP"];   			
	opt.ITP_step             = root["RunOptions"]["ITP_step"]; 				
	opt.RTE_step             = root["RunOptions"]["RTE_step"];
	opt.Q                    = root["RunOptions"]["Q"];
	opt.RTE_only			 = root["RunOptions"]["RTE_only"];
	// opt.workingfile			 = root["RunOptions"]["workingfile"]
	cfg.lookupValue("RunOptions.workingfile",opt.workingfile);
	// opt.name

	opt.startgrid[0]         = root["RunOptions"]["gaussian"];
	opt.startgrid[1]         = root["RunOptions"]["vortices"];
	opt.startgrid[2]		 = root["RunOptions"]["potential"];

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
	opt.scale_factor = complex<double>(0,0); //Scale factor
	opt.t_abs = complex<double>(0,0); //Absolute time 
	opt.name       = "run";


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

#endif // INIT_BH3_H__