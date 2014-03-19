#ifndef INIT_BH3_H__
#define INIT_BH3_H__

#include <libconfig.h++>
#include <string>
#include <cstring>

using namespace libconfig;

void printInitVar(Options &opt)
	{
		std::cout.setf(std::ios::boolalpha);
		std::cout 	<< "Initial Values of this run:" << endl
					<< "Used Configfile: \"" << opt.config << endl
					<< "Gridsize in x direction: " << opt.grid[1] << "\t" << "omega_x = " << opt.omega_x.real() << endl
			  		<< "Gridsize in y direction: " << opt.grid[2] << "\t" << "omega_y = " << opt.omega_y.real() << endl
			  		<< "Expansion Factor: " << opt.exp_factor.real() << "\t" << "Number of particles: " << opt.N << "\t" << "Interaction constant g: " << opt.g << endl
			  		<< "Initial Packet: " << opt.startgrid[0] << "\t" << "Vortices added: " << opt.startgrid[1] << endl
			  		<< "RTE is having potential: " << opt.startgrid[2] << endl
					<< "Runtime of the ITP1-Step: " << opt.n_it_ITP1 << endl
					<< "Runtime of the ITP2-Step: " << opt.n_it_ITP2 << endl
					<< "Runtime of the RTE-Step: " << opt.n_it_RTE << endl << endl;
	}


int init_bh3(int argc, char** argv, Options &opt, vector<double> &snapshot_times)
{


	libconfig::Config cfg;

	string sConfig = opt.config;
	// const char* configfile = argv[1]; /// I would really like to use opt.config here.. ffs. WHY NOT? Something with the scope I bet.


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
	opt.times                = root["RunOptions"]["times"]; 	    			
	opt.ITP_step             = root["RunOptions"]["ITP_step"]; 				
	opt.RTE_step             = root["RunOptions"]["RTE_step"];
	opt.Q                    = root["RunOptions"]["Q"];	
	// opt.name 				 = cfg.lookup("RunOptions.name").c_str();

	opt.startgrid[0]         = root["RunOptions"]["startgrid0"];
	opt.startgrid[1]         = root["RunOptions"]["startgrid1"];
	opt.startgrid[2]		 = root["RunOptions"]["startgrid2"];

	double exp_factor        = root["RunOptions"]["exp_factor"];
	double omega_x_realValue = root["RunOptions"]["omega_x"];  // cfg.lookup("RunOptions.omega_x");
	double omega_y_realValue = root["RunOptions"]["omega_y"];  // cfg.lookup("RunOptions.omega_y");

	opt.exp_factor           = complex<double>(exp_factor,0); //Expansion factor
	opt.omega_x              = complex<double>(omega_x_realValue,0);
	opt.omega_y              = complex<double>(omega_y_realValue,0);	

	}
	catch(const SettingNotFoundException &nfex)
	{
	cerr << endl <<  "Something is wrong here with your config." << endl << endl;

	return 1;

	}

	char* command;
	sprintf(command,"mkdir %s",opt.workingdirectory.c_str());
    system(command);
    chdir(opt.workingdirectory.c_str());



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
	

	
	snapshot_times.resize(5);
	snapshot_times[0] =5000;
	snapshot_times[1] =10000;
	snapshot_times[2] =15000;
	snapshot_times[3] =20000;
	snapshot_times[4] =25000;
//         snapshot_times[6] =17000;
          // snapshot_times[7] =20000;
          // snapshot_times[8] =30000;
          // snapshot_times[9] =50000;
          // snapshot_times[10]=70000;
          // snapshot_times[11]=100000;
          // snapshot_times[12]=150000;
          // snapshot_times[13]=180000;
          // snapshot_times[14]=220000;
          // snapshot_times[15]=250000;
          // snapshot_times[16]=300000;
          // snapshot_times[17]=450000;
          // snapshot_times[18]=500000;
          // snapshot_times[19]=550000;
          // snapshot_times[20]=600000;
          // snapshot_times[21]=600000;
          // snapshot_times[22]=600000;
          // snapshot_times[23]=650000;
          // snapshot_times[24]=700000;
          // snapshot_times[25]=750000;
          // snapshot_times[26]=800000;
          // snapshot_times[27]=840000;
          // snapshot_times[28]=870000;
          // snapshot_times[29]=900000;
          // snapshot_times[30]=930000;
          // snapshot_times[31]=960000;
          // snapshot_times[32]=1000000;
          // snapshot_times[33]=1100000;
          // snapshot_times[34]=1200000;
          // snapshot_times[35]=1300000;
          // snapshot_times[36]=1400000;
          // snapshot_times[37]=1500000;
          // snapshot_times[38]=1600000;
          // snapshot_times[39]=1700000;
          // snapshot_times[40]=1800000;
          // snapshot_times[41]=1900000;
      
		

        /*		for(int i = 0; i < snapshot_times.size(); i++)///////////////////////was passiert hier/????????////====warum die snapshot zeit veraendern???========?????
                {
                int time_res = 1000;
                snapshot_times[i] = (int) (500.*pow(10,((double)i/10.)));
                    //cout<<snapshot_times[i]<<endl;
                        //snapshot_times[i] = i*time_res + 25000;
                        }*/

     return 0;
	
}

#endif // INIT_BH3_H__