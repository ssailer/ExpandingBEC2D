#ifndef INIT_BH3_H__
#define INIT_BH3_H__

#include <libconfig.h++>
#include <string>
#include <cstring>

using namespace libconfig;

 inline double vortex(int a, int b, int x, int y) //Vortex with phase [0,2*pi)          
{
        if(atan2(b-y,a-x)<0){ return 2*M_PI+atan2(b-y,a-x); } //atan2 is defined from [-pi,pi) so it needs to be changed to [0,2*pi)
	else{ return atan2(b-y,a-x); }        
}
	
inline double gauss(double & x,double & y){return (exp(-x*x/2-y*y/2));} //A simple Gaussian

inline void printInitVar(Options &opt)
	{
		std::cout 	<< endl << "Initial Values of this run:" << endl
					<< "Used Configfile: \"" << opt.config << "\"   Remember, all inputs override the config!" << endl << endl
					<< "Gridsize in x direction: " << opt.grid[1] << "\t" << "omega_x = " << opt.omega_x << endl
			  		<< "Gridsize in y direction: " << opt.grid[2] << "\t" << "omega_y = " << opt.omega_y << endl
			  		<< "Expansion Factor: " << opt.exp_factor << "\t" << "Number of particles: " << opt.N << "\t" << "Interaction constant g: " << opt.g << endl
			  		<< "Initial Packet (1 = yes, 0 = false): " << opt.startgrid[0] << "\t" << "Vortices added (1 = yes, 0 = false): " << opt.startgrid[1] << endl
					<< "Total time of the ITP-Step: " << opt.n_it_ITP << endl
					<< "Total time of the RTE-Step: " << opt.n_it_RTE << endl << endl;
	}


void init_bh3(int argc, char** argv, Options &opt, vector<double> &snapshot_times)
{
	libconfig::Config cfg;

	string sConfig = opt.config;
	const char* configfile = argv[1]; /// I would really like to use opt.config here.. ffs. WHY NOT? Something with the scope I bet.


	  // Read the file. If there is an error, report it and exit.
	  try
	  {
	    cfg.readFile(configfile);
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
	opt.n_it_ITP             = root["RunOptions"]["n_it_ITP"];				
	opt.n_save_ITP           = root["RunOptions"]["n_save_ITP"]; 			
	opt.times                = root["RunOptions"]["times"]; 	    			
	opt.ITP_step             = root["RunOptions"]["ITP_step"]; 				
	opt.RTE_step             = root["RunOptions"]["RTE_step"]; 				
	opt.cV                   = root["RunOptions"]["cV"]; 		   	   		
	opt.rV                   = root["RunOptions"]["rV"]; 			   		
	opt.Q                    = root["RunOptions"]["Q"];	
	// opt.name 				 = cfg.lookup("RunOptions.name").c_str();

	opt.startgrid[0]         = root["RunOptions"]["startgrid0"];
	opt.startgrid[1]         = root["RunOptions"]["startgrid1"];

	double exp_factor        = root["RunOptions"]["exp_factor"];
	double omega_x_realValue = cfg.lookup("RunOptions.omega_x");	// root["RunOptions"]["omega_x"];
	double omega_y_realValue = cfg.lookup("RunOptions.omega_y");	// root["RunOptions"]["omega_y"];

	opt.exp_factor           = complex<double>(exp_factor,0); //Expansion factor
	opt.omega_x              = complex<double>(omega_x_realValue,0);
	opt.omega_y              = complex<double>(omega_y_realValue,0);	

	}
	catch(const SettingNotFoundException &nfex)
	{
	cerr << endl <<  "Something is wrong here with your config, using default settings,\nexcept the values you have give me directly." << endl << endl;

	opt.N          = 1000;
	opt.exp_factor = complex<double>(1.0,0);
	opt.omega_x = complex<double>(150.0,0);
	opt.omega_y = complex<double>(100.0,0);
	opt.startgrid[0] = true; // gaussian packet on
	opt.startgrid[1] = false; // vortices off
	opt.min_x      = 2;
	opt.min_y      = 2;	
	opt.grid[0]    = 1;
	opt.grid[1]    = 600;
	opt.grid[2]    = 600;
	opt.grid[3]    = 1;
	opt.g          = 1;
	opt.cV         = 2;
	opt.rV         = 2;
	opt.Q          = 1;
	opt.n_it_RTE   = 101;
	opt.n_save_RTE = 50;
	opt.n_it_ITP   = 10000;
	opt.n_save_ITP = 1000;
	opt.times      = 1;
	
	opt.ITP_step   = 0.000001; //Time-step for the ITP (0.000001)
	opt.RTE_step   = 0.00001; //Time-step for the RTE (0.00001)

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
	
}

#endif // INIT_BH3_H__