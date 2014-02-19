/**************************************************************************
Title: Simulating the Expansion of Turbulent Bose-Einstein Condensates (2D) 
Author: Bartholomew Andrews
Last Update: 22/07/13
Website: www.bartholomewandrews.com
**************************************************************************/
#include <boost/program_options.hpp>
#include <iostream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <complex>
#include <complexgrid.h>
#include <exp_RK4_tools.h>
#include <bh3defaultgrid.h>
#include <omp.h>
#include <libconfig.h++>
// #include <vortexcoordinates.h>


using namespace std;
using namespace libconfig;

namespace // namespace for program options
{ 
  const size_t ERROR_IN_COMMAND_LINE = 1; 
  const size_t SUCCESS = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
 
}

// const int vortex_start=8000; //ITP iterations before the phase disturbances are added (8000<n_it_ITP) 
Options opt;
vector<double> snapshot_times;



// Helper functions

 double vortex(int a, int b, int x, int y) //Vortex with phase [0,2*pi)          
{
        if(atan2(b-y,a-x)<0){ return 2*M_PI+atan2(b-y,a-x); } //atan2 is defined from [-pi,pi) so it needs to be changed to [0,2*pi)
	else{ return atan2(b-y,a-x); }        
}

void init_bh3(int argc, char** argv, Options &opt, vector<double> &snapshot_times);	
	
inline double gauss(double & x,double & y){return (exp(-x*x-y*y));} //A simple Gaussian

inline void printInitVar()
	{
		std::cout 	<< "Initial Values of this run:"<< endl
					<< "Gridsize in x direction: " << opt.grid[1] << "\t" << "omega_x = " << opt.omega_x << endl
			  		<< "Gridsize in y direction: " << opt.grid[2] << "\t" << "omega_y = " << opt.omega_y << endl
			  		<< "Expansion Factor: " << opt.exp_factor << "\t" << "Number of particles: " << opt.N << "\t" << "Interaction constant: " << opt.g << endl
					<< "Total time of the ITP-Step: " << opt.n_it_ITP << endl
					<< "Total time of the RTE-Step: " << opt.n_it_RTE << endl << endl;
	}


//>>>>>main program<<<<< 

int main( int argc, char** argv) 
{	

	// Initialize all option variables

	init_bh3(argc, argv, opt, snapshot_times);
	

	// Beginning of the options block
try 
{ 
    /** Define and parse the program options 
     */ 
    namespace po = boost::program_options; 
    po::options_description desc("Options"); 
    desc.add_options() 
      ("help,h", "Print help messages.") 
      ("xgrid,x",po::value<int>(&opt.grid[1]),"Gridsize in x direction.")
      ("ygrid,y",po::value<int>(&opt.grid[2]),"Gridsize in y direction.")
      ("gauss",po::value<bool>(&opt.startgrid[0]),"Initial Grid has gaussian form.")
      ("vortices",po::value<bool>(&opt.startgrid[1]),"Add Vortices to the grid.")
      ("itp",po::value<int>(&opt.n_it_ITP),"Total runtime of the ITP-Step.")
      ("rte",po::value<int>(&opt.n_it_RTE),"Total runtime of the RTE-Step.")
      ("number,N",po::value<double>(&opt.N),"Number of particles.")
      ("expansion,e",po::value<complex<double> >(&opt.exp_factor),"Expansion Factor");
      ("interaction,g",po::value<double> (&opt.g),"Interaction Constant");

	po::positional_options_description positionalOptions; 
	positionalOptions.add("xgrid", 1); 
	positionalOptions.add("ygrid", 1);   
 
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
        		  << "an expanding coordinate system.The implemented algorithm to solve the GPE" << endl
        		  << "is a 4-th order Runge-Kutta Integration because of this." << std::endl 
                  << desc << std::endl; 
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

    // print the initial values of the run to the console

    printInitVar(); 
		
	// Initialize the needed grid object and run object

	ComplexGrid* startgrid = new ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	RK4* run = new RK4(startgrid,opt);

	// if the given value is true, initialize the startgrid with a gaussian distribution

	if(opt.startgrid[0]==true) 
	{
	for(int i=0;i<opt.grid[1];i++)
    {
        for(int j=0;j<opt.grid[2];j++)
        {
                        
            double xfactor;
            xfactor = gauss(run->x_axis[i],run->y_axis[j]);           
            complex<double> factor (xfactor,0);             
            run->pPsi->at(0,i,j,0) = factor;                       
        }   
    }
	}


	// set the datafile identifier name and save the initial grid

    opt.name = "INIT";
	run->save_2D(run->pPsi,opt);
	
	// runtime checker
	double start,end[2];
	start = omp_get_wtime();

	//====> Imaginary Time Propagation (ITP)
	opt.name = "ITP1";
	opt.n_it_ITP = 8001;
	run->itpToTime(opt);

//////////// VORTICES ////////////////

// 	//*****Vortex Initial Coordinates*****

// //Spacing (6% across and 3% up relative to the lattice grid size)
// const int across=6,up=3;
// int half_across=across/2;

// //Central Vortex (at the origin)
// int x_1=opt.grid[1]/2;
// int y_1=opt.grid[2]/2; 

// //Ring 1 
// // int x_2=(100-across)*opt.grid[1]/200,x_3=(100+across)*opt.grid[1]/200,x_4=(100+half_across)*opt.grid[1]/200,x_5=(100-half_across)*opt.grid[1]/200,x_6=(100-half_across)*opt.grid[1]/200,x_7=(100+half_across)*opt.grid[1]/200;
// // int y_2=opt.grid[2]/2,y_3=opt.grid[2]/2,y_4=(100+up)*opt.grid[2]/200,y_5=(100+up)*opt.grid[2]/200,y_6=(100-up)*opt.grid[2]/200,y_7=(100-up)*opt.grid[2]/200; 

// //Ring 2
// // int x_8=(100-2*across)*opt.grid[1]/200,x_9=(100-across-half_across)*opt.grid[1]/200,x_10=(100-across)*opt.grid[1]/200,x_11=opt.grid[1]/2,x_12=(100+across)*opt.grid[1]/200,x_13=(100+across+half_across)*opt.grid[1]/200,x_14=(100+2*across)*opt.grid[1]/200,x_15=(100+across+half_across)*opt.grid[1]/200,x_16=(100+across)*opt.grid[1]/200,x_17=opt.grid[1]/2,x_18=(100-across)*opt.grid[1]/200,x_19=(100-across-half_across)*opt.grid[1]/200;
// // int y_8=opt.grid[2]/2,y_9=(100+up)*opt.grid[2]/200,y_10=(100+2*up)*opt.grid[2]/200,y_11=(100+2*up)*opt.grid[2]/200,y_12=(100+2*up)*opt.grid[1]/200,y_13=(100+up)*opt.grid[2]/200,y_14=opt.grid[2]/2,y_15=(100-up)*opt.grid[2]/200,y_16=(100-2*up)*opt.grid[2]/200,y_17=(100-2*up)*opt.grid[2]/200,y_18=(100-2*up)*opt.grid[2]/200,y_19=(100-up)*opt.grid[2]/200;

// //Ring 3
// // int x_20=(100-3*across)*opt.grid[1]/200,x_21=(100-2*across-half_across)*opt.grid[1]/200,x_22=(100-2*across)*opt.grid[1]/200,x_23=(100-across-half_across)*opt.grid[1]/200,x_24=(100-half_across)*opt.grid[1]/200,x_25=(100+half_across)*opt.grid[1]/200,x_26=(100+across+half_across)*opt.grid[1]/200,x_27=(100+2*across)*opt.grid[1]/200,x_28=(100+2*across+half_across)*opt.grid[1]/200,x_29=(100+3*across)*opt.grid[1]/200,x_30=(100+2*across+half_across)*opt.grid[1]/200,x_31=(100+2*across)*opt.grid[1]/200,x_32=(100+across+half_across)*opt.grid[1]/200,x_33=(100+half_across)*opt.grid[1]/200,x_34=(100-half_across)*opt.grid[1]/200,x_35=(100-across-half_across)*opt.grid[1]/200,x_36=(100-2*across)*opt.grid[1]/200,x_37=(100-2*across-half_across)*opt.grid[1]/200;
// // int y_20=opt.grid[2]/2,y_21=(100+up)*opt.grid[2]/200,y_22=(100+2*up)*opt.grid[2]/200,y_23=(100+3*up)*opt.grid[2]/200,y_24=(100+3*up)*opt.grid[2]/200,y_25=(100+3*up)*opt.grid[2]/200,y_26=(100+3*up)*opt.grid[2]/200,y_27=(100+2*up)*opt.grid[2]/200,y_28=(100+up)*opt.grid[2]/200,y_29=opt.grid[2]/2,y_30=(100-up)*opt.grid[2]/200,y_31=(100-2*up)*opt.grid[2]/200,y_32=(100-3*up)*opt.grid[2]/200,y_33=(100-3*up)*opt.grid[2]/200,y_34=(100-3*up)*opt.grid[2]/200,y_35=(100-3*up)*opt.grid[2]/200,y_36=(100-2*up)*opt.grid[2]/200,y_37=(100-up)*opt.grid[2]/200;  

// //Ring 4
// // int x_38=(100-4*across)*opt.grid[1]/200,x_39=(100-3*across-half_across)*opt.grid[1]/200,x_40=(100-3*across)*opt.grid[1]/200,x_41=(100-2*across-half_across)*opt.grid[1]/200,x_42=(100-2*across)*opt.grid[1]/200,x_43=(100-across)*opt.grid[1]/200,x_44=opt.grid[1]/2,x_45=(100+across)*opt.grid[1]/200,x_46=(100+2*across)*opt.grid[1]/200,x_47=(100+2*across+half_across)*opt.grid[1]/200,x_48=(100+3*across)*opt.grid[1]/200,x_49=(100+3*across+half_across)*opt.grid[1]/200,x_50=(100+4*across)*opt.grid[1]/200,x_51=(100+3*across+half_across)*opt.grid[1]/200,x_52=(100+3*across)*opt.grid[1]/200,x_53=(100+2*across+half_across)*opt.grid[1]/200,x_54=(100+2*across)*opt.grid[1]/200,x_55=(100+across)*opt.grid[1]/200,x_56=opt.grid[1]/2,x_57=(100-across)*opt.grid[1]/200,x_58=(100-2*across)*opt.grid[1]/200,x_59=(100-2*across-half_across)*opt.grid[1]/200,x_60=(100-3*across)*opt.grid[1]/200,x_61=(100-3*across-half_across)*opt.grid[1]/200;
// // int y_38=opt.grid[2]/2,y_39=(100+up)*opt.grid[2]/200,y_40=(100+2*up)*opt.grid[2]/200,y_41=(100+3*up)*opt.grid[2]/200,y_42=(100+4*up)*opt.grid[2]/200,y_43=(100+4*up)*opt.grid[2]/200,y_44=(100+4*up)*opt.grid[2]/200,y_45=(100+4*up)*opt.grid[2]/200,y_46=(100+4*up)*opt.grid[2]/200,y_47=(100+3*up)*opt.grid[2]/200,y_48=(100+2*up)*opt.grid[2]/200,y_49=(100+up)*opt.grid[2]/200,y_50=opt.grid[2]/2,y_51=(100-up)*opt.grid[2]/200,y_52=(100-2*up)*opt.grid[2]/200,y_53=(100-3*up)*opt.grid[2]/200,y_54=(100-4*up)*opt.grid[2]/200,y_55=(100-4*up)*opt.grid[2]/200,y_56=(100-4*up)*opt.grid[2]/200,y_57=(100-4*up)*opt.grid[2]/200,y_58=(100-4*up)*opt.grid[2]/200,y_59=(100-3*up)*opt.grid[2]/200,y_60=(100-2*up)*opt.grid[2]/200,y_61=(100-up)*opt.grid[2]/200;




// 	double Phase[opt.grid[1]][opt.grid[2]];	
// 	ComplexGrid c = ComplexGrid(opt.grid[0],opt.grid[2],opt.grid[2],opt.grid[3]);

//    	for(int i=0;i<opt.grid[1];i++)
//  	{
//       		for(int j=0;j<opt.grid[2];j++)
//        		{
//  		        Phase[i][j] = run->phase_save(run->pPsi,i,j)+vortex(i,j,x_1,y_1);
// 			  /*Ring 1*/
// 				// +vortex(i,j,x_2,y_2)+vortex(i,j,x_3,y_3)+vortex(i,j,x_4,y_4)+vortex(i,j,x_5,y_5)+vortex(i,j,x_6,y_6)+vortex(i,j,x_7,y_7)
//  			  /*Ring 2*/
//  				// +vortex(i,j,x_8,y_8)+vortex(i,j,x_9,y_9)+vortex(i,j,x_10,y_10)+vortex(i,j,x_11,y_11)+vortex(i,j,x_12,y_12)+vortex(i,j,x_13,y_13)+vortex(i,j,x_14,y_14)+vortex(i,j,x_15,y_15)+vortex(i,j,x_16,y_16)+vortex(i,j,x_17,y_17)+vortex(i,j,x_18,y_18)+vortex(i,j,x_19,y_19)
//  			  /*Ring 3*/
//  				// +vortex(i,j,x_20,y_20)+vortex(i,j,x_21,y_21)+vortex(i,j,x_22,y_22)+vortex(i,j,x_23,y_23)+vortex(i,j,x_24,y_24)+vortex(i,j,x_25,y_25)+vortex(i,j,x_26,y_26)+vortex(i,j,x_27,y_27)+vortex(i,j,x_28,y_28)+vortex(i,j,x_29,y_29)+vortex(i,j,x_30,y_30)+vortex(i,j,x_31,y_31)+vortex(i,j,x_32,y_32)+vortex(i,j,x_33,y_33)+vortex(i,j,x_34,y_34)+vortex(i,j,x_35,y_35)+vortex(i,j,x_36,y_36)+vortex(i,j,x_37,y_37) 
//  			  /*Ring 4*/
//  			  	// +vortex(i,j,x_38,y_38)+vortex(i,j,x_39,y_39)+vortex(i,j,x_40,y_40)+vortex(i,j,x_41,y_41)+vortex(i,j,x_42,y_42)+vortex(i,j,x_43,y_43)+vortex(i,j,x_44,y_44)+vortex(i,j,x_45,y_45)+vortex(i,j,x_46,y_46)+vortex(i,j,x_47,y_47)+vortex(i,j,x_48,y_48)+vortex(i,j,x_49,y_49)+vortex(i,j,x_50,y_50)+vortex(i,j,x_51,y_51)+vortex(i,j,x_52,y_52)+vortex(i,j,x_53,y_53)+vortex(i,j,x_54,y_54)+vortex(i,j,x_55,y_55)+vortex(i,j,x_56,y_56)+vortex(i,j,x_57,y_57)+vortex(i,j,x_58,y_58)+vortex(i,j,x_59,y_59)+vortex(i,j,x_60,y_60)+vortex(i,j,x_61,y_61);

//  			 // compute psi by using the initial psi^2 and adding the phase 

//        		}
//     	}

//     for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ c(0,i,j,0)=polar(abs(run->pPsi->at(0,i,j,0)),Phase[i][i]);} }

//     for(int i=0;i<opt.grid[1];i++){for(int j=0;j<opt.grid[2];j++){ run->pPsi->at(0,i,j,0) = c(0,i,j,0); } }


	// if the given value is true, add vortices to the startgrid
	if(opt.startgrid[1]==true)
    {
    run->pPsi = create_Vortex_start_Grid3(run->pPsi,opt,opt.cV*opt.rV,opt.cV,opt.rV,opt.Q);
   	cout << "Vortices added." << endl;
   	opt.name = "VORT";
   	run->save_2D(run->pPsi,opt);
	}

////// END VORTICES //////////


    opt.name = "ITP2";
	opt.n_it_ITP = 2001;
	run->itpToTime(opt);

	end[0] = omp_get_wtime();
	cout << "ITP took " << end[0] - start << " seconds." << endl << endl;


	//====> Real Time Expansion (RTE)
	opt.name = "RTE";
	run->rteToTime(opt);

	end[1] = omp_get_wtime();
	cout << "RTE took " << end[1] - start << " seconds." << endl;
	cout << "Run finished." << endl;

	// Everything finished here, plots and cleanup remaining	

    delete startgrid;
	delete run;
 
 
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

void init_bh3(int argc, char** argv, Options &opt, vector<double> &snapshot_times)
{
	libconfig::Config cfg;

	  // Read the file. If there is an error, report it and exit.
	  try
	  {
	    cfg.readFile("run.cfg");
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

	 opt.N          = root["RunOptions"]["N"];
	 opt.min_x      = root["RunOptions"]["min_x"]; 					
	 opt.min_y      = root["RunOptions"]["min_y"]; 					
	 opt.grid[0]    = root["RunOptions"]["grid0"];				
	 opt.grid[1]    = root["RunOptions"]["grid1"];				
	 opt.grid[2]    = root["RunOptions"]["grid2"];	   			
	 opt.grid[3]    = root["RunOptions"]["grid3"];				
	 opt.g          = root["RunOptions"]["g"]; 						
	 opt.n_it_RTE   = root["RunOptions"]["n_it_RTE"]; 				
	 opt.n_save_RTE = root["RunOptions"]["n_save_RTE"]; 			
	 opt.n_it_ITP   = root["RunOptions"]["n_it_ITP"];				
	 opt.n_save_ITP = root["RunOptions"]["n_save_ITP"]; 			
	 // auto test1       = root["opt"]["name"];
	 // cout << root["opt"]["name"]; << endl;
	 opt.times      = root["RunOptions"]["times"]; 	    			
	 opt.ITP_step   = root["RunOptions"]["ITP_step"]; 				
	 opt.RTE_step   = root["RunOptions"]["RTE_step"]; 				
	 opt.cV         = root["RunOptions"]["cV"]; 		   	   		
	 opt.rV         = root["RunOptions"]["rV"]; 			   		
	 opt.Q          = root["RunOptions"]["Q"];
	 // opt.startgrid[0] = root["RunOptions"]["startgrid0"];
  	 // opt.startgrid[1] = root["RunOptions"]["startgrid1"];

  	   try
	  {
	    string fufufu2 = root["RunOptions"]["name"].c_str();
	    cout << "Store name: " << fufufu2 << endl << endl;
	  }
	  catch(const SettingNotFoundException &nfex)
	  {
	    cerr << "No 'name' setting in configuration file." << endl;
	  }

	  try
	  {
	  	int fufufu3 = cfg.lookup("cV");
	  	cout << "test read of cV: " << fufufu3 << endl << endl;
	  }
	  catch(const SettingNotFoundException &nfex)
	  {
	  	cerr << "No cV setting found in config file." << endl;
	  }



    // Set Parameters manually (default values)
	//opt.timestepsize = 0.2;
	//opt.delta_t.resize(0);   
    //opt.delta_t[0]=0.2;
    //opt.delta_t[1]=0.4;
 
	// opt.N = 1000;  //normed for 512*512 N=64*50000
	
	// opt.min_x = 2;
	// opt.min_y = 2;
	
	// opt.grid[0] = 1;
	// opt.grid[1] = 600;
	// opt.grid[2] = 600;
	// opt.grid[3] = 1;
	// opt.U = 3e-5;
	
	// opt.g.resize(1);
	// opt.g = 1;

	// opt.cV = 2;
    // opt.rV = 2;
    // opt.Q = 1;
	
	opt.omega_x = complex<double>(100,0);
	opt.omega_y = complex<double>(150,0);
	opt.scale_factor = complex<double>(0,0); //Scale factor
	opt.t_abs = complex<double>(0,0); //Absolute time 
	opt.exp_factor = complex<double>(1.2,0); //Expansion factor
	// opt.n_it_RTE = 101;
	// opt.n_save_RTE = 50;
	// opt.n_it_ITP = 10000;
	// opt.n_save_ITP = 1000;
	// opt.times = 1;
	opt.name = "run"; // Must be an Integer
	opt.startgrid[0] = true; // gaussian packet
    opt.startgrid[1] = false; // add vortices

	// 	opt.klength[0] = 2.0;
	// 	opt.klength[1] = 2.0;
	// 	opt.klength[2] = 2.0;
	
	// opt.ITP_step=0.000001; //Time-step for the ITP (0.000001)
	// opt.RTE_step=0.00001; //Time-step for the RTE (0.00001)
	
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

