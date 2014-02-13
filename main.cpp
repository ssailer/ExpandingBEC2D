/**************************************************************************
Title: Simulating the Expansion of Turbulent Bose-Einstein Condensates (2D) 
Author: Bartholomew Andrews
Last Update: 22/07/13
Website: www.bartholomewandrews.com
**************************************************************************/
#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <cmath>
#include <complex>
#include <complexgrid.h>
#include <exp_RK4_tools.h>
#include <bh3defaultgrid.h>
#include <omp.h>
// #include <2dexpan.h>



using namespace std;


namespace // namespace for program options
{ 
  const size_t ERROR_IN_COMMAND_LINE = 1; 
  const size_t SUCCESS = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
 
}

//*****Parameter Initialisation*****

// const double ITP_step=0.000001; //Time-step for the ITP (0.000001)
// const double RTE_step=0.00001; //Time-step for the RTE (0.00001)
// const int n_x=128,n_y=128; //Lattice size (1500x1500)  lo

const int vortex_start=8000; //ITP iterations before the phase disturbances are added (8000<n_it_ITP) 
// const int n_save_ITP=5000; //Save ITP after every n_save_ITP iterations (initial state is auto saved)
// const int n_it_ITP=15000; //Number of iterations for ITP (10000)
// const int n_save_RTE=500; //Save RTE after every n_save_RTE iterations - intial state is auto saved (500)
// const int n_it_RTE=2501; //Number of iterations for RTE (2501) 
// const int name=1; //Name of output (must be an integer)

//*****Variable Declarations*****

// int counter_ITP=0, counter_RTE=0; //Initialise the variables for the percentage loading

// double phase[opt.grid[1]][opt.grid[2]]; //Sum of all phases (used in the add_vortex() function) 

//*****Vortex Initial Coordinates*****
/*
//Spacing (6% across and 3% up relative to the lattice grid size)
const int across=6,up=3;
int half_across=across/2;

//Central Vortex (at the origin)
int x_1=n_x/2;
int y_1=n_y/2; 

//Ring 1 
int x_2=(100-across)*n_x/200,x_3=(100+across)*n_x/200,x_4=(100+half_across)*n_x/200,x_5=(100-half_across)*n_x/200,x_6=(100-half_across)*n_x/200,x_7=(100+half_across)*n_x/200;
int y_2=n_y/2,y_3=n_y/2,y_4=(100+up)*n_y/200,y_5=(100+up)*n_y/200,y_6=(100-up)*n_y/200,y_7=(100-up)*n_y/200; 

//Ring 2
int x_8=(100-2*across)*n_x/200,x_9=(100-across-half_across)*n_x/200,x_10=(100-across)*n_x/200,x_11=n_x/2,x_12=(100+across)*n_x/200,x_13=(100+across+half_across)*n_x/200,x_14=(100+2*across)*n_x/200,x_15=(100+across+half_across)*n_x/200,x_16=(100+across)*n_x/200,x_17=n_x/2,x_18=(100-across)*n_x/200,x_19=(100-across-half_across)*n_x/200;
int y_8=n_y/2,y_9=(100+up)*n_y/200,y_10=(100+2*up)*n_y/200,y_11=(100+2*up)*n_y/200,y_12=(100+2*up)*n_x/200,y_13=(100+up)*n_y/200,y_14=n_y/2,y_15=(100-up)*n_y/200,y_16=(100-2*up)*n_y/200,y_17=(100-2*up)*n_y/200,y_18=(100-2*up)*n_y/200,y_19=(100-up)*n_y/200;

//Ring 3
int x_20=(100-3*across)*n_x/200,x_21=(100-2*across-half_across)*n_x/200,x_22=(100-2*across)*n_x/200,x_23=(100-across-half_across)*n_x/200,x_24=(100-half_across)*n_x/200,x_25=(100+half_across)*n_x/200,x_26=(100+across+half_across)*n_x/200,x_27=(100+2*across)*n_x/200,x_28=(100+2*across+half_across)*n_x/200,x_29=(100+3*across)*n_x/200,x_30=(100+2*across+half_across)*n_x/200,x_31=(100+2*across)*n_x/200,x_32=(100+across+half_across)*n_x/200,x_33=(100+half_across)*n_x/200,x_34=(100-half_across)*n_x/200,x_35=(100-across-half_across)*n_x/200,x_36=(100-2*across)*n_x/200,x_37=(100-2*across-half_across)*n_x/200;
int y_20=n_y/2,y_21=(100+up)*n_y/200,y_22=(100+2*up)*n_y/200,y_23=(100+3*up)*n_y/200,y_24=(100+3*up)*n_y/200,y_25=(100+3*up)*n_y/200,y_26=(100+3*up)*n_y/200,y_27=(100+2*up)*n_y/200,y_28=(100+up)*n_y/200,y_29=n_y/2,y_30=(100-up)*n_y/200,y_31=(100-2*up)*n_y/200,y_32=(100-3*up)*n_y/200,y_33=(100-3*up)*n_y/200,y_34=(100-3*up)*n_y/200,y_35=(100-3*up)*n_y/200,y_36=(100-2*up)*n_y/200,y_37=(100-up)*n_y/200;  

//Ring 4
int x_38=(100-4*across)*n_x/200,x_39=(100-3*across-half_across)*n_x/200,x_40=(100-3*across)*n_x/200,x_41=(100-2*across-half_across)*n_x/200,x_42=(100-2*across)*n_x/200,x_43=(100-across)*n_x/200,x_44=n_x/2,x_45=(100+across)*n_x/200,x_46=(100+2*across)*n_x/200,x_47=(100+2*across+half_across)*n_x/200,x_48=(100+3*across)*n_x/200,x_49=(100+3*across+half_across)*n_x/200,x_50=(100+4*across)*n_x/200,x_51=(100+3*across+half_across)*n_x/200,x_52=(100+3*across)*n_x/200,x_53=(100+2*across+half_across)*n_x/200,x_54=(100+2*across)*n_x/200,x_55=(100+across)*n_x/200,x_56=n_x/2,x_57=(100-across)*n_x/200,x_58=(100-2*across)*n_x/200,x_59=(100-2*across-half_across)*n_x/200,x_60=(100-3*across)*n_x/200,x_61=(100-3*across-half_across)*n_x/200;
int y_38=n_y/2,y_39=(100+up)*n_y/200,y_40=(100+2*up)*n_y/200,y_41=(100+3*up)*n_y/200,y_42=(100+4*up)*n_y/200,y_43=(100+4*up)*n_y/200,y_44=(100+4*up)*n_y/200,y_45=(100+4*up)*n_y/200,y_46=(100+4*up)*n_y/200,y_47=(100+3*up)*n_y/200,y_48=(100+2*up)*n_y/200,y_49=(100+up)*n_y/200,y_50=n_y/2,y_51=(100-up)*n_y/200,y_52=(100-2*up)*n_y/200,y_53=(100-3*up)*n_y/200,y_54=(100-4*up)*n_y/200,y_55=(100-4*up)*n_y/200,y_56=(100-4*up)*n_y/200,y_57=(100-4*up)*n_y/200,y_58=(100-4*up)*n_y/200,y_59=(100-3*up)*n_y/200,y_60=(100-2*up)*n_y/200,y_61=(100-up)*n_y/200;*/


//*****Function Definitions*****

// double gauss(double x,double y){return (exp(-x*x-y*y));} //A simple Gaussian


	

void init_bh3(int argc, char** argv, Options &opt, vector<double> &snapshot_times);	
	
//>>>>>Main Program<<<<< 
inline double gauss(double & x,double & y){return (exp(-x*x-y*y));} //A simple Gaussian

int main( int argc, char** argv) 
{	
	//////////////// ALL NEED VARIABLES HERE ///////

	Options opt;
	vector<double> snapshot_times;
	init_bh3(argc, argv, opt, snapshot_times);
	
	///////////////// PROGRAM OPTIONS ///////////

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
      // ("name,n",po::value<string>(&opt.name),"Used in naming the data files.") // opt.name is the wrong variable for this purpose, CHANGE THIS!!
      ("itp",po::value<int>(&opt.n_it_ITP),"Total runtime of the ITP-Step.")
      ("rte",po::value<int>(&opt.n_it_RTE),"Total runtime of the RTE-Step.")
      ("expansion,e",po::value<complex<double> >(&opt.exp_factor),"Expansion Factor");
      //("add", "additional options") 
      //("like", "this"); 

	po::positional_options_description positionalOptions; 
	positionalOptions.add("xgrid", 1); 
	positionalOptions.add("ygrid", 1);   
 
    po::variables_map vm; 
    try 
    { 

		po::store(po::command_line_parser(argc, argv).options(desc)
					.positional(positionalOptions).run(), 
          			vm); 
      //po::store(po::parse_command_line(argc, argv, desc),  
      //          vm); // can throw 
 
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

////////////////////////////////////////////////
/////////////// PROGRAM MAIN ///////////////////
////////////////////////////////////////////////

    std::cout << "Initial Values of this run:"<< endl
			  << "Gridsize in x direction: " << opt.grid[1] << "\t" << "omega_x = " << opt.omega_x << endl
			  << "Gridsize in y direction: " << opt.grid[2] << "\t" << "omega_y = " << opt.omega_y << endl
			  << "Expansion Factor: " << opt.exp_factor << endl
			  << "Total time of the ITP-Step: " << opt.n_it_ITP << endl
			  << "Total time of the RTE-Step: " << opt.n_it_RTE << endl << endl;


	
	int cV = 2;
    int rV = 2;
    int Q = 1;
	
	ComplexGrid* startgrid = new ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	RK4* run = new RK4(startgrid,opt);
	// cout << "Run initialized with Grid: "<<opt.grid[1]<<"x"<<opt.grid[2]<< endl;

	for(int i=0;i<opt.grid[1];i++) //Initialise the wavefunction as gaussian by default with default_gridsize
    {
        for(int j=0;j<opt.grid[2];j++)
        {
                        
            double xfactor;
            xfactor = gauss(run->x_axis[i],run->y_axis[j]);           
            complex<double> factor (xfactor,0);             
            run->pPsi->at(0,i,j,0) = factor;                       
        }   
    }

    startgrid = create_Vortex_start_Grid3(startgrid,opt,4,cV,rV,Q);

    opt.name = "INIT";
	run->save_2D(run->pPsi,opt);
	
	double start,end[2];
	start = omp_get_wtime();

	//====> Imaginary Time Propagation (ITP)
	run->itpToTime(opt);

	end[0] = omp_get_wtime();
	cout << "ITP took " << end[0] - start << " seconds." << endl << endl;

	//====> Real Time Expansion (RTE)
	run->rteToTime(opt);

	end[1] = omp_get_wtime();
	cout << "RTE took " << end[1] - start << " seconds." << endl;
	cout << "Run finished." << endl;

	
	//Python Plot
/*
	
		for(int i=0;i<snapshot_times.size(); i++)
  {

	ofstream fs;
	fs.open((dir+string("Spectrum.py")).c_str(), ios_base::trunc | ios_base::out);


	  

        fs << "#!/usr/bin/python" << endl;
	fs << "# -*- coding: utf-8 -*-" << endl;

	fs << "from matplotlib import use" << endl;
	fs << "use('Agg')" << endl;
	fs << "from matplotlib import rc" << endl;
        fs << "import scipy, pylab" << endl;
        fs << "import scipy.optimize" << endl;
	fs << "import sys" << endl;
        fs << "import matplotlib.colors " << endl;
	fs << "import pylab as p "<< endl;
	fs << "import numpy as np" << endl;
	fs << "import matplotlib.pyplot as plt" << endl;
        fs << "from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, MultipleLocator" << endl;
        fs << "from matplotlib.ticker import FixedFormatter" << endl;
	fs << "from mpl_toolkits.axes_grid1 import make_axes_locatable" << endl;
	fs << "import math as m" << endl;
	
	fs << "def main():" << endl;
        
	
	fs << "\tpath = '" << dirname << "Spectrum" << "_Path_"  << snapshot_times[i] <<"time"<<".png'" << endl;	
	fs << "\tfig = plt.figure(figsize=(8.3,5.7), dpi=100)" << endl;	
	fs << "\tdata = np.loadtxt('"<< (dir + string("radial_avgs.dat")).c_str() << "', delimiter='\\t', unpack=True, usecols = (0,1,2,3,4,5))" << endl;
	fs << "\ttime=" << snapshot_times[i] << endl;
        fs << "\tdim ="  << meansFD[i].ares_1.number.size() << endl;
	fs << "\tsnapshot=" << i << endl;

        fs << "\tstart=" << ((meansFD[i].ares_1.number.size())*i)+1 << endl;
	fs << "\tend=" << (meansFD[i].ares_1.number.size())*(i+1) << endl;

        fs << "\tstart2=" << ((meansFD[i].ares_1.k.size())*i)+1 << endl;
	fs << "\tend2=" << (meansFD[i].ares_1.k.size())*(i+1) << endl;

        fs << "\tHealing_k=" << sqrt((opt.N/(opt.grid[0]*opt.grid[1]*opt.grid[2]))*opt.U) << endl;
	

	  fs << "\tfig.subplots_adjust(hspace = 0.10, wspace = 0.40, right  = 0.85, left=0.15, bottom = 0.10, top = 0.90)" << endl;

      	  fs << "\txfit=np.arange(0.03,0.25,0.01)" << endl;
	  fs << "\tyfit=(pow(xfit,(-6))*0.1)  "<< endl;

	  //fs << "\txfit2=np.arange(0.25,2.8,0.01)" << endl;
	  //fs << "\tyfit2=(pow(xfit2,(-2))*20)  "<< endl;

	  // fs << "\txfit3=np.arange(0.03,2.8,0.01)" << endl;
	  // fs << "\tyfit3=(pow(xfit3,(-4))*20)  "<< endl;
       

          //fs << "\txfit4=np.arange(0.03,0.3,0.01)" << endl;
	  //fs << "\tyfit4=(pow(xfit4,(-4.6666666666666666666666666666666666666666))*1)  "<< endl;
       


	//fs << "\tn1=data[2]" << endl;
        //fs << "\tn2=data[1]" << endl;


          fs << "\tn1=data[1]" << endl;
	  fs << "\tn2=data[2]" << endl;
        
	  fs << "\tn_q=data[3]" << endl;
	  fs << "\tn_i=data[4]" << endl;
	  fs << "\tn_c=data[5]" << endl;

      





       
	  fs << "\tn4=n2[start:end]" << endl;
	  fs << "\tn3=n1[start2:end2]" << endl;

	  fs << "\tn_qneu=n_q[start:end]" << endl;
	  fs << "\tn_ineu=n_i[start:end]" << endl;
	  fs << "\tn_cneu=n_c[start:end]" << endl;	
	
	   fs << "\tax = fig.add_subplot(111)" << endl;	
	   fs << "\tim=plt.plot(n3, n4,'b.',n3,n_ineu,'y.',n3, n_qneu,'r.',n3,n_cneu,'g.')" << endl;

	  

	   fs << "\tplt.setp(im, 'markersize', 3)" << endl; 
	   fs << "\tplt.xlim(0.02,3.2)" << endl;


           
	   fs << "\tax.set_title('Time: $"<< snapshot_times[i]<<","<<"Particles:" <<meansFD[i].ares_1.particle_count/opt.N<<"$ ') " << endl;
	   fs << "\tax.set_yscale('log')" << endl;
	   fs << "\tax.set_xscale('log')" << endl; 
           
           fs << "\tax.set_xticks([0.1, 0.3, 1, 2, 3, Healing_k])" << endl;
           fs << "\tax.set_xticklabels(['0.1', '0.3','1', '2', '3','H'])" << endl;
	   fs << "\tax.set_xlabel('$k$')" << endl;
	   fs << "\tax.set_ylabel('$n(k)$')" << endl;
	   fs << "\tax.yaxis.label.set_size(15) " << endl;
	   fs << "\tax.xaxis.label.set_size(15) " << endl;
     
	fs << "\tfig.savefig(path, dpi=100)" << endl;
        fs << "\tsys.exit(3)" << endl;
        fs << "if __name__=='__main__':" << endl;

        fs << "\tmain()" << endl;

  

        fs.close();
	
       system((string("python ") + dir + string("Spectrum.py")).c_str());
       }	
       */
    delete startgrid;
	delete run;

////////////////////////////////////////////////
/////////////// PROGRAM END ////////////////////
////////////////////////////////////////////////

 
 
  } 
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
        // Parameter setzen
	//opt.timestepsize = 0.2;
	//opt.delta_t.resize(0);   
    //opt.delta_t[0]=0.2;
    //opt.delta_t[1]=0.4;

	opt.N = 1000;  //normed for 512*512 N=64*50000
	
	opt.min_x = 4;
	opt.min_y = 4;
	
	opt.grid[0] = 1;
	opt.grid[1] = 256;
	opt.grid[2] = 256;
	opt.grid[3] = 1;
// 	opt.U = 3e-5;
	
// 	opt.g.resize(1);
	opt.g = 15;
	
	opt.omega_x = complex<double>(100,0);
	opt.omega_y = complex<double>(150,0);
	opt.scale_factor = complex<double>(0,0); //Scale factor
	opt.t_abs = complex<double>(0,0); //Absolute time 
	opt.exp_factor = complex<double>(2.0,0); //Expansion factor
	opt.n_it_RTE = 1001;
	opt.n_it_ITP = 10001;
	opt.n_save_RTE = 1000;
	opt.n_save_ITP = 10000;
	opt.times = 1;
	opt.name = "run"; // Must be an Integer
    	
// 	opt.klength[0] = 2.0;
// 	opt.klength[1] = 2.0;
// 	opt.klength[2] = 2.0;
	
	opt.ITP_step=0.000001; //Time-step for the ITP (0.000001)
	opt.RTE_step=0.00001; //Time-step for the RTE (0.00001)
	
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
          // snapshot_times[41]=1900000;*/
      
		

        /*		for(int i = 0; i < snapshot_times.size(); i++)///////////////////////was passiert hier/????????////====warum die snapshot zeit veraendern???========?????
                {
                int time_res = 1000;
                snapshot_times[i] = (int) (500.*pow(10,((double)i/10.)));
                    //cout<<snapshot_times[i]<<endl;
                        //snapshot_times[i] = i*time_res + 25000;
                        }*/
	
}

