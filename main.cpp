/**************************************************************************
Title: Simulating the Expansion of Turbulent Bose-Einstein Condensates (2D) 
Author: Bartholomew Andrews
Last Update: 22/07/13
Website: www.bartholomewandrews.com
**************************************************************************/

#include <iostream>
#include <cmath>
#include <complex>
#include <complexgrid.h>
#include <exp_RK4_tools.h>
// #include <2dexpan.h>



using namespace std;

//*****Parameter Initialisation*****

// const double ITP_step=0.000001; //Time-step for the ITP (0.000001)
// const double RTE_step=0.00001; //Time-step for the RTE (0.00001)
// const int n_x=128,n_y=128; //Lattice size (1500x1500)  lo

const int vortex_start=8000; //ITP iterations before the phase disturbances are added (8000<n_it_ITP) 
const int n_save_ITP=1000; //Save ITP after every n_save_ITP iterations (initial state is auto saved)
const int n_it_ITP=10000; //Number of iterations for ITP (10000)
const int n_save_RTE=100; //Save RTE after every n_save_RTE iterations - intial state is auto saved (500)
const int n_it_RTE=101; //Number of iterations for RTE (2501) 
// const int name=1; //Name of output (must be an integer)

//*****Variable Declarations*****

int counter_ITP=0, counter_RTE=0; //Initialise the variables for the percentage loading

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

	Options opt;
	vector<double> snapshot_times;
	init_bh3(argc, argv, opt, snapshot_times);
// 	complex<double> h_x(((2.*min_x)/opt.grid[1]),0),h_y(((2.*min_y)/opt.grid[2]),0); //Lattice constants (chosen for symmetric initial axis ranges)
// 	
// 	
// 	
// 	ComplexGrid * pPsi;
// 	pPsi = new ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
// 		
// 	
// 		

// 	
// 	ComplexGridPtr pPsiCopy(pPsi);
// 	
// 	  for(i=0;i<opt.grid[1];i++)
//         {
// 	        for(j=0;j<opt.grid[2];j++)
// 	        {
// 		        save_obdm(x_axis[i]*real(lambda_x(t_abs)),y_axis[j]*real(lambda_y(t_abs)),norm(pPsiCopy->at(0,i,j,0)),phase_save(i,j));
// 	        }
//                 blank_line();
//         }
//         blank_line();
	
	
	// initialize_psi(n_x,n_y,h_x,h_y);
	cout << "Starting Grid now" << endl;
	ComplexGrid* startgrid;
	startgrid = new ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	RK4* run;
	run = new RK4(startgrid,opt);
	cout << "Initializing Grid now" << endl;
// 	delete startgrid;
	for(int i=0;i<opt.grid[1];i++) //Initialise the wavefunction
	{
		for(int j=0;j<opt.grid[2];j++)
		{
// 			double x;
// 			double y;
// 			x = run->x_axis[i];
// 			cout << run->x_axis[i] << "   ";
// 			y = run->y_axis[j];
// 			cout << run->y_axis[j] << "   ";
			double xfactor;
			xfactor = gauss(run->x_axis[i],run->y_axis[j]);
			
			complex<double> factor (xfactor,0);
// 			cout << factor << "   ";
			
		        run->pPsi->at(0,i,j,0) = factor;
// 			cout << run->pPsi->at(0,i,j,0) << endl;;
// 			run->pPhase->at(0,i,j,0) = (0,0); //!!!! Check if this should be double, and not complexdouble !!!!
			 
		}	
	}
	cout << "Grid initialized" << endl;
	complex<double> ausgabe;
	ausgabe = run->pPsi->at(0,2,3,0);
	cout << ausgabe << endl;
	
	//====> Imaginary Time Propagation (ITP)
	
	for(int k=0;k<n_it_ITP;k++)	//Time-evolution given by the 4th-order Runge-Kutta iteration			
	{
		run->ITP(run->pPsi,opt);
	


// 		if(k==vortex_start){add_vortex();} //Add the vortex near the end of the imaginary time propagation 
		
	       	if(k>0 && k%n_save_ITP==0){run->save_2D(run->pPsi,opt);} //New block every multiple of n_save_ITP
  
		counter_ITP+=1; //Loading counter for ITP
   		if(counter_ITP%(n_it_ITP/100)==0){cout<<"ITP "<<(counter_ITP/(n_it_ITP/100))<<"%"<<endl;}
	}
	
	//====> Real Time Expansion (RTE)

	for(int k=0;k<n_it_RTE;k++)	//Time-evolution given by the 4th-order Runge-Kutta iteration			
	{
		run->RTE(run->pPsi,opt);

		if(k>=0 && k%n_save_RTE==0){run->save_2D(run->pPsi,opt);} //New block every multiple of n_save_RTE
  
		opt.t_abs.real()+=opt.RTE_step; //Increment absolute time by the time-step t_RTE

		counter_RTE+=1; //Loading counter for RTE
   		if(counter_RTE%(n_it_RTE/100)==0){cout<<"RTE "<<(counter_RTE/(n_it_RTE/100))<<"%"<<endl;}
	}
	delete startgrid;
	delete run;

// 	 //Close the file
	return 0;
}

void init_bh3(int argc, char** argv, Options &opt, vector<double> &snapshot_times)
{
        // Parameter setzen
// 	opt.timestepsize = 0.2;
// 	opt.delta_t.resize(0);   
        //opt.delta_t[0]=0.2;
        //opt.delta_t[1]=0.4;

	opt.N = 1000;  //normed for 512*512 N=64*50000
	
	opt.min_x = 4;
	opt.min_y = 4;
	
	opt.grid[0] = 1;
	opt.grid[1] = 32;
	opt.grid[2] = 32;
	opt.grid[3] = 1;
// 	opt.U = 3e-5;
	
// 	opt.g.resize(1);
	opt.g = 15;
	
	opt.omega_x = (100,0);
	opt.omega_y = (150,0);
	opt.scale_factor = (0,0); //Scale factor
	opt.t_abs = (0,0); //Absolute time 
	opt.exp_factor = (1.5,0); //Expansion factor
	opt.name = 1; // Must be an Integer
    	
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

