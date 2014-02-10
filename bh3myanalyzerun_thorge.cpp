#include <sys/stat.h>
#include <sys/types.h>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <engine.h>

#include <stdio.h>
#include <stdlib.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include "engine.h"



#include <phasegrid.h>
#include <averageclass.h>
#include <bh3cudapropagator.h>
#include <complexgrid.h>
#include <bh3defaultgrid.h>
#include <bh3binaryfile.h>
#include <gauss_random.h>
#include <bh3observables.h>
#include <bh3observablesFD.h>
#include <bh3save.h>
#include <omp.h>
#include <wrapped_cuda_functions.h>

#define FRAME_DURATION 10
#define WAIT_AFTER 200
#define NUM_SINGLE_PLOTS 999999
#define PNG_WIDTH 640 
#define PNG_HEIGHT 480
#define FIELD_WIDTH 6
#define PATHS 1
using namespace std;

void init_bh3(int argc, char** argv, PathOptions &opt, vector<double> &snapshot_times);
void plot(const string &dirname, const PathOptions& opt, vector<double> &snapshot_times, vector <Bh3EvaluationFD::Averages> &meansFD, int hh, int cV, int rV, int Q);
void plot3(ComplexGrid *xx, ComplexGrid *yy, const string &dirname,const PathOptions& opt,const vector<ComplexGrid> &grid1,double time);


////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv) 
{
  	


        int hh=0;
        int cV=2;
        int rV=2;
        int Q =3;

	PathOptions opt;
	vector<double> snapshot_times;
	
	string rm = "rm ";
	
	init_bh3(argc, argv, opt, snapshot_times);
	ComplexGrid::set_fft_planning_rigorosity(FFTW_MEASURE);
	init_random();
	stringstream dstr;
	dstr << "bh3cuda_TS" << opt.timestepsize
	<< "_G" << opt.grid[0] << "_" << opt.grid[1] << "_" << opt.grid[2]
	<< "_N" << opt.N
	<< "_U" << opt.U
	     << "_KL" << opt.klength[0] << "_" << opt.klength[1] << "_" << opt.klength[2]//////////////////////////////////////////////////////////////////////Wofuer opt.klingth////////////
	<< "_a" << opt.u_12;
	string dirname = dstr.str();
	initialize_evaluation_dir(dirname, "_run", list<string>());
	mkdir((dirname + "/temp").c_str(), 0755);
	

       

	AverageClass<vector<Bh3EvaluationFD::Averages> >avFD;
	
	for(int k = 0; k < PATHS; k++)
	{
	  if(init_cuda_device())//////////////////////////////////////was wird hier gemacht///////////////////////////////////////////////////////////////////////////////////////////////
		{
			unsigned int t;
			time_t timer = time(NULL);
			
			ComplexGrid *start1;
			ComplexGrid *start2;
			Bh3Propagator *cp;
			
			double tau_q = 5000.;  ///Wofuer zur Hoelle??????????????????////////////////////////////////////////////////////////////////////////////////////////////////////////
		
                         
						
			start1 = create_Vortex_start_Grid2(opt,4,cV,rV,Q);
			start2 = create_Vortex_start_Grid2(opt,4,cV,rV,Q);
                        //start1 = create_Default_Start_Grid(opt,1);
			//start2 = create_Default_Start_Grid(opt,1);
			//cout<<start1->at(1,0,0)<<endl;
			// start1 = create_flowing_Start_Grid(opt, 1.0, 1.0);
			//start2 = create_flowing_Start_Grid(opt, -1.0, 1.0);

			cp = new Bh3CudaPropagator(opt, *start1, *start2, tau_q);

			cp->set_Omega(0.0);
			delete start1;
			delete start2;
			
			vector<Bh3EvaluationFD::Averages> path_resFD(snapshot_times.size());
                       	//vector<Bh3Evaluation> path_resFD(snapshot_times.size(), Bh3Evaluation(opt,Bh3Evaluation::RSpace));
			//vector<Bh3Evaluation::Averages> path_resFD(snapshot_times.size());


			for(int j = 0; j < snapshot_times.size(); j++)
			{
				#pragma omp parallel
				{
					#pragma omp master
					{
						//for(int i = 0; i < WAIT_AFTER; i++)
						//{    
					  
					               double time = snapshot_times[j];
							cp->propagateToTime(time);
					 std::cout<< "hier kommt die Maus2" <<std::endl;

							stringstream filestr;
							filestr << dirname << "/temp/Bh3RunSnapshot" << j << ".bin";
							string filename = filestr.str();
							Bh3BinaryFile *bf = new Bh3BinaryFile(filename, opt, Bh3BinaryFile::out);////////////////////////////////////////////was hat das
                                                                                                                    /////////was hat das alles auf sich mit den Stringstreams zum schreiben
							///////////////////////////////////////////////////////////////////////////////////////////////////////////////// lesen und der omp auf sich     
                                                                                                                            ////////////////////////////////////////////was hat das
                                                                                                                        ////////////////////////////////////////////was hat das
                                                                                                                         ////////////////////////////////////////////was hat das
							bf->append_snapshot(time, cp->getRGrids(1));
							bf->append_snapshot(cp->get_Omega(), cp->getRGrids(2));
							delete bf;
							cout << "propagating...." << omp_get_thread_num() << endl;
							
							#pragma omp task default(shared) firstprivate(j,filename)
							{

								vector<ComplexGrid> r1;
								vector<ComplexGrid> r2;
                                                                
                                                      		//ComplexGrid *r3;
                                                                //ComplexGrid *r4;	
					                    
								
								PathOptions options;
								//double time;
								double Omega;
								Bh3BinaryFile *temp_file = new Bh3BinaryFile(filename, options, Bh3BinaryFile::in);/////obwohl filename eine Refrenz und auch opti    
                                                        

								options = temp_file->get_options();

                                                                
 
								temp_file->get_snapshot(time, r1);
								temp_file->get_snapshot(Omega, r2);
								
								delete temp_file;
								system((rm + filename).c_str());

								cout << "evaluating ..." << omp_get_thread_num() << endl;

								Bh3EvaluationFD ev(options, true);///////////weiter

								ev.setTime(time);/////////////////////////////////////////////gibt es eine Funktion die die Zeit?????????????
								ev.setOmega(Omega);
								ev.setData(r1,r2,Bh3EvaluationFD::RSpace);

                                                                                                                 				                                                
							        ev.calc_radial_averages(); ////Wichtig//////////FD kacke/////////////////////////////////////////    
							       

							        //r3=ev.gv2;
								//r4=ev.sv2;
								//r5=ev.mv2;
				      
                                                                //plot3(r3, r4, dirname, opt,r1,time);
							
								// output as png
								//plot(dirname, options, j*WAIT_AFTER + i, r1, r2, time);
		
	
							}	// task
						//}	// for i
	
					}

				}
						cout<<"hallddddddo"<<endl;
			}

			avFD.average(path_resFD);  ///Wichtig/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			       					      					  
			delete cp;
			
			cout << "Path CUDA took " << time(NULL) - timer << "seconds" << endl;

		}		// if init_cuda_device
		else
		{
			cout << "Could not initialize Cuda device!" << endl;
		}
		
		if(k%1==0)
		{
                  hh++;
		  vector<Bh3EvaluationFD::Averages> meansFD = avFD.av();
		  plot(dirname, opt, snapshot_times, meansFD, hh, rV, cV, Q);
		  //void save_energies(const string &dirname, int p, const vector<Bh3Evaluation::Averages> *av, const vector<Bh3Evaluation::Averages> *std_dev = NULL)
     		  cout<< "run " << k << " done" << endl;

		}
	}
	
	vector <Bh3EvaluationFD::Averages> meansFD = avFD.av();
       

 
	//plot(dirname, opt, snapshot_times, meansFD);
	stringstream rm_command;                         ///////////////////////////////////////////////////////////////////////////////////////////////////was passiert hier mit dem rm_command
	rm_command << "rm -r " << dirname << "/temp";
	system(rm_command.str().c_str());
	
	return 0;
}

void init_bh3(int argc, char** argv, PathOptions &opt, vector<double> &snapshot_times)
{
	// Parameter setzen
	opt.timestepsize = 0.2;
	opt.delta_t.resize(2);             ////////////////////wird dafuer verwednet um vor einem Zeitschritt gewisse Aufnahmen zu machen um Abletiungen zu berechnen des halb vector<complexGrid>
	opt.delta_t[0]=0.2;
        opt.delta_t[1]=0.4;

	opt.N=1598029824.0;  //normed for 512*512 N=64*50000 128
	
	opt.grid[0] = 128;
	opt.grid[1] = 128;
	opt.grid[2] = 1;
	opt.U = 0.000025;
	opt.u_12 = 0;		//interaction constants relative to U
	opt.u_22 = 1.0;
	
	opt.omega_l = 0.0; //0.0000025; //longitudinal trap frequenzy, absolute and squared
	
	
	opt.klength[0] = 2.0;
	opt.klength[1] = 2.0;
	opt.klength[2] = 2.0;
			

	snapshot_times.resize(2);
       snapshot_times[0] =100;
       snapshot_times[1] =500;
//snapshot_times[2] =5000;
//snapshot_times[3] =7000;
       /*snapshot_times[4] =10000;
       snapshot_times[5] =12000;
       snapshot_times[6] =17000;
       snapshot_times[7] =20000;
       snapshot_times[8] =30000;
       snapshot_times[9] =50000;
       snapshot_times[10]=70000;
       snapshot_times[11]=100000;
       snapshot_times[12]=150000;
       snapshot_times[13]=180000;
       snapshot_times[14]=220000;
       snapshot_times[15]=250000;
       snapshot_times[16]=300000;
       /*snapshot_times[17]=450000;
       snapshot_times[18]=500000;
       snapshot_times[19]=550000;
       snapshot_times[20]=600000;
       snapshot_times[21]=600000;
       snapshot_times[22]=600000;
       snapshot_times[23]=650000;
       snapshot_times[24]=700000;
       snapshot_times[25]=750000;
       snapshot_times[26]=800000;
       snapshot_times[27]=840000;
       snapshot_times[28]=870000;
       snapshot_times[29]=900000;
       snapshot_times[30]=930000;
       snapshot_times[31]=960000;
       snapshot_times[32]=1000000;
       snapshot_times[33]=1100000;
       snapshot_times[34]=1200000;
       snapshot_times[35]=1300000;
       snapshot_times[36]=1400000;
       snapshot_times[37]=1500000;
       snapshot_times[38]=1600000;
       snapshot_times[39]=1700000;
       snapshot_times[40]=1800000;
       snapshot_times[41]=1900000;*/
	
}




void plot3(ComplexGrid *xx, ComplexGrid *yy,const string &dirname, const PathOptions& opt,const vector<ComplexGrid> &grid1,double time)
{
	
 
	stringstream d;
	d << dirname << "/" << "Neu" << "/";
	string dir = d.str();
	system((string("mkdir -p ") + dir).c_str());

double *meineZahl = new double;
   *meineZahl = time;

string s;
ostringstream outStream;
outStream << *meineZahl;
s = outStream.str();


ofstream plotfile;

 plotfile.open((dir + (s+string("time.dat"))).c_str(), ios::out | ios::trunc);

for(int x = 0; x < opt.grid[0]; x++)
	{
		for (int y = 0; y < opt.grid[1]; y++)
		{
			for (int z = 0; z < opt.grid[2]; z++)
			{
			   plotfile<< x <<"\t";
			   plotfile << y <<"\t" ;
			   plotfile <<norm(grid1[0](x,y,z))<<"\t";
			   plotfile<<imag(xx->at(x,y,z))<<"\t" ;
			   plotfile<<imag(yy->at(x,y,z))<<"\t" ;
		           plotfile << endl;
			}
				
		}
	}
       
   plotfile.close();

   delete meineZahl;      
}

void plot(const string &dirname, const PathOptions& opt, vector<double> &snapshot_times,vector <Bh3EvaluationFD::Averages> &meansFD, int hh, int rV, int cV, int Q)
{  
  

	stringstream d;
	d << dirname << "/" << "spectrum" << "/";
	string dir = d.str();
	system((string("mkdir -p ") + dir).c_str());





ofstream plotfileinfo;	

	plotfileinfo.open((dir + string("info.dat")).c_str(), ios::out | ios::trunc);
	
	  plotfileinfo << "Runs:" << "\t" <<  hh << endl;
	  plotfileinfo << "Numbers of Vortices:"<< "\t" << rV <<"\t" << "x"<< "\t" << cV << "\t" << "=" << "\t" << rV*cV << endl;
          plotfileinfo << "Quatization:" << "\t"<< Q << endl;
          plotfileinfo << "Noise:" << "\t"<< "None" << endl;
	  plotfileinfo <<endl;	

     plotfileinfo << endl << endl;	
        

	
   
	plotfileinfo.close();


	
	
		ofstream plotfile;
	

	plotfile.open((dir + string("radial_avgs.dat")).c_str(), ios::out | ios::trunc);
	
	for (int i = 0; i < snapshot_times.size(); i++)
	{		
	  for (int r = 0; r < meansFD[i].ares_1.number.size(); r++)             
		{
		  plotfile << r <<"\t"<< meansFD[i].ares_1.k[r] <<"\t" << (meansFD[i].ares_1.number[r]+meansFD[i].ares_2.number[r])/2 <<"\t";
                  plotfile << (meansFD[i].ares_1.n_qkineticq[r]+meansFD[i].ares_2.n_qkineticq[r])/2 << "\t" << (meansFD[i].ares_1.n_ikinetick[r]+meansFD[i].ares_2.n_ikinetick[r])/2  << "\t";
	          plotfile << (meansFD[i].ares_1.n_ckinetick[r]+meansFD[i].ares_2.n_ckinetick[r])/2 << "\t";
                  plotfile << (meansFD[i].ares_1.n_qkineticq_wo_phase[r]+meansFD[i].ares_2.n_qkineticq_wo_phase[r])/2 << "\t";
                  plotfile << (meansFD[i].ares_1.n_ikinetick_wo_phase[r]+meansFD[i].ares_2.n_ikinetick_wo_phase[r])/2<<"\t";
	          plotfile << (meansFD[i].ares_1.n_ckinetick_wo_phase[r]+meansFD[i].ares_2.n_ckinetick_wo_phase[r])/2<<"\t";
	          plotfile << (meansFD[i].ares_1. ikinetick[r]+meansFD[i].ares_2. ikinetick[r])/2 << "\t" << (meansFD[i].ares_1.ckinetick[r]+meansFD[i].ares_2.ckinetick[r])/2 <<  "\t";
	          plotfile << (meansFD[i].ares_1.kinetick[r]+meansFD[i].ares_2.kinetick[r])/2 <<  "\t" << (meansFD[i].ares_1.pressure[r]+meansFD[i].ares_2.pressure[r])/2 <<  "\t";
	          plotfile << (meansFD[i].ares_1.ikinetick_wo_phase[r]+meansFD[i].ares_2.ikinetick_wo_phase[r])/2 <<  "\t" ;
		  plotfile << (meansFD[i].ares_1.ckinetick_wo_phase[r]+meansFD[i].ares_2.ckinetick_wo_phase[r])/2 << "\t";
                  plotfile << (meansFD[i].ares_1.kinetick_wo_phase[r]+meansFD[i].ares_2.kinetick_wo_phase[r])/2 << "\t" ;
	          plotfile << (meansFD[i].ares_1.pressure_wo_phase[r]+meansFD[i].ares_1.pressure_wo_phase[r])/2 <<"\t";
		  plotfile << endl;
		}
		plotfile << endl << endl;
	}
	
	plotfile.close();






		ofstream plotfile2;
		
   plotfile2.open((dir + string("g2_av.dat")).c_str(), ios::out | ios::trunc);
	
	       for (int i = 0; i < snapshot_times.size(); i++)
	{	    
	              
	 for (int l = 0; l < meansFD[i].ares_1.g2_av.size(); l++) 
	   {	
	     plotfile2 << meansFD[i].ares_1.r[l] << "\t";

             plotfile2 << (((meansFD[i].ares_1.g2_av[l]*sqrt(2)*4)/((meansFD[i].ares_1.nm1/(opt.grid[0]*opt.grid[1]*opt.grid[2]))*(meansFD[i].ares_1.n1/(opt.grid[0]*opt.grid[1]*opt.grid[2]))*(opt.                         grid[0]*opt.grid[1]*opt.grid[2])))+(((meansFD[i].ares_2.g2_av[l]*sqrt(2)*4)/((meansFD[i].ares_2.nm1/(opt.grid[0]*opt.grid[1]*opt.grid[2]))*(meansFD[i].ares_2.n1/(opt.grid[0]*opt.grid[1]*opt.grid[2]))*(opt.grid[0]*opt.grid[1]*opt.grid[2])))))/2 << "\t";

             plotfile2 << (((meansFD[i].ares_1.g2_vv[l]*sqrt(2)*4)/((meansFD[i].ares_1.nm1/(opt.grid[0]*opt.grid[1]*opt.grid[2]))*(meansFD[i].ares_1.n1/(opt.grid[0]*opt.grid[1]*opt.grid[2]))*(opt.                         grid[0]*opt.grid[1]*opt.grid[2])))+(((meansFD[i].ares_2.g2_vv[l]*sqrt(2)*4)/((meansFD[i].ares_2.nm1/(opt.grid[0]*opt.grid[1]*opt.grid[2]))*(meansFD[i].ares_2.n1/(opt.grid[0]*opt.grid[1]*opt.grid[2]))*(opt.grid[0]*opt.grid[1]*opt.grid[2])))))/2 << "\t";


	     plotfile2 << endl; 
	   }
         plotfile2 << endl << endl;
	}
	
	plotfile2.close();





	ofstream plotfile3;
	
	plotfile3.open((dir + string("Vortex_evo.dat")).c_str(), ios::out | ios::trunc);
	
	for (int i = 0; i < snapshot_times.size(); i++)
	{	
            plotfile3 <<  snapshot_times[i] <<"\t";	
	    plotfile3 <<  (meansFD[i].ares_1.nm1+meansFD[i].ares_2.nm1)/2  << "\t"<< (meansFD[i].ares_1.n1+meansFD[i].ares_2.n1)/2<< "\t" ;
            plotfile3 << (meansFD[i].ares_1.num_vortices+meansFD[i].ares_2.num_vortices)/2 << "\t";
	    plotfile3 << (meansFD[i].ares_1.pair_distance_all+meansFD[i].ares_2.pair_distance_all)/2  << "\t";
            plotfile3 <<(meansFD[i].ares_1.pair_distance_nonzero+meansFD[i].ares_2.pair_distance_nonzero)/2 << "\t" ;
	    plotfile3 << (((meansFD[i].ares_1.nm1 + meansFD[i].ares_1.n1)/(opt.grid[0]*opt.grid[1]*opt.grid[2]))+((meansFD[i].ares_2.nm1 + meansFD[i].ares_2.n1)/(opt.grid[0]*opt.grid[1]*opt.grid[2])))/2<< "\t";
            plotfile3 << endl;
	}
	
	plotfile3.close();






	
	ofstream plotfile4;	

	plotfile4.open((dir + string("Particle_evo.dat")).c_str(), ios::out | ios::trunc);
	
	for (int i = 0; i < snapshot_times.size(); i++)
	{		
	  plotfile4 << meansFD[i].ares_1.time << "\t" << (meansFD[i].ares_1.particle_count+meansFD[i].ares_2.particle_count)/2 << "\t";
          plotfile4 << ((meansFD[i].ares_1.particle_count/opt.N)+(meansFD[i].ares_2.particle_count/opt.N))/2 << endl;
	}
	
	plotfile4.close();






	ofstream plotfile5;

	plotfile5.open((dir + string("current.dat")).c_str(), ios::out | ios::trunc);
	
	for (int i = 0; i < snapshot_times.size(); i++)
	{		
	  for (int r = 0; r < meansFD[i].ares_1.number.size(); r++)             ///////////////////von FD klasse/////////////////////////Tauchen hier die radialen Kreise auf//////////////////////////////////
		{
		  plotfile5 << meansFD[i].ares_1.k[r] << "\t"<< (meansFD[i].ares_1.k_current[r]+meansFD[i].ares_2.k_current[r])/2<< "\t";
                  plotfile5 << (meansFD[i].ares_1.E_current[r]+meansFD[i].ares_2.E_current[r])/2 << "\t";
		  plotfile5 << endl;
		}
		plotfile5 << endl << endl;
	}
	
	plotfile5.close();






ofstream plotfile6;

plotfile6.open((dir + string("g1.dat")).c_str(), ios::out | ios::trunc);
	
	       for (int i = 0; i < snapshot_times.size(); i++)
	{	    
	              
	 for (int l = 0; l < meansFD[i].ares_1.g1.size(); l++) 
	   {	
	     plotfile6 << meansFD[i].ares_1.g1_r[l]  << "\t" << (meansFD[i].ares_1.g1r[l]+meansFD[i].ares_2.g1r[l])/2<< "\t";
	     plotfile6 << endl; 
	   }
         plotfile6 << endl << endl;
	}
	
	plotfile6.close();



ofstream plotfile7;





plotfile7.open((dir + string("energies_int.dat")).c_str(), ios::out | ios::trunc);
	
	       for (int i = 0; i < snapshot_times.size(); i++)
	{	    
	              	
	  plotfile7 << meansFD[i].ares_1.time << "\t";
          plotfile7 << (meansFD[i].ares_1.pressure_int+meansFD[i].ares_2.pressure_int)/2 << "\t" << (meansFD[i].ares_1.interaction_int+meansFD[i].ares_2.interaction_int)/2<< "\t";
          plotfile7 << (meansFD[i].ares_1.ckinetic_int+meansFD[i].ares_2.ckinetic_int)/2<< "\t" << (meansFD[i].ares_1.ikinetic_int+meansFD[i].ares_2.ikinetic_int)/2 << "\t" ;
          plotfile7 << (meansFD[i].ares_1.Ekin+meansFD[i].ares_2.Ekin)/2 << "\t";
	     plotfile7 << endl; 
	   
	}
	
	plotfile7.close();







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
        fs <<  "import matplotlib.colors " << endl;
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


			for(int i=0;i<snapshot_times.size(); i++)
  {

	ofstream fs;
	fs.open((dir+string("Spectrum2.py")).c_str(), ios_base::trunc | ios_base::out);


	  

        fs << "#!/usr/bin/python" << endl;
	fs << "# -*- coding: utf-8 -*-" << endl;

	fs << "from matplotlib import use" << endl;
	fs << "use('Agg')" << endl;
	fs << "from matplotlib import rc" << endl;
        fs << "import scipy, pylab" << endl;
        fs << "import scipy.optimize" << endl;
	fs << "import sys" << endl;
	fs << "import pylab as p "<< endl;
	fs << "import numpy as np" << endl;
	fs << "import matplotlib.pyplot as plt" << endl;
        fs << "from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, MultipleLocator" << endl;
        fs << "from matplotlib.ticker import FixedFormatter" << endl;
	fs << "from mpl_toolkits.axes_grid1 import make_axes_locatable" << endl;
	fs << "import math as m" << endl;
	
	fs << "def main():" << endl;
        
	
	fs << "\tpath = '" << dirname << "g2_av" << "_Path_"  << snapshot_times[i] <<"time"<<".png'" << endl;	
	fs << "\tfig = plt.figure(figsize=(8.3,5.7), dpi=100)" << endl;	
	fs << "\tdata = np.loadtxt('"<< (dir + string("g2_av.dat")).c_str() << "', delimiter='\\t', unpack=True, usecols = (0,1))" << endl;
	fs << "\ttime=" << snapshot_times[i] << endl;
        fs << "\tdim ="  << meansFD[i].ares_1.g2_av.size() << endl;
	fs << "\tsnapshot=" << i << endl;

        fs << "\tstart2=" << ((meansFD[i].ares_1.g2_av.size())*i)+1 << endl;
	fs << "\tend2=" << (meansFD[i].ares_1.g2_av.size())*(i+1) << endl;

        fs << "\tstart1=" << ((meansFD[i].ares_1.r.size())*i)+1 << endl;
	fs << "\tend1=" << (meansFD[i].ares_1.r.size())*(i+1) << endl;


	


	

	fs << "\tfig.subplots_adjust(hspace = 0.10, wspace = 0.40, right  = 0.85, left=0.15, bottom = 0.10, top = 0.90)" << endl;
        fs << "\tplt.xlim(0,200)" << endl;
	//fs << "\tplt.ylim(-1,10)" << endl;
     
          fs << "\tn1=data[0]" << endl;
	  fs << "\tn2=data[1]" << endl;




	  fs << "\tn1neu=n1[start1:end1]" << endl;
	  fs << "\tn2neu=n2[start2:end2]" << endl;
	  
		
	 	
	  fs << "\tax = fig.add_subplot(111)" << endl;	
	  fs << "\tim=plt.plot(n1neu, n2neu, 'ro')" << endl;

	  

	  fs << "\tplt.setp(im, 'markersize', 2)" << endl; 
 

         fs << "\tax.set_title('Time: $"<< snapshot_times[i]<<","<<"Particles:" <<opt.N<<"$ ') " << endl;
	 fs << "\tfig.savefig(path, dpi=100)" << endl;
         fs << "\tsys.exit(3)" << endl;
         fs << "if __name__=='__main__':" << endl;

         fs << "\tmain()" << endl;

  

         fs.close();
	
	system((string("python ") + dir + string("Spectrum2.py")).c_str());
	
        
     
	
	}*/

	
				
  	
	


	}







  
























