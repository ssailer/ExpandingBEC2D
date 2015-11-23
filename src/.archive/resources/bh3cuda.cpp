#include <vector>
#include <string>
#include <cstring>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>
#include <fstream>
#include <omp.h>

#include <bh3cudapropagator.h>
#include <complexgrid.h>
#include <bh3defaultgrid.h>
#include <bh3binaryfile.h>
#include <wrapped_cuda_functions.h>
#include <gauss_random.h>

#define PATHNUMBER 1		// Muss ein Vielfaches von NUM_THREADS sein, damit nicht Dateien ueberschrieben werden
#define NUM_THREADS 1

void init_bh3(int argc, char** argv, PathOptions &opt, vector<double> &snapshot_times);

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv) 
{ 
	PathOptions opt;
	vector<double> snapshot_times;
	string dirname = "bh3cuda";

	init_bh3(argc, argv, opt, snapshot_times);
	initialize_binary_dir(dirname, opt);
	init_random();
	ComplexGrid::set_fft_planning_rigorosity(FFTW_MEASURE);
	
	
	#pragma omp parallel /*for schedule(static,1)*/ num_threads(NUM_THREADS)
	//for(int thread=0; thread < NUM_THREADS; thread++)
	{
		int thread = omp_get_thread_num();
		if(init_cuda_device(0))
		{
				
			for(int p = 1; p <= PATHNUMBER/NUM_THREADS; p++)
			{
				unsigned int t;
			
				stringstream filename;
				filename << dirname << "Bh3" << "Path" << p+PATHNUMBER/NUM_THREADS*thread << ".bin";
				Bh3CPUPropagator *cp, *cp_imag;
			
				time_t timer = time(NULL);
				ComplexGrid *start =create_Vortex_start_Grid2(opt,2,2,1,2);

                opt.timestepsize = 0.015;
                cp_imag = new Bh3CPUPropagator(opt, *start, Bh3CPUPropagator::imag);
                cp_imag -> propagateToTime(opt.timestepsize*2000.);
                cp_imag -> renoise();
        
                *start = cp_imag -> getRGrid()[0];
        
                delete cp_imag;

                opt.timestepsize = 0.2;
                cp = new Bh3CPUPropagator(opt, *start, Bh3CPUPropagator::diss);

                delete start;

				if(!cp->start(snapshot_times, filename.str()))
				{
					cout << "Propagator for path " << p+PATHNUMBER/NUM_THREADS*thread << " returned false! File will be invalid!" << endl;
				}
				delete start;
				delete cp;
				cout << "Path Cuda took " << time(NULL) - timer << "seconds" << endl;

				cudaDeviceSynchronize();
				#pragma omp barrier
				
			}
			cudaThreadExit();
		}
	}
	
	return 0;
}

void init_bh3(int argc, char** argv, PathOptions &opt, vector<double> &snapshot_times)
{
	// Parameter setzen
    opt.timestepsize = 0.2;  // Setzt die Schrittweite zwischen zwei Propagationsschritten
	opt.delta_t.resize(0);
	

	opt.N = 3.2e9;  //normed for 512*512 N=64*50000 4*8*100000000.0

	opt.grid[0] = 1; //number of field components
    opt.grid[1] = 1024;
	opt.grid[2] = 1024;
	opt.grid[3] = 1;

	opt.U = 3e-5;        // 0.00003 gibt die Wechselwirkung zwischen zwischen zwei Teilchen gleicher Art an

    opt.g.resize(1);
    opt.g[0] = 1./30000.; //nonlinear dissipation strength in units of U
	
	opt.klength[0] = 2.0;// Einstellen der Konstante vor dem Sinus fuer die K-werte da im k wert eine Verzerrung drin ist
	opt.klength[1] = 2.0;
	opt.klength[2] = 2.0;
	
	snapshot_times.resize(2);
	snapshot_times[0] = 0;
	snapshot_times[1] = 1000;
	//snapshot_times[2] = 100;
	//snapshot_times[3] = 1000000;
	//snapshot_times[4] = 10000000;
	//snapshot_times[5] = 20000000;
	//snapshot_times[6] = 75000;
	//snapshot_times[7] = 100000;
	//snapshot_times[8] = 500000;
	//snapshot_times[8] = 500000;
	//snapshot_times[9] = 1000000;  
}		
