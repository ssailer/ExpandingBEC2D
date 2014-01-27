#include <sys/stat.h>
#include <sys/types.h>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <averageclass.h>
#include <bh3cudapropagator.h>
#include <complexgrid.h>
#include <bh3defaultgrid.h>
#include <bh3binaryfile.h>
#include <gauss_random.h>
#include <wrapped_cuda_functions.h>

#define FRAME_DURATION 20
#define MAXTIME 40000
#define NUM_SINGLE_PLOTS 9999999
#define FIELD_WIDTH 6

using namespace std;

void init_bh3(int argc, char** argv, PathOptions &opt);
void plot(const string &dirname, int snapshot);

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv) 
{
	PathOptions opt;
	string rm = "rm ";
	
	init_bh3(argc, argv, opt);
	ComplexGrid::set_fft_planning_rigorosity(FFTW_MEASURE);
	init_random();
	stringstream dstr;
	dstr << "dpgmovie_TS"; //<< opt.timestepsize
        //<< "_G" << opt.grid[0] << "_" << opt.grid[1] << "_" << opt.grid[2] << "_" << opt.grid[3]
        //<< "_N" << opt.N
        //<< "_U" << opt.U
        //<< "_KL" << opt.klength[0] << "_" << opt.klength[1] << "_" << opt.klength[2];
	string dirname = dstr.str();
	initialize_binary_dir(dirname, opt);
	mkdir((dirname + "temp").c_str(), 0755);
		
	if(init_cuda_device(0))
	{
		unsigned int t;
		time_t timer = time(NULL);
		
		ComplexGrid *start;
		Bh3CudaPropagator *cp1, *cp2;
				
            //start = create_Vortex_start_Grid2(opt, 16,4,4,4);
            //start = create_Default_Start_Grid(opt,1);
        start = new ComplexGrid (opt.grid[0], opt.grid[1], opt.grid[2], opt.grid[3]);
        

/*        opt.timestepsize = 0.015;
        cp_imag = new Bh3CudaPropagator(opt, *start, Bh3CudaPropagator::imag);
        cp_imag -> propagateToTime(opt.timestepsize*2000.);
        cp_imag -> renoise();
        
        *start = cp_imag -> getRGrid()[0];
        
        delete cp_imag;
        

        opt.timestepsize = 0.2;*/
        
        
/*        cp1 = new Bh3CudaPropagator(opt, *start, Bh3CudaPropagator::drivdiss);
        cp1->propagateToTime(15000.);
        *start = cp1->getRGrid()[0];

        delete cp1;
        
        opt.g[0] = 1./100.;*/
        cp2 = new Bh3CudaPropagator(opt, *start, Bh3CudaPropagator::drivdiss);

		delete start;
		
		for(int j = 0; j < MAXTIME/FRAME_DURATION; j++)
		{
            double time = j*FRAME_DURATION;
						 
            cp2->propagateToTime(time);

            stringstream filestr;
                                                  
            filestr << dirname << "temp/Bh3RunSnapshot" << ".h5";
            string filename = filestr.str();                                                
            Bh3BinaryFile *bf = new Bh3BinaryFile(filename, opt, Bh3BinaryFile::out);
            bf->append_snapshot(time, cp2->getRGrid());

            delete bf;
						 
/*							vector<ComplexGrid> r;
							
							PathOptions options;
							double time;
                            Bh3BinaryFile *temp_file = new Bh3BinaryFile(filename, options, Bh3BinaryFile::in);
							options = temp_file->get_options();
							temp_file->get_snapshot(time, r1);
							 
							delete temp_file;
							system((rm + filename).c_str());

                                /*Bh3EvaluationFD ev(options, false);
                                ev.setData(r1,r2,Bh3EvaluationFD::RSpace);
                                ev.setTime(time);
                                ev.calc_radial_averages();*/
                // output as png
            plot(dirname, j);
        }
        
        delete cp2;
		
        cout << "Path CUDA took " << time(NULL) - timer << "seconds" << endl;
		
    }         // if init_cuda_device
    else
    {
        cout << "Could not initialize Cuda device!" << endl;
    }

	stringstream rm_command;
//    stringstream compression;
    
//    compression << "ffmpeg -r 10 -b:v 1800 -i movie%06d.png movie.mp4";
    rm_command << "rm -r " << dirname << "temp";
	
//    system(compression.str().c_str());

    system(rm_command.str().c_str());
		
	return 0;
}

void init_bh3(int argc, char** argv, PathOptions &opt)
{
	// Parameter setzen
	opt.timestepsize = 0.2;
	opt.delta_t.resize(0);

	opt.N = 3.2e9/4.;  
	
    opt.grid[0] = 1;
    opt.grid[1] = 512;
	opt.grid[2] = 512;
	opt.grid[3] = 1;
	opt.U = 3.e-5;

    opt.g.resize(1);
    opt.g[0] = 1./4.;
    	
	opt.klength[0] = 2.0;
	opt.klength[1] = 2.0;
	opt.klength[2] = 2.0;
	
}

void plot(const string &dirname, int snapshot)
{
	stringstream d;
	d << dirname << "temp/Bh3RunSnapshot" << ".h5";
	string dir = d.str();
		
		// Daten plotten
	ofstream fs;
	fs.open((dirname+string("temp/movie.py")).c_str(), ios_base::trunc | ios_base::out);
	
	
	fs << "#!/usr/bin/python" << endl;
	fs << "# -*- coding: utf-8 -*-" << endl;

	fs << "from matplotlib import use" << endl;
	fs << "use('Agg')" << endl;
	fs << "from matplotlib import rc" << endl;
    fs << "rc('font',**{'family':'serif','serif':['Computer Modern Roman'], 'size':14.0})" << endl;
	fs << "rc('text', usetex=True)" << endl;

	fs << "import sys" << endl;
	fs << "import numpy as np" << endl;
	fs << "import matplotlib.pyplot as plt" << endl;
    fs << "import bh3pynary as b" << endl;
    fs << "from mpl_toolkits.axes_grid1 import make_axes_locatable" << endl;
	fs << "import math as m" << endl;
	
	fs << "def main():" << endl;
	
	fs << "\tpath = '" << dirname << "movie" << setfill('0') << setw(FIELD_WIDTH) << ((snapshot) % NUM_SINGLE_PLOTS) << ".png'" << endl;
	
	fs << "\tfig = plt.figure(figsize=(8.,4.), dpi=100)" << endl;
	fs << "\tfig.subplots_adjust(hspace = 0.10, wspace = 0.40, right  = 0.85, left=0.15, bottom = 0.10, top = 0.90)" << endl;
	
    fs << "\tdat = b.bh3pynary('" << dir << "')" << endl;

	fs << "\tpsi, time = dat.get_cgrid(0)" << endl;        
    
        //magic
    fs << "\ta = np.mean(np.absolute(psi[0])**2)" << endl;
    fs << "\tprint(a)" << endl;
    
	fs << "\tax = fig.add_subplot(1,2,1)" << endl;
	fs << "\tim = ax.imshow(np.absolute(psi[0])**2, interpolation='gaussian', origin='lower left', extent=[0,dat.dim[0],0,dat.dim[1]], vmin=0, vmax=1.2*dat.N/(dat.dim[0]*dat.dim[1]), label = 's', cmap = 'jet')" << endl;
	fs << "\tdivider = make_axes_locatable(ax)" << endl;
	fs << "\tcax = divider.append_axes('right', size='8%', pad=0.08)" << endl;

	fs << "\tcbar = fig.colorbar(im, cax=cax)" << endl;
	fs << "\tcbar.ax.set_ylabel(r'Density', rotation='vertical', fontsize = 16)" << endl;		

	fs << "\tax.set_ylabel(r'Grid Position $y$', fontsize= 16)" << endl;
	fs << "\tax.yaxis.set_label_coords(-0.25, 0.5)" << endl;
	fs << "\tax.set_xlabel(r'Grid position $x$', fontsize = 16)" << endl;
	fs << "\tax.xaxis.set_label_coords(1.3, -0.17)" << endl;
	
	fs << "\tax.set_title(r'Time: %06d' %int(time), size = 18, position = (1.35,1.15))" << endl; 
	
	fs << "\tax = fig.add_subplot(1,2,2)" << endl;
	fs << "\tim = ax.imshow(np.angle(psi[0]),interpolation='gaussian', origin='lower left', extent=[0,dat.dim[0],0,dat.dim[1]], label = 's', cmap = 'spectral')" << endl;
	fs << "\tdivider = make_axes_locatable(ax)" << endl;
	fs << "\tcax = divider.append_axes('right', size='8%', pad=0.08)" << endl;

	fs << "\tcbar = fig.colorbar(im, cax=cax)" << endl;
	fs << "\tcbar.ax.set_ylabel(r'Phase', rotation='vertical', fontsize = 16)" << endl;

	fs << "\tax.set_yticklabels([])" << endl; 

	fs << "\tfig.savefig(path, dpi=100)" << endl;

	fs << "\tsys.exit(3)" << endl;

	fs << "if __name__=='__main__':" << endl;
	fs << "\tmain()" << endl;

	fs.close();
	
	system((string("python ") + dirname + string("temp/movie.py")).c_str());
	
	//Dateien loeschen!
//	system((string("rm -r ") + dirname + string("temp")).c_str());
}
