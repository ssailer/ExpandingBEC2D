#define EIGEN_VECTORIZE
#define EIGEN_NO_DEBUG

#include <EXP2D_tools.h>

// void saveDataToHDF5(ComplexGrid* &g, Options &opt){ 
// 	PathOptions pathopt;
	
// 	optToPath(opt,pathopt);
	
// 	double time = opt.n_it_ITP1 + opt.n_it_ITP2;
	
// 	Bh3BinaryFile *bf = new Bh3BinaryFile(opt.workingfile, pathopt, Bh3BinaryFile::out);
	
// 	vector<ComplexGrid> vectork(1);
// 	vectork[0] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
// 	for(int i = 0; i < opt.grid[1];i++){for(int j = 0; j < opt.grid[2]; j++){ vectork.at(0).at(0,i,j,0) = g->at(0,i,j,0) ;}}
// 	bf->append_snapshot(time, vectork);
	
// 	delete bf;
// }

// void readDataFromHDF5(ComplexGrid* &g,Options &opt){
// 	try{
// 	PathOptions pathopt;

//   	double time = opt.n_it_ITP1 + opt.n_it_ITP2;

// 	Bh3BinaryFile *bf = new Bh3BinaryFile(opt.workingfile, pathopt, Bh3BinaryFile::in);

// 	pathopt = bf->get_options();

// 	pathToOpt(pathopt,opt);

// 	vector<ComplexGrid> vectork(1);
// 	vectork[0] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

// 	bf->get_snapshot(time, vectork,0);

// 	for(int i = 0; i < opt.grid[1];i++){for(int j = 0; j < opt.grid[2]; j++){ g->at(0,i,j,0) = vectork.at(0).at(0,i,j,0) ;}}
	
// 	delete bf;
// 	}

// 	catch(std::exception& e) 
// 	{ 
//   		std::cerr << "Reading from HDF5 File failed, whaaaat?: " 
//         		  << e.what() << ", application will now exit" << std::endl; 
// 	} 
// }

void noiseTheGrid(ComplexGrid &g){
   GaussRandom r (get_seed());
   double rvalue;
   for(int i = 0;i < g.width();i++){
    for(int j = 0; j < g.height();j++){
        rvalue = real(g(0,i,j,0)) * 0.1;
        g(0,i,j,0) += r.gauss_random(0.0,rvalue);
    }
   }
}

void pathToOpt(PathOptions &pathopt,Options &opt){
	opt.RTE_step	          = pathopt.timestepsize;
	opt.N                     = pathopt.N;
	opt.grid[0]               = pathopt.grid[0];
	opt.grid[1]               = pathopt.grid[1];
	opt.grid[2]               = pathopt.grid[2];
	opt.grid[3]               = pathopt.grid[3];
	opt.g                     = pathopt.U;
	opt.klength[0]            = pathopt.klength[0];
	opt.klength[1]            = pathopt.klength[1];
	opt.klength[2]            = pathopt.klength[2];
	opt.min_x                 = pathopt.g[0];
	opt.min_y                 = pathopt.g[1];
	opt.exp_factor    		  = complex<double>(pathopt.g[2],0.0);
	opt.omega_x       		  = complex<double>(pathopt.g[3],0.0);
	opt.omega_y       		  = complex<double>(pathopt.g[4],0.0);
	opt.dispersion_x  		  = complex<double>(pathopt.g[5],0.0);
	opt.dispersion_y  		  = complex<double>(pathopt.g[6],0.0);
	opt.t_abs         		  = complex<double>(pathopt.g[7],0.0);
	opt.samplesize    		  = (int)pathopt.g[8];
	opt.runmode 			  = std::to_string(pathopt.g[9] - 10000);
	opt.n_it_RTE              = (int)pathopt.g[10];
	opt.scale_factor          = pathopt.g[11];
};

void optToPath(Options &opt,PathOptions &pathopt){
	pathopt.timestepsize 		= opt.RTE_step;	
	pathopt.N            		= opt.N;
	pathopt.grid[0]      		= opt.grid[0];
	pathopt.grid[1]      		= opt.grid[1];
	pathopt.grid[2]      		= opt.grid[2];
	pathopt.grid[3]      		= opt.grid[3];
	pathopt.U            		= opt.g;
	pathopt.klength[0]   		= opt.klength[0];
	pathopt.klength[1]   		= opt.klength[1];
	pathopt.klength[2]   		= opt.klength[2];
	pathopt.delta_t.resize(0);
	pathopt.g.resize(12);
	pathopt.g[0]         		= opt.min_x;
	pathopt.g[1]         		= opt.min_y;
	pathopt.g[2]         		= opt.exp_factor.real();
	pathopt.g[3]         		= opt.omega_x.real();
	pathopt.g[4]         		= opt.omega_y.real();
	pathopt.g[5]         		= opt.dispersion_x.real();
	pathopt.g[6]         		= opt.dispersion_y.real();
	pathopt.g[7]         		= opt.t_abs.real();
	pathopt.g[8]         		= (double)opt.samplesize;
	pathopt.g[9]         		= 10000 + atof(opt.runmode.c_str());
	pathopt.g[10]        		= opt.n_it_RTE;
	pathopt.g[11]        		= opt.scale_factor.real();
};
