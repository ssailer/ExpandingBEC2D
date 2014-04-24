#define EIGEN_VECTORIZE
#define EIGEN_NO_DEBUG

#include <EXP2D_tools.h>

void saveDataToHDF5(ComplexGrid* &g, Options &opt)
{ 

  PathOptions pathopt;

  	// useless to me
  pathopt.timestepsize = opt.RTE_step;
  for(int i = 0; i<3; i++){pathopt.klength[i] = 2.0;}
  pathopt.delta_t.resize(1);
  pathopt.delta_t[0] = 1.0;
  	// still useless to me

  pathopt.U = opt.g;
  pathopt.N = opt.N;
  for(int i = 0; i<4;i++){ pathopt.grid[i] = opt.grid[i]; }
  pathopt.g.resize(3);
  pathopt.g[0] = real(opt.omega_x);
  pathopt.g[1] = real(opt.omega_y);
  pathopt.g[2] = real(opt.t_abs);

  double time = opt.n_it_ITP1 + opt.n_it_ITP2;

  Bh3BinaryFile *bf = new Bh3BinaryFile(opt.workingfile, pathopt, Bh3BinaryFile::out);

	vector<ComplexGrid> vectork(1);
	vectork[0] = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

	for(int i = 0; i < opt.grid[1];i++){for(int j = 0; j < opt.grid[2]; j++){ vectork.at(0).at(0,i,j,0) = g->at(0,i,j,0) ;}}

	bf->append_snapshot(time, vectork);

  delete bf;
}

void readDataFromHDF5(ComplexGrid* &g,Options &opt)
{
	try{
	PathOptions pathopt;

  	double time = opt.n_it_ITP1 + opt.n_it_ITP2;

	Bh3BinaryFile *bf = new Bh3BinaryFile(opt.workingfile, pathopt, Bh3BinaryFile::in);

	pathopt = bf->get_options();

	opt.g = pathopt.U;
	opt.N = pathopt.N;
	for(int i = 0; i<4; i++){ pathopt.grid[i] = opt.grid[i]; }
	pathopt.g.resize(3);
	opt.omega_x = complex<double>(pathopt.g[0],0.0);
	opt.omega_y = complex<double>(pathopt.g[1],0.0);
	opt.t_abs   = complex<double>(pathopt.g[2],0.0);

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
