#ifndef BH3DEFAULTGRID_H__
#define BH3DEFAULTGRID_H__

#include <bh3binaryfile.h>
#include <exp_RK4_tools.h>
int mypow2(int x, int y); // Computes x^y
inline double vortex(int b, int y, int a, int x) //Vortex with phase [0,2*pi)          
{
        if(atan2(b-y,a-x)<0){ return 2*M_PI+atan2(b-y,a-x); } //atan2 is defined from [-pi,pi) so it needs to be changed to [0,2*pi)
    else{ return atan2(b-y,a-x); }        
}
ComplexGrid *set_grid_to_gaussian(ComplexGrid* &g, Options &opt, std::vector<double> &x, std::vector<double> & y, double & sigma_x, double & sigma_y);
ComplexGrid *add_vortex_to_grid(ComplexGrid* &g, Options &opt,int sigma[2]);
ComplexGrid *create_noise_Start_Grid(ComplexGrid* &g,const Options &opt);


ComplexGrid *create_Vortex_start_Grid3(ComplexGrid* &g,const Options &opt,int Vortexnumber, int rows_y, int columns_x, int Q) ;
ComplexGrid *create_Vortex_start_Grid2(const PathOptions &opt,int Vortexnumber, int rows_y, int columns_x, int Q) ;
ComplexGrid *create_Vortex_start_Grid(const PathOptions &opt,int Vortexnumber);// Vortex Startgitter
ComplexGrid *create_no_noise_Start_Grid(const PathOptions &opt, int d = -1);
ComplexGrid *create_Default_Start_Grid(const PathOptions &opt, int d = -1);
ComplexGrid *create_test2_Start_Grid(const PathOptions &opt, int d = -1);
ComplexGrid *create_oszi_test_Start_Grid(const PathOptions &opt, int d = -1); //n =2 harm oszi eigenstate
ComplexGrid *create_vortex(const PathOptions &opt, int d = -1); //bright soliton solution, so only for U < 0!
ComplexGrid *create_Energy_Start_Grid(const PathOptions &opt, int d);
ComplexGrid *create_inverse_Start_Grid(const ComplexGrid &start);

#endif // BH3DEFAULTGRID_H__
