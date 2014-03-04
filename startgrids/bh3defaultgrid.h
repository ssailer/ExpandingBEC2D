#ifndef BH3DEFAULTGRID_H__
#define BH3DEFAULTGRID_H__

#include <bh3binaryfile.h>
#include <exp_RK4_tools.h>
int mypow2(int x, int y); // Computes x^y
ComplexGrid *add_vortex_to_grid(ComplexGrid* &g, Options &opt,int sigma[2]);
ComplexGrid *create_Vortex_start_Grid3(ComplexGrid* &g,const Options &opt,int Vortexnumber, int rows_y, int columns_x, int Q) ;
ComplexGrid *create_Vortex_start_Grid2(const PathOptions &opt,int Vortexnumber, int rows_y, int columns_x, int Q) ;
ComplexGrid *create_Vortex_start_Grid(const PathOptions &opt,int Vortexnumber);// Vortex Startgitter
ComplexGrid *create_no_noise_Start_Grid(const PathOptions &opt, int d = -1);
ComplexGrid *create_Default_Start_Grid(const PathOptions &opt, int d = -1);
ComplexGrid *create_noise_Start_Grid(const PathOptions &opt);
ComplexGrid *create_test2_Start_Grid(const PathOptions &opt, int d = -1);
ComplexGrid *create_oszi_test_Start_Grid(const PathOptions &opt, int d = -1); //n =2 harm oszi eigenstate
ComplexGrid *create_vortex(const PathOptions &opt, int d = -1); //bright soliton solution, so only for U < 0!
ComplexGrid *create_Energy_Start_Grid(const PathOptions &opt, int d);
ComplexGrid *create_inverse_Start_Grid(const ComplexGrid &start);

#endif // BH3DEFAULTGRID_H__
