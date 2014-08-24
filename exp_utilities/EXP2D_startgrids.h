#ifndef EXP2D_STARTGRIDS_H__
#define EXP2D_STARTGRIDS_H__

#include <EXP2D_MatrixData.h>
#include <EXP2D_tools.h>
#include <complexgrid.h>

// int mypow2(int x, int y); // Computes x^y

// inline double vortex(int b, int y, int a, int x) //Vortex with phase [0,2*pi)          
// {
//         if(atan2(b-y,a-x)<0){ return 2*M_PI+atan2(b-y,a-x); } //atan2 is defined from [-pi,pi) so it needs to be changed to [0,2*pi)
//     else{ return atan2(b-y,a-x); }        
// }

void setGridToDoubleGaussian(MatrixData* &data, Options opt)
{
    double sigma_x = opt.min_x/4;
    double sigma_y = opt.min_y/4;
    double h_x = 2.*opt.min_x/opt.grid[1];
    double h_y = 2.*opt.min_y/opt.grid[2];
    vector<double> x(opt.grid[1]);
    vector<double> y(opt.grid[2]);

    for(int i=0;i<opt.grid[1];i++){x[i]=-opt.min_x+i*h_x;}
    for(int j=0;j<opt.grid[2];j++){y[j]=-opt.min_y+j*h_y;}

    complex<double> value;
    for(int i=0; i < opt.grid[1]; i++){
        for(int j=0; j < opt.grid[2]; j++){
            value = complex<double>((opt.N/(4 * opt.min_x * opt.min_y)) * exp( -((x[i] - opt.min_x/2) * (x[i] - opt.min_x/2))/(2.*sigma_x*sigma_x) - (y[j] * y[j])/(2.*sigma_y*sigma_y) ), 0.0 );
            value += complex<double>((opt.N/(4 * opt.min_x * opt.min_y)) * exp( -((x[i] + opt.min_x/2) * (x[i] + opt.min_x/2))/(2.*sigma_x*sigma_x) - (y[j] * y[j])/(2.*sigma_y*sigma_y) ), 0.0 );
            for(int m = 0; m < data->wavefunction.size(); m++)
                data->wavefunction[m](i,j) = value;
        }
    }
};

void setGridToGaussian(MatrixData* &data, Options opt)
{
    double sigma_x = opt.min_x/4;
    double sigma_y = opt.min_y/4;
    double h_x = 2.*opt.min_x/opt.grid[1];
    double h_y = 2.*opt.min_y/opt.grid[2];
    vector<double> x(opt.grid[1]);
    vector<double> y(opt.grid[2]);

    for(int i=0;i<opt.grid[1];i++){x[i]=-opt.min_x+i*h_x;}
    for(int j=0;j<opt.grid[2];j++){y[j]=-opt.min_y+j*h_y;}

    complex<double> value;
    for(int i=0; i < opt.grid[1]; i++){
        for(int j=0; j < opt.grid[2]; j++){
            value = complex<double>((opt.N/(4 * opt.min_x * opt.min_y )) * exp( -(x[i] * x[i])/(2.*sigma_x*sigma_x) - (y[j] * y[j])/(2.*sigma_y*sigma_y) ), 0.0 );
                data->wavefunction[0](i,j) = value;
        }
    }
};

int addVortices(MatrixData* &g, Options opt){

int x_jump = opt.grid[1] / 8;
int y_jump = opt.grid[2] / 8;
int windingnumber = 1;

ComplexGrid grid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);

vector<Coordinate<int32_t>> c;

for(int y = y_jump; y < opt.grid[2]; y += y_jump*2){
    for(int x = x_jump; x < opt.grid[1]; x += x_jump){
        if(abs2(g->wavefunction[0](x,y)) >= 2){
            c.push_back(grid.make_coord(x,y,0));
        }
    }
}
for(int y = y_jump*2; y < opt.grid[2]; y += y_jump*2){
    for(int x = x_jump/2; x < opt.grid[1]; x += x_jump){
        if(abs2(g->wavefunction[0](x,y)) >= 2){
            c.push_back(grid.make_coord(x,y,0));
        }
    }
}

for(int i = 0; i < c.size(); i++){
    for(int y = 0; y < opt.grid[2]; y++){
        for(int x = 0; x < opt.grid[1]; x++){   
            g->wavefunction[0](x,y) *= polar(1.0, (windingnumber /* * mypow2(-1,i+1)*/ )*vortex( y,c[i].y(),x,c[i].x() )) ;
        }
    }
    // g->wavefunction[0](c) complex<double>(0.0,0.0);
}
opt.vortexnumber += c.size() * windingnumber;
    return opt.vortexnumber;
}

#endif // EXP2D_STARTGRIDS_H__