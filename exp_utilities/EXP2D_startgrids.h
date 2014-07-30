#ifndef EXP2D_STARTGRIDS_H__
#define EXP2D_STARTGRIDS_H__

#include <EXP2D_MatrixData.h>
#include <EXP2D_tools.h>

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

#endif // EXP2D_STARTGRIDS_H__