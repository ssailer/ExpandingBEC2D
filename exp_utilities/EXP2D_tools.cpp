#define EIGEN_VECTORIZE
#define EIGEN_NO_DEBUG

#include <EXP2D_tools.h>

void toDimensionless(Options &opt){
    if(!opt.isDimensionless){
        const double m = 87 * 1.66 * 1.0e-27;
        const double hbar = 1.054 * 1.0e-22;    
        opt.Ag = 2 * opt.min_x / opt.grid[1];
        opt.OmegaG = hbar / ( m * opt.Ag * opt.Ag);
    
        opt.min_x /= opt.Ag;
        opt.min_y /= opt.Ag;
        opt.ITP_step *= opt.OmegaG;
        opt.RTE_step *= opt.OmegaG;
        opt.omega_x *= 2.0 * M_PI / opt.OmegaG;
        opt.omega_y *= 2.0 * M_PI / opt.OmegaG;
        opt.omega_w *= 2.0 * M_PI / opt.OmegaG;
        opt.dispersion_x *= 2.0 * M_PI / opt.OmegaG;
        opt.dispersion_y *= 2.0 * M_PI / opt.OmegaG;
        opt.isDimensionless = true;
    } else {
        cerr << " Trying to convert dimensionless Options to dimensionless Options. Check EXP2D_tools.h" << endl;
    }

}

void fromDimensionless(Options &opt){
    if(opt.isDimensionless){
        opt.min_x *= opt.Ag;
        opt.min_y *= opt.Ag;
        opt.ITP_step /= opt.OmegaG;
        opt.RTE_step /= opt.OmegaG;
        opt.omega_x /= 2.0 * M_PI / opt.OmegaG;
        opt.omega_y /= 2.0 * M_PI / opt.OmegaG;
        opt.omega_w /= 2.0 * M_PI / opt.OmegaG;
        opt.dispersion_x /= 2.0 * M_PI / opt.OmegaG;
        opt.dispersion_y /= 2.0 * M_PI / opt.OmegaG;
        opt.isDimensionless = false;
    } else {
        cerr << " Trying to convert non-dimensionless Options to non-dimensionless Options. Check EXP2D_tools.h" << endl;
    }
}

