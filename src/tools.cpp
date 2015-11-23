#define EIGEN_VECTORIZE
#define EIGEN_NO_DEBUG

#include <tools.h>

void toDimensionlessUnits(Options &opt){
    if(!opt.isDimensionless){
        const double m = 86.9091835 *  1.660538921 * 1.0e-27;
        const double hbar = 1.0545718 * 1.0e-22;
        const double epsilon = 0.1;
        double a0 = sqrt(hbar / ( m * opt.omega_x.real() * 2.0 * M_PI));
        opt.Ag = a0 / (epsilon * epsilon);
        // opt.Ag = 2 * opt.min_x / opt.grid[1];

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

void toPhysicalUnits(Options &opt){
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

