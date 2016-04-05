#define EIGEN_VECTORIZE
#define EIGEN_NO_DEBUG

#include "tools.h"

void toDimensionlessUnits(Options &opt){
    if(!opt.isDimensionless){
        const double m = 86.9091835 *  1.660538921 * 1.0e-27;
        const double hbar = 1.0545718 * 1.0e-22;
        const double epsilon = 0.1;
        double a0 = sqrt(hbar / ( m * opt.omega_x.real() * 2.0 * M_PI));
        opt.Ag = a0 / (epsilon * epsilon);
        // opt.Ag = 2 * opt.min_x / opt.grid[1];

        opt.OmegaG = hbar / ( m * opt.Ag * opt.Ag);

        cout << currentTime() << " Stabilitycheck in physical units: " << endl;
        double h1 = opt.min_x * 2 / opt.grid[1];
        cout << "\t grid spacing = " << h1 << endl;
        cout << "\t coordinate maximum = " << opt.min_x << endl; 
        double k1 = 2 * h1 * h1 / (M_PI * M_PI);
        cout << "\t timestep < spectral resolution " << opt.RTE_step << " < " << k1 << " ";
        if(opt.RTE_step < k1) cout << "Yes." << endl; else cout << "No." << endl;
        cout << endl;
    
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

        cout << currentTime() << " Stabilitycheck in dimensionless units: " << endl;
        double h = opt.min_x * 2 / opt.grid[1];
        cout << "\t grid spacing = " << h << endl;
        cout << "\t coordinate maximum = " << opt.min_x << endl; 
        double k = 2 * h * h / (M_PI * M_PI);
        cout << "\t timestep < spectral resolution " << opt.RTE_step << " < " << k << " ";
        if(opt.RTE_step < k) cout << "Yes." << endl; else cout << "No." << endl;
        cout << endl;

        cout << currentTime() << " Set timestep to maximum appropriate value: " << k << endl << endl;
        opt.RTE_step = k;

    } else {
        cerr << " Trying to convert dimensionless Options to dimensionless Options. Check tools.h" << endl;
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
        cerr << " Trying to convert non-dimensionless Options to non-dimensionless Options. Check tools.h" << endl;
    }
}

