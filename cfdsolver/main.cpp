#include "../solverlib/solver.hpp"


int main(int argc, char* argv[]) {

    
    int Nx = 400, Ny = 200;

    vector<int> CFL_timesteps = {0, 250, 500, 750, 1000, 1500, 2000, 2500, 3000};
    Vector CFLs = {0.5, 1.0, 1.5, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0};

    double bigG = 1.4;
    double Mach = 20;
    double T = 226.15;
    double p = 900;
    double rho = 0.01388;
    double a = sqrt(bigG * p / rho);
    double u = Mach * a;
    double v = 0.0;

    Vector V_inlet = {rho, u, v, p};
    Vector U_inlet(4, 0.0);

    primtocons(U_inlet.data(), V_inlet.data(), bigG, 2); 

    MPI_Init(&argc, &argv);

    bool modelling_real_gas = false;
    bool using_bilinear_interpolation = false; 
    string filename = "../plotfiles/perf_gas_400x200_ramp.dat";

    Solver2D solver(Nx, Ny, CFLs, CFL_timesteps, U_inlet, modelling_real_gas, using_bilinear_interpolation, filename);
    
    solver.solve(); 
 
    MPI_Finalize(); 

    return 0;

}
