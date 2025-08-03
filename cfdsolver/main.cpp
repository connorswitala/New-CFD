#include "../solverlib/solver.hpp"


int main(int argc, char* argv[]) {

    
    int Nx = 2000, Ny = 1000;

    vector<int> CFL_timesteps = {0, 1000, 2000, 5000, 10000, 15000};
    Vector CFLs = {1.0, 2.0, 3.0, 5.0, 10.0, 5.0, 1.0};

    double bigG = 1.4;
    double Mach = 2.5;
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

    bool modelling_real_gas = true;
    bool using_bilinear_interpolation = false; 
    string filename = "../plotfiles/perf_gas_400x200_cylinder.csv";

    Solver2D solver(Nx, Ny, CFLs, CFL_timesteps, U_inlet, modelling_real_gas, using_bilinear_interpolation, filename);
    
    solver.create_ramp_grid(3, 0.75, 15); 
    
    // solver.create_cylinder_grid(0.1, 0.3, 0.45, 0.0001, pi, pi / 2);

    solver.solve();

 
    MPI_Finalize(); 

    return 0;

}
