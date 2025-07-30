#include <fstream>
#include <iostream>
#include "../solverlib/solver.hpp" 
using namespace std;

int main(int argc, char* argv[]) {
    
    int Nx = 400, Ny = 200;

    vector<int> CFL_timesteps = {0, 250, 500, 750, 1000, 1500, 2000, 2500, 3000};
    Vector CFLs = {0.1, 1.0, 1.5, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0};

    double Mach = 20;
    double a = 329.799;
    double u = Mach * a;

    Vector V_inlet = {0.000977525, u, 0.0, 75.95};
    Vector U_inlet(4, 0.0);
    double bigG = 1.4;

    primtocons(U_inlet.data(), V_inlet.data(), bigG, 2); 

    MPI_Init(&argc, &argv);

    bool modelling_real_gas = true;
    bool using_bilinear_interpolation = true; 
    string filename = "../plotfiles/perf_gas_1000x1000_ramp.dat";

    Solver2D solver(Nx, Ny, CFLs, CFL_timesteps, U_inlet, modelling_real_gas, using_bilinear_interpolation, filename);

 
    MPI_Finalize(); 

    return 0;

}