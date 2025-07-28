#include "../solverlib/solver.hpp"


int main(int argc, char* argv[]) {


    int Nx = 10, Ny = 10;
    double CFL = 1.0; 

    Vector V_inlet = {0.088, 500, 0.0, 5474};
    Vector U_inlet(4, 0.0);
    double bigG = 1.4;

    primtocons(U_inlet.data(), V_inlet.data(), bigG, 2); 

    MPI_Init(&argc, &argv);

    bool modelling_real_gas = true;
    bool using_bilinear_interpolation = false; 
    string filename = "../plotfiles/perf_gas_ramp.dat";

    Solver2D solver(Nx, Ny, CFL, U_inlet, modelling_real_gas, using_bilinear_interpolation, filename);
    solver.solve(); 
 
    MPI_Finalize(); 

    return 0;

}