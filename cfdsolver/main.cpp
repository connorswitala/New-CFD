#include "../solverlib/solver.hpp"


int main(int argc, char* argv[]) {

    int Nx = 200, Ny = 100;
    double CFL = 1.0; 


    Vector V_inlet = {0.01, 4000, 0.0, 10000};
    Vector U_inlet(4, 0.0);

    primtocons(U_inlet.data(), V_inlet.data(), 2);

    MPI_Init(&argc, &argv);

    bool modelling_real_gas = false;
    bool using_bilinear_interpolation = false; 
    Solver2D solver(Nx, Ny, CFL, U_inlet, modelling_real_gas, using_bilinear_interpolation);

    solver.solve(); 
 
    MPI_Finalize(); 

    return 0;

}