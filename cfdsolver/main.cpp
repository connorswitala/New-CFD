#include "../solverlib/solver.hpp"


int main(int argc, char* argv[]) {

    int Nx = 100, Ny = 50;
    double CFL = 1.0; 


    Vector V_inlet = {0.01, 1000, 0.0, 10000};
    Vector U_inlet(4, 0.0);

    primtocons(U_inlet.data(), V_inlet.data(), 2);

    MPI_Init(&argc, &argv);

    Solver2D solver(Nx, Ny, CFL, U_inlet);

    solver.solve(); 
 
    MPI_Finalize(); 

}