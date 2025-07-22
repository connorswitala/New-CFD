#include "../solverlib/solver.hpp"


int main(int argc, char* argv[]) {

    int Nx = 400, Ny = 200;
    double CFL = 2.0; 


    Vector V_inlet = {0.01, 4000, 0.0, 10000};
    Vector U_inlet(4, 0.0);

    primtocons(U_inlet.data(), V_inlet.data(), 2);

    MPI_Init(&argc, &argv);

    Solver2D solver(Nx, Ny, CFL, U_inlet);

    solver.solve(); 
 
    MPI_Finalize(); 

    return 0;

}