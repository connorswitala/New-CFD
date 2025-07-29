#include "../solverlib/solver.hpp"


int main(int argc, char* argv[]) {

    
    int Nx = 200, Ny = 100;
    double CFL = 1.0; 

    double Mach = 20;
    double a = 329.799;
    double u = Mach * a;

    Vector V_inlet = {0.000977525, u, 0.0, 75.95};
    Vector U_inlet(4, 0.0);
    double bigG = 1.4;

    primtocons(U_inlet.data(), V_inlet.data(), bigG, 2); 

    MPI_Init(&argc, &argv);

    bool modelling_real_gas = true;
    bool using_bilinear_interpolation = false; 
    string filename = "../plotfiles/real_gas_coupled_ramp.csv";

    Solver2D solver(Nx, Ny, CFL, U_inlet, modelling_real_gas, using_bilinear_interpolation, filename);
    
    solver.solve(); 
 
    MPI_Finalize(); 

    return 0;

}
