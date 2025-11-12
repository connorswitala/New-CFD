#include "../solverlib/solver.hpp"


int main(int argc, char* argv[]) {

    Config cfg; // Type for holding solver configuration variables.

    cfg.real_gas = false;   // Boolean to check if I am using chemical equilibrium solver
    cfg.interp = false;     // Boolean to check if I am bilinear interpolation for chem eqbm. False if above if false.
    cfg.filename = "../plotfiles/perf_gas_400x200_cylinder.csv"; // .dat/.csv filename I want to store data to
    
    cfg.Nx = 400; // number of x-cells
    cfg.Ny = 200; // number of y-cells

    cfg.CFL_timesteps = {0, 200, 400, 800, 1000, 1200}; // Timestep list for when to change CFLs
    cfg.CFLs = {1.0, 2.0, 3.0, 5.0, 10.0, 5.0, 1.0};    // CFLs


    // Flow parameters
    double bigG = 1.4;
    double Mach = 2.5;
    double T = 226.15;
    double p = 900;
    double rho = 0.01388;
    double a = sqrt(bigG * p / rho);
    double u = Mach * a;
    double v = 0.0;

    Vector V_inlet = {rho, u, v, p};
    cfg.U_inlet = Vector(4, 0.0);

    primtocons(cfg.U_inlet.data(), V_inlet.data(), bigG, 2); // Convert primitives to conserved variables

    MPI_Init(&argc, &argv); // Initialize mpi

    Solver2D solver(cfg); // Call Solver constructor
    
    // Create grid type
    solver.create_ramp_grid(3, 0.75, 15);
    //solver.create_cylinder_grid(0.1, 0.3, 0.45, 0.0001, pi, pi / 2);

    // Run the CFD simulation.
    solver.solve();
 
    MPI_Finalize(); 

    return 0;

}
