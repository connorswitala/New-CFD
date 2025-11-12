#include "../solverlib/solver.hpp"


int main(int argc, char* argv[]) {

    Config cfg;

    cfg.real_gas = false;
    cfg.interp = false; 
    cfg.filename = "../plotfiles/perf_gas_400x200_cylinder.csv";
    
    cfg.Nx = 400;
    cfg.Ny = 200;

    cfg.CFL_timesteps = {0, 200, 400, 800, 1000, 1200};
    cfg.CFLs = {1.0, 2.0, 3.0, 5.0, 10.0, 5.0, 1.0};

    double bigG = 1.4;
    double Mach = 2.5;
    double T = 226.15;
    double p = 900;
    double rho = 0.01388;
    double a = sqrt(bigG * p / rho);
    double u = Mach * a;
    double v = 0.0;

    Vector V_inlet = {rho, u, v, p};

    primtocons(cfg.U_inlet.data(), V_inlet.data(), bigG, 2); 

    MPI_Init(&argc, &argv);



    Solver2D solver(cfg);
    
    solver.create_ramp_grid(3, 0.75, 15); 
    //solver.create_cylinder_grid(0.1, 0.3, 0.45, 0.0001, pi, pi / 2);

    solver.solve();

 
    MPI_Finalize(); 

    return 0;

}
