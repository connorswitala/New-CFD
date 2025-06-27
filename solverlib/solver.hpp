#include "../linalglib/linalg.hpp"
#include <mpi.h>
#include <fstream>

// Constant values for number of variables
constexpr double perfgam = 1.4;
constexpr double perfR = 287.0;

void primtocons(double* U, const double* V, const int dimensions);
void constoprim(const double* U, double* V, const int dimensions);


inline double computeInternalEnergy(const double* U, int dimensions) {
    double udotu = 0.0;
    for (int i = 0; i < dimensions; ++i) {
        udotu += U[i + 1] * U[i + 1]; 
    }

    return U[dimensions + 1] / U[0] - 0.5 / (U[0] * U[0]) * udotu;   
}
inline double computePressure(const double* U, int dimensions) {
    double udotu = 0.0;
    for (int i = 0; i < dimensions; ++i) {
        udotu += U[i + 1] * U[i + 1]; 
    }
    return (U[dimensions + 1] - 0.5 / (U[0] * U[0]) * udotu) * (perfgam - 1);
}


class SodSolver1D { 
private:

    int n_vel = 1;
    int n = 3;

    int Nx, Nx_local;
    double CFL; 
    double dt, dx, L, t;   
    string save_filename; 

    Vector U, U_gathered, Flux, F_plus, F_minus, A_plus, A_minus, A_rho, V, V1, V2, Q, W, int1, int2, int3, UL, UR;   


public: 

    int rank, size;
    
    SodSolver1D(int Nx, double CFL);   

    void solve();
    void exchange_ghost_cells(); 
    void compute_dt();
    void compute_fluxes(); 
    void update_U();
    void write_U_to_csv(string filename);

};