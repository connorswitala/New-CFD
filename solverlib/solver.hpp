#pragma once

#include "../linalglib/linalg.hpp"
#include "../gibbslib/gibbs.hpp" 
#include <mpi.h>
#include <fstream>
#include <iomanip> 
#include <cmath> 
#include <memory>
#include <algorithm>  // for std::upper_bound

// Constant values for number of variables
constexpr double perfgam = 1.4;
constexpr double perfR = 287.0;

// This enum class is for setting boundary conditions types
enum class BCType {
	IsothermalWall,
	AdiabaticWall,
	Inlet,
	Outlet,
	Symmetry,
	Undefined
};
// This struct contains the boundary conditions types for each side of the grid (left, right, bottom, top) 
struct BCMap {

	BCType left;
	BCType right;
	BCType bottom;
	BCType top;


	BCMap(BCType left, BCType right, BCType bottom, BCType top) : left(left), right(right), bottom(bottom), top(top) {}
    BCMap() : left(BCType::Undefined), right(BCType::Undefined), bottom(BCType::Undefined), top(BCType::Undefined) {}
};

struct Point {
    double x;
    double y;

    Point(double x_val = 0.0, double y_val = 0.0) : x(x_val), y(y_val) {};
};

void primtocons(double* U, const double* V, double gam, const int dimensions);
void constoprim(const double* U, double* V, double gam, const int dimensions);

ThermoEntry operator+(const ThermoEntry& A, const ThermoEntry& B);
ThermoEntry operator*(const double& s, const ThermoEntry& A);

double NewtonMethod(double max_dist, int n_points, double d_min);

inline double computeInternalEnergy(const double* U, int n_vel) {
    double udotu = 0.0;
    for (int i = 0; i < n_vel; ++i) 
        udotu += U[i + 1] * U[i + 1]; 

    return U[n_vel + 1] / U[0] - 0.5 / (U[0] * U[0]) * udotu;   
}
inline double computePressure(const double* U, const double gam, int n_vel) {
    double udotu = 0.0;
    for (int i = 0; i < n_vel; ++i) 
        udotu += U[i + 1] * U[i + 1]; 
    
    return (U[n_vel + 1] - 0.5 / U[0] * udotu) * (gam - 1);
}
inline double compute_total_enthalpy(const double* V, const double gam, int n_vel) {
    double udotu = 0.0;
    for (int i = 0; i < n_vel; ++i) 
        udotu += V[i + 1] * V[i + 1];
    
    return (V[n_vel + 1] / (gam - 1) + 0.5 * V[0] * udotu + V[n_vel + 1]) / V[0];
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

class Solver2D {
private: 

    int n_vel = 2; // Number of momentum being solved
    int n = 4;     // Number of conserved variables
    int n_ghosts = 2;
    int buffers = 2 * n_ghosts;   

    int Nx, Ny, Nx_local; 
    int num_cells, num_ifaces, num_jfaces, num_loc_cells, num_gathered_cells; 
    double CFL, dt, t, outer_res, inner_res;
    Vector CFL_vec;
    vector<int> CFL_timesteps; 
    string filename;

    // CFD data structures
    Vector local_Nx, U, U_inlet, U_gathered, dU_old, dU_new, iEL, iER, jEB, jET, iFlux, jFlux, iPlus_A, iMinus_A, 
        jPlus_A, jMinus_A, irho_A, jrho_A, irho_flux, jrho_flux, V, V1, V2, Q, W, int1, int2, int3, UL;

    Vector old_rhos, old_es;
    
    // Thermochemical table data structures
    bool real_gas, using_table;
    vector<ThermoEntry> thermochemical_table;
    vector<ThermoEntry> cell_thermo; 
    Vector internal_energies, densities; 

    Chemistry chem;  
    BCMap BCs; 

    // Intermediate math vectors.
    Vector Vi, Vii, Vj, Vjj, Vp, Vm, F_plus, F_minus, A, B, C, F, v, g, alpha, rv1, rv2, rm0, rm1, rm2, I;

    // Grid data structures per rank
    Vector xCenter, yCenter, Volume, iFxNorm, iFyNorm, jFxNorm, jFyNorm, iArea, jArea;
    
    // Entire grid data structures only on rank 0
    Vector x_vertices, y_vertices, x_cellCenters, y_cellCenters,
		iface_xNormals,	iface_yNormals, jface_xNormals, jface_yNormals,
		iAreas,	jAreas, cellVolumes;

public:

    int rank, size;

    Solver2D(int Nx, int Ny, Vector CFL_vec, vector<int> CFL_timesteps, Vector U_inlet, bool real_gas, bool using_table, string filename); 
    
    void initialize(); 

    void solve();



    void exchange_U_ghost_cells();
    void exchange_dU_ghost_cells(); 
    void form_inviscid_boundary_U(); 
    Vector get_U_values(BCType type, double* U_inner, double x_norm, double y_norm);
    void form_inviscid_boundary_E();
    Vector get_E_values(BCType type, double x_norm, double y_norm);

    void compute_dt();
    void compute_inner_res_perf();
    void compute_outer_res_perf();
    void compute_inner_res_real();
    void compute_outer_res_real();


    void compute_ifluxes();
    void compute_jfluxes(); 
    void compute_rho_fluxes(); 

    void perf_line_relaxation(); 
    void relax_left_line_perf();
    void relax_inner_lines_perf(int i);
    void relax_right_line_perf(); 

    void real_line_relaxation(); 
    void relax_left_line_real();
    void relax_inner_lines_real(int i);
    void relax_right_line_real(); 

    void update_U();
    void explicit_update();
    void finalize(); 

    void print_by_rank(Vector Vec, int nx, int ny, int nvars, string name); 
    void create_ramp_grid(double L, double inlet_height, double ramp_angle);
    void create_cylinder_grid(double Cylinder_Radius, double R1, double R2, double d_min, double theta_left, double theta_right); 

    void writeTecplotDat();
    void writeParaviewCSV(); 

    void initialize_chemistry(); 
    void get_real_chemistry();
    void get_perf_chemistry();
    void load_thermochemical_table(); 
    ThermoEntry bilinear_interpolate(double rho, double e);

    bool is_near_inlet_state(const double* U) const;


};