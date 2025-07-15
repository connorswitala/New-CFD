#include "../linalglib/linalg.hpp"
#include <mpi.h>
#include <fstream>
#include <iomanip> 
#include <cmath> 
#include <memory>

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

};

struct Point {
    double x;
    double y;

    Point(double x_val = 0.0, double y_val = 0.0) : x(x_val), y(y_val) {};
};


void primtocons(double* U, const double* V, const int dimensions);
void constoprim(const double* U, double* V, const int dimensions);


inline double computeInternalEnergy(const double* U, int n_vel) {
    double udotu = 0.0;
    for (int i = 0; i < n_vel; ++i) {
        udotu += U[i + 1] * U[i + 1]; 
    }

    return U[n_vel + 1] / U[0] - 0.5 / (U[0] * U[0]) * udotu;   
}
inline double computePressure(const double* U, int n_vel) {
    double udotu = 0.0;
    for (int i = 0; i < n_vel; ++i) {
        udotu += U[i + 1] * U[i + 1]; 
    }
    return (U[n_vel + 1] - 0.5 / (U[0] * U[0]) * udotu) * (perfgam - 1);
}
inline double compute_total_enthalpy(const double* U, int n_vel) {
    return (U[n_vel + 1] + computePressure(&U[0], n_vel)) / U[0]; 
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

    int n_vel = 2;
    int n = 4;

    int Nx, Ny, Nx_local, N_cells, N_local;
    int num_cells, num_ifaces, num_jfaces, num_loc_cells, num_gathered_cells; 
    double CFL, dt, t, outer_res, inner_res;
    string save_filename;

    // CFD data structures
    Vector U, U_inlet, U_gathered, dU_old, dU_new, iEL, iER, jEB, jET, iFlux, jFlux, iPlus_A, iMinus_A, 
        jPlus_A, jMinus_A, irho_A, jrho_A, V, V1, V2, Q, W, int1, int2, int3, UL, UR;
    Vector local_Nx;
    
    // Grid data structures
    Vector xCenter, yCenter, Volume, iFxNorm, iFyNorm, jFxNorm, jFyNorm, iArea, jArea;
    
    Vector x_vertices, y_vertices, x_cellCenters, y_cellCenters,
		iface_xNormals,	iface_yNormals, jface_xNormals, jface_yNormals,
		iAreas,	jAreas,
		cellVolumes;

    Vector A, B, C, F, v, g, alpha, result_vector1, result_vector2, result_matrix, I; 
    
    BCMap BCs; 

public:

    int rank, size;

    Solver2D(int Nx, int Ny, double CFL, Vector U_inlet); 

    void solve();

    void exchange_ghost_cells();
    void form_inviscid_boundary_U(); 
    Vector get_U_values(BCType type, double* U_inner, double x_norm, double y_norm);
    void form_inviscid_boundary_E();
    Vector get_E_values(BCType type, double x_norm, double y_norm);

    void compute_dt();
    void compute_inner_res();
    void compute_outer_res();

    void compute_ifluxes();
    void compute_jfluxes(); 
    void compute_rho_fluxes(); 
    void relax_left_line();
    void relax_inner_lines();
    void relax_right_line(); 
    void update_U();

    void create_ramp_grid(double L, double inlet_height, double ramp_angle);

};