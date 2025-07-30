#include "solver.hpp"


// Conversion functions
void primtocons(double* U, const double* V, const double gam, const int dimensions) {
    
    U[0] = V[0]; // Set density

    double udotu = 0.0; // Set u dot u = 0.0

    // This loop computes each momentum and total udotu for n dimensions
    for (int i = 0; i < dimensions; ++i) {
        U[i + 1] = V[0] * V[i + 1]; 
        udotu += V[i + 1] * V[i + 1];  
    }

    U[dimensions + 1] = V[dimensions + 1] / (gam - 1) + 0.5 * V[0] * udotu; // Set total energy 
}
void constoprim(const double* U, double* V, const double gam, const int dimensions) {
    
    V[0] = U[0]; // Set density

    double udotu = 0.0; // Set u dot u = 0.0;

    // This loop computes each velocit;y and total udotu for n dimensions
    for (int i = 0; i < dimensions; ++i) {
        V[i + 1] = U[i + 1] / U[0];
        udotu += U[i + 1] * U[i + 1];
    }

    V[dimensions + 1] = (U[dimensions + 1] - 0.5 / U[0] * udotu) * (gam - 1); // Set pressure
}

ThermoEntry operator+(const ThermoEntry& A, const ThermoEntry& B) {
	ThermoEntry result;
	result.rho = A.rho + B.rho;
	result.e = A.e + B.e;
	result.p = A.p + B.p;
	result.T = A.T + B.T;
	result.R = A.R + B.R;
	result.cv = A.cv + B.cv;
	result.gamma = A.gamma + B.gamma;
	result.dpdrho = A.dpdrho + B.dpdrho;
	result.dpde = A.dpde + B.dpde;
    result.a = A.a + B.a;
	return result;
}
ThermoEntry operator*(const double& s, const ThermoEntry& A) {
	ThermoEntry result;
	result.rho = s * A.rho;
	result.e = s * A.e;
	result.p = s * A.p;
	result.T = s * A.T;
	result.R = s * A.R;
	result.cv = s * A.cv;
	result.gamma = s * A.gamma;
	result.dpdrho = s * A.dpdrho;
	result.dpde = s * A.dpde;
    result.a = s * A.a;
	return result;
}

double NewtonMethod(double max_dist, int n_points, double d_min) {
	double k = 1, k_new = 1 / 2, ratio = fabs(k - k_new);
	double func, func_prime;

	while (ratio >= 0.00000000001) {
		func = d_min - max_dist * (exp(k / (n_points - 1)) - 1) / (exp(k) - 1);
		func_prime = -max_dist * (((1 / (n_points - 1) * exp(k / (n_points - 1))) * (exp(k) - 1) - (exp(k / (n_points - 1)) - 1) * exp(k)) / ((exp(k) - 1) * (exp(k) - 1)));
		k_new = k - func / func_prime;
		ratio = fabs(k - k_new);
		k = k_new;
	}

	return k;
}


/////////////////////////////////////////////////
//////// Sod Shock Tube Solver Functions ////////
/////////////////////////////////////////////////
SodSolver1D::SodSolver1D(int Nx, double CFL) : Nx(Nx), CFL(CFL) {


    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    Nx_local = Nx / size;

    L = 2.0;
    dx = L / (Nx - 1);
    
    t = 0.0;
    dt = 1e-5;

    U_gathered = Vector((Nx + 2) * n, 0.0);
    UL = Vector(n, 0.0);
    UR = Vector(n, 0.0);

    if (rank == 0) {
        
        Vector VL = {1.0, 0.0, 1.0};
        Vector VR = {0.125, 0.0, 0.1};
        primtocons(UL.data(), VL.data(), perfgam, n_vel);
        primtocons(UR.data(), VR.data(), perfgam, n_vel);

        for (int i = 0; i < Nx + 2; ++i) {
            for (int k = 0; k < n; ++k) {
                U_gathered[i * n + k] = (i <= Nx / 2) ? UL[k] : UR[k];
            }
        }
    }

    MPI_Bcast(UL.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(UR.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
    U = Vector((Nx_local + 2) * n, 0.0); 
    Flux = Vector((Nx_local + 1) * n, 0.0);
    F_plus = Vector((Nx_local + 1) * n, 0.0);
    F_minus = Vector((Nx_local + 1) * n, 0.0);
    A_plus = Vector((Nx_local + 1) * n * n, 0.0);
    A_minus = Vector((Nx_local + 1) * n * n, 0.0);  

    V = Vector(n, 0.0);
    V1 = Vector(n, 0.0);
    V2 = Vector(n, 0.0);
    Q = Vector(n, 0.0);
    W = Vector(n, 0.0);
    int1 = Vector(n * n, 0.0);
    int2 = Vector(n * n, 0.0);
    int3 = Vector(n * n, 0.0);

    MPI_Scatter(U_gathered.data() + n, Nx_local * n, MPI_DOUBLE,
                U.data() + n, Nx_local * n, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
}
void SodSolver1D::exchange_ghost_cells() {
    MPI_Status status_left, status_right;

    // Exchange with left neighbor
    if (rank > 0) {
        MPI_Sendrecv(U.data() + n, n, MPI_DOUBLE, rank - 1, 0,
                     U.data(), n, MPI_DOUBLE, rank - 1, 1,
                     MPI_COMM_WORLD, &status_left);
    }

    // Exchange with right neighbor
    if (rank < size - 1) {
        MPI_Sendrecv(U.data() + Nx_local * n, n, MPI_DOUBLE, rank + 1, 1,
                     U.data() + (Nx_local + 1) * n, n, MPI_DOUBLE, rank + 1, 0,
                     MPI_COMM_WORLD, &status_right);
    }

    if (rank == 0) {
        for (int k = 0; k < n; ++k) {
            U[k] = UL[k];  // ghost cell at index 0 * n + k
        }
    }

    // Set right ghost cell (index Nx_local+1) on last rank to UR (boundary right state)
    if (rank == size - 1) {
        int ghost_idx = (Nx_local + 1) * n;
        for (int k = 0; k < n; ++k) {
            U[ghost_idx + k] = UR[k];
        }
    }
}
void SodSolver1D::compute_dt() { 

    double local_min_dt = std::numeric_limits<double>::max();

    for (int i = 0; i < Nx_local; ++i) { // skip ghost cells
        int idx = i * 3;
        double rho = U[idx + 0];
        double mom = U[idx + 1];
        double E   = U[idx + 2];

        double u = mom / rho;
        double e = E / rho - 0.5 * u * u;
        double p = (perfgam - 1) * rho * e;
        double a = std::sqrt(perfgam * p / rho);

        double local_speed = std::abs(u) + a;

        if (local_speed > 1e-8) {
            double local_dt = CFL * dx / local_speed;
            local_min_dt = std::min(local_min_dt, local_dt);
        }
    }

    // Get global minimum dt across all ranks
    MPI_Allreduce(&local_min_dt, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);


}
void SodSolver1D::write_U_to_csv(string filename) { 


    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "FAILED to open file: " << filename << '\n';
    }
  
    file << "rho,u,p\n";
    Vector V(3, 0.0); 

    for (int i = 0; i < Nx + 2; ++i) {
        int idx = i * 3;
        constoprim(&U_gathered[idx], V.data(), perfgam, n_vel);
        file << V[0] << ',' << V[1] << ',' << V[2] << '\n';
    }

    file.close();
}
void SodSolver1D::compute_fluxes() {

    for (int i = 0; i <= Nx_local; ++i) {
            int idx = i * n;
            int iidx = (i + 1) * n;
            int aidx = i * n * n;

            // Positive flux calculation
            constoprim(&U[idx], V.data(), perfgam, n_vel);
            double rho = V[0];
            double u = V[1];
            double p = V[2];
            double a = sqrt(perfgam * p / rho);

            double lp = 0.5 * (u + a + fabs(u + a));
            double lm = 0.5 * (u - a + fabs(u - a));
            double l = 0.5 * (u + fabs(u));
            double lt = 0.5 * (lp - lm);
            double lc = 0.5 * (lp + lm - 2 * l);
   
            fill(int1.begin(), int1.end(), 0.0);
            for (int k = 0; k < n; ++k) int1[k * n + k] = l;

            V1 = {lc / (a * a), (u * lc + a * lt) / (a * a), ((U[idx + 2] + computePressure(&U[idx], perfgam, n_vel)) / U[idx] * lc + a * u * lt) / (a * a)};
            Q = {0.5 * u * u * (perfgam - 1), -u * (perfgam - 1), (perfgam - 1)};
            V2 = {lt / a, u * lt / a + lc, (U[idx + 2] + computePressure(&U[idx], perfgam,  n_vel)) / U[idx] * lt / a + u * lc};
            W = {-u, 1, 0};

            outer_product(V1.data(), Q.data(), int2.data(), n);
            outer_product(V2.data(), W.data(), int3.data(), n);

            for (int k = 0; k < n * n; ++k) 
                A_plus[aidx + k] = int1[k] + int2[k] + int3[k];

            matvec_mult(&A_plus[aidx], &U[idx], &F_plus[idx], n);

            // Negative flux calculation
            constoprim(&U[iidx], V.data(), perfgam, n_vel);
            rho = V[0]; 
            u = V[1]; 
            p = V[2];
            a = sqrt(perfgam * p / rho);

            lp = 0.5 * (u + a - fabs(u + a));
            lm = 0.5 * (u - a - fabs(u - a));
            l = 0.5 * (u - fabs(u));
            lt = 0.5 * (lp - lm);
            lc = 0.5 * (lp + lm - 2 * l);

            fill(int1.begin(), int1.end(), 0.0);
            for (int k = 0; k < n; ++k) int1[k * n + k] = l;

            V1 = {lc / (a * a), (u * lc + a * lt) / (a * a), ((U[iidx + 2] + computePressure(&U[iidx], perfgam, n_vel)) / U[iidx] * lc + a * u * lt) / (a * a)};
            Q = {0.5 * u * u * (perfgam - 1), -u * (perfgam - 1), (perfgam - 1)};
            V2 = {lt / a, u * lt / a + lc, (U[iidx + 2] + computePressure(&U[iidx], perfgam, n_vel)) / U[iidx] * lt / a + u * lc};
            W = {-u, 1, 0};

            outer_product(V1.data(), Q.data(), int2.data(), n);
            outer_product(V2.data(), W.data(), int3.data(), n);

            for (int k = 0; k < n * n; ++k) A_minus[aidx + k] = int1[k] + int2[k] + int3[k];

            matvec_mult(&A_minus[aidx], &U[iidx], &F_minus[idx], n);

            // Total flux calculation
            for (int k = 0; k < n; ++k) 
                Flux[idx + k] = F_plus[idx + k] + F_minus[idx + k];
    }
}
void SodSolver1D::update_U() {
    for (int i = 0; i < Nx_local; ++i) {
            int idx = (i + 1) * n;
            int fidx = i * n;
            int fiidx = (i + 1) * n;
            for (int k = 0; k < n; ++k) {
                U[idx + k] += -dt / dx * (-Flux[fidx + k] + Flux[fiidx + k]);
            }
    }
}
void SodSolver1D::solve() {
    
    int counter = 0;

    exchange_ghost_cells();      
   
    while (t <= 0.5) {
        compute_dt();    
        exchange_ghost_cells();      
        compute_fluxes();      
        update_U();

        // std::cout << "Rank " << rank << " U = [";
        // for (int i = 0; i < (Nx_local + 2) * n; ++i) {
        //     std::cout << U[i];
        //     if (i != (Nx_local + 2) * n - 1) std::cout << ", ";
        // }
        // std::cout << "]\n";
        // std::cout.flush();

        t += dt; 

        if (counter % 100 == 0 && rank == 0) {
            cout << "t = " << t << "\tdt = " << dt << std::endl;
        }
        counter++;

    }

    MPI_Gather(U.data() + n, Nx_local * n, MPI_DOUBLE,
           U_gathered.data(), Nx_local * n, MPI_DOUBLE,
           0, MPI_COMM_WORLD);


    if (rank == 0) {
        write_U_to_csv("sod_shock_mpi.csv");
        cout << "Program finished!" << endl;
    }
}


////////////////////////////////////
//////// Solver2D Functions ////////
////////////////////////////////////

Solver2D::Solver2D(int Nx, int Ny, Vector CFL_vec, vector<int> CFL_timesteps, Vector U_inlet, bool real_gas, bool using_table, string filename) : 
    Nx(Nx), Ny(Ny), CFL_vec(CFL_vec), CFL_timesteps(CFL_timesteps), U_inlet(U_inlet), real_gas(real_gas), using_table(using_table), filename(filename) {
    
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    t = 0.0; 
    inner_res = 1.0;
    outer_res = 1.0;
    CFL = 1.0;

    /** 
     *      local_Nx(size) creates a vector that holds how many columns each rank gets. This is important 
     * since we want to use any number of processes for any number of cells in the x-direction. It will have size
     * Nx / size (+1 if r < remainder of Nx_global / size).
     * 
     *      sendcounts(size) creates a vector that holds how many elements to send to each rank. Since it is the total 
     * number of cells being sent for a block, it will have size local_Nx[r] * Ny.
     * 
     *      displacements(size) is the starting index in the sendbuf for rank r.
     */


    // GRID CREATION //

    // create_ramp_grid(3, 0.75, 20); 
    create_cylinder_grid(0.1, 0.3, 0.45, 0.0001, 3 * pi / 2, pi / 2);    

    U_gathered = Vector( num_gathered_cells * n, 0.0); 
    U = Vector((Nx_local + 2) * (Ny + 2) * n, 0.0);  
    dU_old = Vector((Nx_local + 2) * (Ny + 2) * n, 0.0); 
    dU_new = Vector((Nx_local + 2) * (Ny + 2) * n, 0.0); 

    // Important Flux Vectors
    iFlux = Vector( num_ifaces * n, 0.0);
    jFlux = Vector( num_jfaces * n, 0.0); 
    irho_flux = Vector( num_ifaces * n, 0.0);
    jrho_flux = Vector( num_jfaces * n, 0.0); 

    // Important Jacobians
    iPlus_A = Vector( num_ifaces * n * n, 0.0);
    iMinus_A = Vector( num_ifaces * n * n, 0.0);
    jPlus_A = Vector( num_jfaces * n * n, 0.0);
    jMinus_A = Vector( num_jfaces * n * n, 0.0);
    irho_A = Vector( num_ifaces * n * n, 0.0);
    jrho_A = Vector( num_jfaces * n * n, 0.0); 

    // For flux calculations:
    Vi = Vector(n); 
    Vii = Vector(n);
    Vj = Vector(n);
    Vjj = Vector(n);
    Vp = Vector(n);
    Vm = Vector(n);
    F_plus = Vector(n);
    F_minus = Vector(n);

    // For line relaxation
    A = Vector(n * n, 0.0); 
    B = Vector(n * n, 0.0);
    C = Vector(n * n, 0.0); 
    F = Vector(n, 0.0); 

    alpha = Vector(n * n, 0.0);
    v = Vector(Ny * n, 0.0);
    g = Vector(Ny * n * n, 0.0); 

    // Intermediate vectors and matrices;
    rm0 = Vector(n * n, 0.0); 
    rm1 = Vector(n * n, 0.0);
    rm2 = Vector(n * n, 0.0);
    rv1 = Vector(n, 0.0);
    rv2 = Vector(n, 0.0); 
    I = {1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 1, 0,
         0, 0, 0, 1};  

    // These hold the E matrices that never change for the boundary conditions. The first
    // indice is for left or bottom, while the second is for right and top. 
    iEL = Vector(Ny * n * n, 0.0);
    iER = Vector(Ny * n * n, 0.0); 
    jEB = Vector(Nx_local * n * n, 0.0);  
    jET = Vector(Nx_local * n * n, 0.0);   

    // Intermediate matrices for calculations
    V = Vector(n, 0.0);
    V1 = Vector(n, 0.0);
    V2 = Vector(n, 0.0);
    Q = Vector(n, 0.0);
    W = Vector(n, 0.0);
    int1 = Vector(n * n, 0.0);
    int2 = Vector(n * n, 0.0);
    int3 = Vector(n * n, 0.0);

    // This sets all cells to the inflow condition
    for (int i = 0; i < Nx_local + 2; ++i) {
        for (int j = 0; j < Ny + 2; ++j) {
            for (int k = 0; k < n; ++k) {
                U[(i * (Ny + 2) + j) * n + k] = U_inlet[k];
            }
        }
    }  

    /**
     *  Chemistry initialization
     */

    cell_thermo = vector<ThermoEntry>((Nx_local * 2) * (Ny + 2));

    // This checks to see if you are solving a perfect gas EOS, or using Gibbs minimization.
    if (real_gas) {

        if (using_table) {

            densities = Vector(500, 0.0);
            internal_energies = Vector(500, 0.0); 

            for (int i = 0; i < 500; ++i) {
                internal_energies[i] = 717 * 600 + (5e7 - 717 * 600) / (500 - 1) * i;
                densities[i] = 5e-5 + (1 - 5e-5)/(500 - 1) * i;
            }

            thermochemical_table = vector<ThermoEntry>(500 * 500); 
            load_thermochemical_table(); 
        }
        else if (using_table == false) {
            Chemistry chem;
        }
    }

    initialize_chemistry();

    if (rank == 0)
        cout << "Chemistry initialized" << endl;

    form_inviscid_boundary_E();  
    
}


/** 
 * These functions take care of all boundary condition related problems. 'exchance_ghost_cells()' calls everything
 * necessary. It first trades information from the buffer regions of the local U blocks. Then it updates U on the
 * global boundaries. The E matrices are taken care of at the constructor since they only need to be made once. 
 */

void Solver2D::form_inviscid_boundary_U() {
    Vector Uholder(n, 0.0); 

    if (rank == 0) {
        for (int j = 0; j < Ny; ++j) {
            // Interior cell (1, j)
            int inner_idx  = (1 * (Ny + 2) + j + 1) * n;
            int ghost_idx  = (0 * (Ny + 2) + j + 1) * n;

            Uholder = get_U_values(BCs.left, &U[inner_idx], iFxNorm[j], iFyNorm[j]);

            for (int k = 0; k < n; ++k) 
                U[ghost_idx + k] = Uholder[k];    
        }
    }

    if (rank == size - 1) {
        for (int j = 0; j < Ny; ++j) {
            // Interior cell (Nx_local, j)
            int inner_idx = (Nx_local * (Ny + 2) + j + 1) * n;
            int ghost_idx = ((Nx_local + 1) * (Ny + 2) + j + 1) * n;

            Uholder = get_U_values(BCs.right, &U[inner_idx],
                                iFxNorm[Nx_local * Ny + j],
                                iFyNorm[Nx_local * Ny + j]);

            for (int k = 0; k < n; ++k) 
                U[ghost_idx + k] = Uholder[k];
        }
    }

    for (int i = 0; i < Nx_local; ++i) {

        // Bottom ghost cell (j = 0), reference interior cell at j = 1
        int bottom_interior_idx = ((i + 1) * (Ny + 2) + 1) * n;
        int bottom_ghost_idx    = ((i + 1) * (Ny + 2) + 0) * n;

        Uholder = get_U_values(BCs.bottom, &U[bottom_interior_idx], 
                            jFxNorm[i * (Ny + 1) + 0], jFyNorm[i * (Ny + 1) + 0]);

        for (int k = 0; k < n; ++k) 
            U[bottom_ghost_idx + k] = Uholder[k]; 

        // Top ghost cell (j = Ny + 1), reference interior cell at j = Ny
        int top_interior_idx = ((i + 1) * (Ny + 2) + Ny) * n;
        int top_ghost_idx    = ((i + 1) * (Ny + 2) + Ny + 1) * n;

        Uholder = get_U_values(BCs.top, &U[top_interior_idx], 
                            jFxNorm[i * (Ny + 1) + Ny], jFyNorm[i * (Ny + 1) + Ny]);

        for (int k = 0; k < n; ++k) 
            U[top_ghost_idx + k] = Uholder[k];
    }

}
Vector Solver2D::get_U_values(BCType type, double* U_inner, double x_norm, double y_norm) {
    
    Vector ghost(n , 0.0); 
    double u = 0.0, v = 0.0; 

    switch (type) {

    case BCType::Inlet:
        return U_inlet; 

    case BCType::Outlet:
        for (int k = 0; k < n; ++k) 
            ghost[k] = U_inner[k];
        return ghost; 

    case BCType::IsothermalWall:

		u = U_inner[1] / U_inner[0];
		v = U_inner[2] / U_inner[0];

		ghost[0] = U_inner[0];
		ghost[1] = U_inner[0] * (u - 2 * (u * x_norm + v * y_norm) * x_norm);
		ghost[2] = U_inner[0] * (v - 2 * (u * x_norm + v * y_norm) * y_norm);
		ghost[3] = U_inner[3];

        return ghost;
        

    case BCType::AdiabaticWall:

		u = U_inner[1] / U_inner[0];
		v = U_inner[2] / U_inner[0];

		ghost[0] = U_inner[0];
		ghost[1] = U_inner[0] * (u - 2 * (u * x_norm + v * y_norm) * x_norm);
		ghost[2] = U_inner[0] * (v - 2 * (u * x_norm + v * y_norm) * y_norm);
		ghost[3] = U_inner[3];

        return ghost;

    case BCType::Symmetry:

		u = U_inner[1] / U_inner[0];
		v = U_inner[2] / U_inner[0];

		ghost[0] = U_inner[0];
		ghost[1] = U_inner[0] * (u - 2 * (u * x_norm + v * y_norm) * x_norm);
		ghost[2] = U_inner[0] * (v - 2 * (u * x_norm + v * y_norm) * y_norm);
		ghost[3] = U_inner[3];

        return ghost;

    default:
        throw invalid_argument("Unknown boundary condition type.");
    }
}
void Solver2D::form_inviscid_boundary_E() {

    Vector Eholder(n * n, 0.0);

    // LEFT boundary (only rank 0)
    if (rank == 0) {        
        for (int j = 0; j < Ny; ++j) {
            Eholder = get_E_values(BCs.left, iFxNorm[j], iFyNorm[j]);
            for (int k = 0; k < n * n; ++k)  
                iEL[j * n * n + k] = Eholder[k];
        }
    }

    // RIGHT boundary (only final rank)
    if (rank == size - 1) {
        for (int j = 0; j < Ny; ++j) {
            int idx = Nx_local * Ny + j;
            Eholder = get_E_values(BCs.right, iFxNorm[idx], iFyNorm[idx]);
            for (int k = 0; k < n * n; ++k) 
                iER[j * n * n + k] = Eholder[k];

        }
    }

    // TOP and BOTTOM boundaries (all ranks)
    for (int i = 0; i < Nx_local; ++i) {
        int bot_idx = i * (Ny + 1);
        int top_idx = bot_idx + Ny;
 
        Eholder = get_E_values(BCs.bottom, jFxNorm[bot_idx], jFyNorm[bot_idx]); 
        for (int k = 0; k < n * n; ++k) 
            jEB[i * n * n + k] = Eholder[k];


        Eholder = get_E_values(BCs.top, jFxNorm[top_idx], jFyNorm[top_idx]);
        for (int k = 0; k < n * n; ++k) 
            jET[i * n * n + k] = Eholder[k];
     
    }
}
Vector Solver2D::get_E_values(BCType type, double x_norm, double y_norm) {

    Vector holder =     {1, 0, 0, 0,
                        0, 1 - 2 * x_norm * x_norm, -2 * x_norm * y_norm, 0,
                        0, -2 * x_norm * y_norm, 1 - 2 * y_norm * y_norm, 0,
                        0, 0, 0, 1};

    switch (type) {
        
    case BCType::Inlet: 
        holder = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};     
        return holder; 
    
    case BCType::Outlet:
        return I; 

    case BCType::IsothermalWall:
        return holder;

    case BCType::AdiabaticWall:
        return holder; 

    case BCType::Symmetry:
        return holder; 

    default:
        throw invalid_argument("Unknown boundary condition type.");        
    }
}
void Solver2D::exchange_U_ghost_cells() {
    MPI_Status status_left, status_right;

    int col_size = (Ny + 2) * n; // Exchange all variables for entire column of Ny cells. 

    Vector send_leftU(col_size), recv_leftU(col_size);
    Vector send_rightU(col_size), recv_rightU(col_size);

    /** This nested for loop creates vectors for each rank that just holds the innermost
     * real cells of the rank's chunk. These will be used to put into the ghost cell column of
     * the neighboring rank. 
     */


    for (int j = 0; j < Ny + 2; ++j) {
        for (int k = 0; k < n; ++k) {

            int local_left = (1 * (Ny + 2) + j) * n + k;
            int local_right = (Nx_local * (Ny + 2) + j ) * n + k;

            send_leftU[j * n + k] = U[local_left];
            send_rightU[j * n + k] = U[local_right];
        }
    }

    /** This part sends and receives the left or right column depending on the if condition.  
     * If rank > 0, then it send the left column since rank 0 doesnt send. If rank < size - 1, 
     * it sends the right column since rank <size - 1> doesn't have a right column to send. 
     * MPI_Sendrecv looks like this:
     * 
     * MPI_Sendrecv(    <data being send>
     *                  <number of data points being sent>
     *                  <MPI data type being sent>
     *                  <destination of sent data>
     *                  <send tag(just needs to be different)>
     *                  <place to receive data>
     *                  <number of data points being received>
     *                  <MPI data type>
     *                  <source of received data>
     *                  <receive tag>
     *                  <communication>
     *                  <MPI_status *status>
     *                                              )
    */

    if (rank > 0) {
        MPI_Sendrecv(send_leftU.data(), col_size, MPI_DOUBLE, rank - 1, 0,
                     recv_leftU.data(), col_size, MPI_DOUBLE, rank - 1, 1,
                     MPI_COMM_WORLD, &status_left);
    }

    if (rank < size - 1) {
        MPI_Sendrecv(send_rightU.data(), col_size, MPI_DOUBLE, rank + 1, 1,
                     recv_rightU.data(), col_size, MPI_DOUBLE, rank + 1, 0,
                     MPI_COMM_WORLD, &status_right);
    }


    /** This nested for loop takes the send and receive buffs and puts the data
     * from it into the ghost cell spots for each rank. 
     */

    for (int j = 0; j < Ny + 2; ++j) {
        for (int k = 0; k < n; ++k) {

            if (rank > 0) {
                int ghost_left = (0 * (Ny + 2) + j) * n + k;
                U[ghost_left] = recv_leftU[j * n + k];
            }

            if (rank < size - 1) {
                int ghost_right = ( (Nx_local + 1) * (Ny + 2) + j) * n + k;
                U[ghost_right] = recv_rightU[j * n + k];
            }

        }
    }

    form_inviscid_boundary_U(); 
}
void Solver2D::exchange_dU_ghost_cells() {
     MPI_Status status_left, status_right;

    int col_size = (Ny + 2) * n; // Exchange all variables for entire column of Ny cells. 

    Vector send_leftdU(col_size), recv_leftdU(col_size);
    Vector send_rightdU(col_size), recv_rightdU(col_size);

     for (int j = 0; j < Ny + 2; ++j) {
        for (int k = 0; k < n; ++k) {

            int local_left = (1 * (Ny + 2) + j) * n + k;
            int local_right = (Nx_local * (Ny + 2) + j ) * n + k;

            send_leftdU[j * n + k] = dU_old[local_left];
            send_rightdU[j * n + k] = dU_old[local_right];
        }
    }


    if (rank > 0) {
        MPI_Sendrecv(send_leftdU.data(), col_size, MPI_DOUBLE, rank - 1, 0,
                     recv_leftdU.data(), col_size, MPI_DOUBLE, rank - 1, 1,
                     MPI_COMM_WORLD, &status_left);
    }

    if (rank < size - 1) {
        MPI_Sendrecv(send_rightdU.data(), col_size, MPI_DOUBLE, rank + 1, 1,
                     recv_rightdU.data(), col_size, MPI_DOUBLE, rank + 1, 0,
                     MPI_COMM_WORLD, &status_right); 
    }


    /** This nested for loop takes the send and receive buffs and puts the data
     * from it into the ghost cell spots for each rank. 
     */

    for (int j = 0; j < Ny + 2; ++j) {
        for (int k = 0; k < n; ++k) {

            if (rank > 0) {
                int ghost_left =  j * n + k;
                dU_old[ghost_left] = recv_leftdU[j * n + k];
            }

            if (rank < size - 1) {
                int ghost_right = ((Nx_local + 1) * (Ny + 2) + j) * n + k;
                dU_old[ghost_right] = recv_rightdU[j * n + k];
            }

        }
    }

}

void Solver2D::solve() {

    double start_time = MPI_Wtime();

    if (rank == 0) cout << "Starting solve!" << endl;
    int counter = 0;   
    int num = 0; 
    CFL = 1.0;

    while (outer_res > 1e-6) {

        if (counter == CFL_timesteps[num]) {
            CFL = CFL_vec[num];
            num++;
        }

        exchange_U_ghost_cells(); 

        if (real_gas) {
            get_real_chemistry();
        }
        else {
            get_perf_chemistry();
        }

        compute_dt();

        compute_ifluxes();
        compute_jfluxes();

        if (real_gas) {
            real_line_relaxation(); 
        }
        else {
            perf_line_relaxation();
        }

        update_U(); 

        if (real_gas) {
            compute_outer_res_real(); 
        }
        else {
            compute_outer_res_perf();
        }

        counter++;

        if (counter % 50 == 0) {
            double end_time = MPI_Wtime();
            if (rank == 0) { 
                cout << "Iteration: " << counter
                << "\t Inner residual: " << fixed << scientific << setprecision(4) << inner_res
                << "\t Outer residual: " << fixed << scientific << setprecision(4) << outer_res
                << "\tdt: " << fixed << scientific << setprecision(5) << dt 
                << "\tCFL: " << fixed << scientific << setprecision(3) << CFL 
                << "\tElapsed time: " << end_time - start_time << " seconds." << endl;
            }
        }        

        if (counter % 500 == 0)
            finalize();

    }    

    double end_time = MPI_Wtime();

    if (rank == 0) {
        cout << "Solving finished!" << endl;
        cout << "Elapsed time: " << end_time - start_time << " seconds\n";
    }

    finalize(); 

}

void Solver2D::create_ramp_grid(double L, double inlet_height, double ramp_angle) {

    x_vertices = Vector( (Nx + 1) * (Ny + 1), 0.0),   
    y_vertices = Vector( (Nx + 1) * (Ny + 1), 0.0),
    x_cellCenters = Vector(Nx * Ny, 0.0),
    y_cellCenters = Vector(Nx * Ny, 0.0),
    iface_xNormals = Vector((Nx + 1) * Ny, 0.0),
    iface_yNormals = Vector((Nx + 1) * Ny, 0.0),
    jface_xNormals = Vector(Nx * (Ny + 1), 0.0),
    jface_yNormals = Vector(Nx * (Ny + 1), 0.0),
    iAreas = Vector((Nx + 1) * Ny, 0.0),
    jAreas = Vector(Nx * (Ny + 1), 0.0),
    cellVolumes = Vector(Nx * Ny, 0.0);

    BCs = BCMap(BCType::Inlet, BCType::Outlet, BCType::Symmetry, BCType::Symmetry);

    if (rank == 0) {
        int i,j;
        double dx = L / Nx;
        double theta_rad = ramp_angle * 3.141592653 / 180.0;
        double slope = tan(theta_rad);
        double ramp_start_x = L / 3.0;
        double ramp_end_x = 2.0 * L / 3.0;
        double ramp_length = ramp_end_x - ramp_start_x;
        double ramp_end_y = slope * ramp_length;  // height at end of ramp

        for (int i = 0; i <= Nx; ++i) {
            for (int j = 0; j <= Ny; ++j) {

                int x_index = i * (Ny + 1) + j; 
                double x = i * dx;

                x_vertices[x_index] = x;  
                if ( (x < ramp_start_x) && (x + dx > ramp_start_x) ) ramp_start_x = x;
                if ( (x < ramp_end_x) && (x + dx > ramp_end_x) ) {
                    ramp_end_x = x;
                    ramp_length = ramp_end_x - ramp_start_x;
                    ramp_end_y = slope * ramp_length;
                }
            }
        }        

        for (int i = 0; i <= Nx; ++i) {

            double x = i * dx;
            double y_bot, y_top;

            y_top = inlet_height;

            if (x <= ramp_start_x) {
                y_bot = 0.0;
            } else if (x <= ramp_end_x) {
                double x_rel = x - ramp_start_x;
                y_bot = slope * x_rel;
            } else {
                y_bot = ramp_end_y;
            }

            for (int j = 0; j <= Ny; ++j) {

                int y_index = i * (Ny + 1) + j;
                double frac = static_cast<double>(j) / Ny;
                
                y_vertices[y_index] = y_bot + frac * (y_top - y_bot);
                
            }
        }

        // Edge vectors
        Point AB, BC, CD, DA;

        // Calculates cell centers and volumes
        for (i = 0; i < Nx; ++i) {
            for (j = 0; j < Ny; ++j) {

                int ij = i * (Ny + 1) + j;
                int iij = (i + 1) * (Ny + 1) + j;   
                int ijj = i * (Ny + 1) + j + 1; 
                int iijj = (i + 1) * (Ny + 1) + j + 1;   

                DA = { x_vertices[ij] - x_vertices[ijj], y_vertices[ij] - y_vertices[ijj] }; 
                AB = { x_vertices[iij] - x_vertices[ij], y_vertices[iij] - y_vertices[ij] };
                BC = { x_vertices[iijj] - x_vertices[iij], y_vertices[iijj] - y_vertices[iij] };
                CD = { x_vertices[ijj] - x_vertices[iijj], y_vertices[ijj] - y_vertices[iijj] }; 

                int cell_ij = i * Ny + j;   

                x_cellCenters[cell_ij] = (x_vertices[ij] + x_vertices[iij] + x_vertices[iijj] + x_vertices[ijj]) / 4;
                y_cellCenters[cell_ij] =	(y_vertices[ij] + y_vertices[iij] + y_vertices[iijj] + y_vertices[ijj]) / 4;

                cellVolumes[cell_ij] = 0.5 * fabs(DA.x * AB.y - DA.y * AB.x) + 0.5 * fabs(BC.x * CD.y - BC.y * CD.x);

                // cout << "xcenter: " << fixed << setprecision(3) << x_cellCenters[cell_ij] 
                // << "\tycenter: " << fixed << setprecision(3) << y_cellCenters[cell_ij] 
                // << "\tcell volume: " << fixed << setprecision(3) << cellVolumes[cell_ij] << endl; 
            }
        }

        // Calculates geometries for i-faces
        for (i = 0; i <= Nx; ++i) {
            for (j = 0; j < Ny; ++j) {

                int ij = i * (Ny + 1) + j;
                int ijj = i * (Ny + 1) + j + 1;   

                int face_ij = i * Ny + j;  

                AB = { x_vertices[ijj] - x_vertices[ij], y_vertices[ijj] - y_vertices[ij] };

                iAreas[face_ij] = sqrt(AB.x * AB.x + AB.y * AB.y); 

                iface_xNormals[face_ij] = AB.y / fabs(iAreas[face_ij]);
                iface_yNormals[face_ij] = AB.x / fabs(iAreas[face_ij]);

                // cout << "i area: " << fixed << setprecision(3) << iAreas[face_ij] 
                // << "\ti-x norm: " << fixed << setprecision(3) << iface_xNormals[face_ij] 
                // << "\ti-y norm: " << fixed << setprecision(3) << iface_yNormals[face_ij] << endl;
            }
        }

        // Calculates geometries for j-faces
        for (i = 0; i < Nx; ++i) {
            for (j = 0; j <= Ny; ++j) {

                int ij = i * (Ny + 1) + j;
                int iij = (i + 1) * (Ny + 1) + j;  

                CD = { x_vertices[iij] - x_vertices[ij], y_vertices[iij] - y_vertices[ij] };

                int face_ij = i * (Ny + 1) + j;

                jAreas[face_ij] = sqrt(CD.x * CD.x + CD.y * CD.y);

                jface_xNormals[face_ij] = -CD.y / fabs(jAreas[face_ij]);
                jface_yNormals[face_ij] = CD.x / fabs(jAreas[face_ij]);


                // cout << "j area: " << fixed << setprecision(3) <<  jAreas[face_ij] 
                // << "\tj-x norm: " << fixed << setprecision(3) << jface_xNormals[face_ij] 
                // << "\tj-y norm: " << fixed << setprecision(3) << jface_yNormals[face_ij] << endl;
            }
        }
    } 

    local_Nx = Vector(size);
    vector<int> cc_sendcounts(size), cc_displacements(size), 
            if_sendcounts(size), if_displacements(size) , 
            jf_sendcounts(size), jf_displacements(size);

    int base = Nx / size;
    int rem = Nx % size;

    for (int r = 0; r < size; ++r) {
        local_Nx[r] = base + (r < rem ? 1 : 0);
    }


    int cc_offset = 0;
    if (rank == 0) {
        for (int r = 0; r < size; ++r) {
            cc_displacements[r] = cc_offset;
            cc_sendcounts[r] = local_Nx[r] * Ny;
            cc_offset += cc_sendcounts[r];
        }
    }

    // Broadcast each rank's local_Nx[rank] so they know their size
    Nx_local = local_Nx[rank]; 
 
    num_ifaces = (Nx_local + 1) * Ny; 
    num_jfaces = (Ny + 1) * Nx_local; 
    num_loc_cells = Nx_local * Ny; 
    num_gathered_cells = Nx * Ny;

    // Scatter cell-centered data
    xCenter = Vector(num_loc_cells);  
    MPI_Scatterv(x_cellCenters.data(), cc_sendcounts.data(), cc_displacements.data(), MPI_DOUBLE,
                 xCenter.data(), num_loc_cells, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    // Scatter yCenters from grid into new vectors
    yCenter = Vector(num_loc_cells); 
    MPI_Scatterv(y_cellCenters.data(), cc_sendcounts.data(), cc_displacements.data(), MPI_DOUBLE,
                 yCenter.data(), num_loc_cells, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    // Scatter volumes form grid into new vectors
    Volume = Vector(num_loc_cells); 
    MPI_Scatterv(cellVolumes.data(), cc_sendcounts.data(), cc_displacements.data(), MPI_DOUBLE,
                 Volume.data(), num_loc_cells, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 


    // For i-faces now
    if (rank == 0) {
        int offset = 0;
        for (int r = 0; r < size; ++r) {
            if_sendcounts[r] = (local_Nx[r] + 1) * Ny;
            if_displacements[r] = offset;
            offset += if_sendcounts[r] - Ny;
        }
    }

    // Scatter i-face normals and areas
    iFxNorm = Vector(num_ifaces, 0.0);
    MPI_Scatterv(iface_xNormals.data(), if_sendcounts.data(), if_displacements.data(), MPI_DOUBLE,
                 iFxNorm.data(), num_ifaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    iFyNorm = Vector(num_ifaces, 0.0);
    MPI_Scatterv(iface_yNormals.data(), if_sendcounts.data(), if_displacements.data(), MPI_DOUBLE,
                 iFyNorm.data(), num_ifaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    iArea = Vector(num_ifaces, 0.0);
    MPI_Scatterv(iAreas.data(), if_sendcounts.data(), if_displacements.data(), MPI_DOUBLE,
                 iArea.data(), num_ifaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);              


    // For j-faces now
    if (rank == 0) {
        int jface_offset = 0;
        for (int r = 0; r < size; ++r) {
            jf_sendcounts[r]    = local_Nx[r] * (Ny + 1); // includes overlap
            jf_displacements[r] = jface_offset;
            jface_offset     += local_Nx[r] * (Ny + 1);
        }
    }

    // Scatter j-face normals and areas
    jFxNorm = Vector(num_jfaces, 0.0);
    MPI_Scatterv(jface_xNormals.data(), jf_sendcounts.data(), jf_displacements.data(), MPI_DOUBLE,
                 jFxNorm.data(), num_jfaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    jFyNorm = Vector(num_jfaces, 0.0);
    MPI_Scatterv(jface_yNormals.data(), jf_sendcounts.data(), jf_displacements.data(), MPI_DOUBLE,
                 jFyNorm.data(), num_jfaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    jArea = Vector(num_jfaces, 0.0);
    MPI_Scatterv(jAreas.data(), jf_sendcounts.data(), jf_displacements.data(), MPI_DOUBLE,
                 jArea.data(), num_jfaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

} 
void Solver2D::create_cylinder_grid(double Cylinder_Radius, double R1, double R2, double d_min, double theta_left, double theta_right) {
    
    x_vertices = Vector( (Nx + 1) * (Ny + 1), 0.0),   
    y_vertices = Vector( (Nx + 1) * (Ny + 1), 0.0),
    x_cellCenters = Vector(Nx * Ny, 0.0),
    y_cellCenters = Vector(Nx * Ny, 0.0),
    iface_xNormals = Vector((Nx + 1) * Ny, 0.0),
    iface_yNormals = Vector((Nx + 1) * Ny, 0.0),
    jface_xNormals = Vector(Nx * (Ny + 1), 0.0),
    jface_yNormals = Vector(Nx * (Ny + 1), 0.0),
    iAreas = Vector((Nx + 1) * Ny, 0.0),
    jAreas = Vector(Nx * (Ny + 1), 0.0),
    cellVolumes = Vector(Nx * Ny, 0.0);

    BCs = BCMap(BCType::Outlet, BCType::Outlet, BCType::Symmetry, BCType::Inlet); 

    if (rank == 0) {

        const int Ntheta = Nx + 1, Nr = Ny + 1;
        double R_max, k1; 

        Vector theta(Ntheta, 0.0);
        Vector r(Ntheta * Nr, 0.0); 
        double dtheta = (theta_right - theta_left) / (Ntheta - 1); 

        for (int i = 0; i < Ntheta; ++i) {
            r[i * Nr] = Cylinder_Radius;
            theta[i] = theta_left + i * dtheta;
        }

        for (int i = 0; i < Ntheta; ++i) {
            R_max = R1 + (R2 - R1) * cos(theta[i]);
            k1 = NewtonMethod(R_max, Nr, d_min);

            for (int j = 0; j < Nr; ++j) {
                r[i * Nr + j] = r[i * Nr] + R_max * ((exp(k1 * j / (Nr - 1)) - 1) / (exp(k1) - 1)); 
            }
        }

        for (int i = 0; i < Nx + 1; ++i) {
            for (int j = 0; j < Ny + 1; ++j) {
                x_vertices[i * (Ny + 1) + j] = r[i * (Ny + 1) + j] * cos(theta[i]);
                y_vertices[i * (Ny + 1) + j] = r[i * (Ny + 1) + j] * sin(theta[i]);
            }
        }
    

        Point AB, BC, CD, DA;
        int i, j;

        // Calculates cell centers and volumes
        for (i = 0; i < Nx; ++i) {
            for (j = 0; j < Ny; ++j) {

                int ij = i * (Ny + 1) + j;
                int iij = (i + 1) * (Ny + 1) + j;   
                int ijj = i * (Ny + 1) + j + 1; 
                int iijj = (i + 1) * (Ny + 1) + j + 1;   

                DA = { x_vertices[ij] - x_vertices[ijj], y_vertices[ij] - y_vertices[ijj] }; 
                AB = { x_vertices[iij] - x_vertices[ij], y_vertices[iij] - y_vertices[ij] };
                BC = { x_vertices[iijj] - x_vertices[iij], y_vertices[iijj] - y_vertices[iij] };
                CD = { x_vertices[ijj] - x_vertices[iijj], y_vertices[ijj] - y_vertices[iijj] }; 

                int cell_ij = i * Ny + j;   

                x_cellCenters[cell_ij] = (x_vertices[ij] + x_vertices[iij] + x_vertices[iijj] + x_vertices[ijj]) / 4;
                y_cellCenters[cell_ij] =	(y_vertices[ij] + y_vertices[iij] + y_vertices[iijj] + y_vertices[ijj]) / 4;

                cellVolumes[cell_ij] = 0.5 * fabs(DA.x * AB.y - DA.y * AB.x) + 0.5 * fabs(BC.x * CD.y - BC.y * CD.x);

                // cout << "xcenter: " << fixed << setprecision(3) << x_cellCenters[cell_ij] 
                // << "\tycenter: " << fixed << setprecision(3) << y_cellCenters[cell_ij] 
                // << "\tcell volume: " << fixed << setprecision(3) << cellVolumes[cell_ij] << endl; 
            }
        }

        // Calculates geometries for i-faces
        for (i = 0; i <= Nx; ++i) {
            for (j = 0; j < Ny; ++j) {

                int ij = i * (Ny + 1) + j;
                int ijj = i * (Ny + 1) + j + 1;   

                int face_ij = i * Ny + j;  

                AB = { x_vertices[ijj] - x_vertices[ij], y_vertices[ijj] - y_vertices[ij] };

                iAreas[face_ij] = sqrt(AB.x * AB.x + AB.y * AB.y); 

                iface_xNormals[face_ij] = AB.y / fabs(iAreas[face_ij]);
                iface_yNormals[face_ij] = -AB.x / fabs(iAreas[face_ij]);

                // cout << "i area: " << fixed << setprecision(3) << iAreas[face_ij] 
                // << "\ti-x norm: " << fixed << setprecision(3) << iface_xNormals[face_ij] 
                // << "\ti-y norm: " << fixed << setprecision(3) << iface_yNormals[face_ij] << endl;
            }
        }

        // Calculates geometries for j-faces
        for (i = 0; i < Nx; ++i) {
            for (j = 0; j <= Ny; ++j) {

                int ij = i * (Ny + 1) + j;
                int iij = (i + 1) * (Ny + 1) + j;  

                CD = { x_vertices[iij] - x_vertices[ij], y_vertices[iij] - y_vertices[ij] };

                int face_ij = i * (Ny + 1) + j;

                jAreas[face_ij] = sqrt(CD.x * CD.x + CD.y * CD.y);

                jface_xNormals[face_ij] = -CD.y / fabs(jAreas[face_ij]);
                jface_yNormals[face_ij] = CD.x / fabs(jAreas[face_ij]);


                // cout << "j area: " << fixed << setprecision(3) <<  jAreas[face_ij] 
                // << "\tj-x norm: " << fixed << setprecision(3) << jface_xNormals[face_ij] 
                // << "\tj-y norm: " << fixed << setprecision(3) << jface_yNormals[face_ij] << endl;
            }
        }
    }

    local_Nx = Vector(size);
    vector<int> cc_sendcounts(size), cc_displacements(size), 
            if_sendcounts(size), if_displacements(size) , 
            jf_sendcounts(size), jf_displacements(size);

    int base = Nx / size;
    int rem = Nx % size;

    for (int r = 0; r < size; ++r) {
        local_Nx[r] = base + (r < rem ? 1 : 0);
    }


    int cc_offset = 0;
    if (rank == 0) {
        for (int r = 0; r < size; ++r) {
            cc_displacements[r] = cc_offset;
            cc_sendcounts[r] = local_Nx[r] * Ny;
            cc_offset += cc_sendcounts[r];
        }
    }

    // Broadcast each rank's local_Nx[rank] so they know their size
    Nx_local = local_Nx[rank]; 
 
    num_ifaces = (Nx_local + 1) * Ny; 
    num_jfaces = (Ny + 1) * Nx_local; 
    num_loc_cells = Nx_local * Ny; 
    num_gathered_cells = Nx * Ny;

    // Scatter cell-centered data
    xCenter = Vector(num_loc_cells);  
    MPI_Scatterv(x_cellCenters.data(), cc_sendcounts.data(), cc_displacements.data(), MPI_DOUBLE,
                 xCenter.data(), num_loc_cells, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    // Scatter yCenters from grid into new vectors
    yCenter = Vector(num_loc_cells); 
    MPI_Scatterv(y_cellCenters.data(), cc_sendcounts.data(), cc_displacements.data(), MPI_DOUBLE,
                 yCenter.data(), num_loc_cells, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    // Scatter volumes form grid into new vectors
    Volume = Vector(num_loc_cells); 
    MPI_Scatterv(cellVolumes.data(), cc_sendcounts.data(), cc_displacements.data(), MPI_DOUBLE,
                 Volume.data(), num_loc_cells, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 


    // For i-faces now
    if (rank == 0) {
        int offset = 0;
        for (int r = 0; r < size; ++r) {
            if_sendcounts[r] = (local_Nx[r] + 1) * Ny;
            if_displacements[r] = offset;
            offset += if_sendcounts[r] - Ny;
        }
    }

    // Scatter i-face normals and areas
    iFxNorm = Vector(num_ifaces, 0.0);
    MPI_Scatterv(iface_xNormals.data(), if_sendcounts.data(), if_displacements.data(), MPI_DOUBLE,
                 iFxNorm.data(), num_ifaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    iFyNorm = Vector(num_ifaces, 0.0);
    MPI_Scatterv(iface_yNormals.data(), if_sendcounts.data(), if_displacements.data(), MPI_DOUBLE,
                 iFyNorm.data(), num_ifaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    iArea = Vector(num_ifaces, 0.0);
    MPI_Scatterv(iAreas.data(), if_sendcounts.data(), if_displacements.data(), MPI_DOUBLE,
                 iArea.data(), num_ifaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);              


    // For j-faces now
    if (rank == 0) {
        int jface_offset = 0;
        for (int r = 0; r < size; ++r) {
            jf_sendcounts[r]    = local_Nx[r] * (Ny + 1); // includes overlap
            jf_displacements[r] = jface_offset;
            jface_offset     += local_Nx[r] * (Ny + 1);
        }
    }

    // Scatter j-face normals and areas
    jFxNorm = Vector(num_jfaces, 0.0);
    MPI_Scatterv(jface_xNormals.data(), jf_sendcounts.data(), jf_displacements.data(), MPI_DOUBLE,
                 jFxNorm.data(), num_jfaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    jFyNorm = Vector(num_jfaces, 0.0);
    MPI_Scatterv(jface_yNormals.data(), jf_sendcounts.data(), jf_displacements.data(), MPI_DOUBLE,
                 jFyNorm.data(), num_jfaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    jArea = Vector(num_jfaces, 0.0);
    MPI_Scatterv(jAreas.data(), jf_sendcounts.data(), jf_displacements.data(), MPI_DOUBLE,
                 jArea.data(), num_jfaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

}


void Solver2D::compute_ifluxes() {

    int ci, cii, fi; 
    double weight, pi, pii, dp, nx, ny, rho, u, v, p, ho, uprime, a, lp, lm, l, lt, lc;  
    ThermoEntry CTi, CTii, CTp, CTm;
    double g = 5.72;
    double pe, pp, up_half, p_half, rho_half, dpdrho_half; 
    Vector U_half(4), Ui(4), Uii(4);  
    // cout << "I FLUXES:" << endl; 

    for (int i = 0; i < Nx_local + 1; ++i) {
        for (int j = 0; j < Ny; ++j) { 
            
            ci = i * (Ny + 2) + j + 1;
            cii = (i + 1) * (Ny + 2) + j + 1; 
            fi = i * Ny + j;

            nx = iFxNorm[fi];
            ny = iFyNorm[fi]; 

            for (int k = 0; k < n; ++k) {
                Ui[k] = U[ci * n + k];
                Uii[k] = U[cii * n + k];
                U_half[k] = 0.5 * (Ui[k] + Uii[k]);
            }

            CTi = cell_thermo[ci];
            CTii = cell_thermo[cii];

            constoprim(Ui.data(), Vi.data(), CTi.gamma, n_vel);
            constoprim(Uii.data(), Vii.data(), CTii.gamma, n_vel);  

            pi = CTi.p;
            pii = CTii.p; 
            dp = fabs(pii - pi)/ min(pi, pii);
            weight = 1 - 0.5 / ((g * g * dp * dp) + 1);

            CTp = weight * CTi + (1 - weight) * CTii;
            CTm = (1 - weight) * CTi + weight * CTii;

            for (int k = 0; k < n; ++k) {
                Vp[k] = weight * Vi[k] + (1 - weight) * Vii[k];
                Vm[k] = (1 - weight) * Vi[k] + weight * Vii[k];
            }

            // Compute i_rhoA 
            if (real_gas) {
                up_half = (0.5 * Ui[1] / Ui[0] + 0.5 * Uii[1] / Uii[0]) * nx + (0.5 * Ui[2] / Ui[0] + 0.5 * Uii[2] / Uii[0]) * ny;
                p_half = 0.5 * (pi + pii); 
                rho_half = 0.5 * (Ui[0] + Uii[0]);
                dpdrho_half = 0.5 * (CTi.dpdrho + CTii.dpdrho);
                
                irho_A[fi * n * n + 4] = (dpdrho_half - p_half / rho_half) * nx;
                irho_A[fi * n * n + 8] = (dpdrho_half - p_half / rho_half) * ny;
                irho_A[fi * n * n + 12] = (dpdrho_half - p_half / rho_half) * up_half;

                matvec_mult(&irho_A[fi * n * n], U_half.data(), &irho_flux[fi * n], n);
            }

            // Positive flux calculation

            rho = Vp[0];
            u = Vp[1];
            v = Vp[2];
            V[3] = CTp.p;
            p = CTp.p;  
            a = CTp.a;   
            ho = compute_total_enthalpy(Vp.data(), CTp.gamma, n_vel);
            pp = CTp.dpdrho - CTp.e / rho * CTp.dpde + 0.5 * (Vp[1] * Vp[1] + Vp[2] * Vp[2]) / rho * CTp.dpde; 
            pe = 1 / rho * CTp.dpde; 

            uprime = u * nx + v * ny; 

            lp = 0.5 * (uprime + a + fabs(uprime + a));
            lm = 0.5 * (uprime - a + fabs(uprime - a));
            l = 0.5 * (uprime + fabs(uprime));
            lt = 0.5 * (lp - lm);
            lc = 0.5 * (lp + lm - 2 * l);

            fill(int1.begin(), int1.end(), 0.0);
            for (int k = 0; k < n; ++k) int1[k * n + k] = l;

            V1 = {lc / (a * a), (u * lc + a * nx * lt)/(a * a), (v * lc + a * ny * lt)/(a * a), (ho * lc + a * uprime * lt) / (a * a)};
            Q = {pp, -u * pe, -v * pe, pe};

            V2 = {lt / a, u * lt / a + nx * lc, v * lt / a + ny * lc, ho * lt / a + uprime * lc};
            W = {-uprime, nx, ny, 0}; 

            outer_product(V1.data(), Q.data(), int2.data(), n);
            outer_product(V2.data(), W.data(), int3.data(), n); 

            for (int k = 0; k < n * n; ++k) 
                iPlus_A[fi * n * n + k]  = int1[k] + int2[k] + int3[k];

            matvec_mult(&iPlus_A[fi * n * n], Ui.data(), F_plus.data(), n);  

             // Negative flux calculation

            rho = Vm[0];
            u = Vm[1];
            v = Vm[2];
            Vm[3] = CTm.p;
            p = CTm.p; 
            a = CTm.a;   
            ho = compute_total_enthalpy(Vm.data(), CTm.gamma, n_vel);
            pp = CTm.dpdrho - CTm.e / rho * CTm.dpde + 0.5 * (Vm[1] * Vm[1] + Vm[2] * Vm[2]) / rho * CTm.dpde; 
            pe = 1 / rho * CTm.dpde; 

            uprime = u * nx + v * ny; 

            lp = 0.5 * (uprime + a - fabs(uprime + a));
            lm = 0.5 * (uprime - a - fabs(uprime - a));
            l = 0.5 * (uprime - fabs(uprime));
            lt = 0.5 * (lp - lm);
            lc = 0.5 * (lp + lm - 2 * l);

            fill(int1.begin(), int1.end(), 0.0);
            for (int k = 0; k < n; ++k) int1[k * n + k] = l;

            V1 = {lc / (a * a), (u * lc + a * nx * lt)/(a * a), (v * lc + a * ny * lt)/(a * a), (ho * lc + a * uprime * lt) / (a * a)};
            Q = {pp, -u * pe, -v * pe, pe};

            V2 = {lt / a, u * lt / a + nx * lc, v * lt / a + ny * lc, ho * lt / a + uprime * lc};
            W = {-uprime, nx, ny, 0}; 

            outer_product(V1.data(), Q.data(), int2.data(), n);
            outer_product(V2.data(), W.data(), int3.data(), n); 

            for (int k = 0; k < n * n; ++k) 
                iMinus_A[fi * n * n + k]  = int1[k] + int2[k] + int3[k]; 

            matvec_mult(&iMinus_A[fi * n * n], Uii.data(), F_minus.data(), n);  

            // cout << i << ", " << j << ": ";
            for (int k = 0; k < n; ++k) {
                iFlux[fi * n + k] = F_plus[k] + F_minus[k]; 
                // cout << iFlux[fi * n + k] << " ";
            }
            // cout << "\t nx = " << nx << " ny = " << ny << endl;

        }
    }
}
void Solver2D::compute_jfluxes() {

    int cj, cjj, fj; 
    double weight, pj, pjj, dp, nx, ny, rho, u, v, p, ho, uprime, a, lp, lm, l, lt, lc;  
    Vector Uj(4), Ujj(4), U_half(4);
    ThermoEntry CTj, CTjj, CTp, CTm;
    double pe, pp, up_half, p_half, rho_half, dpdrho_half; 
    double g = 5.72;
    // cout << "J FLUXES: " << endl; 

    for (int i = 0; i < Nx_local; ++i) {
        for (int j = 0; j < Ny + 1; ++j) { 
            
            cj = (i + 1) * (Ny + 2) + j;
            cjj = (i + 1) * (Ny + 2) + j + 1; 
            fj = i * (Ny + 1) + j;

            nx = jFxNorm[fj];
            ny = jFyNorm[fj];  

            for (int k = 0; k < n; ++k) {
                Uj[k] = U[cj * n + k];
                Ujj[k] = U[cjj * n + k];
                U_half[k] = 0.5 * (Uj[k] + Ujj[k]);
            }

            CTj = cell_thermo[cj];
            CTjj = cell_thermo[cjj];

            constoprim(&U[cj * n], Vj.data(), CTj.gamma, n_vel);
            constoprim(&U[cjj * n], Vjj.data(), CTjj.gamma, n_vel);   

            pj = CTj.p;
            pjj = CTjj.p; 
            dp = fabs(pjj - pj)/ min(pj, pjj); 
            weight = 1 - 0.5 / ((g * g * dp * dp) + 1);

            CTp = weight * CTj + (1 - weight) * CTjj;
            CTm = (1 - weight) * CTj + weight * CTjj;

            for (int k = 0; k < n; ++k) {
                Vp[k] = weight * Vj[k] + (1 - weight) * Vjj[k];
                Vm[k] = (1 - weight) * Vj[k] + weight * Vjj[k];
            }

            // Compute i_rhoA 
            if (real_gas) {

                up_half = (0.5 * Uj[1] / Uj[0] + 0.5 * Ujj[1] / Ujj[0]) * nx + (0.5 * Uj[2] / Uj[0] + 0.5 * Ujj[2] / Ujj[0]) * ny;
                p_half = 0.5 * (pj + pjj); 
                rho_half = 0.5 * (Uj[0] + Ujj[0]);
                dpdrho_half = 0.5 * (CTj.dpdrho + CTjj.dpdrho); 
                
                jrho_A[fj * n * n + 4] = (dpdrho_half - p_half / rho_half) * nx;
                jrho_A[fj * n * n + 8] = (dpdrho_half - p_half / rho_half) * ny;
                jrho_A[fj * n * n + 12] = (dpdrho_half - p_half / rho_half) * up_half;

                matvec_mult(&jrho_A[fj * n * n], U_half.data(), &jrho_flux[fj * n], n);
            }

            // Positive flux calculation

            rho = Vp[0];
            u = Vp[1];
            v = Vp[2];
            V[3] = CTp.p;
            p = Vp[3]; 
            a = CTp.a;   
            ho = compute_total_enthalpy(Vp.data(), CTp.gamma, n_vel);
            pp = CTp.dpdrho - CTp.e / rho * CTp.dpde + 0.5 * (Vp[1] * Vp[1] + Vp[2] * Vp[2]) / rho * CTp.dpde; 
            pe = 1 / rho * CTp.dpde; 

            uprime = u * nx + v * ny; 

            lp = 0.5 * (uprime + a + fabs(uprime + a));
            lm = 0.5 * (uprime - a + fabs(uprime - a));
            l = 0.5 * (uprime + fabs(uprime));
            lt = 0.5 * (lp - lm);
            lc = 0.5 * (lp + lm - 2 * l);

            fill(int1.begin(), int1.end(), 0.0);
            for (int k = 0; k < n; ++k) int1[k * n + k] = l;

            V1 = {lc / (a * a), (u * lc + a * nx * lt)/(a * a), (v * lc + a * ny * lt)/(a * a), (ho * lc + a * uprime * lt) / (a * a)};
            Q = {pp, -u * pe, -v * pe, pe};

            V2 = {lt / a, u * lt / a + nx * lc, v * lt / a + ny * lc, ho * lt / a + uprime * lc};
            W = {-uprime, nx, ny, 0}; 

            outer_product(V1.data(), Q.data(), int2.data(), n);
            outer_product(V2.data(), W.data(), int3.data(), n); 

            for (int k = 0; k < n * n; ++k) 
                jPlus_A[fj * n * n + k]  = int1[k] + int2[k] + int3[k];

            matvec_mult(&jPlus_A[fj * n * n], Uj.data(), F_plus.data(), n);  

             // Negative flux calculation

            rho = Vm[0];
            u = Vm[1];
            v = Vm[2];
            V[3] = CTm.p;
            p = Vm[3]; 
            a = CTm.a;   
            ho = compute_total_enthalpy(Vm.data(), CTm.gamma, n_vel);
            pp = CTm.dpdrho - CTm.e / rho * CTm.dpde + 0.5 * (Vm[1] * Vm[1] + Vm[2] * Vm[2]) / rho * CTm.dpde; 
            pe = 1 / rho * CTm.dpde; 

            uprime = u * nx + v * ny; 

            lp = 0.5 * (uprime + a - fabs(uprime + a));
            lm = 0.5 * (uprime - a - fabs(uprime - a));
            l = 0.5 * (uprime - fabs(uprime));
            lt = 0.5 * (lp - lm);
            lc = 0.5 * (lp + lm - 2 * l);

            fill(int1.begin(), int1.end(), 0.0);
            for (int k = 0; k < n; ++k) int1[k * n + k] = l;

            V1 = {lc / (a * a), (u * lc + a * nx * lt)/(a * a), (v * lc + a * ny * lt)/(a * a), (ho * lc + a * uprime * lt) / (a * a)};
            Q = {pp, -u * pe, -v * pe, pe};

            V2 = {lt / a, u * lt / a + nx * lc, v * lt / a + ny * lc, ho * lt / a + uprime * lc};
            W = {-uprime, nx, ny, 0}; 

            outer_product(V1.data(), Q.data(), int2.data(), n);
            outer_product(V2.data(), W.data(), int3.data(), n); 

            for (int k = 0; k < n * n; ++k) 
                jMinus_A[fj * n * n + k]  = int1[k] + int2[k] + int3[k]; 

            matvec_mult(&jMinus_A[fj * n * n], Ujj.data(), F_minus.data(), n);  

            // cout << i << ", " << j << ": ";
            for (int k = 0; k < n; ++k) {
                jFlux[fj * n + k] = F_plus[k] + F_minus[k];  
                // cout << jFlux[fj * n + k] << " ";
            }
            // cout << "\t nx = " << nx << " ny = " << ny << endl; 
        }
    }
}

void Solver2D::perf_line_relaxation() {
    int in_counter = 0;

    while (in_counter < 8) {

        relax_left_line_perf();  
        
        if (rank > 0 && rank < size - 1) {
            for (int i = 0; i < Nx_local; ++i) {
                relax_inner_lines_perf(i);
            }
        }

        relax_right_line_perf();

        for (int i = 0; i < Nx_local + 2; ++i) {
            for (int j = 0; j < Ny + 2; ++j) {
                int idx = (i * (Ny + 2) + j) * n;
                for (int k = 0; k < n; ++k) {
                    dU_old[idx + k] = dU_new[idx + k];
                }
            }
        }

        exchange_dU_ghost_cells(); 

        compute_inner_res_perf();

        in_counter++;
    }
}
void Solver2D::relax_left_line_perf() {

    if (rank == 0) {

        // ====================== Top cell ====================== // 
        int i = 0, j = Ny - 1;
        int LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
            BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
            CC = i * Ny + j;
           
        matmat_mult(&iEL[j * n * n], &iPlus_A[LF * n * n], rm1.data(), n);
        matmat_mult(&jET[i * n * n], &jMinus_A[TF * n * n], rm2.data(), n);
        matvec_mult(&iMinus_A[RF * n * n], &dU_old[((i + 2) * (Ny + 2) + j + 1) * n], rv1.data(), n); 
    
        for (int k = 0; k < n * n; ++k) {
            A[k] = 
            Volume[CC] / dt * I[k]          // Time-derivative part
            - (iMinus_A[LF * n * n + k] + rm1[k] ) * iArea[LF] // Left-face contribution
            + iPlus_A[RF * n * n + k] * iArea[RF]   // Right-face contribution
            - jMinus_A[BF * n * n + k] * jArea[BF]  // Bottom-face contribution
            + (jPlus_A[TF * n * n + k] + rm2[k]) * jArea[TF]; // Top-face contribution 


            C[k] = -jPlus_A[BF * n * n + k] * jArea[BF];
        }

        for (int k = 0; k < n; ++k) 
            F[k] = iFlux[LF * n + k] * iArea[LF]
                - iFlux[RF * n + k] * iArea[RF]
                + jFlux[BF * n + k] * jArea[BF]
                - jFlux[TF * n + k] * jArea[TF]
                - rv1[k] * iArea[RF];
        
        alpha = A;
        matrix_divide(alpha.data(), F.data(), &v[j * n], n, 1);  //n = rows of x, m = columns of x
        matrix_divide(alpha.data(), C.data(), &g[j * n * n], n, n); 

        // ======================= Inner cells ======================== //
        for (int j = Ny - 2; j > 0; --j) {

            LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
            BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
            CC = i * Ny + j;

            matmat_mult(&iEL[j * n * n], &iPlus_A[LF * n * n], rm1.data(), n);
            matvec_mult(&iMinus_A[RF * n * n], &dU_old[((i + 2) * (Ny + 2) + j + 1) * n], rv1.data(), n);

            for (int k =  0; k < n * n; ++k) {
                B[k] = jMinus_A[TF * n * n + k] * jArea[TF]; 

                A[k] = 
                    Volume[CC] / dt * I[k]          // Time-derivative part
                    - (iMinus_A[LF * n * n + k] + rm1[k] ) * iArea[LF] // Left-face contribution
                    + iPlus_A[RF * n * n + k] * iArea[RF]   // Right-face contribution
                    - jMinus_A[BF * n * n + k] * jArea[BF]  // Bottom-face contribution
                    + jPlus_A[TF * n * n + k] * jArea[TF]; // Top-face contribution

                C[k] = -jPlus_A[BF * n * n + k] * jArea[BF];
            }

            for (int k = 0; k < n; ++k) 
                F[k] = iFlux[LF * n + k] * iArea[LF]
                    - iFlux[RF * n + k] * iArea[RF]
                    + jFlux[BF * n + k] * jArea[BF]
                    - jFlux[TF * n + k] * jArea[TF]
                    - rv1[k] * iArea[RF];

            // Form alpha
            matmat_mult(B.data(), &g[(j + 1) * n * n], rm1.data(), n);
            for (int k = 0; k < n * n; ++k) 
                alpha[k] = A[k] - rm1[k];
        
            // Form g
            matrix_divide(alpha.data(), C.data(), &g[j * n * n], n, n); 

            // Form v
            matvec_mult(B.data(), &v[(j + 1) * n], rv1.data(), n); 
            for (int k = 0; k < n; ++k) 
                rv2[k] = F[k] - rv1[k];        
            matrix_divide(alpha.data(), rv2.data(), &v[j * n], n, 1);
        }
 
        // ==================== Bottom cell ===================== //
        j = 0;

        LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
        BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
        CC = i * Ny + j;

        matmat_mult(&iEL[j * n * n], &iPlus_A[LF * n * n], rm1.data(), n);
        matmat_mult(&jEB[i * n * n], &jPlus_A[BF * n * n], rm2.data(), n);
        matvec_mult(&iMinus_A[RF * n * n], &dU_old[((i + 2) * (Ny + 2) + j + 1) * n], rv1.data(), n);

        for (int k = 0; k < n * n; ++k) {

            B[k] = jMinus_A[TF * n * n + k] * jArea[TF]; 

            A[k] = 
            Volume[CC] / dt * I[k]          // Time-derivative part
            - (iMinus_A[LF * n * n + k] + rm1[k] ) * iArea[LF] // Left-face contribution
            + iPlus_A[RF * n * n + k] * iArea[RF]   // Right-face contribution
            - (jMinus_A[BF * n * n + k] + rm2[k]) * jArea[BF]  // Bottom-face contribution
            + jPlus_A[TF * n * n + k] * jArea[TF]; // Top-face contribution
        }

        for (int k = 0; k < n; ++k) 
            F[k] = iFlux[LF * n + k] * iArea[LF]
                - iFlux[RF * n + k] * iArea[RF]
                + jFlux[BF * n + k] * jArea[BF]
                - jFlux[TF * n + k] * jArea[TF]
                - rv1[k] * iArea[RF];
        
        // Form alpha
        matmat_mult(B.data(), &g[(j + 1) * n * n], rm1.data(), n);
        for (int k = 0; k < n * n; ++k) 
            alpha[k] = A[k] - rm1[k];

        // Form v
        matvec_mult(B.data(), &v[(j + 1) * n], rv1.data(), n);
        for (int k = 0; k < n; ++k) 
            rv2[k] = F[k] - rv1[k];
        matrix_divide(alpha.data(), rv2.data(), &v[j * n], n, 1); 

        // ====================== Solve for dU_new ====================== //
        for (int k = 0; k < n; ++k)
            dU_new[((i + 1) * (Ny + 2) + j + 1) * n + k] = v[0 * n + k];

        for (int j = 1; j < Ny; ++j) {
            matvec_mult(&g[j * n * n], &dU_new[((i + 1) * (Ny + 2) + j) * n], rv1.data(), n);
            for (int k = 0; k < n; ++k)
                dU_new[((i + 1) * (Ny + 2) + j + 1) * n + k] = v[j * n + k] - rv1[k];
            
        }


        for (int i = 1; i < Nx_local - 1; ++i) {
            relax_inner_lines_perf(i); 
        }
    }
}
void Solver2D::relax_inner_lines_perf(int i) { 

    // ====================== Top cell ====================== // 
    int j = Ny - 1;
    int LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
        BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
        CC = i * Ny + j;

    matmat_mult(&jET[i * n * n], &jMinus_A[TF * n * n], rm1.data(), n);
    matvec_mult(&iPlus_A[LF * n * n], &dU_old[(i * (Ny + 2) + j + 1) * n], rv1.data(), n);
    matvec_mult(&iMinus_A[RF * n * n], &dU_old[((i + 2) * (Ny + 2) + j + 1) * n], rv2.data(), n); 

    for (int k = 0; k < n * n; ++k) {
        A[k] = 
        Volume[CC] / dt * I[k]          // Time-derivative part
        - iMinus_A[LF * n * n + k] * iArea[LF] // Left-face contribution
        + iPlus_A[RF * n * n + k] * iArea[RF]   // Right-face contribution
        - jMinus_A[BF * n * n + k] * jArea[BF]  // Bottom-face contribution
        + (jPlus_A[TF * n * n + k] + rm1[k]) * jArea[TF]; // Top-face contribution

        C[k] = -jPlus_A[BF * n * n + k] * jArea[BF];
    }

    for (int k = 0; k < n; ++k) 
        F[k] = iFlux[LF * n + k] * iArea[LF]
            - iFlux[RF * n + k] * iArea[RF]
            + jFlux[BF * n + k] * jArea[BF]
            - jFlux[TF * n + k] * jArea[TF]
            + rv1[k] * iArea[LF]
            - rv2[k] * iArea[RF];
    
    alpha = A;
    matrix_divide(alpha.data(), F.data(), &v[j * n], n, 1);  //n = rows of x, m = columns of x
    matrix_divide(alpha.data(), C.data(), &g[j * n * n], n, n); 

    // ======================= Inner cells ======================== //
    for (int j = Ny - 2; j > 0; --j) {

        LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
        BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
        CC = i * Ny + j;

        matvec_mult(&iPlus_A[LF * n * n], &dU_old[(i * (Ny + 2) + j + 1) * n], rv1.data(), n);
        matvec_mult(&iMinus_A[RF * n * n], &dU_old[((i + 2) * (Ny + 2) + j + 1) * n], rv2.data(), n); 

        for (int k =  0; k < n * n; ++k) {
            B[k] = jMinus_A[TF * n * n + k] * jArea[TF]; 

            A[k] = 
                Volume[CC] / dt * I[k]          // Time-derivative part
                - iMinus_A[LF * n * n + k] * iArea[LF] // Left-face contribution
                + iPlus_A[RF * n * n + k] * iArea[RF]   // Right-face contribution
                - jMinus_A[BF * n * n + k] * jArea[BF]  // Bottom-face contribution
                + jPlus_A[TF * n * n + k] * jArea[TF]; // Top-face contribution

            C[k] = -jPlus_A[BF * n * n + k] * jArea[BF];
        }

        for (int k = 0; k < n; ++k) 
            F[k] = iFlux[LF * n + k] * iArea[LF]
                - iFlux[RF * n + k] * iArea[RF]
                + jFlux[BF * n + k] * jArea[BF]
                - jFlux[TF * n + k] * jArea[TF]
                + rv1[k] * iArea[LF]
                - rv2[k] * iArea[RF];

        // Form alpha
        matmat_mult(B.data(), &g[(j + 1) * n * n], rm1.data(), n);
        for (int k = 0; k < n * n; ++k) 
            alpha[k] = A[k] - rm1[k];
    
        // Form g
        matrix_divide(alpha.data(), C.data(), &g[j * n * n], n, n); 

        // Form v
        matvec_mult(B.data(), &v[(j + 1) * n], rv1.data(), n); 
        for (int k = 0; k < n; ++k) 
            rv2[k] = F[k] - rv1[k];        
        matrix_divide(alpha.data(), rv2.data(), &v[j * n], n, 1);
    }
    // ==================== Bottom cell ===================== //
    j = 0;

    LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
    BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
    CC = i * Ny + j;

    matmat_mult(&jEB[i * n * n], &jPlus_A[BF * n * n], rm1.data(), n);
    matvec_mult(&iPlus_A[LF * n * n], &dU_old[(i * (Ny + 2) + j + 1) * n], rv1.data(), n);
    matvec_mult(&iMinus_A[RF * n * n], &dU_old[((i + 2) * (Ny + 2) + j + 1) * n], rv2.data(), n); 

    for (int k = 0; k < n * n; ++k) {

        B[k] = jMinus_A[TF * n * n + k] * jArea[TF]; 

        A[k] = 
        Volume[CC] / dt * I[k]          // Time-derivative part
        - iMinus_A[LF * n * n + k] * iArea[LF] // Left-face contribution
        + iPlus_A[RF * n * n + k] * iArea[RF]   // Right-face contribution
        - (jMinus_A[BF * n * n + k] + rm1[k]) * jArea[BF]  // Bottom-face contribution
        + jPlus_A[TF * n * n + k] * jArea[TF]; // Top-face contribution
    }

    for (int k = 0; k < n; ++k) 
        F[k] = iFlux[LF * n + k] * iArea[LF]
            - iFlux[RF * n + k] * iArea[RF]
            + jFlux[BF * n + k] * jArea[BF]
            - jFlux[TF * n + k] * jArea[TF]
            + rv1[k] * iArea[LF]
            - rv2[k] * iArea[RF];
    
    // Form alpha
    matmat_mult(B.data(), &g[(j + 1) * n * n], rm1.data(), n);
    for (int k = 0; k < n * n; ++k) 
        alpha[k] = A[k] - rm1[k];

    // Form v
    matvec_mult(B.data(), &v[(j + 1) * n], rv1.data(), n);
    for (int k = 0; k < n; ++k) 
        rv2[k] = F[k] - rv1[k];
    matrix_divide(alpha.data(), rv2.data(), &v[j * n], n, 1); 

    // ====================== Solve for dU_new ====================== //
    for (int k = 0; k < n; ++k)
        dU_new[((i + 1) * (Ny + 2) + j + 1) * n + k] = v[0 * n + k];

    for (int j = 1; j < Ny; ++j) {
        matvec_mult(&g[j * n * n], &dU_new[((i + 1) * (Ny + 2) + j) * n], rv1.data(), n);
        for (int k = 0; k < n; ++k) 
            dU_new[((i + 1) * (Ny + 2) + j + 1) * n + k] = v[j * n + k] - rv1[k];
    }
}
void Solver2D::relax_right_line_perf() {


    if (rank == size - 1) {

        for (int i = 0; i < Nx_local - 1; ++i) {
            relax_inner_lines_perf(i); 
        }

        // ====================== Top cell ====================== // 
        int i = Nx_local - 1, j = Ny - 1; 
        int LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
            BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
            CC = i * Ny + j;

        matmat_mult(&iER[j * n * n], &iMinus_A[RF * n * n], rm1.data(), n);
        matmat_mult(&jET[i * n * n], &jMinus_A[TF * n * n], rm2.data(), n);
        matvec_mult(&iPlus_A[LF * n * n], &dU_old[(i * (Ny + 2) + j + 1) * n], rv1.data(), n);
    

        for (int k = 0; k < n * n; ++k) {
            A[k] = 
            Volume[CC] / dt * I[k]          // Time-derivative part
            - iMinus_A[LF * n * n + k] * iArea[LF] // Left-face contribution
            + (iPlus_A[RF * n * n + k] + rm1[k]) * iArea[RF]   // Right-face contribution
            - jMinus_A[BF * n * n + k] * jArea[BF]  // Bottom-face contribution
            + (jPlus_A[TF * n * n + k] + rm2[k]) * jArea[TF]; // Top-face contribution

            C[k] = -jPlus_A[BF * n * n + k] * jArea[BF];
        }

        for (int k = 0; k < n; ++k) 
            F[k] = iFlux[LF * n + k] * iArea[LF]
                - iFlux[RF * n + k] * iArea[RF]
                + jFlux[BF * n + k] * jArea[BF]
                - jFlux[TF * n + k] * jArea[TF]
                + rv1[k] * iArea[LF];
        
        alpha = A;
        matrix_divide(alpha.data(), F.data(), &v[j * n], n, 1);  //n = rows of x, m = columns of x
        matrix_divide(alpha.data(), C.data(), &g[j * n * n], n, n); 


        // ======================= Inner cells ======================== //
        for (int j = Ny - 2; j > 0; --j) {

            LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
            BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
            CC = i * Ny + j;

            matmat_mult(&iER[j * n * n], &iMinus_A[RF * n * n], rm1.data(), n);
            matvec_mult(&iPlus_A[LF * n * n], &dU_old[(i * (Ny + 2) + j + 1) * n], rv1.data(), n);            

            for (int k =  0; k < n * n; ++k) {
                B[k] = jMinus_A[TF * n * n + k] * jArea[TF]; 

                A[k] = 
                    Volume[CC] / dt * I[k]          // Time-derivative part
                    - iMinus_A[LF * n * n + k] * iArea[LF] // Left-face contribution
                    + (iPlus_A[RF * n * n + k] + rm1[k]) * iArea[RF]   // Right-face contribution
                    - jMinus_A[BF * n * n + k] * jArea[BF]  // Bottom-face contribution
                    + jPlus_A[TF * n * n + k] * jArea[TF]; // Top-face contribution

                C[k] = -jPlus_A[BF * n * n + k] * jArea[BF];
            }

            for (int k = 0; k < n; ++k) 
                F[k] = iFlux[LF * n + k] * iArea[LF]
                    - iFlux[RF * n + k] * iArea[RF]
                    + jFlux[BF * n + k] * jArea[BF]
                    - jFlux[TF * n + k] * jArea[TF]
                    + rv1[k] * iArea[LF];

            // Form alpha
            matmat_mult(B.data(), &g[(j + 1) * n * n], rm1.data(), n);
            for (int k = 0; k < n * n; ++k) 
                alpha[k] = A[k] - rm1[k];
        
            // Form g
            matrix_divide(alpha.data(), C.data(), &g[j * n * n], n, n); 

            // Form v
            matvec_mult(B.data(), &v[(j + 1) * n], rv1.data(), n); 
            for (int k = 0; k < n; ++k) 
                rv2[k] = F[k] - rv1[k];        
            matrix_divide(alpha.data(), rv2.data(), &v[j * n], n, 1);
        }

        // ==================== Bottom cell ===================== //
        j = 0;

        LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
        BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
        CC = i * Ny + j;

        matmat_mult(&iER[j * n * n], &iMinus_A[RF * n * n], rm1.data(), n);
        matmat_mult(&jEB[i * n * n], &jPlus_A[BF * n * n], rm2.data(), n);
        matvec_mult(&iPlus_A[LF * n * n], &dU_old[(i * (Ny + 2) + j + 1) * n], rv1.data(), n);

        for (int k = 0; k < n * n; ++k) {

            B[k] = jMinus_A[TF * n * n + k] * jArea[TF]; 

            A[k] = 
            Volume[CC] / dt * I[k]          // Time-derivative part
            - iMinus_A[LF * n * n + k] * iArea[LF] // Left-face contribution
            + (iPlus_A[RF * n * n + k] + rm1[k]) * iArea[RF]   // Right-face contribution
            - (jMinus_A[BF * n * n + k] + rm2[k]) * jArea[BF]  // Bottom-face contribution
            + jPlus_A[TF * n * n + k] * jArea[TF]; // Top-face contribution
        }

        for (int k = 0; k < n; ++k) 
            F[k] = iFlux[LF * n + k] * iArea[LF]
                - iFlux[RF * n + k] * iArea[RF]
                + jFlux[BF * n + k] * jArea[BF]
                - jFlux[TF * n + k] * jArea[TF]
                + rv1[k] * iArea[LF];
        
        // Form alpha
        matmat_mult(B.data(), &g[(j + 1) * n * n], rm1.data(), n);
        for (int k = 0; k < n * n; ++k) 
            alpha[k] = A[k] - rm1[k];

        // Form v
        matvec_mult(B.data(), &v[(j + 1) * n], rv1.data(), n);
        for (int k = 0; k < n; ++k) 
            rv2[k] = F[k] - rv1[k];
        matrix_divide(alpha.data(), rv2.data(), &v[j * n], n, 1); 

        // ====================== Solve for dU_new ====================== //
        for (int k = 0; k < n; ++k)
            dU_new[((i + 1) * (Ny + 2) + j + 1) * n + k] = v[0 * n + k];

        for (int j = 1; j < Ny; ++j) {
            matvec_mult(&g[j * n * n], &dU_new[((i + 1) * (Ny + 2) + j) * n], rv1.data(), n);
            for (int k = 0; k < n; ++k) 
                dU_new[((i + 1) * (Ny + 2) + j + 1) * n + k] = v[j * n + k] - rv1[k];
        }        
    }
}

void Solver2D::real_line_relaxation() {
    int in_counter = 0;

    while (in_counter < 8) {

        relax_left_line_real();  
        
        if (rank > 0 && rank < size - 1) {
            for (int i = 0; i < Nx_local; ++i) {
                relax_inner_lines_real(i);
            }
        }

        relax_right_line_real();

        for (int i = 0; i < Nx_local + 2; ++i) {
            for (int j = 0; j < Ny + 2; ++j) {
                int idx = (i * (Ny + 2) + j) * n;
                for (int k = 0; k < n; ++k) {
                    dU_old[idx + k] = dU_new[idx + k];
                }
            }
        }

        exchange_dU_ghost_cells(); 
        compute_inner_res_real();
        in_counter++;
    }
}
void Solver2D::relax_left_line_real() {

    if (rank == 0) {

        // ====================== Top cell ====================== // 
        int i = 0, j = Ny - 1;
        int LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
            BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
            CC = i * Ny + j;
         
        // Left face
        for (int k = 0; k < n * n; ++ k) 
            rm0[k] = iPlus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k];        

        matmat_mult(&iEL[j * n * n], rm0.data(), rm1.data(), n);

        // Top face
        for (int k = 0; k < n * n; ++k) 
            rm0[k] = jMinus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k];
        
        matmat_mult(&jET[i * n * n], rm0.data(), rm2.data(), n);

        // Right face
        for (int k = 0; k < n * n; ++k) 
            rm0[k] = iMinus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k];

        matvec_mult(rm0.data(), &dU_old[((i + 2) * (Ny + 2) + j + 1) * n], rv1.data(), n); 
    
        for (int k = 0; k < n * n; ++k) {
            A[k] = 
            Volume[CC] / dt * I[k]          // Time-derivative part
            - (iMinus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k] + rm1[k] ) * iArea[LF] // Left-face contribution
            + (iPlus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k] ) * iArea[RF]   // Right-face contribution
            - (jMinus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k] ) * jArea[BF]  // Bottom-face contribution
            + (jPlus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k] + rm2[k] ) * jArea[TF]; // Top-face contribution 


            C[k] = -(jPlus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k] ) * jArea[BF];
        }

        for (int k = 0; k < n; ++k) 
            F[k] = (iFlux[LF * n + k] - irho_flux[LF * n + k] ) * iArea[LF]
                - (iFlux[RF * n + k] - irho_flux[RF * n + k] ) * iArea[RF]
                + (jFlux[BF * n + k] - jrho_flux[BF * n + k] ) * jArea[BF]
                - (jFlux[TF * n + k] - jrho_flux[TF * n + k] ) * jArea[TF]
                - rv1[k] * iArea[RF];
        
        alpha = A;
        matrix_divide(alpha.data(), F.data(), &v[j * n], n, 1);  //n = rows of x, m = columns of x
        matrix_divide(alpha.data(), C.data(), &g[j * n * n], n, n); 

        // ======================= Inner cells ======================== //
        for (int j = Ny - 2; j > 0; --j) {

            LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
            BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
            CC = i * Ny + j;

            // Left face
            for (int k = 0; k < n * n; ++ k) 
                rm0[k] = iPlus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k];            

            matmat_mult(&iEL[j * n * n], rm0.data(), rm1.data(), n);

            // Right face
            for (int k = 0; k < n * n; ++k) 
                rm0[k] = iMinus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k];

            matvec_mult(rm0.data(), &dU_old[((i + 2) * (Ny + 2) + j + 1) * n], rv1.data(), n); 

            for (int k =  0; k < n * n; ++k) {
                B[k] = (jMinus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k] ) * jArea[TF]; 

                A[k] = 
                    Volume[CC] / dt * I[k]          // Time-derivative part
                    - (iMinus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k] + rm1[k] ) * iArea[LF] // Left-face contribution
                    + (iPlus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k] ) * iArea[RF]   // Right-face contribution
                    - (jMinus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k] ) * jArea[BF]  // Bottom-face contribution
                    + (jPlus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k] ) * jArea[TF]; // Top-face contribution

                C[k] = -(jPlus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k] ) * jArea[BF];
            }

            for (int k = 0; k < n; ++k) 
                F[k] = (iFlux[LF * n + k] - irho_flux[LF * n + k] ) * iArea[LF]
                    - (iFlux[RF * n + k]  - irho_flux[RF * n + k] ) * iArea[RF]
                    + (jFlux[BF * n + k] - jrho_flux[BF * n + k] ) * jArea[BF]
                    - (jFlux[TF * n + k] - jrho_flux[TF * n + k] ) * jArea[TF]
                    - rv1[k] * iArea[RF];

            // Form alpha
            matmat_mult(B.data(), &g[(j + 1) * n * n], rm1.data(), n);
            for (int k = 0; k < n * n; ++k) 
                alpha[k] = A[k] - rm1[k];
        
            // Form g
            matrix_divide(alpha.data(), C.data(), &g[j * n * n], n, n); 

            // Form v
            matvec_mult(B.data(), &v[(j + 1) * n], rv1.data(), n); 
            for (int k = 0; k < n; ++k) 
                rv2[k] = F[k] - rv1[k];        
            matrix_divide(alpha.data(), rv2.data(), &v[j * n], n, 1);
        }
 
        // ==================== Bottom cell ===================== //
        j = 0;

        LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
        BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
        CC = i * Ny + j;

        // Left face
        for (int k = 0; k < n * n; ++ k) 
            rm0[k] = iPlus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k];            

        matmat_mult(&iEL[j * n * n], rm0.data(), rm1.data(), n);

        // Bottom face
        for (int k = 0; k < n * n; ++ k) 
            rm0[k] = jPlus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k];            

        matmat_mult(&jEB[i * n * n], rm0.data(), rm2.data(), n);

        // Right face
        for (int k = 0; k < n * n; ++k) 
            rm0[k] = iMinus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k];

        matvec_mult(rm0.data(), &dU_old[((i + 2) * (Ny + 2) + j + 1) * n], rv1.data(), n); 

        for (int k = 0; k < n * n; ++k) {

            B[k] = (jMinus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k] ) * jArea[TF]; 

            A[k] = 
            Volume[CC] / dt * I[k]          // Time-derivative part
            - (iMinus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k] + rm1[k] ) * iArea[LF] // Left-face contribution
            + (iPlus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k] ) * iArea[RF]   // Right-face contribution
            - (jMinus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k] + rm2[k]) * jArea[BF]  // Bottom-face contribution
            + (jPlus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k] ) * jArea[TF]; // Top-face contribution
        }

        for (int k = 0; k < n; ++k) 
            F[k] = (iFlux[LF * n + k] - irho_flux[LF * n + k] ) * iArea[LF]
                - (iFlux[RF * n + k] - irho_flux[RF * n + k] ) * iArea[RF]
                + (jFlux[BF * n + k] - jrho_flux[BF * n + k] ) * jArea[BF]
                - (jFlux[TF * n + k] - jrho_flux[TF * n + k] ) * jArea[TF]
                - rv1[k] * iArea[RF];
        
        // Form alpha
        matmat_mult(B.data(), &g[(j + 1) * n * n], rm1.data(), n);
        for (int k = 0; k < n * n; ++k) 
            alpha[k] = A[k] - rm1[k];

        // Form v
        matvec_mult(B.data(), &v[(j + 1) * n], rv1.data(), n);
        for (int k = 0; k < n; ++k) 
            rv2[k] = F[k] - rv1[k];
        matrix_divide(alpha.data(), rv2.data(), &v[j * n], n, 1); 

        // ====================== Solve for dU_new ====================== //
        for (int k = 0; k < n; ++k)
            dU_new[((i + 1) * (Ny + 2) + j + 1) * n + k] = v[0 * n + k];

        for (int j = 1; j < Ny; ++j) {
            matvec_mult(&g[j * n * n], &dU_new[((i + 1) * (Ny + 2) + j) * n], rv1.data(), n);
            for (int k = 0; k < n; ++k)
                dU_new[((i + 1) * (Ny + 2) + j + 1) * n + k] = v[j * n + k] - rv1[k];
            
        }


        for (int i = 1; i < Nx_local - 1; ++i) {
            relax_inner_lines_real(i); 
        }
    }
}
void Solver2D::relax_inner_lines_real(int i) { 

    // ====================== Top cell ====================== // 
    int j = Ny - 1;
    int LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
        BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
        CC = i * Ny + j;
    
    // Top face
    for (int k = 0; k < n * n; ++k) 
        rm0[k] = jMinus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k];

    matmat_mult(&jET[i * n * n], rm0.data(), rm1.data(), n);

    // Left face
    for (int k = 0; k < n * n; ++k) 
        rm0[k] = iPlus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k]; 

    matvec_mult(rm0.data(), &dU_old[(i * (Ny + 2) + j + 1) * n], rv1.data(), n);

    // Right face
    for (int k = 0; k < n * n; ++k) 
        rm0[k] = iMinus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k];  

    matvec_mult(rm0.data(), &dU_old[((i + 2) * (Ny + 2) + j + 1) * n], rv2.data(), n); 

    for (int k = 0; k < n * n; ++k) {
        A[k] = 
        Volume[CC] / dt * I[k]          // Time-derivative part
        - (iMinus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k] ) * iArea[LF] // Left-face contribution
        + (iPlus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k] ) * iArea[RF]   // Right-face contribution
        - (jMinus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k] ) * jArea[BF]  // Bottom-face contribution
        + (jPlus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k] + rm1[k]) * jArea[TF]; // Top-face contribution

        C[k] = -(jPlus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k] ) * jArea[BF];
    }

    for (int k = 0; k < n; ++k) 
        F[k] = (iFlux[LF * n + k] - irho_flux[LF * n + k] ) * iArea[LF]
            - (iFlux[RF * n + k] - irho_flux[RF * n + k] ) * iArea[RF]
            + (jFlux[BF * n + k] - jrho_flux[BF * n + k] ) * jArea[BF]
            - (jFlux[TF * n + k] - jrho_flux[TF * n + k] ) * jArea[TF]
            + rv1[k] * iArea[LF]
            - rv2[k] * iArea[RF];
    
    alpha = A;
    matrix_divide(alpha.data(), F.data(), &v[j * n], n, 1);  //n = rows of x, m = columns of x
    matrix_divide(alpha.data(), C.data(), &g[j * n * n], n, n); 

    // ======================= Inner cells ======================== //
    for (int j = Ny - 2; j > 0; --j) {

        LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
        BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
        CC = i * Ny + j;

        // Left face
        for (int k = 0; k < n * n; ++k) 
            rm0[k] = iPlus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k]; 

        matvec_mult(rm0.data(), &dU_old[(i * (Ny + 2) + j + 1) * n], rv1.data(), n);

        // Right face
        for (int k = 0; k < n * n; ++k) 
            rm0[k] = iMinus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k]; 

        matvec_mult(rm0.data(), &dU_old[((i + 2) * (Ny + 2) + j + 1) * n], rv2.data(), n); 

        for (int k =  0; k < n * n; ++k) {
            B[k] = (jMinus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k] ) * jArea[TF]; 

            A[k] = 
                Volume[CC] / dt * I[k]          // Time-derivative part
                - (iMinus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k] ) * iArea[LF] // Left-face contribution
                + (iPlus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k] ) * iArea[RF]   // Right-face contribution
                - (jMinus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k] ) * jArea[BF]  // Bottom-face contribution
                + (jPlus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k] ) * jArea[TF]; // Top-face contribution

            C[k] = -(jPlus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k] ) * jArea[BF];
        }

        for (int k = 0; k < n; ++k) 
            F[k] = (iFlux[LF * n + k] - irho_flux[LF * n + k] ) * iArea[LF]
                - (iFlux[RF * n + k] - irho_flux[RF * n + k] ) * iArea[RF]
                + (jFlux[BF * n + k] - jrho_flux[BF * n + k] ) * jArea[BF]
                - (jFlux[TF * n + k] - jrho_flux[TF * n + k] ) * jArea[TF]
                + rv1[k] * iArea[LF]
                - rv2[k] * iArea[RF];

        // Form alpha
        matmat_mult(B.data(), &g[(j + 1) * n * n], rm1.data(), n);
        for (int k = 0; k < n * n; ++k) 
            alpha[k] = A[k] - rm1[k];
    
        // Form g
        matrix_divide(alpha.data(), C.data(), &g[j * n * n], n, n); 

        // Form v
        matvec_mult(B.data(), &v[(j + 1) * n], rv1.data(), n); 
        for (int k = 0; k < n; ++k) 
            rv2[k] = F[k] - rv1[k];        
        matrix_divide(alpha.data(), rv2.data(), &v[j * n], n, 1);
    }
    // ==================== Bottom cell ===================== //
    j = 0;

    LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
    BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
    CC = i * Ny + j;

    // Bottom face
    for (int k = 0; k < n * n; ++k) 
        rm0[k] = jPlus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k];

    matmat_mult(&jEB[i * n * n], rm0.data(), rm1.data(), n);

    // Left face
    for (int k = 0; k < n * n; ++k) 
        rm0[k] = iPlus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k]; 

    matvec_mult(rm0.data(), &dU_old[(i * (Ny + 2) + j + 1) * n], rv1.data(), n);

    // Right face
    for (int k = 0; k < n * n; ++k) 
        rm0[k] = iMinus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k]; 

    matvec_mult(rm0.data(), &dU_old[((i + 2) * (Ny + 2) + j + 1) * n], rv2.data(), n); 


    for (int k = 0; k < n * n; ++k) {

        B[k] = (jMinus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k] ) * jArea[TF]; 

        A[k] = 
        Volume[CC] / dt * I[k]          // Time-derivative part
        - (iMinus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k] ) * iArea[LF] // Left-face contribution
        + (iPlus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k] ) * iArea[RF]   // Right-face contribution
        - (jMinus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k] + rm1[k] ) * jArea[BF]  // Bottom-face contribution
        + (jPlus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k] )* jArea[TF]; // Top-face contribution
    }

    for (int k = 0; k < n; ++k) 
        F[k] = (iFlux[LF * n + k] - irho_flux[LF * n + k] ) * iArea[LF]
            - (iFlux[RF * n + k] - irho_flux[RF * n + k] ) * iArea[RF]
            + (jFlux[BF * n + k] - jrho_flux[BF * n + k] ) * jArea[BF]
            - (jFlux[TF * n + k] - jrho_flux[TF * n + k] ) * jArea[TF]
            + rv1[k] * iArea[LF]
            - rv2[k] * iArea[RF];
    
    // Form alpha
    matmat_mult(B.data(), &g[(j + 1) * n * n], rm1.data(), n);
    for (int k = 0; k < n * n; ++k) 
        alpha[k] = A[k] - rm1[k];

    // Form v
    matvec_mult(B.data(), &v[(j + 1) * n], rv1.data(), n);
    for (int k = 0; k < n; ++k) 
        rv2[k] = F[k] - rv1[k];
    matrix_divide(alpha.data(), rv2.data(), &v[j * n], n, 1); 

    // ====================== Solve for dU_new ====================== //
    for (int k = 0; k < n; ++k)
        dU_new[((i + 1) * (Ny + 2) + j + 1) * n + k] = v[0 * n + k];

    for (int j = 1; j < Ny; ++j) {
        matvec_mult(&g[j * n * n], &dU_new[((i + 1) * (Ny + 2) + j) * n], rv1.data(), n);
        for (int k = 0; k < n; ++k) 
            dU_new[((i + 1) * (Ny + 2) + j + 1) * n + k] = v[j * n + k] - rv1[k];
    }
}
void Solver2D::relax_right_line_real() {

    if (rank == size - 1) {

        for (int i = 0; i < Nx_local - 1; ++i) {
            relax_inner_lines_real(i); 
        }

        // ====================== Top cell ====================== // 
        int i = Nx_local - 1, j = Ny - 1; 
        int LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
            BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
            CC = i * Ny + j;

        // Right face
        for (int k = 0; k < n * n; ++k) 
            rm0[k] = iMinus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k];
        
        matmat_mult(&iER[j * n * n], rm0.data(), rm1.data(), n);

        // Left face
        for (int k = 0; k < n * n; ++k) 
            rm0[k] = iPlus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k];

        matvec_mult(rm0.data(), &dU_old[(i * (Ny + 2) + j + 1) * n], rv1.data(), n);

        // Top face
        for (int k = 0; k < n * n; ++k) 
            rm0[k] = jMinus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k];

        matmat_mult(&jET[i * n * n], rm0.data(), rm2.data(), n);
    

        for (int k = 0; k < n * n; ++k) {
            A[k] = 
            Volume[CC] / dt * I[k]          // Time-derivative part
            - (iMinus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k] ) * iArea[LF] // Left-face contribution
            + (iPlus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k] + rm1[k]) * iArea[RF]   // Right-face contribution
            - (jMinus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k] ) * jArea[BF]  // Bottom-face contribution
            + (jPlus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k] + rm2[k]) * jArea[TF]; // Top-face contribution

            C[k] = -(jPlus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k] ) * jArea[BF];
        }

        for (int k = 0; k < n; ++k) 
            F[k] = (iFlux[LF * n + k] - irho_flux[LF * n + k] ) * iArea[LF]
                - (iFlux[RF * n + k] - irho_flux[RF * n + k] ) * iArea[RF]
                + (jFlux[BF * n + k] - jrho_flux[BF * n + k] ) * jArea[BF]
                - (jFlux[TF * n + k] - jrho_flux[TF * n + k] ) * jArea[TF]
                + rv1[k] * iArea[LF];
        
        alpha = A;
        matrix_divide(alpha.data(), F.data(), &v[j * n], n, 1);  //n = rows of x, m = columns of x
        matrix_divide(alpha.data(), C.data(), &g[j * n * n], n, n); 


        // ======================= Inner cells ======================== //
        for (int j = Ny - 2; j > 0; --j) {

            LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
            BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
            CC = i * Ny + j;

            // Right face
            for (int k = 0; k < n * n; ++k) 
                rm0[k] = iMinus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k];
            
            matmat_mult(&iER[j * n * n], rm0.data(), rm1.data(), n);

            // Left face
            for (int k = 0; k < n * n; ++k) 
                rm0[k] = iPlus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k];

            matvec_mult(rm0.data(), &dU_old[(i * (Ny + 2) + j + 1) * n], rv1.data(), n);


            for (int k =  0; k < n * n; ++k) {
                B[k] = (jMinus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k] ) * jArea[TF]; 

                A[k] = 
                    Volume[CC] / dt * I[k]          // Time-derivative part
                    - (iMinus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k] ) * iArea[LF] // Left-face contribution
                    + (iPlus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k] + rm1[k]) * iArea[RF]   // Right-face contribution
                    - (jMinus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k] ) * jArea[BF]  // Bottom-face contribution
                    + (jPlus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k] ) * jArea[TF]; // Top-face contribution

                C[k] = -(jPlus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k] ) * jArea[BF];
            }

            for (int k = 0; k < n; ++k) 
                F[k] = (iFlux[LF * n + k] - irho_flux[LF * n + k] ) * iArea[LF]
                    - (iFlux[RF * n + k] - irho_flux[RF * n + k] ) * iArea[RF]
                    + (jFlux[BF * n + k] - jrho_flux[BF * n + k] ) * jArea[BF]
                    - (jFlux[TF * n + k] - jrho_flux[TF * n + k] ) * jArea[TF]
                    + rv1[k] * iArea[LF];

            // Form alpha
            matmat_mult(B.data(), &g[(j + 1) * n * n], rm1.data(), n);
            for (int k = 0; k < n * n; ++k) 
                alpha[k] = A[k] - rm1[k];
        
            // Form g
            matrix_divide(alpha.data(), C.data(), &g[j * n * n], n, n); 

            // Form v
            matvec_mult(B.data(), &v[(j + 1) * n], rv1.data(), n); 
            for (int k = 0; k < n; ++k) 
                rv2[k] = F[k] - rv1[k];        
            matrix_divide(alpha.data(), rv2.data(), &v[j * n], n, 1);
        }

        // ==================== Bottom cell ===================== //
        j = 0;

        LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
        BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
        CC = i * Ny + j;

        // Right face
        for (int k = 0; k < n * n; ++k) 
            rm0[k] = iMinus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k];
        
        matmat_mult(&iER[j * n * n], rm0.data(), rm1.data(), n);

        // Left face
        for (int k = 0; k < n * n; ++k) 
            rm0[k] = iPlus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k];

        matvec_mult(rm0.data(), &dU_old[(i * (Ny + 2) + j + 1) * n], rv1.data(), n);

        // Bottom face
        for (int k = 0; k < n * n; ++k)
            rm0[k] = jPlus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k];
      
        matmat_mult(&jEB[i * n * n], rm0.data(), rm2.data(), n); 

        for (int k = 0; k < n * n; ++k) {

            B[k] = (jMinus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k] ) * jArea[TF]; 

            A[k] = 
            Volume[CC] / dt * I[k]          // Time-derivative part
            - (iMinus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k] ) * iArea[LF] // Left-face contribution
            + (iPlus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k] + rm1[k]) * iArea[RF]   // Right-face contribution
            - (jMinus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k] + rm2[k]) * jArea[BF]  // Bottom-face contribution
            + (jPlus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k] ) * jArea[TF]; // Top-face contribution
        }

        for (int k = 0; k < n; ++k) 
            F[k] = (iFlux[LF * n + k] - irho_flux[LF * n + k] ) * iArea[LF]
                - (iFlux[RF * n + k] - irho_flux[RF * n + k] ) * iArea[RF]
                + (jFlux[BF * n + k] - jrho_flux[BF * n + k] ) * jArea[BF]
                - (jFlux[TF * n + k] - jrho_flux[TF * n + k] ) * jArea[TF]
                + rv1[k] * iArea[LF];
        
        // Form alpha
        matmat_mult(B.data(), &g[(j + 1) * n * n], rm1.data(), n);
        for (int k = 0; k < n * n; ++k) 
            alpha[k] = A[k] - rm1[k];

        // Form v
        matvec_mult(B.data(), &v[(j + 1) * n], rv1.data(), n);
        for (int k = 0; k < n; ++k) 
            rv2[k] = F[k] - rv1[k];
        matrix_divide(alpha.data(), rv2.data(), &v[j * n], n, 1); 

        // ====================== Solve for dU_new ====================== //
        for (int k = 0; k < n; ++k)
            dU_new[((i + 1) * (Ny + 2) + j + 1) * n + k] = v[0 * n + k];

        for (int j = 1; j < Ny; ++j) {
            matvec_mult(&g[j * n * n], &dU_new[((i + 1) * (Ny + 2) + j) * n], rv1.data(), n);
            for (int k = 0; k < n; ++k) 
                dU_new[((i + 1) * (Ny + 2) + j + 1) * n + k] = v[j * n + k] - rv1[k];
        }        
    }
}

void Solver2D::update_U() {

    for (int i = 1; i < Nx_local + 1; ++i) {
        for (int j = 1; j < Ny + 1; ++j) {
            int idx = (i * (Ny + 2) + j) * n;
            for (int k = 0; k < n; ++k) {
                U[idx + k] += dU_old[idx + k];
            }
        }
    }

}

void Solver2D::compute_dt() {
    Vector V(4); 
    double c, max_new = 0.0, max_old = 0.0, sidesum, l, r, b, t;

    for (int i = 0; i < Nx_local; ++i) {
        for (int j = 0; j < Ny; ++j) {
            constoprim(&U[((i + 1) * (Ny + 2) + j + 1) * n], V.data(), cell_thermo[(i + 1) * (Ny + 2) + j + 1].gamma, n_vel);
            c = cell_thermo[(i + 1) * (Ny + 2) + j + 1].a;  

            l = (fabs(V[1] * iFxNorm[i * Ny + j] + V[2] * iFyNorm[i * Ny + j]) + c) * iArea[i * Ny + j];   
			r = (fabs(V[1] * iFxNorm[(i + 1) * Ny + j] + V[2] * iFyNorm[(i + 1) * Ny + j]) + c) * iArea[(i + 1) * Ny + j];   
            b = (fabs(V[1] * jFxNorm[i * (Ny + 1) + j] + V[2] * jFyNorm[i * (Ny + 1) + j]) + c) * jArea[i * (Ny + 1) + j];   
			t = (fabs(V[1] * jFxNorm[i * (Ny + 1) + j + 1] + V[2] * jFyNorm[i * (Ny + 1) + j + 1]) + c) * jArea[i * (Ny + 1) + j + 1];   

            sidesum = l + r + b + t;

            max_new = sidesum/(2 * Volume[i * Ny + j]);
            if (max_new >= max_old) max_old = max_new;

        }
    }

    double dt_local = max_old; 

    double dt_global;

    MPI_Allreduce(&dt_local, &dt_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); 

    dt = CFL/dt_global; 

}
void Solver2D::compute_inner_res_perf() {

    Vector X(n, 0.0), Y(n, 0.0), Z(n, 0.0), ID = {1, 0, 0, 0}; 

    inner_res = 0.0; 
    double inner_res_local = 0.0;
    double F, res = 0.0;

    for (int i = 2; i < Nx_local; ++i) {
        for (int j = 2; j < Ny; ++j) {

            int LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
            BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
            CC = i * Ny + j;

            for (int k =  0; k < n ; ++k) {
                Y[k] = jMinus_A[TF * n * n + k] * jArea[TF]; 

                X[k] = 
                    Volume[CC] / dt * ID[k]          // Time-derivative part
                    - iMinus_A[LF * n * n + k] * iArea[LF] // Left-face contribution
                    + iPlus_A[RF * n * n + k] * iArea[RF]   // Right-face contribution
                    - jMinus_A[BF * n * n + k] * jArea[BF]  // Bottom-face contribution
                    + jPlus_A[TF * n * n + k] * jArea[TF]; // Top-face contribution

                Z[k] = -jPlus_A[BF * n * n + k] * jArea[BF];
            }

            
            F = iFlux[LF * n] * iArea[LF]
                - iFlux[RF * n] * iArea[RF]
                + jFlux[BF * n] * jArea[BF]
                - jFlux[TF * n] * jArea[TF]
                + dot_product(&iPlus_A[LF * n * n], &dU_old[(i * (Ny + 2) + j + 1) * n], n) * iArea[LF]
                - dot_product(&iMinus_A[RF * n * n], &dU_old[((i + 2) * (Ny + 2) + j + 1) * n], n) * iArea[RF];

            res = dot_product(Y.data(), &dU_old[((i + 1) * (Ny + 2) + j + 2) * n], n)
                    + dot_product(X.data(), &dU_old[((i + 1) * (Ny + 2) + j + 1) * n], n) 
                    + dot_product(Z.data(), &dU_old[((i + 1) * (Ny + 2) + j) * n], n) - F; 

            inner_res_local += res * res;
        }
    }

    MPI_Allreduce(&inner_res_local, &inner_res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    inner_res = sqrt(inner_res); 
}
void Solver2D::compute_inner_res_real() {

    Vector X(n, 0.0), Y(n, 0.0), Z(n, 0.0), ID = {1, 0, 0, 0}; 

    inner_res = 0.0; 
    double inner_res_local = 0.0;
    double F, res = 0.0;
    Vector vv1(4), vv2(4); 

    for (int i = 2; i < Nx_local; ++i) {
        for (int j = 2; j < Ny; ++j) {

            int LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
            BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
            CC = i * Ny + j;

            // Left face
            for (int k = 0; k < n; ++k) 
                vv1[k] = iPlus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k]; 

            // Right face
            for (int k = 0; k < n; ++k) 
                vv2[k] = iMinus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k]; 


            for (int k =  0; k < n; ++k) {
                Y[k] = (jMinus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k] ) * jArea[TF]; 

                X[k] = 
                    Volume[CC] / dt * I[k]          // Time-derivative part
                    - (iMinus_A[LF * n * n + k] - 0.5 * irho_A[LF * n * n + k] ) * iArea[LF] // Left-face contribution
                    + (iPlus_A[RF * n * n + k] - 0.5 * irho_A[RF * n * n + k] ) * iArea[RF]   // Right-face contribution
                    - (jMinus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k] ) * jArea[BF]  // Bottom-face contribution
                    + (jPlus_A[TF * n * n + k] - 0.5 * jrho_A[TF * n * n + k] ) * jArea[TF]; // Top-face contribution

                Z[k] = -(jPlus_A[BF * n * n + k] - 0.5 * jrho_A[BF * n * n + k] ) * jArea[BF];
            }

        
            F = (iFlux[LF * n] - irho_flux[LF * n] ) * iArea[LF]
                - (iFlux[RF * n] - irho_flux[RF * n] ) * iArea[RF]
                + (jFlux[BF * n] - jrho_flux[BF * n] ) * jArea[BF]
                - (jFlux[TF * n] - jrho_flux[TF * n] ) * jArea[TF]
                + dot_product(vv1.data(), &dU_old[(i * (Ny + 2) + j + 1) * n], n) * iArea[LF]
                - dot_product(vv2.data(), &dU_old[((i + 2) * (Ny + 2) + j + 1) * n], n) * iArea[RF];


            res = dot_product(Y.data(), &dU_old[((i + 1) * (Ny + 2) + j + 2) * n], n)
                    + dot_product(X.data(), &dU_old[((i + 1) * (Ny + 2) + j + 1) * n], n) 
                    + dot_product(Z.data(), &dU_old[((i + 1) * (Ny + 2) + j) * n], n) - F; 

            inner_res_local += res * res;
        }
    }

    MPI_Allreduce(&inner_res_local, &inner_res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    inner_res = sqrt(inner_res); 
}

void Solver2D::compute_outer_res_perf() {
    outer_res = 0.0;
    double outer_res_local = 0.0;
    double intres;

    for (int i = 0; i < Nx_local; ++i) {
        for (int j = 0; j < Ny; ++j) {

            intres = (- iFlux[(i * Ny + j) * n] * iArea[i * Ny + j]
                        + iFlux[((i + 1) * Ny + j) * n] * iArea[(i + 1) * Ny + j]
                        - jFlux[(i * (Ny + 1) + j) * n] * jArea[(i * (Ny + 1) + j)]
                        + jFlux[(i * (Ny + 1) + j + 1) * n] * jArea[i * (Ny + 1) + j + 1] ) / Volume[i * Ny + j];

            outer_res_local += intres * intres;      
        }
    }
    MPI_Allreduce(&outer_res_local, &outer_res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    outer_res = sqrt(outer_res); 
}
void Solver2D::compute_outer_res_real() {
    outer_res = 0.0;
    double outer_res_local = 0.0;
    double intres;

    for (int i = 0; i < Nx_local; ++i) {
        for (int j = 0; j < Ny; ++j) {

            intres = (- (iFlux[(i * Ny + j) * n] - irho_flux[(i * Ny + j) * n] ) * iArea[i * Ny + j]
                        + (iFlux[((i + 1) * Ny + j) * n] - irho_flux[((i + 1) * Ny + j) * n] ) * iArea[(i + 1) * Ny + j]
                        - (jFlux[(i * (Ny + 1) + j) * n] - jrho_flux[(i * (Ny + 1) + j) * n] ) * jArea[(i * (Ny + 1) + j)]
                        + (jFlux[(i * (Ny + 1) + j + 1) * n] - jrho_flux[(i * (Ny + 1) + j + 1) * n] ) * jArea[i * (Ny + 1) + j + 1] ) / Volume[i * Ny + j];

            outer_res_local += intres * intres;      
        }
    }

    MPI_Allreduce(&outer_res_local, &outer_res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    outer_res = sqrt(outer_res); 
}

void Solver2D::explicit_update() {

    // cout << endl << "dU_old: " << endl;
    for (int i = 0; i < Nx_local; ++i) {
        for (int j = 0; j < Ny; ++j) {

            int LF = (i * Ny + j), RF = ((i + 1) * Ny + j), 
                BF = (i * (Ny + 1) + j), TF = (i * (Ny + 1) + j + 1),
                CC = i * Ny + j;

                // cout << i + 1 << ", " << j + 1 << ": ";
                for (int k = 0; k < n; ++k) {
                    dU_old[((i + 1) * (Ny + 2) + j + 1) * n + k] = - dt / Volume[CC] * (- iFlux[LF * n + k] * iArea[LF]
                                                                 + iFlux[RF * n + k] * iArea[RF]
                                                                 - jFlux[BF * n + k] * jArea[BF]
                                                                 + jFlux[TF * n + k] * jArea[TF] );

                    // cout << dU_old[((i + 1) * (Ny + 2) + j + 1) * n + k] << " ";
                }
            // cout << endl;
        }
    }
}

void Solver2D::print_by_rank(Vector Vec, int nx, int ny, int nvars, string name){

    for (int r = 0; r < size; ++r) {
        MPI_Barrier(MPI_COMM_WORLD);  // Synchronize before each rank prints
        if (rank == r) {
            std::cout << "========== Rank " << rank << " ==========" << std::endl;
            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) { 
                    int cell_idx = (i * ny + j) * nvars;  // assuming i-fastest
                    std::cout << name << " [" << i << ", " << j << "]: ";
                    for (int k = 0; k < nvars; ++k) {
                        std::cout << Vec[cell_idx + k] << " ";
                    }
                    std::cout << "\n";
                }
            }
            std::cout << std::flush;
        }
    }

}
void Solver2D::finalize() {

    Vector U_sendbuf(Nx_local * Ny * n, 0.0);

    for (int i = 0; i < Nx_local; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < n; ++k) {
                U_sendbuf[(i * Ny + j) * n + k] =
                    U[((i + 1) * (Ny + 2) + (j + 1)) * n + k];
            }
        }
    }

    vector<int> recvcounts(size);
    vector<int> displs(size);

    int offset = 0;
    for (int r = 0; r < size; ++r) {
        recvcounts[r] = local_Nx[r] * Ny * n; // number of elements each rank sends
        displs[r] = offset;
        offset += recvcounts[r];
    }

    MPI_Gatherv(
        U_sendbuf.data(), Nx_local * Ny * n, MPI_DOUBLE,
        U_gathered.data(), recvcounts.data(), displs.data(), MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );        

    writeTecplotDat(); 
}

void Solver2D::load_thermochemical_table() {
    const int N = 500;

    string thermofilename = "../gibbslib/thermochemical_table.csv";
    ifstream file(thermofilename); 

    if (!file.is_open()) {
        cerr << "Error: could not open " << thermofilename << "\n";
        return;
    }

    string line;
    getline(file, line); // skip header

    // Resize table and coordinate vectors
    int count = 0;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (!getline(file, line)) {
                cerr << "Error: unexpected end of file at entry " << count << "\n";
                return;
            }

            stringstream ss(line);
            string token;
            ThermoEntry entry;

            auto read = [&](double& val) {
                if (!getline(ss, token, ',')) {
                    cerr << "Error: malformed line at entry " << count << "\n";
                    return false;
                }
                val = stod(token);
                return true;
            };

            if (!read(entry.rho)    || !read(entry.e)    || !read(entry.p) ||
                !read(entry.T)      || !read(entry.R)    || !read(entry.cv) ||
                !read(entry.gamma)  || !read(entry.dpdrho) ||!read(entry.dpde) || !getline(ss, token)) {
                cerr << "Error: incomplete or bad format at entry " << count << "\n";
                return;
            }

            entry.a = stod(token); // last field, no comma

            thermochemical_table[i * N + j] = entry;

            ++count;
        }
    }

    if (rank == 0) 
        cout << "Successfully loaded " << count << " thermochemical entries.\n";
}
ThermoEntry Solver2D::bilinear_interpolate(double rho, double e) {
    constexpr int N = 500;

    // Return default entry if out of bounds
    if (rho < 5e-5 || rho > 1.0 ||
        e   < 717 * 600 || e > 5e7) {
        cout << "Warning: interpolation input (rho=" << rho << ", e=" << e << ") out of bounds.\n";
        return ThermoEntry{};
    }

    // Locate bounding indices
    auto it_rho = upper_bound(densities.begin(), densities.end(), rho);
    int i = max(0, int(it_rho - densities.begin()) - 1);
    i = min(i, N - 2);

    auto it_e = upper_bound(internal_energies.begin(), internal_energies.end(), e);
    int j = max(0, int(it_e - internal_energies.begin()) - 1);
    j = min(j, N - 2);

    double x1 = densities[i],     x2 = densities[i+1];
    double y1 = internal_energies[j], y2 = internal_energies[j+1];

    double x = rho;
    double y = e;

    const ThermoEntry& Q11 = thermochemical_table[i * N + j];
    const ThermoEntry& Q21 = thermochemical_table[(i+1) * N + j];
    const ThermoEntry& Q12 = thermochemical_table[i * N + (j+1)];
    const ThermoEntry& Q22 = thermochemical_table[(i+1) * N + (j+1)];

    auto bilinear = [&](double f11, double f21, double f12, double f22) -> double {
        double dx = x2 - x1;
        double dy = y2 - y1;

        if (dx == 0 || dy == 0)
            return f11; // fallback: degenerate grid spacing

        double t = (x - x1) / dx;
        double u = (y - y1) / dy;

        return (1 - t) * (1 - u) * f11 +
               t       * (1 - u) * f21 +
               (1 - t) * u       * f12 +
               t       * u       * f22;
    };

    ThermoEntry result;
    result.rho    = rho;
    result.e      = e;
    result.p      = bilinear(Q11.p, Q21.p, Q12.p, Q22.p);
    result.T      = bilinear(Q11.T, Q21.T, Q12.T, Q22.T);
    result.R      = bilinear(Q11.R, Q21.R, Q12.R, Q22.R);
    result.cv     = bilinear(Q11.cv, Q21.cv, Q12.cv, Q22.cv);
    result.gamma  = bilinear(Q11.gamma, Q21.gamma, Q12.gamma, Q22.gamma);
    result.dpdrho = bilinear(Q11.dpdrho, Q21.dpdrho, Q12.dpdrho, Q22.dpdrho);
    result.dpde   = bilinear(Q11.dpde, Q21.dpde, Q12.dpde, Q22.dpde);
    result.a      = bilinear(Q11.a, Q21.a, Q12.a, Q22.a);

    return result;
}
void Solver2D::initialize_chemistry() {

    if (real_gas) {

        double density, energy; 

        for (int i = 0; i < Nx_local + 2; ++i) {
            for (int j = 0; j < Ny + 2; ++j) {

                int idx = i * (Ny + 2) + j;

                density = U[idx * n];
                energy =  computeInternalEnergy(&U[idx * n], n_vel);
                
                if (energy < 717 * 650) {

                    cell_thermo[idx].rho = U[idx * n];
                    cell_thermo[idx].e = computeInternalEnergy(&U[idx * n], n_vel);
                    cell_thermo[idx].p = computePressure(&U[idx * n], perfgam, n_vel);
                    cell_thermo[idx].R = 287.0;
                    cell_thermo[idx].T = cell_thermo[idx].p / (cell_thermo[idx].rho * cell_thermo[idx].R);
                    cell_thermo[idx].cv = 717.0;
                    cell_thermo[idx].gamma = perfgam;
                    cell_thermo[idx].dpdrho = (perfgam - 1) * cell_thermo[idx].e;
                    cell_thermo[idx].dpde = (perfgam - 1) * cell_thermo[idx].rho;
                    cell_thermo[idx].a = sqrt(perfgam * cell_thermo[idx].p / cell_thermo[idx].rho);
                }
                else {
                    if (using_table) {         
                        cell_thermo[idx] = bilinear_interpolate(density, energy);                     
                    }
                    else {                    
                        cell_thermo[idx] = chem.compute_equilibrium_thermodynamic_variables(density, energy); 
                    }
                }


            }
        }
    }       
    else {

        for (int i = 0; i < Nx_local + 2; ++i) {
            for (int j = 0; j < Ny + 2; ++j) {
                int idx = i * (Ny + 2) + j;

                cell_thermo[idx].rho = U[idx * n];
                cell_thermo[idx].e = computeInternalEnergy(&U[idx * n], n_vel);
                cell_thermo[idx].p = computePressure(&U[idx * n], perfgam, n_vel);
                cell_thermo[idx].R = 287.0;
                cell_thermo[idx].T = cell_thermo[idx].p / (cell_thermo[idx].rho * cell_thermo[idx].R);
                cell_thermo[idx].cv = 717.0;
                cell_thermo[idx].gamma = perfgam;
                cell_thermo[idx].dpdrho = (perfgam - 1) * cell_thermo[idx].e;
                cell_thermo[idx].dpde = (perfgam - 1) * cell_thermo[idx].rho;
                cell_thermo[idx].a = sqrt(perfgam * cell_thermo[idx].p / cell_thermo[idx].rho);

            }
        }
    }
}
void Solver2D::get_real_chemistry() {

    double density, energy; 

    for (int i = 0; i < Nx_local + 2; ++i) {
        for (int j = 0; j < Ny + 2; ++j) {
            
            int idx = i * (Ny + 2) + j;
            density = U[idx * n];
            energy =  computeInternalEnergy(&U[idx * n], n_vel);

            if (energy < 717 * 650) {
                cell_thermo[idx].rho = U[idx * n];
                cell_thermo[idx].e = computeInternalEnergy(&U[idx * n], n_vel);
                cell_thermo[idx].p = computePressure(&U[idx * n], perfgam, n_vel);
                cell_thermo[idx].R = 287.0;
                cell_thermo[idx].T = cell_thermo[idx].p / (cell_thermo[idx].rho * cell_thermo[idx].R);
                cell_thermo[idx].cv = 717.0;
                cell_thermo[idx].gamma = perfgam;
                cell_thermo[idx].dpdrho = (perfgam - 1) * cell_thermo[idx].e;
                cell_thermo[idx].dpde = (perfgam - 1) * cell_thermo[idx].rho;
                cell_thermo[idx].a = sqrt(perfgam * cell_thermo[idx].p / cell_thermo[idx].rho);
            }
            else {
                if (using_table) {         
                    cell_thermo[idx] = bilinear_interpolate(density, energy);                     
                }
                else {                    
                    cell_thermo[idx] = chem.compute_equilibrium_thermodynamic_variables(density, energy); 
                }
            }
        }
    }

}
void Solver2D::get_perf_chemistry() {

    for (int i = 0; i < Nx_local + 2; ++i) {
        for (int j = 0; j < Ny + 2; ++j) {
            int idx = i * (Ny + 2) + j;

            cell_thermo[idx].rho = U[idx * n];
            cell_thermo[idx].e = computeInternalEnergy(&U[idx * n], n_vel);
            cell_thermo[idx].p = computePressure(&U[idx * n], perfgam, n_vel);
            cell_thermo[idx].R = 287.0;
            cell_thermo[idx].cv = 717.0;
            cell_thermo[idx].gamma = perfgam;
            cell_thermo[idx].dpdrho = (perfgam - 1) * cell_thermo[idx].e;
            cell_thermo[idx].dpde = (perfgam - 1) * cell_thermo[idx].rho;
            cell_thermo[idx].a = sqrt(perfgam * cell_thermo[idx].p / cell_thermo[idx].rho);

        }
    }
}


void Solver2D::writeTecplotDat() {
    if (rank != 0) return;

    ofstream file(filename);
    file << "VARIABLES = \"x\", \"y\", \"density\", \"u-vel\", \"v-vel\", \"pressure\", \"temperature\", \"a\", \"e\" \n";
    file << "ZONE T=\"Flow Field\", I=" << Nx + 1 << ", J=" << Ny + 1 << ", F=BLOCK\n";
    file << "VARLOCATION=([3-9]=CELLCENTERED)\n";

    vector<ThermoEntry> thermo(Nx * Ny); 

    if (real_gas) {
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                int idx = i * Ny + j;

                double density = U_gathered[idx * n];
                double energy =  computeInternalEnergy(&U_gathered[idx * n], n_vel);

                if (energy < 717 * 600) {

                    thermo[idx].rho = U_gathered[idx * n];
                    thermo[idx].e = computeInternalEnergy(&U_gathered[idx * n], n_vel);
                    thermo[idx].p = computePressure(&U_gathered[idx * n], perfgam, n_vel);
                    thermo[idx].R = 287.0;
                    thermo[idx].T = thermo[idx].p / (thermo[idx].rho * thermo[idx].R);
                    thermo[idx].cv = 717.0;
                    thermo[idx].gamma = perfgam;
                    thermo[idx].dpdrho = (perfgam - 1) * thermo[idx].e;
                    thermo[idx].dpde = (perfgam - 1) * thermo[idx].rho;
                    thermo[idx].a = sqrt(perfgam * thermo[idx].p / thermo[idx].rho);
                }
                else {
                    if (using_table == false) {
                        thermo[idx] = chem.compute_equilibrium_thermodynamic_variables(density, energy);
                    }
                    else {
                        thermo[idx] = bilinear_interpolate(density, energy); 
                    }
                }
            }
        }
    }
    else {

        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {

                int idx = i * Ny + j;
                
                thermo[idx].rho = U_gathered[idx * n];
                thermo[idx].e = computeInternalEnergy(&U_gathered[idx * n], n_vel);
                thermo[idx].p = computePressure(&U_gathered[idx * n], perfgam, n_vel);
                thermo[idx].R = 287.0;
                thermo[idx].T = thermo[idx].p / (thermo[idx].rho * thermo[idx].R);
                thermo[idx].cv = 717.0;
                thermo[idx].gamma = perfgam;
                thermo[idx].dpdrho = (perfgam - 1) * thermo[idx].e;
                thermo[idx].dpde = (perfgam - 1) * thermo[idx].rho;
                thermo[idx].a = sqrt(perfgam * thermo[idx].p / thermo[idx].rho);
            }
        }

    }

    // Write X vertices (node-centered, i-fastest)
    for (int j = 0; j < Ny + 1; ++j) {
        for (int i = 0; i < Nx + 1; ++i) {
            int idx = i * (Ny + 1) + j;
            file << x_vertices[idx] << "\n";
        }
    }

    // Write Y vertices (node-centered, i-fastest)
    for (int j = 0; j < Ny + 1; ++j) {
        for (int i = 0; i < Nx + 1; ++i) {
            int idx = i * (Ny + 1) + j;
            file << y_vertices[idx] << "\n";
        }
    }

    // Allocate temporary variable to hold primitive variables
    Vector V(n, 0.0);

    // Create buffers for each cell-centered variable
    Vector density(Nx * Ny), uvel(Nx * Ny), vvel(Nx * Ny), pressure(Nx * Ny), temperature(Nx * Ny);

    // Extract and store all primitive variables first (i-fastest)
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idx = i * Ny + j;
            constoprim(&U_gathered[idx * n], V.data(), thermo[idx].gamma, n_vel);
            density[idx]  = V[0];
            uvel[idx]     = V[1];
            vvel[idx]     = V[2];
            pressure[idx] = thermo[idx].p;
        }
    }

    // Write each variable block (i-fastest order)
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idx = i * Ny + j;
            file << density[idx] << "\n";
        }
    }

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idx = i * Ny + j;
            file << uvel[idx] << "\n";
        }
    }

    for (int j = 0; j < Ny; ++j) {  
        for (int i = 0; i < Nx; ++i) {
            int idx = i * Ny + j;
            file << vvel[idx] << "\n";
        }
    }

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idx = i * Ny + j;
            file << pressure[idx] << "\n";
        }
    }

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idx = i * Ny + j;
            file << thermo[idx].T << "\n";
        }
    }

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idx = i * Ny + j;
            file << thermo[idx].a << "\n";
        }
    }

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idx = i * Ny + j;
            file << thermo[idx].e << "\n"; 
        }
    }


    cout << "File created: " << filename << std::endl;
}
void Solver2D::writeParaviewCSV() {
    if (rank != 0) return;

    cout << "Writing file..." << endl;
    ofstream file(filename);

    file << "x,y,density,u-vel,v-vel,pressure,temperature,a,e\n";

    vector<ThermoEntry> thermo(Nx * Ny); 

    if (real_gas) {

        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {

                int idx = i * Ny + j;
                double density = U_gathered[idx * n];
                double energy = computeInternalEnergy(&U_gathered[idx * n], n_vel);


                if (energy < 717 * 650) {

                    thermo[idx].rho   = U_gathered[idx * n];
                    thermo[idx].e     = computeInternalEnergy(&U_gathered[idx * n], n_vel);
                    thermo[idx].p     = computePressure(&U_gathered[idx * n], perfgam, n_vel);
                    thermo[idx].R     = 287.0;
                    thermo[idx].T = thermo[idx].p / (thermo[idx].rho * thermo[idx].R);
                    thermo[idx].cv    = 717.0;
                    thermo[idx].gamma = perfgam;
                    thermo[idx].dpdrho = (perfgam - 1) * thermo[idx].e;
                    thermo[idx].dpde   = (perfgam - 1) * thermo[idx].rho;
                    thermo[idx].a     = sqrt(perfgam * thermo[idx].p / thermo[idx].rho);

                }
                else {
                    if (!using_table)
                        thermo[idx] = chem.compute_equilibrium_thermodynamic_variables(density, energy);
                    else
                        thermo[idx] = bilinear_interpolate(density, energy);
                }
            }
        }
    } 
    else {
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                int idx = i * Ny + j;

                thermo[idx].rho   = U_gathered[idx * n];
                thermo[idx].e     = computeInternalEnergy(&U_gathered[idx * n], n_vel);
                thermo[idx].p     = computePressure(&U_gathered[idx * n], perfgam, n_vel);
                thermo[idx].R     = 287.0;
                thermo[idx].T = thermo[idx].p / (thermo[idx].rho * thermo[idx].R);
                thermo[idx].cv    = 717.0;
                thermo[idx].gamma = perfgam;
                thermo[idx].dpdrho = (perfgam - 1) * thermo[idx].e;
                thermo[idx].dpde   = (perfgam - 1) * thermo[idx].rho;
                thermo[idx].a     = sqrt(perfgam * thermo[idx].p / thermo[idx].rho);
            }
        }
    }

    Vector V(n, 0.0);
    Vector density(Nx * Ny), uvel(Nx * Ny), vvel(Nx * Ny), pressure(Nx * Ny), temperature(Nx * Ny);

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idx = i * Ny + j;
            constoprim(&U_gathered[idx * n], V.data(), thermo[idx].gamma, n_vel);
            density[idx]  = V[0];
            uvel[idx]     = V[1];
            vvel[idx]     = V[2];
            pressure[idx] = thermo[idx].p;
            temperature[idx] = thermo[idx].T;
        }
    }

    // Write all variables as a flat CSV table
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idx_cell = i * Ny + j;
            int idx_vert = (i + 1) * (Ny + 1) + (j + 1); // center of cell

            // For accurate cell-centered plotting, use the average of four surrounding nodes
            // Here we approximate using lower-left corner of cell (or change this)
            double x = 0.25 * (x_vertices[i * (Ny + 1) + j] +
                               x_vertices[(i + 1) * (Ny + 1) + j] +
                               x_vertices[i * (Ny + 1) + (j + 1)] +
                               x_vertices[(i + 1) * (Ny + 1) + (j + 1)]);

            double y = 0.25 * (y_vertices[i * (Ny + 1) + j] +
                               y_vertices[(i + 1) * (Ny + 1) + j] +
                               y_vertices[i * (Ny + 1) + (j + 1)] +
                               y_vertices[(i + 1) * (Ny + 1) + (j + 1)]);

            file << x << "," << y << "," << density[idx_cell] << "," << uvel[idx_cell] << "," 
                 << vvel[idx_cell] << "," << pressure[idx_cell] << "," 
                 << thermo[idx_cell].T << "," << thermo[idx_cell].a << "," 
                 << thermo[idx_cell].e << "\n";
        }
    }

    cout << "CSV file written for ParaView: " << filename << endl;
}
