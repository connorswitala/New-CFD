#include "solver.hpp"


// Conversion functions
void primtocons(double* U, const double* V, const int dimensions) {
    
    U[0] = V[0]; // Set density

    double udotu = 0.0; // Set u dot u = 0.0

    // This loop computes each momentum and total udotu for n dimensions
    for (int i = 0; i < dimensions; ++i) {
        U[i + 1] = V[0] * V[i + 1]; 
        udotu += V[i + 1] * V[i + 1];  
    }

    U[dimensions + 1] = V[dimensions + 1] / (perfgam - 1) + 0.5 * V[0] * udotu; // Set total energy 
}
void constoprim(const double* U, double* V, const int dimensions) {
    
    V[0] = U[0]; // Set density

    double udotu = 0.0; // Set u dot u = 0.0;

    // This loop computes each velocit;y and total udotu for n dimensions
    for (int i = 0; i < dimensions; ++i) {
        V[i + 1] = U[i + 1] / U[0];
        udotu += U[i + 1] * U[i + 1];
    }

    V[dimensions + 1] = (U[dimensions + 1] - 0.5 / U[0] * udotu) * (perfgam - 1); // Set pressure
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
        primtocons(UL.data(), VL.data(), n_vel);
        primtocons(UR.data(), VR.data(), n_vel);

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
        constoprim(&U_gathered[idx], V.data(), n_vel);
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
            constoprim(&U[idx], V.data(), n_vel);
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

            V1 = {lc / (a * a), (u * lc + a * lt) / (a * a), ((U[idx + 2] + computePressure(&U[idx], n_vel)) / U[idx] * lc + a * u * lt) / (a * a)};
            Q = {0.5 * u * u * (perfgam - 1), -u * (perfgam - 1), (perfgam - 1)};
            V2 = {lt / a, u * lt / a + lc, (U[idx + 2] + computePressure(&U[idx], n_vel)) / U[idx] * lt / a + u * lc};
            W = {-u, 1, 0};

            outer_product(V1.data(), Q.data(), int2.data(), n);
            outer_product(V2.data(), W.data(), int3.data(), n);

            for (int k = 0; k < n * n; ++k) A_plus[aidx + k] = int1[k] + int2[k] + int3[k];
            matvec_mult(&A_plus[aidx], &U[idx], &F_plus[idx], n);

            // Negative flux calculation
            constoprim(&U[iidx], V.data(), n_vel);
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

            V1 = {lc / (a * a), (u * lc + a * lt) / (a * a), ((U[iidx + 2] + computePressure(&U[iidx], n_vel)) / U[iidx] * lc + a * u * lt) / (a * a)};
            Q = {0.5 * u * u * (perfgam - 1), -u * (perfgam - 1), (perfgam - 1)};
            V2 = {lt / a, u * lt / a + lc, (U[iidx + 2] + computePressure(&U[iidx], n_vel)) / U[iidx] * lt / a + u * lc};
            W = {-u, 1, 0};

            outer_product(V1.data(), Q.data(), int2.data(), n);
            outer_product(V2.data(), W.data(), int3.data(), n);

            for (int k = 0; k < n * n; ++k) A_minus[aidx + k] = int1[k] + int2[k] + int3[k];

            matvec_mult(&A_minus[aidx], &U[iidx], &F_minus[idx], n);

            // Total flux calculation
            for (int k = 0; k < n; ++k) Flux[idx + k] = F_plus[idx + k] + F_minus[idx + k];
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

Solver2D::Solver2D(int Nx, int Ny, double CFL, Vector U_inlet, BCMap BCs) : Nx(Nx), Ny(Ny), CFL(CFL) {
    
    RampGrid grid(Nx, Ny, 3, 1, 15);  // This initializes the grid on rank 0 to later be scattered using MPI.
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    t = 0.0; 

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
    vector<int> local_Nx(size), cc_sendcounts(size), cc_displacements(size), 
                if_sendcounts(size), if_displacements(size) , 
                jf_sendcounts(size), jf_displacements(size);

    int base = Nx / size;
    int rem = Nx % size; 

    int cc_offset = 0;
    for (int r = 0; r < size; ++r) {
        local_Nx[r]       = base + (r < rem ? 1 : 0);
        cc_displacements[r]  = cc_offset;
        cc_sendcounts[r]     = local_Nx[r] * Ny;
        cc_offset           += local_Nx[r] * Ny;
    }

    int Nx_local = local_Nx[rank]; 
    int num_ifaces = (Nx_local + 1) * Ny; 
    int num_jfaces = (Ny + 1) * Nx_local; 
    int num_loc_cells = Nx_local * Ny; 
    int num_gathered_cells = Nx * Ny;

    U_gathered = Vector( num_gathered_cells * n, 0.0); 
    U = Vector(((Nx_local + 2) * (Ny + 2)) * n, 0.0);  

    // Important Flux Vectors
    iFlux = Vector( num_ifaces * n, 0.0);
    jFlux = Vector( num_jfaces * n, 0.0); 

    // Important Jacobians
    iPlus_A = Vector( num_ifaces * n * n, 0.0);
    iMinus_A = Vector( num_ifaces * n * n, 0.0);
    jPlus_A = Vector( num_jfaces * n * n, 0.0);
    jMinus_A = Vector( num_jfaces * n * n, 0.0);
    irho_A = Vector( num_ifaces * n * n, 0.0);
    jrho_A = Vector( num_jfaces * n * n, 0.0); 

    // These hold the E matrices that never change for the boundary conditions. The first
    // indice is for left or bottom, while the second is for right and top. 
    iE = Vector( 2 * (Nx_local + 1) * n * n, 0.0);
    jE = Vector(2 * (Ny + 1) * n * n, 0.0); 

    // =========== Local geometry data structures for each rank and scattering. ============== 

    // Scatter xCenters from grid into new vectors
    Vector xCenter(num_loc_cells, 0.0);  
    MPI_Scatterv(grid.x_cellCenters.data(), cc_sendcounts.data(), cc_displacements.data(), MPI_DOUBLE,
                 xCenter.data(), num_loc_cells, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    // Scatter yCenters from grid into new vectors
    Vector yCenter(num_loc_cells, 0.0); 
    MPI_Scatterv(grid.y_cellCenters.data(), cc_sendcounts.data(), cc_displacements.data(), MPI_DOUBLE,
                 yCenter.data(), num_loc_cells, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    // Scatter volumes form grid into new vectors
    Vector Volume(num_loc_cells, 0.0); 
    MPI_Scatterv(grid.cellVolumes.data(), cc_sendcounts.data(), cc_displacements.data(), MPI_DOUBLE,
                 Volume.data(), num_loc_cells, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    // For i-faces now
    int iface_offset = 0;
    for (int r = 0; r < size; ++r) {
        if_sendcounts[r]    = (local_Nx[r] + 1) * Ny; // includes overlap
        if_displacements[r] = iface_offset;
        iface_offset     += (local_Nx[r] + 1) * Ny;
    }

    // Scatter i-face normals and areas
    Vector iFxNorm(num_ifaces, 0.0);
    MPI_Scatterv(grid.iface_xNormals.data(), if_sendcounts.data(), if_displacements.data(), MPI_DOUBLE,
                 iFxNorm.data(), num_ifaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    Vector iFyNorm(num_ifaces, 0.0);
    MPI_Scatterv(grid.iface_yNormals.data(), if_sendcounts.data(), if_displacements.data(), MPI_DOUBLE,
                 iFyNorm.data(), num_ifaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    Vector iArea(num_ifaces, 0.0);
    MPI_Scatterv(grid.iAreas.data(), if_sendcounts.data(), if_displacements.data(), MPI_DOUBLE,
                 iArea.data(), num_ifaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);              

    // For j-faces now
    int jface_offset = 0;
    for (int r = 0; r < size; ++r) {
        jf_sendcounts[r]    = local_Nx[r] * (Ny + 1); // includes overlap
        jf_displacements[r] = jface_offset;
        jface_offset     += local_Nx[r] * (Ny + 1);
    }

    // Scatter j-face normals and areas
    Vector jFxNorm(num_jfaces, 0.0);
    MPI_Scatterv(grid.jface_xNormals.data(), jf_sendcounts.data(), jf_displacements.data(), MPI_DOUBLE,
                 jFxNorm.data(), num_jfaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    Vector jFyNorm(num_jfaces, 0.0);
    MPI_Scatterv(grid.jface_yNormals.data(), jf_sendcounts.data(), jf_displacements.data(), MPI_DOUBLE,
                 jFyNorm.data(), num_jfaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 

    Vector jArea(num_jfaces, 0.0);
    MPI_Scatterv(grid.jAreas.data(), jf_sendcounts.data(), jf_displacements.data(), MPI_DOUBLE,
                 jArea.data(), num_jfaces, MPI_DOUBLE,
                 0, MPI_COMM_WORLD); 


    // Intermediate matrices for calculations
    V = Vector(n, 0.0);
    V1 = Vector(n, 0.0);
    V2 = Vector(n, 0.0);
    Q = Vector(n, 0.0);
    W = Vector(n, 0.0);
    int1 = Vector(n * n, 0.0);
    int2 = Vector(n * n, 0.0);
    int3 = Vector(n * n, 0.0);

    // Only rank 0 initializes the full domain with U_inlet initial conditions
    if (rank == 0) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                for (int k = 0; k < n; ++k) {
                    int idx = (i * Ny + j) * n + k; // flattened index: cell-major, then state variable
                    U_gathered[idx] = U_inlet[k];
                }
            }
        }
    }

    // Prepare receive buffer on each rank for local cells (without ghosts)
    Vector recv_buf(Nx_local * Ny * n, 0.0);
    vector<int> u_sendcounts(size), u_displacements(size);

    // Example setup (adjust if you already have this from previous code):
    int u_offset = 0; 
    for (int r = 0; r < size; ++r) {
        u_sendcounts[r] = local_Nx[r] * Ny * n;
        u_displacements[r] = u_offset; // offset in terms of doubles, not cells!
        u_offset += local_Nx[r] * Ny * n;
    }

    // Scatter the U data across ranks (full cells only, no ghosts)
    MPI_Scatterv(U_gathered.data(), u_sendcounts.data(), u_displacements.data(), MPI_DOUBLE,
                recv_buf.data(), Nx_local * Ny * n, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    // Copy scattered data into local U array including ghost cells (offset by 1)
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx_local; ++i) {
            for (int k = 0; k < n; ++k) {
                int idx_recv = (i * Ny + j) * n + k; // index in recv_buf
                int idx_local = ((i + 1) * (Ny + 2) + (j + 1)) * n + k; // index in local U with ghost cells
                U[idx_local] = recv_buf[idx_recv];
            }
        }
    } 


}

Vector Solver2D::inviscid_boundary_U(BCType type, const double* U, double x_norm, double y_norm) {
    Vector ghost(n * n, 0.0); 

    return ghost; 

}

Vector Solver2D::inviscid_boundary_E(BCType type, double x_norm, double y_norm) {

    Vector E(n * n, 0.0); 

    switch (type) {
        
    case BCType::Inlet: 
        return E;
    
    case BCType::Outlet:
        return identity(n); 

    case BCType::IsothermalWall:
        E = {1, 0, 0, 0,
            0, 1 - 2 * x_norm * x_norm, -2 * x_norm * y_norm, 0,
            0, 2 - x_norm * y_norm, 1 - 2 * y_norm * y_norm, 0,
            0, 0, 0, 1};
        return E;

    case BCType::AdiabaticWall:
        E = {1, 0, 0, 0,
            0, 1 - 2 * x_norm * x_norm, -2 * x_norm * y_norm, 0,
            0, 2 - x_norm * y_norm, 1 - 2 * y_norm * y_norm, 0,
            0, 0, 0, 1};
        return E;

    case BCType::Symmetry:
        E = {1, 0, 0, 0,
                0, 1 - 2 * x_norm * x_norm, -2 * x_norm * y_norm, 0,
                0, 2 - x_norm * y_norm, 1 - 2 * y_norm * y_norm, 0,
                0, 0, 0, 1};
        return E;

    default:
        throw invalid_argument("Unknown boundary condition type.");        
    }
}


void Solver2D::exchange_ghost_cells() {
    MPI_Status status_left, status_right;

    int col_size = Ny * n; // Exchange all variables for entire column of Ny cells. 

    Vector send_left(col_size), recv_left(col_size);
    Vector send_right(col_size), recv_right(col_size);

    /** This nested for loop creates vectors for each rank that just holds the innermost
     * real cells of the rank's chunk. These will be used to put into the ghost cell column of
     * the neighboring rank. 
     */


    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < n; ++k) {

            int local_left = (1 * (Ny + 2) + (j + 1)) * n + k;
            int local_right = (Nx_local * (Ny + 2) + (j + 1)) * n + k;

            send_left[j * n + k] = U[local_left];
            send_right[j * n + k] = U[local_right];
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
        MPI_Sendrecv(send_left.data(), col_size, MPI_DOUBLE, rank - 1, 0,
                     recv_left.data(), col_size, MPI_DOUBLE, rank - 1, 1,
                     MPI_COMM_WORLD, &status_left);
    }

    if (rank < size - 1) {
        MPI_Sendrecv(send_right.data(), col_size, MPI_DOUBLE, rank + 1, 1,
                     recv_right.data(), col_size, MPI_DOUBLE, rank + 1, 0,
                     MPI_COMM_WORLD, &status_right);
    }


    /** This nested for loop takes the send and receive buffs and puts the data
     * from it into the ghost cell spots for each rank. 
     */

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < n; ++k) {

            if (rank > 0) {
                int ghost_left = (0 * (Ny + 2) + (j + 1)) * n + k;
                U[ghost_left] = recv_left[j * n + k];
            }

            if (rank < size - 1) {
                int ghost_right = (Nx_local * (Ny + 2) + (j + 1)) * n + k;
                U[ghost_right] = recv_right[j * n + k];
            }

        }
    }


}




int Solver2D::cell_index(int i, int j, int k) {
    return k * (Nx_local + 2) * (Ny + 2)
         + j * (Nx_local + 2)
         + i;
}
int Solver2D::x_index(int i, int j, int k) {
    return k * (Nx_local + 1) * Ny + j * (Nx_local + 1) + i;
}
int Solver2D::y_index(int i, int j, int k) {
    return k * Nx_local * (Ny + 1) + j * Nx_local + i;
}
void Solver2D::compute_fluxes() {




}









