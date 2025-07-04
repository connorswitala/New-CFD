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

// SodSolver1D Constructor
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