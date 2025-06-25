#include "../linalglib/linalg.hpp"
#include <fstream>
#include <string>
#include <mpi.h>

#include <cerrno>  // for errno
#include <cstring> // for strerror

constexpr double gam = 1.4;
constexpr double R = 287.0;


void primtocons(double* U, const double* V) {
    U[0] = V[0];
    U[1] = V[0] * V[1]; 
    U[2] = V[2] / (gam - 1) + 0.5 * V[0] * (V[1] * V[1]); 
}
void constoprim(double* V, const double* U) {
    V[0] = U[0];
    V[1] = U[1] / U[0];
    V[2] = (U[2] - U[0]/2 * (V[1] * V[1])) * (gam - 1);
}
inline double computeInternalEnergy(const double* U) {
    return U[2] / U[0] - 0.5 * (U[1] * U[1] / (U[0] * U[0]));
}
inline double computePressure(const double* U) {
    return (U[2] - 0.5 * U[0] * (U[1] * U[1] / (U[0] * U[0]))) * (gam - 1); 
}
void write_U_to_csv(const Vector& U, int Nx, const string filename) {


    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "FAILED to open file: " << filename << '\n';
        std::cerr << "errno = " << errno << " (" << std::strerror(errno) << ")\n";
    }
  
    file << "rho,u,p\n";
    Vector V(3, 0.0); 

    for (int i = 0; i < Nx + 2; ++i) {
        int idx = i * 3; 
        constoprim(V.data(), &U[idx]);
        double rho = V[0]; 
        double u = V[1];
        double p = V[2];
        file << rho << ',' << u << ',' << p << '\n';
    }

    file.close();
}
double compute_max_dt(const std::vector<double>& U, int Nx, double dx, double cfl) {
    double min_dt = std::numeric_limits<double>::max();

    for (int i = 1; i <= Nx; ++i) { // skip ghost cells
        int idx = i * 3;
        double rho = U[idx + 0];
        double mom = U[idx + 1];
        double E   = U[idx + 2];

        double u = mom / rho;
        double e = E / rho - 0.5 * u * u;
        double p = (gam - 1) * rho * e;
        double a = std::sqrt(gam * p / rho);

        double local_speed = std::abs(u) + a;

        if (local_speed > 1e-8) { // avoid division by zero
            double dt = dx / local_speed;
            min_dt = std::min(min_dt, dt);
        }
    }

    return cfl * min_dt;
}


int main(int argc, char* argv[]) {
    
    MPI_Init(&argc, &argv);

    double start_time = MPI_Wtime(); 
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int counter = 0;
    int Nx_total = 24000;
    const int n = 3;
    const double CFL = 0.9;

    if (argc >= 2) {
        Nx_total = std::stoi(argv[1]);  // get number of cells from command line
        if (rank == 0) cout << "Nx = " << Nx_total << endl;
    }

    int Nx_local = Nx_total / size;
    int start = rank * Nx_local;
    int end = (rank == size - 1) ? Nx_total : (rank + 1) * Nx_local;

    Vector U((Nx_local + 2) * n, 0.0);
    Vector Flux((Nx_local + 1) * n, 0.0);
    Vector F_plus((Nx_local + 1) * n, 0.0);
    Vector F_minus((Nx_local + 1) * n, 0.0);
    Vector A_plus((Nx_local + 1) * n * n, 0.0);
    Vector A_minus((Nx_local + 1) * n * n, 0.0);

    Vector V(n, 0.0), V1(n), V2(n), Q(n), W(n);
    Vector int1(n * n), int2(n * n), int3(n * n);
    Vector UL(n), UR(n);

    Vector VL = {1.0, 0.0, 1.0};
    Vector VR = {0.125, 0.0, 0.1};
    primtocons(UL.data(), VL.data());
    primtocons(UR.data(), VR.data());

    double L = 2.0;
    double S = L / (Nx_total - 1);
    double t = 0.0;
    double dt = 1e-4;

    for (int i = 0; i <= Nx_local + 1; ++i) {
        for (int k = 0; k < n; ++k) {
            U[i * n + k] = (start + i <= Nx_total / 2) ? UL[k] : UR[k];
        }
    }

    while (t <= 0.5) {


        for (int i = 0; i <= Nx_local; ++i) {
            int idx = i * n;
            int iidx = (i + 1) * n;
            int aidx = i * n * n;

            // Positive flux calculation
            constoprim(V.data(), &U[idx]);
            double rho = V[0];
            double u = V[1];
            double p = V[2];
            double a = sqrt(gam * p / rho);

            double lp = 0.5 * (u + a + fabs(u + a));
            double lm = 0.5 * (u - a + fabs(u - a));
            double l = 0.5 * (u + fabs(u));
            double lt = 0.5 * (lp - lm);
            double lc = 0.5 * (lp + lm - 2 * l);

            for (int k = 0; k < n; ++k) int1[k * n + k] = l;

            V1 = {lc / (a * a), (u * lc + a * lt) / (a * a), ((U[idx + 2] + computePressure(&U[idx])) / U[idx] * lc + a * u * lt) / (a * a)};
            Q = {0.5 * u * u * (gam - 1), -u * (gam - 1), (gam - 1)};
            V2 = {lt / a, u * lt / a + lc, (U[idx + 2] + computePressure(&U[idx])) / U[idx] * lt / a + u * lc};
            W = {-u, 1, 0};

            outer_product(V1.data(), Q.data(), int2.data(), n);
            outer_product(V2.data(), W.data(), int3.data(), n);

            for (int k = 0; k < n * n; ++k) A_plus[aidx + k] = int1[k] + int2[k] + int3[k];
            matvec_mult(&A_plus[aidx], &U[idx], &F_plus[idx], n);

            // Negative flux calculation
            constoprim(V.data(), &U[iidx]);
            rho = V[0]; u = V[1]; p = V[2];
            a = sqrt(gam * p / rho);
            lp = 0.5 * (u + a - fabs(u + a));
            lm = 0.5 * (u - a - fabs(u - a));
            l = 0.5 * (u - fabs(u));
            lt = 0.5 * (lp - lm);
            lc = 0.5 * (lp + lm - 2 * l);

            for (int k = 0; k < n; ++k) int1[k * n + k] = l;
            V1 = {lc / (a * a), (u * lc + a * lt) / (a * a), ((U[iidx + 2] + computePressure(&U[iidx])) / U[iidx] * lc + a * u * lt) / (a * a)};
            Q = {0.5 * u * u * (gam - 1), -u * (gam - 1), (gam - 1)};
            V2 = {lt / a, u * lt / a + lc, (U[iidx + 2] + computePressure(&U[iidx])) / U[iidx] * lt / a + u * lc};
            W = {-u, 1, 0};

            outer_product(V1.data(), Q.data(), int2.data(), n);
            outer_product(V2.data(), W.data(), int3.data(), n);
            for (int k = 0; k < n * n; ++k) A_minus[aidx + k] = int1[k] + int2[k] + int3[k];
            matvec_mult(&A_minus[aidx], &U[iidx], &F_minus[idx], n);

            // Total flux calculation
            for (int k = 0; k < n; ++k) Flux[idx + k] = F_plus[idx + k] + F_minus[idx + k];
        }

        // Update U
        for (int i = 0; i <= Nx_local; ++i) {
            int idx = (i + 1) * n;
            int fidx = i * n;
            int fiidx = (i + 1) * n;
            for (int k = 0; k < n; ++k) {
                U[idx + k] += -dt / S * (-Flux[fidx + k] + Flux[fiidx + k]);
            }
        }

        t += dt;
        if (counter % 100 == 0 && rank == 0) {
            std::cout << "t = " << t << "\tdt = " << dt << std::endl;
        }
        counter++;
        dt = compute_max_dt(U, Nx_local, S, CFL);
    }

    // Gather all into U_gathered
    vector<double> U_gathered;
    if (rank == 0) U_gathered.resize((Nx_total + 2) * n);

    MPI_Gather(U.data() + n, Nx_local * n, MPI_DOUBLE, U_gathered.data() + n, Nx_local * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        write_U_to_csv(U_gathered, Nx_total, "sod_shock_mpi.csv");
        std::cout << "Program finished!" << std::endl;
    }

    double end_time = MPI_Wtime();    // End timer
    double elapsed = end_time - start_time;

    if (rank == 0) {
        std::cout << "Elapsed time: " << elapsed << " seconds\n";
    }

    MPI_Finalize();
    return 0;
}


