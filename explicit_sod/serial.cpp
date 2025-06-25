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



int main(int argc, char** argv) {

    // MPI_Init(&argc, &argv);
    // int rank, size;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &size);

    int counter = 0;
    const int Nx = 4000;
    const int n = 3;
    const double CFL = 0.9; 

    Vector A_plus( (Nx + 1) * n * n, 0.0);
    Vector A_minus( (Nx + 1) * n * n, 0.0);
    Vector Flux((Nx + 1) * n, 0.0);
    Vector F_plus((Nx + 1) * n, 0.0);
    Vector F_minus((Nx + 1) * n, 0.0); 
    Vector U((Nx + 2) * n, 0.0); 

    double L = 2.0;
    double S = L / (Nx - 1); 

    Vector V(n, 0.0);
    Vector VL = {1.0, 0.0, 1.0};
    Vector VR = {0.125, 0.0, 0.1};

    Vector UL(n, 0.0), UR(n, 0.0);
    
    primtocons(UL.data(), VL.data());
    primtocons(UR.data(), VR.data());


    double t = 0.0;

    for (int i = 0; i <= Nx + 1; ++i) {
        for (int k = 0; k < n; ++k) {
            U[i * n + k] = (i <= (Nx + 2) / 2) ? UL[k] : UR[k];
        }
    }

    double dt = compute_max_dt(U, Nx, S, 1.0);;
    Vector V1(n, 0.0), V2(n, 0.0), Q(n, 0.0), W(n, 0.0);  
    Vector int1(n * n, 0.0), int2(n * n, 0.0), int3(n * n, 0.0), int4(n * n, 0.0);   
    double rho, u, p, a, lp, lm, l, lc, lt; 


    while (t <= 0.5) {

        for (int k = 0; k < n; ++k) {
            U[k] = UL[k];              
            U[(Nx + 1) * n + k] = UR[k]; 
        }


        for (int i = 0; i <= Nx; ++i) {
            int idx = i * n; 
            int iidx = (i + 1) * n;
            int aidx = i * n * n;

            constoprim(V.data(), &U[idx]);

            rho = V[0];
            u = V[1];
            p = V[2]; 
            a = sqrt(gam * p / rho); 
            

            // Compute positive flux Jacobian
            lp = 0.5 * (u + a + fabs(u + a));
            lm = 0.5 * (u - a + fabs(u - a));
            l = 0.5 * (u + fabs(u)); 
            lt = 0.5 * (lp - lm);
            lc = 0.5 * (lp + lm - 2 * l);

            for (int k = 0; k < n; ++k) {
                int1[k * n + k] = l;
            }
            
            V1 = {lc / (a * a), (u * lc + a * lt)/(a * a), ( (U[idx + 2] + computePressure(&U[idx]))/U[idx] * lc + a * u * lt)/(a * a) };
            Q = {0.5 * u * u * (gam - 1), -u * (gam - 1), (gam - 1)};

            V2 = {lt / a, u * lt / a + lc, (U[idx + 2] + computePressure(&U[idx]))/U[idx] * lt / a + u * lc};
            W = {-u, 1, 0};

            outer_product(V1.data(), Q.data(), int2.data(), n);   
            outer_product(V2.data(), W.data(), int3.data(), n); 

            for (int k = 0; k < n * n; ++k) {
                A_plus[aidx + k] = int2[k] + int3[k] + int1[k];
            }

            matvec_mult(&A_plus[aidx], &U[idx], &F_plus[idx], n);

           
            // Compute negative flux Jacobian
            constoprim(V.data(), &U[iidx]);

            rho = V[0];
            u = V[1];
            p = V[2]; 
            a = sqrt(gam * p / rho); 

            lp = 0.5 * (u + a - fabs(u + a));
            lm = 0.5 * (u - a - fabs(u - a));
            l = 0.5 * (u - fabs(u)); 
            lt = 0.5 * (lp - lm);
            lc = 0.5 * (lp + lm - 2 * l);

            for (int k = 0; k < n; ++k) {
                int1[k * n + k] = l;
            }
            
            V1 = {lc / (a * a), (u * lc + a * lt)/(a * a), ( (U[iidx + 2] + computePressure(&U[iidx]))/U[iidx] * lc + a * u * lt)/(a * a) };
            Q = {0.5 * u * u * (gam - 1), -u * (gam - 1), (gam - 1)};

            V2 = {lt / a, u * lt / a + lc, (U[iidx + 2] + computePressure(&U[iidx]))/U[iidx] * lt / a + u * lc};
            W = {-u, 1, 0};

            outer_product(V1.data(), Q.data(), int2.data(), n);   
            outer_product(V2.data(), W.data(), int3.data(), n); 

            for (int k = 0; k < n * n; ++k) {
                A_minus[aidx + k] = int2[k] + int3[k] + int1[k];
            }

            matvec_mult(&A_minus[aidx], &U[iidx], &F_minus[idx], n);
            
            for (int k = 0; k < n; ++k) {
                Flux[idx + k] = F_plus[idx + k] + F_minus[idx + k];
            }
        }


        for (int i = 0; i <= Nx; ++i) {
            int fidx = i * n; 
            int fiidx = (i + 1) * n;
            int uidx = (i + 1) * n;
            
            for (int k = 0; k < n; ++k) {
                U[uidx + k] += -dt / S * (-Flux[fidx + k] + Flux[fiidx + k]); 
            }

        }


        t += dt;
        if (counter % 1000 == 0) {
            cout << "t = " << t << "\tdt = " << dt << endl;
        }        
        counter++;

        dt = compute_max_dt(U, Nx, S, CFL);

    }

    string filename = "sod_shock.csv";   
    write_U_to_csv(U, Nx, filename);
    cout << "Program finished!" << endl;

    return 0;

}