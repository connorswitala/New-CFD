#include "../linalglib/linalg.hpp"

constexpr double gam = 1.4;
constexpr double R = 287.0;


void primtocons(double* U, const double* V) {
    U[0] = V[0];
    U[1] = V[0] * V[1]; 
    U[2] = V[2] / (gam - 1) + V[0] / 2 * (V[1] * V[1]); 
}
void constoprim(double* V, const double* U) {
    V[0] = U[0];
    V[1] = U[1] / U[0];
    V[2] = (U[2] - U[0]/2 * (V[1] * V[1])) * (gam - 1);
}
inline double computeInternalEnergy(const double* U) {
    return U[2] / U[0] - 0.5 * (U[1] * U[1] / (U[0] * U[0]));
}

int main() {

    
    const int Nx = 100;
    const int n = 3;

    Vector A_plus( (Nx + 1) * n * n, 0.0);
    Vector A_minus( (Nx + 1) * n, 0.0);
    Vector Flux((Nx + 1) * n, 0.0);
    Vector F_plus((Nx + 1) * n, 0.0);
    Vector F_minus((Nx + 1) * n, 0.0); 
    Vector U((Nx + 2) * n, 0.0); 

    double L = 1.0;
    double S = L / (Nx - 1); 

    Vector V(n, 0.0);
    Vector VL = {1.0, 0.0, 1.0};
    Vector VR = {0.125, 0.0, 0.1};

    Vector UL(n, 0.0), UR(n, 0.0);
    
    primtocons(UL.data(), VL.data());
    primtocons(UR.data(), VR.data());


    double t = 0.0;

    for (int i = 0; i < Nx; ++i) {

        if (i < Nx / 2) {        
            for (int k = 0; k < n; ++k) {
                U[i * n + k] = UL[k]; 
            }
        }
        else {
            for (int k = 0; k < n; ++k) {
                U[i * n + k] = UR[k];
            }
        }
    }

    const double dt = 1e-5;
    Vector V1(n, 0.0), V2(n, 0.0), Q(n, 0.0), W(n, 0.0);  
    Vector int1(n * n, 0.0), int2(n * n, 0.0), int3(n * n, 0.0), int4(n * n, 0.0);   
    double rho, u, p, a, lp, lm, l, lc, lt, e; 

    cout << "INITIALIZE IS FINE" << endl; 

    while (t <= 0.5) {

        for (int k = 0; k < n; ++k) {
            U[k] = UL[k];              
            U[(Nx + 1) * n + k] = UR[k]; 
        }


        for (int i = 0; i <= Nx; ++i) {
            int idx = i * n; 
            int iidx = (i + 1) * n;
            int aidx = i * n * n;
            int aiidx = (i + 1) * n * n;

            constoprim(V.data(), &U[idx]);

            rho = V[0];
            u = V[1];
            p = V[2]; 
            a = sqrt(gam * p / rho); 
            e = computeInternalEnergy(&U[idx]); 

            // Compute positive flux Jacobian
            lp = 0.5 * (u + a + fabs(u + a));
            lm = 0.5 * (u - a + fabs(u - a));
            l = 0.5 * (u + fabs(u)); 
            lt = 0.5 * (lp - lm);
            lc = 0.5 * (lp + lm - 2 * l);

            for (int k = 0; k < n; ++k) {
                int1[k * n + k] = l;
            }
            
            V1 = {lc / (a * a), (u * lc + a * lt)/(a * a), (U[idx + 2]/U[idx] * lc + a * u * lt)/(a * a) };
            Q = {(gam - 1) * U[idx], -u * (gam - 1), (gam - 1)};

            V2 = {lt / a, u * lt / a + lc, U[idx + 2]/U[idx] * lt / a + u * lc};
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
            e = computeInternalEnergy(&U[iidx]); 

            lp = 0.5 * (u + a - fabs(u + a));
            lm = 0.5 * (u - a - fabs(u - a));
            l = 0.5 * (u - fabs(u)); 
            lt = 0.5 * (lp - lm);
            lc = 0.5 * (lp + lm - 2 * l);

            for (int k = 0; k < n; ++k) {
                int1[k * n + k] = l;
            }
            
            V1 = {lc / (a * a), (u * lc + a * lt)/(a * a), (U[iidx + 2]/U[iidx] * lc + a * u * lt)/(a * a) };
            Q = {(gam - 1) * U[iidx], -u * (gam - 1), (gam - 1)};

            V2 = {lt / a, u * lt / a + lc, U[iidx + 2]/U[iidx] * lt / a + u * lc};
            W = {-u, 1, 0};

            outer_product(V1.data(), Q.data(), int2.data(), n);   
            outer_product(V2.data(), W.data(), int3.data(), n); 

            for (int k = 0; k < n * n; ++k) {
                A_minus[aiidx + k] = int2[k] + int3[k] + int1[k];
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
                U[uidx + k] += dt / S * (Flux[fidx + k] - Flux[fiidx + k]); 
            }

        }


        t += dt;
    

    }

    cout << "PROGRAM FINISHED" << endl;
    return 0;

}