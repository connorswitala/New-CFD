#include "../linalglib/linalg.hpp"

int main() {

    const int n = 4, Nx = 100, Ny = 100;

    Vector A(Nx * Ny * n * n, 0.0);
    Vector U(Nx * Ny * n, 0.0);
    Vector F(Nx * Ny * n, 0.0); 

    int i = 42, j = 17;
    int cell = i * Ny + j; 

    for (int k = 0; k < n; ++k) {
        cout << F[cell * n + k] << endl;
    }


    for (int k = 0; k < 16; ++k) { 
        A[cell * n * n + k] = k; 
    }

    // Fill vector U at (i,j)
    for (int k = 0; k < n; ++k) {
        U[cell * n + k] = 1.0 * k;  // example values
    }

    matvec_mult(&A[cell * n * n], &U[cell * n], &F[cell * n], n);

    for (int k = 0; k < n; ++k) {
        cout << F[cell * n + k] << endl;
    }

    return 0;
}