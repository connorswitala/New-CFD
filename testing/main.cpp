#include <fstream>
#include <iostream>
#include "../linalglib/linalg.hpp"
using namespace std;



int main() {
    
    const int n = 4, m = 4;
    double A[n * n] = { 4, 3, 0, 0, 3, 4, -1, 0, 0, -1, 4, 2, 0, 0, 2, 4 };
    double B[n * m] = { 5, 6, -1, -2, 17, 14, 19, 1, 0, 0, 2, 4, 15, 20, 20, 3};
    double X[n * m];

    matmat_mult(A, B, X, n);

    print_matrix(X, n, m);



    return 0;
}