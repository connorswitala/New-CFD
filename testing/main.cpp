#include <fstream>
#include <iostream>
#include "../linalglib/linalg.hpp"
using namespace std;

int main() {
    
    const int n = 4, m = 1;
    double A[n * n] = { 4, 3, 0, 0, 3, 4, -1, 0, 0, -1, 4, 2, 0, 0, 2, 4 };
    double B[n * m] = { 24, 30, -24, -24 };
    double X[n * m];

    matrix_divide(A, B, X, n, m);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            cout << X[j * n + i] << " ";
        }
        cout << endl;
    }


    return 0;
}