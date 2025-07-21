#pragma once

#include <iostream> 
#include <vector>
#include <cassert>
#include <stdexcept>
#include <cmath> 
#include <iomanip>

using namespace std;
constexpr double pi = 3.141592653;

typedef vector<double> Vector;
inline Vector identity(int n) {
    Vector res(n * n, 0.0); 
    for (int i = 0; i < n; ++i) {
        res[i * n + i] = 1.0;
    }
    return res;
}

double dot_product(const double* vec1, const double* vec2, int n);
void outer_product(const double* vec1, const double* vec2, double* res, int n);

void vecvec_add(const double* vec1, const double* vec2, double* res, int n);
void vecvec_sub(const double* vec1, const double* vec2, double* res, int n);


void matvec_mult(const double* mat, const double* vec, double* res, int n);
void matmat_add(const double* mat1, const double* mat2, double* res, int n);
void matmat_sub(const double* mat1, const double* mat2, double* res, int n);
void matrix_divide(double* A, const double* B, double* X, int n, int m);
void matmat_mult(const double* mat1, const double* mat2, double* res, int n); 
void print_matrix(const double* M, int rows, int cols);