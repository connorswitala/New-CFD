
#include "linalg.hpp"

// Matrix vector arithmetic
void matvec_mult(const double* mat, const double* vec, double* res, int n) {
    for (int i = 0; i < n; ++i) {
        res[i] = 0.0;
        for (int j = 0; j < n; ++j) {
            res[i] += mat[i * n + j] * vec[j];
        }
    }
}

// Pure vector arithmetic
double dot_product(const double* vec1, const double* vec2, int n) {
    double result = 0.0;
    for (int i = 0; i < n; ++i) {
        result += vec1[i] * vec2[i];
    }
    return result; 
}
void vecvec_add(const double* vec1, const double* vec2, double* res, int n) {
    for (int i = 0; i < n; ++i) {
        res[i] = vec1[i] + vec2[i];
    }
}
void outer_product(const double* vec1, const double* vec2, double* res, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            res[i * n + j] = vec1[i] * vec2[j];
        }
    }
}
void vecvec_sub(const double* vec1, const double* vec2, double* res, int n) {
    for (int i = 0; i < n; ++i) {
        res[i] = vec1[i] + vec2[i];
    }
}


// Pure matrix arithmetic
void matmat_add(const double* mat1, const double* mat2, double* res, int n) {

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            res[i * n + j] = mat1[i * n + j] + mat2[i * n + j];
        }
    }

}
void matmat_sub(const double* mat1, const double* mat2, double* res, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            res[i * n + j] = mat1[i * n + j] - mat2[i * n + j];
        }
    }

}
void matmat_mult(const double* mat1, const double* mat2, double* res, int n) {  
    for (int i = 0; i < n; ++i) {           // row of A
        for (int j = 0; j < n; ++j) {       // column of B
            res[i * n + j] = 0.0;
            for (int k = 0; k < n; ++k) {   // shared index
                res[i * n + j] += mat1[i * n + k] * mat2[k * n + j];  
            }
        }
    }
}
