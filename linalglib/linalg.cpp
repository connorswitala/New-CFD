
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

void matrix_divide(double* A, const double* B, double* X, int n, int m) {
    // Copy A because LU is in-place and we don't want to destroy the original
    std::vector<double> LU(A, A + n * n);

    // In-place LU Decomposition (no pivoting)
    for (int k = 0; k < n; ++k) {
        for (int i = k + 1; i < n; ++i) {
            LU[i * n + k] /= LU[k * n + k];
            for (int j = k + 1; j < n; ++j) {
                LU[i * n + j] -= LU[i * n + k] * LU[k * n + j];
            }
        }
    }

    std::vector<double> y(n, 0.0);

    for (int col = 0; col < m; ++col) {
        // Forward Substitution: Ly = b
        for (int i = 0; i < n; ++i) {
            y[i] = B[i * m + col];  // row-major access
            for (int j = 0; j < i; ++j)
                y[i] -= LU[i * n + j] * y[j];
        }

        // Backward Substitution: Ux = y
        for (int i = n - 1; i >= 0; --i) {
            double sum = y[i];
            for (int j = i + 1; j < n; ++j)
                sum -= LU[i * n + j] * X[j * m + col];  // row-major
            X[i * m + col] = sum / LU[i * n + i];
        }
    }
}

void print_matrix(const double* M, int rows, int cols) {
    std::cout << std::fixed << std::setprecision(4);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << std::setw(10) << M[i * cols + j] << " ";
        }
        std::cout << "\n";
    }
}