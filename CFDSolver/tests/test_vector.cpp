#include "linearAlgebra/Vector.hpp"
#include "linearAlgebra/ConjugateGradient.hpp"
#include "linearAlgebra/SparseMatrix.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace CFD;

// ══════════════════════════════════════════════
// VECTOR TESTS
// ══════════════════════════════════════════════

void test_construction() {
    Vector v(5, 1.0);
    assert(v.size() == 5);
    for (std::size_t i = 0; i < v.size(); ++i)
        assert(v[i] == 1.0);
    std::cout << "PASSED: Vector construction\n";
}

void test_arithmetic() {
    Vector a(3, 2.0);
    Vector b(3, 3.0);
    Vector c = a + b;
    assert(c[0] == 5.0);
    Vector d = b - a;
    assert(d[0] == 1.0);
    Vector e = a * 3.0;
    assert(e[0] == 6.0);
    std::cout << "PASSED: Vector arithmetic\n";
}

void test_dot_product() {
    Vector a = { 1.0, 2.0, 3.0 };
    Vector b = { 4.0, 5.0, 6.0 };
    double result = a.dot(b);
    assert(result == 32.0);
    std::cout << "PASSED: Vector dot product\n";
}

void test_norm() {
    Vector v = { 3.0, 4.0 };
    double n = v.norm2();
    assert(std::abs(n - 5.0) < 1e-10);
    std::cout << "PASSED: Vector L2 norm\n";
}

// ══════════════════════════════════════════════
// SPARSE MATRIX TESTS
// ══════════════════════════════════════════════

void test_sparse_construction() {
    SparseMatrix A(4, 4);
    assert(A.rows() == 4);
    assert(A.cols() == 4);
    assert(A.nnz() == 0);
    std::cout << "PASSED: SparseMatrix construction\n";
}

void test_sparse_assembly() {
    /**
     * Build this matrix:
     * | 4  0  0  2 |
     * | 0  3  0  0 |
     * | 0  0  7  1 |
     * | 0  0  0  5 |
     */
    SparseMatrix A(4, 4);

    std::vector<std::size_t> rows = { 0, 0, 1, 2, 2, 3 };
    std::vector<std::size_t> cols = { 0, 3, 1, 2, 3, 3 };
    std::vector<double>      vals = { 4, 2, 3, 7, 1, 5 };

    A.setFromTriplets(rows, cols, vals);

    assert(A.nnz() == 6);
    A.printInfo();
    A.printDense();
    std::cout << "PASSED: SparseMatrix assembly\n";
}

void test_sparse_multiply() {
    /**
     * Matrix:
     * | 2  0 |
     * | 0  3 |
     *
     * Vector x = [1, 2]
     * Expected: A*x = [2, 6]
     */
    SparseMatrix A(2, 2);

    std::vector<std::size_t> rows = { 0, 1 };
    std::vector<std::size_t> cols = { 0, 1 };
    std::vector<double>      vals = { 2.0, 3.0 };

    A.setFromTriplets(rows, cols, vals);

    Vector x = { 1.0, 2.0 };
    Vector y = A * x;

    assert(std::abs(y[0] - 2.0) < 1e-10);
    assert(std::abs(y[1] - 6.0) < 1e-10);
    std::cout << "PASSED: SparseMatrix multiply\n";
}

void test_sparse_tridiagonal() {
    /**
     * Classic tridiagonal matrix (like heat diffusion!):
     * |  2 -1  0  0 |
     * | -1  2 -1  0 |
     * |  0 -1  2 -1 |
     * |  0  0 -1  2 |
     *
     * x = [1, 1, 1, 1]
     * Expected: A*x = [1, 0, 0, 1]
     */
    int N = 4;
    SparseMatrix A(N, N);

    std::vector<std::size_t> rows, cols;
    std::vector<double>      vals;

    for (int i = 0; i < N; ++i) {
        // Diagonal
        rows.push_back(i); cols.push_back(i);
        vals.push_back(2.0);
        // Lower diagonal
        if (i > 0) {
            rows.push_back(i); cols.push_back(i - 1);
            vals.push_back(-1.0);
        }
        // Upper diagonal
        if (i < N - 1) {
            rows.push_back(i); cols.push_back(i + 1);
            vals.push_back(-1.0);
        }
    }

    A.setFromTriplets(rows, cols, vals);
    A.printDense();

    Vector x(N, 1.0);
    Vector y = A * x;

    assert(std::abs(y[0] - 1.0) < 1e-10);
    assert(std::abs(y[1] - 0.0) < 1e-10);
    assert(std::abs(y[2] - 0.0) < 1e-10);
    assert(std::abs(y[3] - 1.0) < 1e-10);
    std::cout << "PASSED: SparseMatrix tridiagonal\n";
}

void test_diagonal_extraction() {
    SparseMatrix A(3, 3);

    std::vector<std::size_t> rows = { 0, 1, 2 };
    std::vector<std::size_t> cols = { 0, 1, 2 };
    std::vector<double>      vals = { 5.0, 8.0, 3.0 };

    A.setFromTriplets(rows, cols, vals);
    Vector d = A.diagonal();

    assert(std::abs(d[0] - 5.0) < 1e-10);
    assert(std::abs(d[1] - 8.0) < 1e-10);
    assert(std::abs(d[2] - 3.0) < 1e-10);
    std::cout << "PASSED: SparseMatrix diagonal extraction\n";
}
// ══════════════════════════════════════════════
// CONJUGATE GRADIENT TESTS
// ══════════════════════════════════════════════

void test_cg_diagonal() {
    /**
     * Simplest possible test — diagonal matrix:
     * | 4  0 |   | x1 |   | 8  |
     * | 0  3 | * | x2 | = | 6  |
     *
     * Exact solution: x = [2, 2]
     */
    SparseMatrix A(2, 2);
    std::vector<std::size_t> rows = { 0, 1 };
    std::vector<std::size_t> cols = { 0, 1 };
    std::vector<double>      vals = { 4.0, 3.0 };
    A.setFromTriplets(rows, cols, vals);

    Vector b = { 8.0, 6.0 };
    Vector x(2, 0.0);

    ConjugateGradient cg(1e-10, 100);
    CGResult result = cg.solve(A, b, x);

    assert(std::abs(x[0] - 2.0) < 1e-8);
    assert(std::abs(x[1] - 2.0) < 1e-8);
    assert(result.converged);
    std::cout << "PASSED: CG diagonal system\n";
}

void test_cg_tridiagonal() {
    /**
     * Tridiagonal system — exactly what heat diffusion produces!
     *
     * |  2 -1  0 |   | x1 |   | 1 |
     * | -1  2 -1 | * | x2 | = | 0 |
     * |  0 -1  2 |   | x3 |   | 1 |
     *
     * Exact solution: x = [1, 1, 1]
     */
    int N = 3;
    SparseMatrix A(N, N);

    std::vector<std::size_t> rows, cols;
    std::vector<double>      vals;

    for (int i = 0; i < N; ++i) {
        rows.push_back(i); cols.push_back(i);
        vals.push_back(2.0);
        if (i > 0) {
            rows.push_back(i); cols.push_back(i - 1);
            vals.push_back(-1.0);
        }
        if (i < N - 1) {
            rows.push_back(i); cols.push_back(i + 1);
            vals.push_back(-1.0);
        }
    }

    A.setFromTriplets(rows, cols, vals);

    Vector b = { 1.0, 0.0, 1.0 };
    Vector x(N, 0.0);

    ConjugateGradient cg(1e-10, 1000);
    CGResult result = cg.solve(A, b, x, true);
    result.print();

    assert(std::abs(x[0] - 1.0) < 1e-8);
    assert(std::abs(x[1] - 1.0) < 1e-8);
    assert(std::abs(x[2] - 1.0) < 1e-8);
    assert(result.converged);
    std::cout << "PASSED: CG tridiagonal system\n";
}

void test_cg_larger_system() {
    /**
     * Larger tridiagonal system — N=50
     * Tests performance and convergence on a
     * realistic sized problem
     *
     * Same pattern as 1D heat diffusion matrix!
     */
    int N = 50;
    SparseMatrix A(N, N);

    std::vector<std::size_t> rows, cols;
    std::vector<double>      vals;

    for (int i = 0; i < N; ++i) {
        rows.push_back(i); cols.push_back(i);
        vals.push_back(2.0);
        if (i > 0) {
            rows.push_back(i); cols.push_back(i - 1);
            vals.push_back(-1.0);
        }
        if (i < N - 1) {
            rows.push_back(i); cols.push_back(i + 1);
            vals.push_back(-1.0);
        }
    }

    A.setFromTriplets(rows, cols, vals);
    A.printInfo();

    // RHS — all ones
    Vector b(N, 1.0);
    Vector x(N, 0.0);

    ConjugateGradient cg(1e-10, 1000);
    CGResult result = cg.solve(A, b, x);
    result.print();

    assert(result.converged);
    std::cout << "PASSED: CG larger system (N=50)\n";
}
// ══════════════════════════════════════════════
// MAIN
// ══════════════════════════════════════════════

int main() {
    std::cout << "\n== CFD Vector Tests ==\n";
    test_construction();
    test_arithmetic();
    test_dot_product();
    test_norm();

    std::cout << "\n== CFD SparseMatrix Tests ==\n";
    test_sparse_construction();
    test_sparse_assembly();
    test_sparse_multiply();
    test_sparse_tridiagonal();
    test_diagonal_extraction();

    std::cout << "\n== CFD ConjugateGradient Tests ==\n";
    test_cg_diagonal();
    test_cg_tridiagonal();
    test_cg_larger_system();

    std::cout << "\nAll tests passed!\n";
    return 0;
}