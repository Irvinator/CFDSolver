#include "linearAlgebra/Vector.hpp"
#include "linearAlgebra/SparseMatrix.hpp"
#include "linearAlgebra/BiCGSTAB.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

    int main()
    {
        try
        {
            std::cout << "==============================\n";
            std::cout << "       BiCGSTAB Test\n";
            std::cout << "==============================\n\n";

            // --------------------------------------------------
            // Matrix:
            //
            // [ 4  1  0 ]
            // [ 1  4  1 ]
            // [ 0  1  4 ]
            //
            // Exact solution:
            //
            // x = [1, 2, 3]
            //
            // Therefore:
            //
            // b = [6, 12, 14]
            // --------------------------------------------------

            CFD::SparseMatrix A(3, 3);

            // COO triplets
            std::vector<std::size_t> rowIdx =
            {
                0, 0,
                1, 1, 1,
                2, 2
            };

            std::vector<std::size_t> colIdx =
            {
                0, 1,
                0, 1, 2,
                1, 2
            };

            std::vector<double> values =
            {
                4.0, 1.0,
                1.0, 4.0, 1.0,
                1.0, 4.0
            };

            // Convert COO triplets to CSR
            A.setFromTriplets(rowIdx, colIdx, values);

            std::cout << "Matrix information:\n";
            A.printInfo();

            std::cout << "\nMatrix:\n";
            A.printDense();

            // --------------------------------------------------
            // RHS
            // --------------------------------------------------

            CFD::Vector b(3);

            b[0] = 6.0;
            b[1] = 12.0;
            b[2] = 14.0;

            // Initial guess
            CFD::Vector x(3, 0.0);

            // --------------------------------------------------
            // Solve
            // --------------------------------------------------

            CFD::BiCGSTAB solver(1e-10, 1000);

            CFD::BiCGSTABResult result =
                solver.solve(A, b, x, true);

            result.print();

            // --------------------------------------------------
            // Expected solution
            // --------------------------------------------------

            const double expected[] =
            {
                1.0,
                2.0,
                3.0
            };

            std::cout << "\nSolution:\n";

            for (std::size_t i = 0; i < x.size(); ++i)
            {
                std::cout
                    << "x[" << i << "] = "
                    << x[i]
                    << "  expected = "
                    << expected[i]
                    << "\n";
            }

            // --------------------------------------------------
            // Verify solution
            // --------------------------------------------------

            const double tolerance = 1e-8;

            for (std::size_t i = 0; i < x.size(); ++i)
            {
                if (std::abs(x[i] - expected[i]) > tolerance)
                {
                    std::cerr
                        << "\nTEST FAILED: incorrect solution.\n";

                    return 1;
                }
            }

            // Check convergence
            if (!result.converged)
            {
                std::cerr
                    << "\nTEST FAILED: solver did not converge.\n";

                return 1;
            }

            // Check the actual residual A*x - b
            CFD::Vector Ax = A.multiply(x);
            CFD::Vector residual = Ax - b;

            double residualNorm =
                std::sqrt(residual.dot(residual));

            std::cout
                << "\nFinal residual norm = "
                << residualNorm
                << "\n";

            if (residualNorm > tolerance)
            {
                std::cerr
                    << "\nTEST FAILED: residual too large.\n";

                return 1;
            }

            std::cout << "\n==============================\n";
            std::cout << "       TEST PASSED!\n";
            std::cout << "==============================\n";

            return 0;
        }
        catch (const std::exception& e)
        {
            std::cerr
                << "\nTEST FAILED: "
                << e.what()
                << "\n";

            return 1;
        }
    }
