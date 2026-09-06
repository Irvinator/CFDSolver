#pragma once

#include "Vector.hpp"
#include "SparseMatrix.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace CFD
{
    struct BiCGSTABResult
    {
        int iterations;
        double finalResidual;
        bool converged;

        void print() const
        {
            std::cout
                << "BiCGSTAB Result:"
                << " | Iterations: " << iterations
                << " | Residual: " << finalResidual
                << " | Converged: "
                << (converged ? "YES" : "NO")
                << "\n";
        }
    };

    class BiCGSTAB
    {
    public:

        explicit BiCGSTAB(
            double tolerance = 1e-8,
            int maxIterations = 10000)
            : tolerance_(tolerance),
            maxIterations_(maxIterations)
        {
        }

        void setTolerance(double tolerance)
        {
            tolerance_ = tolerance;
        }

        void setMaxIterations(int maxIterations)
        {
            maxIterations_ = maxIterations;
        }

        BiCGSTABResult solve(
            const SparseMatrix& A,
            const Vector& b,
            Vector& x,
            bool verbose = false) const
        {
            if (A.rows() != A.cols())
                throw std::runtime_error(
                    "BiCGSTAB: matrix must be square");

            if (A.rows() != b.size())
                throw std::runtime_error(
                    "BiCGSTAB: matrix and RHS size mismatch");

            if (x.size() != b.size())
                x.resize(b.size(), 0.0);

            Vector r = b - A.multiply(x);
            Vector r0 = r;

            double rhoOld = 1.0;
            double alpha = 1.0;
            double omega = 1.0;

            Vector p(b.size(), 0.0);
            Vector v(b.size(), 0.0);

            double residual = std::sqrt(r.dot(r));

            if (residual < tolerance_)
                return { 0, residual, true };

            for (int k = 0; k < maxIterations_; ++k)
            {
                const double rhoNew = r0.dot(r);

                if (std::abs(rhoNew) < 1e-300)
                    throw std::runtime_error(
                        "BiCGSTAB: breakdown");

                const double beta =
                    (rhoNew / rhoOld) * (alpha / omega);

                p = r + (p - v * omega) * beta;

                v = A.multiply(p);

                const double denominator = r0.dot(v);

                if (std::abs(denominator) < 1e-300)
                    throw std::runtime_error(
                        "BiCGSTAB: breakdown");

                alpha = rhoNew / denominator;

                Vector s = r - v * alpha;

                residual = std::sqrt(s.dot(s));

                if (residual < tolerance_)
                {
                    x += p * alpha;

                    return {
                        k + 1,
                        residual,
                        true
                    };
                }

                Vector t = A.multiply(s);

                const double tDotT = t.dot(t);

                if (std::abs(tDotT) < 1e-300)
                    throw std::runtime_error(
                        "BiCGSTAB: breakdown");

                omega = t.dot(s) / tDotT;

                x += p * alpha + s * omega;

                r = s - t * omega;

                residual = std::sqrt(r.dot(r));

                if (verbose)
                {
                    std::cout
                        << "Iteration " << k + 1
                        << " | Residual = "
                        << residual << "\n";
                }

                if (residual < tolerance_)
                {
                    return {
                        k + 1,
                        residual,
                        true
                    };
                }

                if (std::abs(omega) < 1e-300)
                    throw std::runtime_error(
                        "BiCGSTAB: breakdown");

                rhoOld = rhoNew;
            }

            return {
                maxIterations_,
                residual,
                false
            };
        }

    private:

        double tolerance_;
        int maxIterations_;
    };
}
