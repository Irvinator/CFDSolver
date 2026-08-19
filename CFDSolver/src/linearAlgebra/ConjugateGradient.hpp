#pragma once
#include "Vector.hpp"
#include "SparseMatrix.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace CFD {

    /**
     * Struct to hold the result of a CG solve
     */
    struct CGResult {
        int    iterations;    // Number of iterations taken
        double finalResidual; // Final residual norm
        bool   converged;     // Did it converge?

        void print() const {
            std::cout << "CG Result:"
                << " | Iterations: " << iterations
                << " | Residual: " << finalResidual
                << " | Converged: " << (converged ? "YES" : "NO")
                << "\n";
        }
    };

    /**
     * Conjugate Gradient Solver
     *
     * Solves: A * x = b
     * Requirements: A must be symmetric positive definite (SPD)
     *
     * This is satisfied by your FVM diffusion matrices!
     *
     * Algorithm:
     *   r0 = b - A*x0          (initial residual)
     *   p0 = r0                 (initial search direction)
     *
     *   for k = 0, 1, 2, ...
     *     alpha_k = (r_k . r_k) / (p_k . A*p_k)
     *     x_{k+1} = x_k + alpha_k * p_k
     *     r_{k+1} = r_k - alpha_k * A*p_k
     *     beta_k  = (r_{k+1} . r_{k+1}) / (r_k . r_k)
     *     p_{k+1} = r_{k+1} + beta_k * p_k
     */
    class ConjugateGradient {
    public:

        // ── Constructor ───────────────────────────────────────────────
        explicit ConjugateGradient(double tol = 1e-8,
            int    maxIter = 10000)
            : tol_(tol), maxIter_(maxIter) {
        }

        // ── Setters ───────────────────────────────────────────────────
        void setTolerance(double tol) { tol_ = tol; }
        void setMaxIter(int maxIter) { maxIter_ = maxIter; }

        // ── Main solve function ───────────────────────────────────────
        /**
         * Solves A * x = b
         *
         * @param A       — system matrix (must be SPD)
         * @param b       — right hand side vector
         * @param x       — solution vector (modified in place)
         * @param verbose — print convergence history
         * @return CGResult struct with iterations, residual, converged
         */
        CGResult solve(const SparseMatrix& A,
            const Vector& b,
            Vector& x,
            bool verbose = false) const
        {
            // ── Validate inputs ───────────────────────────────────────
            if (A.rows() != A.cols())
                throw std::runtime_error("CG: Matrix must be square");
            if (A.rows() != b.size())
                throw std::runtime_error("CG: A and b size mismatch");

            // ── Initialise x if wrong size ────────────────────────────
            if (x.size() != b.size())
                x.resize(b.size(), 0.0);

            // ── Step 1: r = b - A*x ───────────────────────────────────
            Vector r = b - A.multiply(x);
            Vector p = r;                    // Search direction

            double rDotR = r.dot(r);         // r^T * r

            // Check if already converged
            if (std::sqrt(rDotR) < tol_)
                return { 0, std::sqrt(rDotR), true };

            if (verbose) {
                std::cout << "\nCG Solver starting...\n";
                std::cout << "Initial residual: "
                    << std::sqrt(rDotR) << "\n";
            }

            // ── Main iteration loop ───────────────────────────────────
            for (int k = 0; k < maxIter_; ++k) {

                // Ap = A * p
                Vector Ap = A.multiply(p);

                // alpha = (r.r) / (p . A*p)
                double pAp = p.dot(Ap);
                if (std::abs(pAp) < 1e-300)
                    throw std::runtime_error(
                        "CG: breakdown — pAp is zero");

                double alpha = rDotR / pAp;

                // x = x + alpha * p
                x += p * alpha;

                // r = r - alpha * Ap
                r -= Ap * alpha;

                // Check convergence
                double rDotRNew = r.dot(r);
                double residual = std::sqrt(rDotRNew);

                if (verbose && (k % 10 == 0)) {
                    std::cout << "  Iter " << k
                        << " | residual = " << residual
                        << "\n";
                }

                if (residual < tol_) {
                    if (verbose)
                        std::cout << "Converged at iteration "
                        << k + 1
                        << " | residual = "
                        << residual << "\n";
                    return { k + 1, residual, true };
                }

                // beta = (r_new . r_new) / (r . r)
                double beta = rDotRNew / rDotR;

                // p = r + beta * p
                p = r + p * beta;

                rDotR = rDotRNew;
            }

            // Did not converge
            double finalRes = std::sqrt(rDotR);
            std::cerr << "CG WARNING: did not converge in "
                << maxIter_ << " iterations."
                << " Final residual: " << finalRes << "\n";
            return { maxIter_, finalRes, false };
        }

    private:
        double tol_;
        int    maxIter_;
    };

} // namespace CFD