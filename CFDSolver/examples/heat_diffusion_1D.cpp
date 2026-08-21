#include "linearAlgebra/Vector.hpp"
#include "linearAlgebra/SparseMatrix.hpp"
#include "linearAlgebra/ConjugateGradient.hpp"
#include "fields/ScalarField.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using Clock = std::chrono::high_resolution_clock;

namespace {

    namespace Params {
        constexpr int    N = 100;
        constexpr double L = 1.0;
        constexpr double alpha = 1e-4;
        constexpr double dt = 1.0;
        constexpr double tEnd = 20000.0;
        constexpr double T_hot = 1.0;
        constexpr double T_cold = 0.0;
        constexpr double CG_TOL = 1e-10;
        constexpr int    CG_MAXITER = 10000;
        constexpr int    OUTPUT_EVERY = 750;
        constexpr double STEADY_TOL = 1e-10;
    }

    inline double analyticalSteady(double x) {
        return Params::T_hot - (Params::T_hot - Params::T_cold) * x / Params::L;
    }

    void printHeader(const std::string& title) {
        std::cout << "\n== " << title << ' ';
        int pad = std::max(0, 52 - static_cast<int>(title.size()));
        std::cout << std::string(static_cast<std::size_t>(pad), '=') << "\n";
    }

    void writeCSV(const std::string& filename,
        const std::vector<double>& x,
        const CFD::ScalarField& T,
        double time)
    {
        std::ofstream f(filename);
        if (!f) {
            throw std::runtime_error("Failed to open file: " + filename);
        }

        f << std::fixed << std::setprecision(10);
        f << "x,T,T_steady,error_to_steady,time\n";
        for (int i = 0; i < T.size(); ++i) {
            const double Ta = analyticalSteady(x[i]);
            const double err = std::abs(T[i] - Ta);
            f << x[i] << ','
                << T[i] << ','
                << Ta << ','
                << err << ','
                << time << '\n';
        }
    }

    struct ErrorMetrics {
        double maxErr = 0.0;
        double l2Err = 0.0;
    };

    ErrorMetrics computeSteadyError(const std::vector<double>& x,
        const CFD::ScalarField& T)
    {
        ErrorMetrics e;
        double sumSq = 0.0;

        for (int i = 0; i < T.size(); ++i) {
            const double Ta = analyticalSteady(x[i]);
            const double err = std::abs(T[i] - Ta);
            e.maxErr = std::max(e.maxErr, err);
            sumSq += err * err;
        }

        e.l2Err = std::sqrt(sumSq / static_cast<double>(T.size()));
        return e;
    }

} // namespace

int main() {
    try {
        using namespace CFD;

        const auto tStart = Clock::now();

        std::filesystem::create_directories("output");

        std::cout << "====================================================\n";
        std::cout << "  CFDSolver - 1D Heat Diffusion (Implicit FVM)\n";
        std::cout << "====================================================\n";

        if (Params::N <= 0) {
            throw std::invalid_argument("Params::N must be positive");
        }
        if (Params::L <= 0.0 || Params::alpha <= 0.0 || Params::dt <= 0.0) {
            throw std::invalid_argument("L, alpha, and dt must be positive");
        }
        if (Params::tEnd <= 0.0) {
            throw std::invalid_argument("tEnd must be positive");
        }

        const double dx = Params::L / static_cast<double>(Params::N);
        const double r = Params::alpha * Params::dt / (dx * dx);
        const int nSteps = static_cast<int>(std::ceil(Params::tEnd / Params::dt));
        const double diffusionTimeScale = Params::L * Params::L / Params::alpha;

        std::vector<double> x(Params::N);
        for (int i = 0; i < Params::N; ++i) {
            x[i] = (static_cast<double>(i) + 0.5) * dx;
        }

        printHeader("Setup");
        std::cout << "Cells               : " << Params::N << "\n";
        std::cout << "Domain length       : " << Params::L << "\n";
        std::cout << "dx                  : " << dx << "\n";
        std::cout << "alpha               : " << Params::alpha << "\n";
        std::cout << "dt                  : " << Params::dt << "\n";
        std::cout << "r = alpha*dt/dx^2   : " << r << "\n";
        std::cout << "tEnd                : " << Params::tEnd << "\n";
        std::cout << "nSteps              : " << nSteps << "\n";
        std::cout << "L^2/alpha           : " << diffusionTimeScale << "\n";
        std::cout << "tEnd / (L^2/alpha)  : " << Params::tEnd / diffusionTimeScale << "\n";

        printHeader("Matrix Assembly");

        std::vector<std::size_t> rowIdx;
        std::vector<std::size_t> colIdx;
        std::vector<double> vals;
        rowIdx.reserve(3 * Params::N);
        colIdx.reserve(3 * Params::N);
        vals.reserve(3 * Params::N);

        for (int i = 0; i < Params::N; ++i) {
            if (i == 0) {
                rowIdx.push_back(i);
                colIdx.push_back(i);
                vals.push_back(1.0 + 3.0 * r);

                if (Params::N > 1) {
                    rowIdx.push_back(i);
                    colIdx.push_back(i + 1);
                    vals.push_back(-r);
                }
            }
            else if (i == Params::N - 1) {
                rowIdx.push_back(i);
                colIdx.push_back(i - 1);
                vals.push_back(-r);

                rowIdx.push_back(i);
                colIdx.push_back(i);
                vals.push_back(1.0 + 3.0 * r);
            }
            else {
                rowIdx.push_back(i);
                colIdx.push_back(i - 1);
                vals.push_back(-r);

                rowIdx.push_back(i);
                colIdx.push_back(i);
                vals.push_back(1.0 + 2.0 * r);

                rowIdx.push_back(i);
                colIdx.push_back(i + 1);
                vals.push_back(-r);
            }
        }

        SparseMatrix A(Params::N, Params::N);
        A.setFromTriplets(rowIdx, colIdx, vals);
        A.printInfo();

        ScalarField T(Params::N, 0.0, "Temperature");
        Vector b(Params::N, 0.0);
        Vector Tnew(Params::N, 0.0);
        Vector Told(Params::N, 0.0);

        ConjugateGradient cg(Params::CG_TOL, Params::CG_MAXITER);

        std::vector<double> linearResiduals;
        std::vector<int> cgIters;
        std::vector<double> times;
        std::vector<double> deltaInfHistory;
        linearResiduals.reserve(static_cast<std::size_t>(nSteps));
        cgIters.reserve(static_cast<std::size_t>(nSteps));
        times.reserve(static_cast<std::size_t>(nSteps));
        deltaInfHistory.reserve(static_cast<std::size_t>(nSteps));

        printHeader("Time Integration");

        double time = 0.0;
        bool reachedSteady = false;
        int completedSteps = 0;

        for (int step = 1; step <= nSteps; ++step) {
            for (int i = 0; i < Params::N; ++i) {
                Told[i] = T[i];
                b[i] = Told[i];
            }

            b[0] += 2.0 * r * Params::T_hot;
            b[Params::N - 1] += 2.0 * r * Params::T_cold;

            CGResult res = cg.solve(A, b, Tnew, false);
            if (!res.converged) {
                std::cerr << "Warning: CG did not converge at step " << step << "\n";
            }

            double deltaInf = 0.0;
            for (int i = 0; i < Params::N; ++i) {
                deltaInf = std::max(deltaInf, std::abs(Tnew[i] - Told[i]));
            }

            T.fromVector(Tnew);
            time += Params::dt;
            completedSteps = step;

            linearResiduals.push_back(res.finalResidual);
            cgIters.push_back(res.iterations);
            times.push_back(time);
            deltaInfHistory.push_back(deltaInf);

            if (step % Params::OUTPUT_EVERY == 0 || step == 1) {
                const double maxT = T.max();
                const ErrorMetrics err = computeSteadyError(x, T);

                std::cout << "t=" << std::fixed << std::setprecision(1) << std::setw(9) << time
                    << " | CG iters=" << std::setw(4) << res.iterations
                    << " | linRes=" << std::scientific << std::setprecision(3) << res.finalResidual
                    << " | dTinf=" << deltaInf
                    << " | Tmax=" << std::fixed << std::setprecision(6) << maxT
                    << " | maxErrSteady=" << std::scientific << err.maxErr
                    << "\n";

                writeCSV("output/heat_t" + std::to_string(static_cast<int>(std::round(time))) + ".csv",
                    x, T, time);
            }

            if (step > 5 && deltaInf < Params::STEADY_TOL) {
                std::cout << "\nSteady-state change tolerance reached at t = "
                    << time << "\n";
                reachedSteady = true;
                break;
            }
        }

        printHeader("Validation Against Steady Solution");
        std::cout << std::left
            << std::setw(12) << "x"
            << std::setw(16) << "Numerical"
            << std::setw(16) << "Steady"
            << std::setw(16) << "Abs Error"
            << "\n"
            << std::string(60, '-') << "\n";

        const int stride = std::max(1, Params::N / 10);
        for (int i = 0; i < Params::N; i += stride) {
            const double Ta = analyticalSteady(x[i]);
            const double err = std::abs(T[i] - Ta);
            std::cout << std::fixed << std::setprecision(6)
                << std::left
                << std::setw(12) << x[i]
                << std::setw(16) << T[i]
                << std::setw(16) << Ta
                << std::setw(16) << err
                << '\n';
        }
        if ((Params::N - 1) % stride != 0) {
            const int i = Params::N - 1;
            const double Ta = analyticalSteady(x[i]);
            const double err = std::abs(T[i] - Ta);
            std::cout << std::fixed << std::setprecision(6)
                << std::left
                << std::setw(12) << x[i]
                << std::setw(16) << T[i]
                << std::setw(16) << Ta
                << std::setw(16) << err
                << '\n';
        }

        const ErrorMetrics finalErr = computeSteadyError(x, T);
        const double avgCg = cgIters.empty()
            ? 0.0
            : std::accumulate(cgIters.begin(), cgIters.end(), 0.0) / static_cast<double>(cgIters.size());

        std::cout << std::string(60, '-') << "\n";
        std::cout << "Max error vs steady : " << std::scientific << finalErr.maxErr << "\n";
        std::cout << "L2  error vs steady : " << finalErr.l2Err << "\n";
        std::cout << "Average CG iters    : " << std::fixed << std::setprecision(3) << avgCg << "\n";
        std::cout << "Completed steps     : " << completedSteps << "\n";
        std::cout << "Reached steady      : " << (reachedSteady ? "yes" : "no") << "\n";

        printHeader("Output Files");
        writeCSV("output/heat_final.csv", x, T, time);
        std::cout << "Written: output/heat_final.csv\n";

        {
            std::ofstream f("output/convergence.csv");
            if (!f) {
                throw std::runtime_error("Failed to open file: output/convergence.csv");
            }
            f << "step,time,linear_residual,cg_iterations,deltaT_inf\n";
            f << std::scientific << std::setprecision(8);
            for (std::size_t i = 0; i < times.size(); ++i) {
                f << (i + 1) << ','
                    << times[i] << ','
                    << linearResiduals[i] << ','
                    << cgIters[i] << ','
                    << deltaInfHistory[i] << '\n';
            }
        }
        std::cout << "Written: output/convergence.csv\n";

        const auto tEndClock = Clock::now();
        const double elapsed = std::chrono::duration<double>(tEndClock - tStart).count();

        printHeader("Performance");
        std::cout << "Total runtime [s]   : " << std::fixed << std::setprecision(6) << elapsed << "\n";
        std::cout << "Per completed step  : "
            << (completedSteps > 0 ? 1000.0 * elapsed / static_cast<double>(completedSteps) : 0.0)
            << " ms\n";

        std::cout << "\n====================================================\n";
        if (reachedSteady && finalErr.maxErr < 1e-3) {
            std::cout << "VALIDATION PASSED\n";
        }
        else {
            std::cout << "VALIDATION NOT YET PASSED OR TRANSIENT NOT FINISHED\n";
        }
        std::cout << "Final max error vs steady = " << std::scientific << finalErr.maxErr << "\n";
        std::cout << "====================================================\n";

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << '\n';
        return 1;
    }
}
