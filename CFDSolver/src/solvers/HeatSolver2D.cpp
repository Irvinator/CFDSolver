#include "HeatSolver2D.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace CFD {

    namespace {

        double analyticalLinearX(double x, double width, double T_west, double T_east)
        {
            return T_west * (1.0 - x / width) + T_east * (x / width);
        }

    } // namespace

    HeatSolver2D::HeatSolver2D(const Mesh& mesh,
        double alpha,
        double dt,
        double tEnd)
        : mesh_(mesh), alpha_(alpha), dt_(dt), tEnd_(tEnd),
        T_(mesh.numberOfCells(), 0.0, "Temperature"),
        A_(mesh.numberOfCells(), mesh.numberOfCells())
    {
        validateInputs();

        rx_ = alpha_ * dt_ / (mesh_.getDx() * mesh_.getDx());
        ry_ = alpha_ * dt_ / (mesh_.getDy() * mesh_.getDy());

        printSetup();
        assembleMatrix();
    }

    void HeatSolver2D::validateInputs() const
    {
        if (alpha_ <= 0.0) {
            throw std::invalid_argument("HeatSolver2D: alpha must be positive");
        }
        if (dt_ <= 0.0) {
            throw std::invalid_argument("HeatSolver2D: dt must be positive");
        }
        if (tEnd_ <= 0.0) {
            throw std::invalid_argument("HeatSolver2D: tEnd must be positive");
        }
        if (mesh_.numberOfCells() <= 0) {
            throw std::invalid_argument("HeatSolver2D: mesh must contain cells");
        }
    }

    void HeatSolver2D::setBCs(const BoundaryConditions& bcs)
    {
        bcs_ = bcs;
        std::cout << "Boundary conditions:\n"
            << "  West  = " << bcs_.T_west << "\n"
            << "  East  = " << bcs_.T_east << "\n"
            << "  South = " << bcs_.T_south << "\n"
            << "  North = " << bcs_.T_north << "\n";
    }

    void HeatSolver2D::setIC(double T0)
    {
        T_.fill(T0);
        std::cout << "Initial condition: uniform T = " << T0 << "\n";
    }

    void HeatSolver2D::setIC(const std::function<double(double, double)>& f)
    {
        for (int idx = 0; idx < mesh_.numberOfCells(); ++idx) {
            const Mesh::Cell& c = mesh_.getCell(idx);
            T_[idx] = f(c.x, c.y);
        }
        std::cout << "Initial condition: custom function\n";
    }

    void HeatSolver2D::setOutputFreq(int freq)
    {
        if (freq <= 0) {
            throw std::invalid_argument("HeatSolver2D: output frequency must be positive");
        }
        outputFreq_ = freq;
    }

    void HeatSolver2D::setSteadyTolerance(double tol)
    {
        if (tol <= 0.0) {
            throw std::invalid_argument("HeatSolver2D: steady tolerance must be positive");
        }
        steadyTol_ = tol;
    }

    void HeatSolver2D::enableSteadyStop(bool enable)
    {
        stopAtSteady_ = enable;
    }

    void HeatSolver2D::assembleMatrix()
    {
        const int N = mesh_.numberOfCells();

        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> vals;
        rows.reserve(static_cast<std::size_t>(5 * N));
        cols.reserve(static_cast<std::size_t>(5 * N));
        vals.reserve(static_cast<std::size_t>(5 * N));

        for (int idx = 0; idx < N; ++idx) {
            const Mesh::Cell& c = mesh_.getCell(idx);

            double aP = 1.0;

            if (c.west != -1) {
                rows.push_back(static_cast<std::size_t>(idx));
                cols.push_back(static_cast<std::size_t>(c.west));
                vals.push_back(-rx_);
                aP += rx_;
            }
            else {
                aP += 2.0 * rx_;
            }

            if (c.east != -1) {
                rows.push_back(static_cast<std::size_t>(idx));
                cols.push_back(static_cast<std::size_t>(c.east));
                vals.push_back(-rx_);
                aP += rx_;
            }
            else {
                aP += 2.0 * rx_;
            }

            if (c.south != -1) {
                rows.push_back(static_cast<std::size_t>(idx));
                cols.push_back(static_cast<std::size_t>(c.south));
                vals.push_back(-ry_);
                aP += ry_;
            }
            else {
                aP += 2.0 * ry_;
            }

            if (c.north != -1) {
                rows.push_back(static_cast<std::size_t>(idx));
                cols.push_back(static_cast<std::size_t>(c.north));
                vals.push_back(-ry_);
                aP += ry_;
            }
            else {
                aP += 2.0 * ry_;
            }

            rows.push_back(static_cast<std::size_t>(idx));
            cols.push_back(static_cast<std::size_t>(idx));
            vals.push_back(aP);
        }

        A_.setFromTriplets(rows, cols, vals);

        std::cout << "Matrix assembled: ";
        A_.printInfo();
    }

    Vector HeatSolver2D::buildRHS() const
    {
        const int N = mesh_.numberOfCells();
        Vector b(static_cast<std::size_t>(N), 0.0);

        for (int idx = 0; idx < N; ++idx) {
            const Mesh::Cell& c = mesh_.getCell(idx);
            b[static_cast<std::size_t>(idx)] = T_[idx];

            if (c.westBoundary) {
                b[static_cast<std::size_t>(idx)] += 2.0 * rx_ * bcs_.T_west;
            }
            if (c.eastBoundary) {
                b[static_cast<std::size_t>(idx)] += 2.0 * rx_ * bcs_.T_east;
            }
            if (c.southBoundary) {
                b[static_cast<std::size_t>(idx)] += 2.0 * ry_ * bcs_.T_south;
            }
            if (c.northBoundary) {
                b[static_cast<std::size_t>(idx)] += 2.0 * ry_ * bcs_.T_north;
            }
        }

        return b;
    }

    void HeatSolver2D::step(ConjugateGradient& cg, std::vector<int>& cgIters)
    {
        const Vector b = buildRHS();
        Vector Tnew = T_.toVector();
        const Vector Told = T_.toVector();

        const CGResult res = cg.solve(A_, b, Tnew, false);
        if (!res.converged) {
            throw std::runtime_error("HeatSolver2D: Conjugate Gradient failed to converge");
        }

        lastDeltaInf_ = 0.0;
        for (std::size_t i = 0; i < Tnew.size(); ++i) {
            lastDeltaInf_ = std::max(lastDeltaInf_, std::abs(Tnew[i] - Told[i]));
        }

        T_.fromVector(Tnew);
        time_ += dt_;
        ++step_;
        cgIters.push_back(res.iterations);

        if (stopAtSteady_ && step_ > 5 && lastDeltaInf_ < steadyTol_) {
            reachedSteady_ = true;
        }

        if (step_ % outputFreq_ == 0 || step_ == 1) {
            printProgress(res.iterations, res.finalResidual);
        }
    }

    HeatSolver2D::SolverResult HeatSolver2D::run(const std::string& outputDir)
    {
        using Clock = std::chrono::high_resolution_clock;
        const auto tStart = Clock::now();

        if (!outputDir.empty()) {
            std::filesystem::create_directories(outputDir);
        }

        const int estimatedSteps = static_cast<int>(std::ceil(tEnd_ / dt_));
        ConjugateGradient cg(1e-10, 10000);
        std::vector<int> cgIters;
        cgIters.reserve(static_cast<std::size_t>(estimatedSteps));

        std::cout << "\n== Time Loop =========================================\n"
            << "Estimated steps    : " << estimatedSteps << "\n"
            << "Final time target  : " << tEnd_ << "\n";

        while (!finished()) {
            step(cg, cgIters);

            if (!outputDir.empty() && (step_ % outputFreq_ == 0 || step_ == 1)) {
                writeCSV(outputDir + "/heat2D_t" + std::to_string(static_cast<int>(std::round(time_))) + ".csv");
            }

            if (reachedSteady_) {
                std::cout << "Steady-state change tolerance reached at t = " << time_ << "\n";
                break;
            }
        }

        if (!outputDir.empty()) {
            writeCSV(outputDir + "/heat2D_final.csv");
        }

        const auto tStop = Clock::now();
        const double elapsed = std::chrono::duration<double>(tStop - tStart).count();

        SolverResult result;
        result.runtime = elapsed;
        result.totalSteps = step_;
        result.reachedSteady = reachedSteady_;
        result.finalDeltaInf = lastDeltaInf_;
        result.avgCgIterations = cgIters.empty()
            ? 0.0
            : std::accumulate(cgIters.begin(), cgIters.end(), 0.0) / static_cast<double>(cgIters.size());

        std::cout << "\n== Summary ===========================================\n"
            << "Runtime [s]        : " << elapsed << "\n"
            << "Steps              : " << step_ << "\n"
            << "Reached steady     : " << (reachedSteady_ ? "yes" : "no") << "\n"
            << "Final dTinf        : " << lastDeltaInf_ << "\n"
            << "Average CG iters   : " << result.avgCgIterations << "\n"
            << "T min              : " << T_.min() << "\n"
            << "T max              : " << T_.max() << "\n"
            << "T mean             : " << T_.mean() << "\n";

        return validate(result);
    }

    const ScalarField& HeatSolver2D::T() const
    {
        return T_;
    }

    double HeatSolver2D::time() const
    {
        return time_;
    }

    int HeatSolver2D::steps() const
    {
        return step_;
    }

    bool HeatSolver2D::finished() const
    {
        return time_ >= tEnd_;
    }

    void HeatSolver2D::writeCSV(const std::string& filename) const
    {
        std::ofstream f(filename);
        if (!f) {
            throw std::runtime_error("HeatSolver2D: failed to open output file: " + filename);
        }

        f << std::fixed << std::setprecision(10);
        f << "x,y,T,time\n";

        for (int idx = 0; idx < mesh_.numberOfCells(); ++idx) {
            const Mesh::Cell& c = mesh_.getCell(idx);
            f << c.x << ','
                << c.y << ','
                << T_[idx] << ','
                << time_ << '\n';
        }
    }
    bool HeatSolver2D::canValidateLinearXCase() const
    {
        return false;
    }
    HeatSolver2D::SolverResult HeatSolver2D::validate(SolverResult result) const
    {
        std::cout << "\n== Validation ========================================\n";

        if (!canValidateLinearXCase()) {
            std::cout << "Case               : skipped\n";
            std::cout << "Reason             : no compatible analytical solution "
                "implemented for current 2D Dirichlet boundary conditions\n";
            result.validated = false;
            result.maxError = 0.0;
            result.l2Error = 0.0;
            return result;
        }


        const double W = mesh_.getWidth();
        double maxErr = 0.0;
        double sumSq = 0.0;
        const int N = mesh_.numberOfCells();

        for (int idx = 0; idx < N; ++idx) {
            const Mesh::Cell& c = mesh_.getCell(idx);
            const double Ta = analyticalLinearX(c.x, W, bcs_.T_west, bcs_.T_east);
            const double err = std::abs(T_[idx] - Ta);
            maxErr = std::max(maxErr, err);
            sumSq += err * err;
        }

        result.maxError = maxErr;
        result.l2Error = std::sqrt(sumSq / static_cast<double>(N));
        result.validated = (maxErr < 1e-3);

        std::cout << "\n== Validation ========================================\n"
            << "Case               : linear steady state in x\n"
            << "Max error          : " << result.maxError << "\n"
            << "L2 error           : " << result.l2Error << "\n"
            << "Validation         : " << (result.validated ? "passed" : "failed") << "\n";

        return result;
    }

    void HeatSolver2D::printSetup() const
    {
        std::cout << "\n== HeatSolver2D Setup ================================\n"
            << "Grid               : " << mesh_.getNx() << " x " << mesh_.getNy() << "\n"
            << "Cells              : " << mesh_.numberOfCells() << "\n"
            << "dx                 : " << mesh_.getDx() << "\n"
            << "dy                 : " << mesh_.getDy() << "\n"
            << "alpha              : " << alpha_ << "\n"
            << "dt                 : " << dt_ << "\n"
            << "tEnd               : " << tEnd_ << "\n"
            << "rx                 : " << rx_ << "\n"
            << "ry                 : " << ry_ << "\n";
    }

    void HeatSolver2D::printProgress(int cgIters, double cgResidual) const
    {
        std::cout << "t=" << std::fixed << std::setprecision(4) << std::setw(10) << time_
            << " | step=" << std::setw(6) << step_
            << " | CG iters=" << std::setw(4) << cgIters
            << " | CG res=" << std::scientific << std::setprecision(3) << cgResidual
            << " | Tmin=" << std::fixed << std::setprecision(6) << T_.min()
            << " | Tmax=" << T_.max()
            << " | Tmean=" << T_.mean()
            << "\n";
    }

} // namespace CFD
