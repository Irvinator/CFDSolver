#pragma once

#include "../mesh/mesh2D.h"
#include "../fields/ScalarField.hpp"
#include "../linearAlgebra/SparseMatrix.hpp"
#include "../linearAlgebra/ConjugateGradient.hpp"

#include <filesystem>
#include <functional>
#include <string>
#include <vector>

namespace CFD {

    class HeatSolver2D {
    public:
        struct BoundaryConditions {
            double T_west = 0.0;
            double T_east = 0.0;
            double T_south = 0.0;
            double T_north = 0.0;
        };

        struct SolverResult {
            double maxError = 0.0;
            double l2Error = 0.0;
            double runtime = 0.0;
            int    totalSteps = 0;
            bool   validated = false;
            bool   reachedSteady = false;
            double finalDeltaInf = 0.0;
            double avgCgIterations = 0.0;
        };

        HeatSolver2D(const Mesh& mesh,
            double alpha,
            double dt,
            double tEnd);

        void setBCs(const BoundaryConditions& bcs);
        void setIC(double T0);
        void setIC(const std::function<double(double, double)>& f);
        void setOutputFreq(int freq);
        void setSteadyTolerance(double tol);
        void enableSteadyStop(bool enable);

        SolverResult run(const std::string& outputDir = "");

        const ScalarField& T() const;
        double time() const;
        int steps() const;
        bool finished() const;

        void writeCSV(const std::string& filename) const;

    private:
        const Mesh& mesh_;

        double alpha_ = 0.0;
        double dt_ = 0.0;
        double tEnd_ = 0.0;

        double rx_ = 0.0;
        double ry_ = 0.0;

        double time_ = 0.0;
        int    step_ = 0;
        int    outputFreq_ = 10;

        double steadyTol_ = 1e-10;
        bool   stopAtSteady_ = true;
        bool   reachedSteady_ = false;
        double lastDeltaInf_ = 0.0;

        BoundaryConditions bcs_;
        ScalarField        T_;
        SparseMatrix       A_;

        void validateInputs() const;
        void assembleMatrix();
        Vector buildRHS() const;
        void step(ConjugateGradient& cg, std::vector<int>& cgIters);

        SolverResult validate(SolverResult result) const;
        bool canValidateLinearXCase() const;
        void printSetup() const;
        void printProgress(int cgIters, double cgResidual) const;
    };

} // namespace CFD
