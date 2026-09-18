#pragma once

#include "mesh/mesh2D.h"
#include "BCs/BC.h"
#include "linearAlgebra/SparseMatrix.hpp"
#include "linearAlgebra/Vector.hpp"
#include "fields/StaggeredFields.h"

#include <cstddef>
#include <vector>

namespace CFD
{
    class StaggeredSIMPLE
    {
    private:

        Mesh& mesh;
        StaggeredFields& fields;

        BoundaryCondition& northBC;
        BoundaryCondition& southBC;
        BoundaryCondition& eastBC;
        BoundaryCondition& westBC;

        // ============================================================
        // LINEAR SYSTEMS
        // ============================================================

        SparseMatrix uMatrix;
        Vector uRHS;

        SparseMatrix vMatrix;
        Vector vRHS;

        SparseMatrix pressureMatrix;
        Vector pressureRHS;
        Vector pressureCorrection;

        // ============================================================
        // SOLVER SETTINGS
        // ============================================================

        double relaxationPressure{ 0.3 };
        double relaxationVelocity{ 0.7 };

        double rho{ 1.0 };
        double mu{ 0.01 };

        // Continuity is in mass-flow units. Momentum residuals are
        // dimensionless relative algebraic L1 residuals.
        double convergenceTolerance{ 1.0e-6 };
        double momentumTolerance_{ 1.0e-6 };

        std::size_t maxIterations{ 1000 };

        // ============================================================
        // ITERATION STATE
        // ============================================================

        std::size_t iteration{ 0 };

        bool converged_{ false };
        bool finished_{ false };

        // Legacy alias for continuityResidual_.
        double residual{ 0.0 };

        double continuityResidual_{ 0.0 };
        double uMomentumResidual_{ 0.0 };
        double vMomentumResidual_{ 0.0 };

        // ============================================================
        // SIMPLE VELOCITY CORRECTION COEFFICIENTS
        // ============================================================

        /*
         * dU has one value for every U face:
         *
         * (nx + 1) x ny
         */
        std::vector<double> dU;

        /*
         * dV has one value for every V face:
         *
         * nx x (ny + 1)
         */
        std::vector<double> dV;

        // ============================================================
        // INTERNAL METHODS
        // ============================================================

        void applyBoundaryConditions();

        void assembleUMomentum();
        void solveUMomentum();
        void assembleVMomentum();
        void solveVMomentum();

        void assemblePressureCorrection();
        void solvePressureCorrection();

        void correctPressure();
        void correctVelocities();

        void updateMomentumResiduals();

        double calculateResidual();
        bool   checkConvergence();

    public:

        // ============================================================
        // CONSTRUCTOR
        // ============================================================

        StaggeredSIMPLE(
            Mesh& mesh,
            StaggeredFields& fields,
            BoundaryCondition& northBC,
            BoundaryCondition& southBC,
            BoundaryCondition& eastBC,
            BoundaryCondition& westBC);

        // ============================================================
        // SETTINGS
        // ============================================================

        void setPressureRelaxation(double value);
        void setVelocityRelaxation(double value);

        void setConvergenceTolerance(double value);
        void setMomentumConvergenceTolerance(double value);

        void setMaxIterations(std::size_t value);

        void setDensity(double value);
        void setViscosity(double value);

        // ============================================================
        // GETTERS
        // ============================================================

        double getPressureRelaxation() const;
        double getVelocityRelaxation() const;

        double getConvergenceTolerance() const;
        double getMomentumConvergenceTolerance() const;

        double getDensity() const;
        double getViscosity() const;

        std::size_t getMaxIterations() const;
        std::size_t getIteration() const;

        // ============================================================
        // RESIDUALS
        // ============================================================

        // getResidual() is a backwards-compatible alias for
        // getContinuityResidual().
        double getResidual() const;
        double getContinuityResidual() const;
        double getUMomentumResidual() const;
        double getVMomentumResidual() const;

        // ============================================================
        // ITERATIVE SOLVER INTERFACE
        // ============================================================

        void step();

        bool finished() const;

        bool converged() const;

        // ============================================================
        // FULL SOLVE
        // ============================================================

        void solve();
    };
}
