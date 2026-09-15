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

        // ------------------------------------------------------------
        // Linear systems
        // ------------------------------------------------------------

        SparseMatrix uMatrix;
        Vector uRHS;

        SparseMatrix vMatrix;
        Vector vRHS;

        SparseMatrix pressureMatrix;
        Vector pressureRHS;
        Vector pressureCorrection;

        // ------------------------------------------------------------
        // Solver settings
        // ------------------------------------------------------------

        double relaxationPressure{ 0.3 };
        double relaxationVelocity{ 0.7 };

        double rho{ 1.0 };
        double mu{ 0.01 };

        double convergenceTolerance{ 1.0e-6 };
        double momentumTolerance_{ 1.0e-6 };

        std::size_t maxIterations{ 1000 };

        // ------------------------------------------------------------
        // Iteration state
        // ------------------------------------------------------------

        std::size_t iteration{ 0 };

        bool converged_{ false };
        bool finished_{ false };

        // Legacy residual.
        // This remains equal to the continuity residual.
        double residual{
            0.0
        };

        double continuityResidual_{
            0.0
        };

        double uMomentumResidual_{
            0.0
        };

        double vMomentumResidual_{
            0.0
        };

        // ------------------------------------------------------------
        // SIMPLE velocity correction coefficients
        // ------------------------------------------------------------

        std::vector<double> dU;
        std::vector<double> dV;

        // ------------------------------------------------------------
        // Internal methods
        // ------------------------------------------------------------

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
        bool checkConvergence();

    public:

        // ------------------------------------------------------------
        // Constructor
        // ------------------------------------------------------------

        StaggeredSIMPLE(
            Mesh& mesh,
            StaggeredFields& fields,
            BoundaryCondition& northBC,
            BoundaryCondition& southBC,
            BoundaryCondition& eastBC,
            BoundaryCondition& westBC);

        // ------------------------------------------------------------
        // Settings
        // ------------------------------------------------------------

        void setPressureRelaxation(double value);
        void setVelocityRelaxation(double value);

        void setConvergenceTolerance(double value);
        void setMomentumConvergenceTolerance(double value);

        void setMaxIterations(std::size_t value);

        void setDensity(double value);
        void setViscosity(double value);

        // ------------------------------------------------------------
        // Settings getters
        // ------------------------------------------------------------

        double getPressureRelaxation() const;
        double getVelocityRelaxation() const;

        double getConvergenceTolerance() const;
        double getMomentumConvergenceTolerance() const;

        double getDensity() const;
        double getViscosity() const;

        std::size_t getMaxIterations() const;

        // ------------------------------------------------------------
        // Iteration information
        // ------------------------------------------------------------

        std::size_t getIteration() const;

        // ------------------------------------------------------------
        // Residuals
        // ------------------------------------------------------------

        // Legacy getter.
        // Returns continuity residual.
        double getResidual() const;

        double getContinuityResidual() const;
        double getUMomentumResidual() const;
        double getVMomentumResidual() const;

        // ------------------------------------------------------------
        // Iterative solver interface
        // ------------------------------------------------------------

        /*
         * Performs exactly ONE SIMPLE iteration.
         *
         * This is the function used by the GUI animation thread.
         */
        void step();

        /*
         * Returns true when either:
         *
         * 1. The SIMPLE solution has converged, or
         * 2. Maximum iterations have been reached.
         */
        bool finished() const;

        /*
         * Returns whether the current solution has actually converged.
         */
        bool converged() const;

        // ------------------------------------------------------------
        // Full solve
        // ------------------------------------------------------------

        /*
         * Runs SIMPLE until convergence or maximum iterations.
         *
         * This is retained for non-GUI / normal solver usage.
         */
        void solve();
    };
    

}
