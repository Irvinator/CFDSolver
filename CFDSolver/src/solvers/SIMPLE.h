#pragma once

#include "mesh/mesh2D.h"
#include "fields/Fields.h"
#include "BCs/BC.h"
#include "linearAlgebra/SparseMatrix.hpp"
#include "linearAlgebra/Vector.hpp"
#include "linearAlgebra/BiCGSTAB.h"

#include <cstddef>
#include <vector>

namespace CFD
{
    class SIMPLE
    {
    private:
        // Core model references
        Mesh& mesh;
        Fields& fields;

        // Boundary conditions per side
        BoundaryCondition& northBC;
        BoundaryCondition& southBC;
        BoundaryCondition& eastBC;
        BoundaryCondition& westBC;

        // Momentum systems
        SparseMatrix uMatrix;
        Vector uRHS;

        SparseMatrix vMatrix;
        Vector vRHS;

        // Pressure-correction system
        SparseMatrix pressureMatrix;
        Vector pressureRHS;
        Vector pressureCorrection;

        // SIMPLE parameters
        double relaxationPressure{ 0.3 };
        double relaxationVelocity{ 0.7 };
        double rho{ 1.0 };
        double mu{ 0.01 };
        double convergenceTolerance{ 1.0e-6 };
        std::size_t maxIterations{ 1000 };

        // Solver state
        std::size_t iteration{ 0 };
        double residual{ 0.0 };

        // Face mass flux storage
        std::vector<double> fluxEast;
        std::vector<double> fluxWest;
        std::vector<double> fluxNorth;
        std::vector<double> fluxSouth;

        // Boundary condition enforcement
        void applyBoundaryConditions();

        // Face flux construction
        void calculateFaceFluxes();

        // Momentum equations
        void assembleUMomentum();
        void solveUMomentum();

        void assembleVMomentum();
        void solveVMomentum();

        // Pressure correction
        void assemblePressureCorrection();
        void solvePressureCorrection();

        // SIMPLE corrections
        void correctPressure();
        void correctVelocity();

        // Convergence monitoring
        double calculateResidual();
        bool checkConvergence();

    public:
        SIMPLE(
            Mesh& mesh,
            Fields& fields,
            BoundaryCondition& northBC,
            BoundaryCondition& southBC,
            BoundaryCondition& eastBC,
            BoundaryCondition& westBC);

        // Solver settings
        void setPressureRelaxation(double value);
        void setVelocityRelaxation(double value);
        void setConvergenceTolerance(double value);
        void setMaxIterations(std::size_t value);
        void setDensity(double value);

        // Optional transport-property setter used by momentum assembly
        void setViscosity(double value);

        // Solver settings getters
        double getPressureRelaxation() const;
        double getVelocityRelaxation() const;
        double getConvergenceTolerance() const;
        double getDensity() const;
        double getViscosity() const;
        std::size_t getMaxIterations() const;

        // Solver state getters
        std::size_t getIteration() const;
        double getResidual() const;

        // Main solver
        void solve();
    };
}