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

        SparseMatrix uMatrix;
        Vector uRHS;

        SparseMatrix vMatrix;
        Vector vRHS;

        SparseMatrix pressureMatrix;
        Vector pressureRHS;
        Vector pressureCorrection;

        double relaxationPressure{ 0.3 };
        double relaxationVelocity{ 0.7 };
        double rho{ 1.0 };
        double mu{ 0.01 };
        double convergenceTolerance{ 1.0e-6 };
        std::size_t maxIterations{ 1000 };

        std::size_t iteration{ 0 };
        double residual{ 0.0 };

        std::vector<double> dU;
        std::vector<double> dV;

        void applyBoundaryConditions();

        void assembleUMomentum();
        void solveUMomentum();

        void assembleVMomentum();
        void solveVMomentum();

        void assemblePressureCorrection();
        void solvePressureCorrection();

        void correctPressure();
        void correctVelocities();

        double calculateResidual();
        bool checkConvergence();

    public:
        StaggeredSIMPLE(
            Mesh& mesh,
            StaggeredFields& fields,
            BoundaryCondition& northBC,
            BoundaryCondition& southBC,
            BoundaryCondition& eastBC,
            BoundaryCondition& westBC);

        void setPressureRelaxation(double value);
        void setVelocityRelaxation(double value);
        void setConvergenceTolerance(double value);
        void setMaxIterations(std::size_t value);
        void setDensity(double value);
        void setViscosity(double value);

        double getPressureRelaxation() const;
        double getVelocityRelaxation() const;
        double getConvergenceTolerance() const;
        double getDensity() const;
        double getViscosity() const;
        std::size_t getMaxIterations() const;
        std::size_t getIteration() const;
        double getResidual() const;

        void solve();
    };
}
