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
        Mesh& mesh;
        Fields& fields;

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

        double relaxationPressure{ 0.2 };
        double relaxationVelocity{ 0.3 };
        double rho{ 1.0 };
        double mu{ 0.01 };
        double convergenceTolerance{ 1.0e-6 };
        std::size_t maxIterations{ 1000 };

        std::size_t iteration{ 0 };
        double residual{ 0.0 };

        std::vector<double> fluxEast;
        std::vector<double> fluxWest;
        std::vector<double> fluxNorth;
        std::vector<double> fluxSouth;

        std::vector<double> invDiagU;
        std::vector<double> invDiagV;

        void applyBoundaryConditions();
        void initializeFluxStorage();

        void assembleUMomentum();
        void solveUMomentum();

        void assembleVMomentum();
        void solveVMomentum();

        void buildRhieChowFaceFluxes();

        void assemblePressureCorrection();
        void solvePressureCorrection();

        void correctPressure();
        void correctFaceFluxes();
        void correctVelocity();


        // Convergence
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


        // Main solver
        void solve();
    };
}
