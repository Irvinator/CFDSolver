
#pragma once

#include "mesh/mesh2D.h"
#include "fields/Fields.h"
#include "BC/BC.h"

#include <cstddef>

namespace CFD
{
    class SIMPLE
    {
    private:

        // References to the CFD model
        Mesh& mesh;
        Fields& fields;
        BoundaryCondition& boundaryConditions;


        // -------------------------
        // SIMPLE parameters
        // -------------------------

        double relaxationPressure;
        double relaxationVelocity;

        double convergenceTolerance;

        std::size_t maxIterations;


        // -------------------------
        // Solver functions
        // -------------------------

        void solveUMomentum();

        void solveVMomentum();

        void solvePressureCorrection();

        void correctPressure();

        void correctVelocity();

        double calculateResidual();


    public:

        // Constructor
        SIMPLE(
            Mesh& mesh,
            Fields& fields,
            BoundaryConditions& boundaryConditions
        );


        // -------------------------
        // Solver settings
        // -------------------------

        void setPressureRelaxation(double value);

        void setVelocityRelaxation(double value);

        void setConvergenceTolerance(double value);

        void setMaxIterations(std::size_t value);


        // -------------------------
        // Get solver settings
        // -------------------------

        double getPressureRelaxation() const;

        double getVelocityRelaxation() const;

        double getConvergenceTolerance() const;

        std::size_t getMaxIterations() const;


        // -------------------------
        // Main solver
        // -------------------------

        void solve();
    };
}


