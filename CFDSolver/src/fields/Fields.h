//combines the required pressure and velocity fields into a 
//single vector for output or further processing. 
//pressure corecctor is required for the simple solver algorithm.

#pragma once

#include "ScalarField.hpp"
#include "VectorField.h"

#include "mesh/mesh2D.h"

#include <cstddef>
#include <stdexcept>

namespace CFD {

    class Fields {
    public:

        // Default constructor
        Fields() = default;

        // Construct fields using the number of cells
        explicit Fields(std::size_t nCells)
            : pressure(nCells, 0.0, "pressure"),
            pressureCorrection(nCells, 0.0, "pressureCorrection"),
            velocity(nCells, 0.0, "velocity") {
        }

        // Construct fields directly from a mesh
        explicit Fields(const Mesh& mesh)
            : Fields(mesh.getNx()* mesh.getNy()) {
        }

        // Number of cells
        std::size_t size() const {
            return pressure.size();
        }

        // Check whether fields contain no cells
        bool empty() const {
            return pressure.empty();
        }

        // Set all fields to zero
        void initialise(double pressureValue = 0.0,
            double uValue = 0.0,
            double vValue = 0.0) {

            pressure.fill(pressureValue);
            pressureCorrection.fill(0.0);

            velocity.getx().fill(uValue);
            velocity.gety().fill(vValue);
        }

        // Physical fields
        ScalarField pressure;
        ScalarField pressureCorrection;

        VectorField velocity;
    };

} // namespace CFD
