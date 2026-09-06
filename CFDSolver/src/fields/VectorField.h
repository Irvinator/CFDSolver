#pragma once

#include "ScalarField.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace CFD {

    class VectorField {
    public:

        // Default constructor
        VectorField() = default;

        // Create a vector field with nCells
        // x and y components are initialised to initVal
        explicit VectorField(std::size_t nCells, double initVal = 0.0)
            : x_(nCells, initVal),
            y_(nCells, initVal) {
        }

        // Create a named vector field
        VectorField(
            std::size_t nCells,
            double initVal,
            const std::string& name
        )
            : x_(nCells, initVal, name + "_x"),
            y_(nCells, initVal, name + "_y") {
        }

        // Access x-component
        ScalarField& getx() {
            return x_;
        }

        const ScalarField& getx() const {
            return x_;
        }

        // Access y-component
        ScalarField& gety() {
            return y_;
        }

        const ScalarField& gety() const {
            return y_;
        }

        // Number of cells
        std::size_t getsize() const {
            return x_.size();
        }

        // Check whether field is empty
        bool empty() const {
            return x_.empty();
        }

        // Set both components to the same value
        void fill(double val) {
            x_.fill(val);
            y_.fill(val);
        }

        // Set x and y components independently
        void fill(double xVal, double yVal) {
            x_.fill(xVal);
            y_.fill(yVal);
        }

        // Magnitude of the vector at a particular cell
        double magnitude(std::size_t i) const {
            if (i >= getsize())
                throw std::out_of_range("VectorField::magnitude index out of range");

            return std::sqrt(
                x_[i] * x_[i] +
                y_[i] * y_[i]
            );
        }

        // Maximum velocity magnitude across the field
        double maxMagnitude() const {
            if (empty())
                return 0.0;

            double maxMag = 0.0;

            for (std::size_t i = 0; i < getsize(); ++i) {
                maxMag = std::max(maxMag, magnitude(i));
            }

            return maxMag;
        }

        // Set a particular vector value
        void set(std::size_t i, double xValue, double yValue) {
            if (i >= getsize())
                throw std::out_of_range("VectorField::set index out of range");

            x_[i] = xValue;
            y_[i] = yValue;
        }

        // Get the name of the x-component
        const std::string& xName() const {
            return x_.name();
        }

        // Get the name of the y-component
        const std::string& yName() const {
            return y_.name();
        }

    private:

        // x-component of the vector field
        ScalarField x_;

        // y-component of the vector field
        ScalarField y_;
    };

} // namespace CFD