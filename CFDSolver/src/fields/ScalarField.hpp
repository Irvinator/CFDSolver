#pragma once
#include "../linearAlgebra/Vector.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace CFD {

    class ScalarField {
    public:
        ScalarField() = default;

        explicit ScalarField(std::size_t nCells, double initVal = 0.0)
            : data_(nCells, initVal) {
        }

        ScalarField(std::size_t nCells, double initVal, const std::string& name)
            : data_(nCells, initVal), name_(name) {
        }

        double& operator[](std::size_t i) {
            assert(i < data_.size());
            return data_[i];
        }

        const double& operator[](std::size_t i) const {
            assert(i < data_.size());
            return data_[i];
        }

        std::size_t size() const { return data_.size(); }
        bool empty() const { return data_.empty(); }

        double min() const {
            if (data_.empty())
                throw std::runtime_error("ScalarField::min on empty field");
            return *std::min_element(data_.begin(), data_.end());
        }

        double max() const {
            if (data_.empty())
                throw std::runtime_error("ScalarField::max on empty field");
            return *std::max_element(data_.begin(), data_.end());
        }

        double mean() const {
            if (data_.empty()) return 0.0;
            return std::accumulate(data_.begin(), data_.end(), 0.0) /
                static_cast<double>(data_.size());
        }

        double norm2() const {
            double sum = 0.0;
            for (double v : data_) sum += v * v;
            return std::sqrt(sum);
        }

        ScalarField operator+(const ScalarField& rhs) const {
            checkSize(rhs);
            ScalarField result(size());
            for (std::size_t i = 0; i < size(); ++i)
                result[i] = data_[i] + rhs[i];
            return result;
        }

        ScalarField operator-(const ScalarField& rhs) const {
            checkSize(rhs);
            ScalarField result(size());
            for (std::size_t i = 0; i < size(); ++i)
                result[i] = data_[i] - rhs[i];
            return result;
        }

        ScalarField operator*(double scalar) const {
            ScalarField result(size());
            for (std::size_t i = 0; i < size(); ++i)
                result[i] = data_[i] * scalar;
            return result;
        }

        ScalarField& operator+=(const ScalarField& rhs) {
            checkSize(rhs);
            for (std::size_t i = 0; i < size(); ++i)
                data_[i] += rhs[i];
            return *this;
        }

        ScalarField& operator-=(const ScalarField& rhs) {
            checkSize(rhs);
            for (std::size_t i = 0; i < size(); ++i)
                data_[i] -= rhs[i];
            return *this;
        }

        void fill(double val) {
            std::fill(data_.begin(), data_.end(), val);
        }

        void setName(const std::string& name) { name_ = name; }
        const std::string& name() const { return name_; }

        Vector toVector() const {
            Vector v(size());
            for (std::size_t i = 0; i < size(); ++i)
                v[i] = data_[i];
            return v;
        }

        void fromVector(const Vector& v) {
            if (v.size() != size())
                throw std::runtime_error("ScalarField::fromVector size mismatch");
            for (std::size_t i = 0; i < size(); ++i)
                data_[i] = v[i];
        }

        const std::vector<double>& data() const { return data_; }
        std::vector<double>& data() { return data_; }

        void print() const {
            std::cout << "ScalarField '" << name_ << "' ["
                << size() << " cells]";
            if (!empty()) {
                std::cout << " min=" << min()
                    << " max=" << max()
                    << " mean=" << mean();
            }
            std::cout << "\n";
        }

    private:
        std::vector<double> data_;
        std::string name_ = "unnamed";

        void checkSize(const ScalarField& rhs) const {
            if (size() != rhs.size())
                throw std::runtime_error("ScalarField size mismatch");
        }
    };

} // namespace CFD