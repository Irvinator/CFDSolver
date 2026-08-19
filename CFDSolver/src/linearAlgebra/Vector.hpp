#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <cassert>

namespace CFD {

    class Vector {
    public:
        // ── Constructors ─────────────────────────────────────────────
        Vector() = default;

        explicit Vector(std::size_t size, double val = 0.0) //size_t is an unsigned int to hold largest object
            : data_(size, val) {
        }

        Vector(std::initializer_list<double> list)
            : data_(list) {
        }

        // ── Size ─────────────────────────────────────────────────────
        std::size_t size() const { return data_.size(); }
        bool empty()       const { return data_.empty(); }

        // ── Element access ───────────────────────────────────────────
        double& operator[](std::size_t i) {
            assert(i < data_.size() && "Index out of bounds");
            return data_[i];
        }
        const double& operator[](std::size_t i) const {
            assert(i < data_.size() && "Index out of bounds");
            return data_[i];
        }

        // ── Arithmetic ───────────────────────────────────────────────
        Vector operator+(const Vector& rhs) const {
            checkSize(rhs);
            Vector result(data_.size());
            for (std::size_t i = 0; i < data_.size(); ++i)
                result[i] = data_[i] + rhs[i];
            return result;
        }

        Vector operator-(const Vector& rhs) const {
            checkSize(rhs);
            Vector result(data_.size());
            for (std::size_t i = 0; i < data_.size(); ++i)
                result[i] = data_[i] - rhs[i];
            return result;
        }

        Vector operator*(double scalar) const {
            Vector result(data_.size());
            for (std::size_t i = 0; i < data_.size(); ++i)
                result[i] = data_[i] * scalar;
            return result;
        }

        friend Vector operator*(double scalar, const Vector& v) {
            return v * scalar;
        }

        Vector& operator+=(const Vector& rhs) {
            checkSize(rhs);
            for (std::size_t i = 0; i < data_.size(); ++i)
                data_[i] += rhs[i];
            return *this;
        }

        Vector& operator-=(const Vector& rhs) {
            checkSize(rhs);
            for (std::size_t i = 0; i < data_.size(); ++i)
                data_[i] -= rhs[i];
            return *this;
        }

        // ── Dot product ──────────────────────────────────────────────
        double dot(const Vector& rhs) const {
            checkSize(rhs);
            double sum = 0.0;
            for (std::size_t i = 0; i < data_.size(); ++i)
                sum += data_[i] * rhs[i];
            return sum;
        }

        // ── Norms ────────────────────────────────────────────────────
        double norm2()   const { return std::sqrt(dot(*this)); }
        double normInf() const {
            double maxVal = 0.0;
            for (const auto& v : data_)
                maxVal = std::max(maxVal, std::abs(v));
            return maxVal;
        }

        // ── Utilities ────────────────────────────────────────────────
        void fill(double val) {
            std::fill(data_.begin(), data_.end(), val);
        }

        void resize(std::size_t n, double val = 0.0) {
            data_.assign(n, val);
        }

        void print(const std::string& name = "Vector") const {
            std::cout << name << " [" << data_.size() << "]: ";
            for (const auto& v : data_)
                std::cout << v << " ";
            std::cout << "\n";
        }

        // ── Raw data access ──────────────────────────────────────────
        const std::vector<double>& data() const { return data_; }
        std::vector<double>& data() { return data_; }

    private:
        std::vector<double> data_;

        void checkSize(const Vector& rhs) const {
            if (data_.size() != rhs.size()) {
                throw std::runtime_error("Vector size mismatch");
            }
        }
    };

} // namespace CFD