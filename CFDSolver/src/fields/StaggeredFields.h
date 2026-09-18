#pragma once

#include <vector>
#include <cstddef>

namespace CFD
{
    class StaggeredFields
    {
    private:
        int nx_{ 0 };
        int ny_{ 0 };

        std::vector<double> p_;
        std::vector<double> u_;
        std::vector<double> v_;

    public:
        StaggeredFields() = default;

        StaggeredFields(int nx, int ny)
            : nx_(nx),
            ny_(ny),
            p_(static_cast<std::size_t>(nx* ny), 0.0),
            u_(static_cast<std::size_t>((nx + 1)* ny), 0.0),
            v_(static_cast<std::size_t>(nx* (ny + 1)), 0.0)
        {
        }

        void resize(int nx, int ny)
        {
            nx_ = nx;
            ny_ = ny;
            p_.assign(static_cast<std::size_t>(nx * ny), 0.0);
            u_.assign(static_cast<std::size_t>((nx + 1) * ny), 0.0);
            v_.assign(static_cast<std::size_t>(nx * (ny + 1)), 0.0);
        }

        int nx() const { return nx_; }
        int ny() const { return ny_; }

        int pIndex(int i, int j) const { return j * nx_ + i; }
        int uIndex(int i, int j) const { return j * (nx_ + 1) + i; }
        int vIndex(int i, int j) const { return j * nx_ + i; }

        double& p(int i, int j) { return p_[static_cast<std::size_t>(pIndex(i, j))]; }
        double& u(int i, int j) { return u_[static_cast<std::size_t>(uIndex(i, j))]; }
        double& v(int i, int j) { return v_[static_cast<std::size_t>(vIndex(i, j))]; }

        double p(int i, int j) const { return p_[static_cast<std::size_t>(pIndex(i, j))]; }
        double u(int i, int j) const { return u_[static_cast<std::size_t>(uIndex(i, j))]; }
        double v(int i, int j) const { return v_[static_cast<std::size_t>(vIndex(i, j))]; }

        std::vector<double>& pressureData() { return p_; }
        std::vector<double>& uData() { return u_; }
        std::vector<double>& vData() { return v_; }

        const std::vector<double>& pressureData() const { return p_; }
        const std::vector<double>& uData() const { return u_; }
        const std::vector<double>& vData() const { return v_; }

        void initialise(double pValue = 0.0, double uValue = 0.0, double vValue = 0.0)
        {
            std::fill(p_.begin(), p_.end(), pValue);
            std::fill(u_.begin(), u_.end(), uValue);
            std::fill(v_.begin(), v_.end(), vValue);
        }
    };
}
