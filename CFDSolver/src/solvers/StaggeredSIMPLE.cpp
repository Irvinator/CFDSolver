#include "solvers/StaggeredSIMPLE.h"

#include "linearAlgebra/BiCGSTAB.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace CFD
{
    namespace
    {
        constexpr double SMALL = 1.0e-14;

        double safeDiagonal(double value)
        {
            if (!std::isfinite(value) || std::abs(value) < SMALL)
            {
                return SMALL;
            }
            return value;
        }
    }

    static double getDiagonal(const SparseMatrix& matrix, std::size_t P)
    {
        if (P >= matrix.rows())
        {
            return 0.0;
        }

        const auto& rp = matrix.rowPtr();
        const auto& ci = matrix.colIndices();
        const auto& vv = matrix.values();

        if (rp.size() < matrix.rows() + 1)
        {
            return 0.0;
        }

        if (rp[P] > rp[P + 1] || rp[P + 1] > ci.size() || rp[P + 1] > vv.size())
        {
            return 0.0;
        }

        for (std::size_t k = rp[P]; k < rp[P + 1]; ++k)
        {
            if (ci[k] == P)
            {
                return vv[k];
            }
        }

        return 0.0;
    }

    static double relativeAlgebraicResidual(
        const SparseMatrix& matrix,
        const Vector& solution,
        const Vector& rhs)
    {
        const auto& rowPtr = matrix.rowPtr();
        const auto& colIdx = matrix.colIndices();
        const auto& values = matrix.values();

        if (rowPtr.size() != matrix.rows() + 1 ||
            solution.size() != matrix.cols() ||
            rhs.size() != matrix.rows())
        {
            return std::numeric_limits<double>::infinity();
        }

        double l1Residual = 0.0;
        double l1Rhs = 0.0;
        for (std::size_t row = 0; row < matrix.rows(); ++row)
        {
            double ax = 0.0;
            for (std::size_t k = rowPtr[row]; k < rowPtr[row + 1]; ++k)
            {
                ax += values[k] * solution[colIdx[k]];
            }
            l1Residual += std::abs(ax - rhs[row]);
            l1Rhs += std::abs(rhs[row]);
        }

        // A unit floor keeps an identically zero RHS from producing an
        // undefined relative residual.
        return l1Residual / std::max(l1Rhs, 1.0);
    }

    StaggeredSIMPLE::StaggeredSIMPLE(
        Mesh& mesh,
        StaggeredFields& fields,
        BoundaryCondition& northBC,
        BoundaryCondition& southBC,
        BoundaryCondition& eastBC,
        BoundaryCondition& westBC)
        :
        mesh(mesh),
        fields(fields),
        northBC(northBC),
        southBC(southBC),
        eastBC(eastBC),
        westBC(westBC),
        residual(std::numeric_limits<double>::infinity())
    {
        dU.assign(static_cast<std::size_t>((mesh.getNx() + 1) * mesh.getNy()), 0.0);
        dV.assign(static_cast<std::size_t>(mesh.getNx() * (mesh.getNy() + 1)), 0.0);
    }

    void StaggeredSIMPLE::setPressureRelaxation(double value)
    {
        if (!std::isfinite(value) || value <= 0.0 || value > 1.0)
        {
            throw std::runtime_error("StaggeredSIMPLE: pressure relaxation must be between 0 and 1");
        }
        relaxationPressure = value;
    }

    void StaggeredSIMPLE::setVelocityRelaxation(double value)
    {
        if (!std::isfinite(value) || value <= 0.0 || value > 1.0)
        {
            throw std::runtime_error("StaggeredSIMPLE: velocity relaxation must be between 0 and 1");
        }
        relaxationVelocity = value;
    }

    void StaggeredSIMPLE::setConvergenceTolerance(double value)
    {
        if (!std::isfinite(value) || value <= 0.0)
        {
            throw std::runtime_error("StaggeredSIMPLE: continuity tolerance must be positive");
        }
        convergenceTolerance = value;
    }

    void StaggeredSIMPLE::setMomentumConvergenceTolerance(double value)
    {
        if (!std::isfinite(value) || value <= 0.0)
        {
            throw std::runtime_error("StaggeredSIMPLE: momentum tolerance must be positive");
        }
        momentumTolerance_ = value;
    }

    void StaggeredSIMPLE::setMaxIterations(std::size_t value)
    {
        if (value == 0)
        {
            throw std::runtime_error("StaggeredSIMPLE: maximum iterations must be greater than zero");
        }
        maxIterations = value;
    }

    void StaggeredSIMPLE::setDensity(double value)
    {
        if (!std::isfinite(value) || value <= 0.0)
        {
            throw std::runtime_error("StaggeredSIMPLE: density must be positive");
        }
        rho = value;
    }

    void StaggeredSIMPLE::setViscosity(double value)
    {
        if (!std::isfinite(value) || value <= 0.0)
        {
            throw std::runtime_error("StaggeredSIMPLE: viscosity must be positive");
        }
        mu = value;
    }

    double StaggeredSIMPLE::getPressureRelaxation() const { return relaxationPressure; }
    double StaggeredSIMPLE::getVelocityRelaxation() const { return relaxationVelocity; }
    double StaggeredSIMPLE::getConvergenceTolerance() const { return convergenceTolerance; }
    double StaggeredSIMPLE::getMomentumConvergenceTolerance() const { return momentumTolerance_; }
    double StaggeredSIMPLE::getDensity() const { return rho; }
    double StaggeredSIMPLE::getViscosity() const { return mu; }
    std::size_t StaggeredSIMPLE::getMaxIterations() const { return maxIterations; }
    std::size_t StaggeredSIMPLE::getIteration() const { return iteration; }
    double StaggeredSIMPLE::getResidual() const { return residual; }
    double StaggeredSIMPLE::getContinuityResidual() const { return continuityResidual_; }
    double StaggeredSIMPLE::getUMomentumResidual() const { return uMomentumResidual_; }
    double StaggeredSIMPLE::getVMomentumResidual() const { return vMomentumResidual_; }

    void StaggeredSIMPLE::applyBoundaryConditions()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();

        // Normal velocity components live on their respective boundary faces.
        for (int j = 0; j < ny; ++j)
        {
            if (westBC.hasU()) fields.u(0, j) = westBC.getU();
        }

        for (int i = 0; i < nx; ++i)
        {
            if (southBC.hasV()) fields.v(i, 0) = southBC.getV();
            if (northBC.hasV()) fields.v(i, ny) = northBC.getV();
        }

        // Do not directly overwrite U at j=0 or j=ny-1, or V at i=0 or
        // i=nx-1. Those components are tangential to the boundary and their
        // no-slip Dirichlet values are imposed through half-cell diffusion
        // coefficients and source terms in the momentum equations.
        //
        // Pressure is cell-centred. A prescribed outlet pressure is a face
        // condition for U momentum, not p(nx - 1, j).
    }

    void StaggeredSIMPLE::assembleUMomentum()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();
        const std::size_t nU = static_cast<std::size_t>((nx + 1) * ny);

        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> values;
        uRHS = Vector(nU, 0.0);

        const double dx = mesh.getDx();
        const double dy = mesh.getDy();
        const double Ae = mesh.eastWestFaceArea();
        const double An = mesh.northSouthFaceArea();
        const double Dw = mu * An / dy;

        rows.reserve(5 * nU);
        cols.reserve(5 * nU);
        values.reserve(5 * nU);

        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i <= nx; ++i)
            {
                const std::size_t P = static_cast<std::size_t>(fields.uIndex(i, j));

                // U is normal to west/east boundaries and is therefore stored
                // directly on those faces. At south/north it is tangential.
                const bool fixed = (i == 0 && westBC.hasU());
                if (fixed)
                {
                    rows.push_back(P); cols.push_back(P); values.push_back(1.0);
                    uRHS[P] = westBC.getU();
                    dU[P] = 0.0;
                    continue;
                }

                const double aE = (i < nx) ? mu * Ae / dx : 0.0;
                const double aW = (i > 0) ? mu * Ae / dx : 0.0;
                const double aN = (j < ny - 1) ? Dw : (northBC.hasU() ? 2.0 * Dw : 0.0);
                const double aS = (j > 0) ? Dw : (southBC.hasU() ? 2.0 * Dw : 0.0);

                const double aPNoRelax = aE + aW + aN + aS + SMALL;
                const double aP = safeDiagonal(aPNoRelax / relaxationVelocity);
                dU[P] = Ae / aP;

                const double pW = (i > 0) ? fields.p(i - 1, j) : fields.p(0, j);
                const double pE = (i < nx) ? fields.p(i, j) :
                    (eastBC.hasPressure() ? eastBC.getPressure() : fields.p(nx - 1, j));

                double source = (pW - pE) * Ae +
                    ((1.0 - relaxationVelocity) / relaxationVelocity) *
                    aPNoRelax * fields.u(i, j);

                // For a tangential Dirichlet wall, the U location is dy/2
                // from the wall. The ghost-cell relation gives 2*Dw.
                if (j == 0 && southBC.hasU()) source += 2.0 * Dw * southBC.getU();
                if (j == ny - 1 && northBC.hasU()) source += 2.0 * Dw * northBC.getU();

                rows.push_back(P); cols.push_back(P); values.push_back(aP);
                if (i < nx) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(fields.uIndex(i + 1, j))); values.push_back(-aE); }
                if (i > 0) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(fields.uIndex(i - 1, j))); values.push_back(-aW); }
                if (j < ny - 1) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(fields.uIndex(i, j + 1))); values.push_back(-aN); }
                if (j > 0) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(fields.uIndex(i, j - 1))); values.push_back(-aS); }
                uRHS[P] = source;
            }
        }

        uMatrix = SparseMatrix(nU, nU);
        uMatrix.setFromTriplets(rows, cols, values);
    }

    void StaggeredSIMPLE::solveUMomentum()
    {
        assembleUMomentum();
        const std::size_t nU = static_cast<std::size_t>((mesh.getNx() + 1) * mesh.getNy());
        Vector solution(nU, 0.0);
        for (std::size_t k = 0; k < nU; ++k) solution[k] = fields.uData()[k];

        BiCGSTAB solver(1.0e-8, 10000);
        const BiCGSTABResult result = solver.solve(uMatrix, uRHS, solution, false);
        if (!result.converged)
        {
            throw std::runtime_error("StaggeredSIMPLE: U momentum solver failed to converge");
        }

        for (std::size_t k = 0; k < nU; ++k) fields.uData()[k] = solution[k];
    }

    void StaggeredSIMPLE::assembleVMomentum()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();
        const std::size_t nV = static_cast<std::size_t>(nx * (ny + 1));

        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> values;
        vRHS = Vector(nV, 0.0);

        const double dx = mesh.getDx();
        const double dy = mesh.getDy();
        const double Ae = mesh.eastWestFaceArea();
        const double An = mesh.northSouthFaceArea();
        const double Dv = mu * Ae / dx;

        rows.reserve(5 * nV);
        cols.reserve(5 * nV);
        values.reserve(5 * nV);

        for (int j = 0; j <= ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t P = static_cast<std::size_t>(fields.vIndex(i, j));

                // V is normal to south/north faces. At west/east it is
                // tangential and is handled by the half-cell wall stencil.
                const bool fixed =
                    (j == 0 && southBC.hasV()) ||
                    (j == ny && northBC.hasV());

                if (fixed)
                {
                    rows.push_back(P); cols.push_back(P); values.push_back(1.0);
                    vRHS[P] = (j == 0) ? southBC.getV() : northBC.getV();
                    dV[P] = 0.0;
                    continue;
                }

                const double aE = (i < nx - 1) ? Dv : (eastBC.hasV() ? 2.0 * Dv : 0.0);
                const double aW = (i > 0) ? Dv : (westBC.hasV() ? 2.0 * Dv : 0.0);
                const double aN = (j < ny) ? mu * An / dy : 0.0;
                const double aS = (j > 0) ? mu * An / dy : 0.0;

                const double aPNoRelax = aE + aW + aN + aS + SMALL;
                const double aP = safeDiagonal(aPNoRelax / relaxationVelocity);
                dV[P] = An / aP;

                const double pS = (j > 0) ? fields.p(i, j - 1) : fields.p(i, 0);
                const double pN = (j < ny) ? fields.p(i, j) : fields.p(i, ny - 1);
                double source = (pS - pN) * An +
                    ((1.0 - relaxationVelocity) / relaxationVelocity) *
                    aPNoRelax * fields.v(i, j);

                if (i == 0 && westBC.hasV()) source += 2.0 * Dv * westBC.getV();
                if (i == nx - 1 && eastBC.hasV()) source += 2.0 * Dv * eastBC.getV();

                rows.push_back(P); cols.push_back(P); values.push_back(aP);
                if (i < nx - 1) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(fields.vIndex(i + 1, j))); values.push_back(-aE); }
                if (i > 0) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(fields.vIndex(i - 1, j))); values.push_back(-aW); }
                if (j < ny) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(fields.vIndex(i, j + 1))); values.push_back(-aN); }
                if (j > 0) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(fields.vIndex(i, j - 1))); values.push_back(-aS); }
                vRHS[P] = source;
            }
        }

        vMatrix = SparseMatrix(nV, nV);
        vMatrix.setFromTriplets(rows, cols, values);
    }

    void StaggeredSIMPLE::solveVMomentum()
    {
        assembleVMomentum();
        const std::size_t nV = static_cast<std::size_t>(mesh.getNx() * (mesh.getNy() + 1));
        Vector solution(nV, 0.0);
        for (std::size_t k = 0; k < nV; ++k) solution[k] = fields.vData()[k];

        BiCGSTAB solver(1.0e-8, 10000);
        const BiCGSTABResult result = solver.solve(vMatrix, vRHS, solution, false);
        if (!result.converged)
        {
            throw std::runtime_error("StaggeredSIMPLE: V momentum solver failed to converge");
        }

        for (std::size_t k = 0; k < nV; ++k) fields.vData()[k] = solution[k];
    }

    void StaggeredSIMPLE::assemblePressureCorrection()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();
        const std::size_t nP = static_cast<std::size_t>(nx * ny);

        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> values;
        pressureRHS = Vector(nP, 0.0);

        const double Ae = mesh.eastWestFaceArea();
        const double An = mesh.northSouthFaceArea();
        rows.reserve(5 * nP);
        cols.reserve(5 * nP);
        values.reserve(5 * nP);

        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t P = static_cast<std::size_t>(fields.pIndex(i, j));

                // A fixed east outlet pressure imposes p'_E = 0 and anchors
                // the system. Without it, retain one arbitrary reference.
                const bool isReference = !eastBC.hasPressure() && i == nx - 1 && j == ny - 1;
                if (isReference)
                {
                    rows.push_back(P); cols.push_back(P); values.push_back(1.0);
                    pressureRHS[P] = 0.0;
                    continue;
                }

                double aE = 0.0;
                double aW = 0.0;
                double aN = 0.0;
                double aS = 0.0;

                // dU and dV already contain the velocity-face area. There is
                // deliberately no additional division by dx or dy here.
                if (i < nx - 1)
                    aE = rho * Ae * dU[static_cast<std::size_t>(fields.uIndex(i + 1, j))];
                else if (eastBC.hasPressure())
                    aE = rho * Ae * dU[static_cast<std::size_t>(fields.uIndex(nx, j))];

                if (i > 0)
                    aW = rho * Ae * dU[static_cast<std::size_t>(fields.uIndex(i, j))];
                if (j < ny - 1)
                    aN = rho * An * dV[static_cast<std::size_t>(fields.vIndex(i, j + 1))];
                if (j > 0)
                    aS = rho * An * dV[static_cast<std::size_t>(fields.vIndex(i, j))];

                const double aP = safeDiagonal(aE + aW + aN + aS);
                const double continuity =
                    rho * Ae * fields.u(i + 1, j) - rho * Ae * fields.u(i, j) +
                    rho * An * fields.v(i, j + 1) - rho * An * fields.v(i, j);

                rows.push_back(P); cols.push_back(P); values.push_back(aP);
                if (i < nx - 1 && aE > SMALL) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(fields.pIndex(i + 1, j))); values.push_back(-aE); }
                if (i > 0 && aW > SMALL) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(fields.pIndex(i - 1, j))); values.push_back(-aW); }
                if (j < ny - 1 && aN > SMALL) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(fields.pIndex(i, j + 1))); values.push_back(-aN); }
                if (j > 0 && aS > SMALL) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(fields.pIndex(i, j - 1))); values.push_back(-aS); }
                pressureRHS[P] = -continuity;
            }
        }

        pressureMatrix = SparseMatrix(nP, nP);
        pressureMatrix.setFromTriplets(rows, cols, values);
    }

    void StaggeredSIMPLE::solvePressureCorrection()
    {
        assemblePressureCorrection();
        const std::size_t nP = static_cast<std::size_t>(mesh.getNx() * mesh.getNy());
        pressureCorrection = Vector(nP, 0.0);

        BiCGSTAB solver(1.0e-8, 10000);
        const BiCGSTABResult result = solver.solve(pressureMatrix, pressureRHS, pressureCorrection, false);
        if (!result.converged)
        {
            throw std::runtime_error("StaggeredSIMPLE: pressure correction solver failed to converge");
        }
    }

    void StaggeredSIMPLE::correctPressure()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();

        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                fields.p(i, j) += relaxationPressure *
                    pressureCorrection[static_cast<std::size_t>(fields.pIndex(i, j))];
            }
        }

        // Do not overwrite p(nx - 1, j): the specified outlet pressure is a
        // face value, already used in U momentum and p'_E = 0 treatment.
    }

    void StaggeredSIMPLE::correctVelocities()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();

        for (int j = 0; j < ny; ++j)
        {
            for (int i = 1; i < nx; ++i)
            {
                const std::size_t idx = static_cast<std::size_t>(fields.uIndex(i, j));
                const double pW = pressureCorrection[static_cast<std::size_t>(fields.pIndex(i - 1, j))];
                const double pE = pressureCorrection[static_cast<std::size_t>(fields.pIndex(i, j))];
                fields.u(i, j) += dU[idx] * (pW - pE);
            }
        }

        if (eastBC.hasPressure())
        {
            for (int j = 0; j < ny; ++j)
            {
                const std::size_t idx = static_cast<std::size_t>(fields.uIndex(nx, j));
                const double pW = pressureCorrection[static_cast<std::size_t>(fields.pIndex(nx - 1, j))];
                fields.u(nx, j) += dU[idx] * pW; // p'_E = 0 at outlet face.
            }
        }

        for (int j = 1; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t idx = static_cast<std::size_t>(fields.vIndex(i, j));
                const double pS = pressureCorrection[static_cast<std::size_t>(fields.pIndex(i, j - 1))];
                const double pN = pressureCorrection[static_cast<std::size_t>(fields.pIndex(i, j))];
                fields.v(i, j) += dV[idx] * (pS - pN);
            }
        }

        applyBoundaryConditions();
    }

    void StaggeredSIMPLE::updateMomentumResiduals()
    {
        // Reassemble with the pressure and velocity fields after correction.
        // This tests the actual current SIMPLE iterate, not merely whether
        // the predictor linear systems were solved accurately.
        assembleUMomentum();
        Vector u(static_cast<std::size_t>((mesh.getNx() + 1) * mesh.getNy()), 0.0);
        for (std::size_t k = 0; k < u.size(); ++k) u[k] = fields.uData()[k];
        uMomentumResidual_ = relativeAlgebraicResidual(uMatrix, u, uRHS);

        assembleVMomentum();
        Vector v(static_cast<std::size_t>(mesh.getNx() * (mesh.getNy() + 1)), 0.0);
        for (std::size_t k = 0; k < v.size(); ++k) v[k] = fields.vData()[k];
        vMomentumResidual_ = relativeAlgebraicResidual(vMatrix, v, vRHS);
    }

    double StaggeredSIMPLE::calculateResidual()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();
        const double Ae = mesh.eastWestFaceArea();
        const double An = mesh.northSouthFaceArea();

        double totalResidual = 0.0;
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                const double continuity =
                    rho * Ae * fields.u(i + 1, j) - rho * Ae * fields.u(i, j) +
                    rho * An * fields.v(i, j + 1) - rho * An * fields.v(i, j);
                totalResidual += std::abs(continuity);
            }
        }
        return totalResidual;
    }

    bool StaggeredSIMPLE::checkConvergence()
    {
        continuityResidual_ = calculateResidual();
        residual = continuityResidual_; // Backward-compatible legacy getter.

        return continuityResidual_ < convergenceTolerance &&
            uMomentumResidual_ < momentumTolerance_ &&
            vMomentumResidual_ < momentumTolerance_;
    }

    void StaggeredSIMPLE::solve()
    {
        iteration = 0;
        residual = std::numeric_limits<double>::infinity();
        continuityResidual_ = std::numeric_limits<double>::infinity();
        uMomentumResidual_ = std::numeric_limits<double>::infinity();
        vMomentumResidual_ = std::numeric_limits<double>::infinity();

        for (iteration = 1; iteration <= maxIterations; ++iteration)
        {
            applyBoundaryConditions();
            solveUMomentum();
            solveVMomentum();
            solvePressureCorrection();
            correctPressure();
            correctVelocities();

            updateMomentumResiduals();
            const bool converged = checkConvergence();

            std::cout << "Staggered SIMPLE iteration " << iteration
                << " | continuity = " << continuityResidual_
                << " | U momentum = " << uMomentumResidual_
                << " | V momentum = " << vMomentumResidual_ << "\n";

            if (converged)
            {
                std::cout << "\nStaggered SIMPLE converged after " << iteration
                    << " iterations.\n";
                return;
            }
        }

        std::cout << "\nStaggered SIMPLE reached maximum iterations.\n"
            << "Final continuity residual = " << continuityResidual_ << "\n"
            << "Final U-momentum residual = " << uMomentumResidual_ << "\n"
            << "Final V-momentum residual = " << vMomentumResidual_ << "\n";
    }
}
