#include "solvers/StaggeredSIMPLE.h"

#include "linearAlgebra/BiCGSTAB.h"

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

    static double getDiagonal(
        const SparseMatrix& matrix,
        std::size_t P)
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

        if (rp[P] > rp[P + 1] ||
            rp[P + 1] > ci.size() ||
            rp[P + 1] > vv.size())
        {
            return 0.0;
        }

        for (std::size_t k = rp[P];
            k < rp[P + 1];
            ++k)
        {
            if (ci[k] == P)
            {
                return vv[k];
            }
        }

        return 0.0;
    }


    // ================================================================
    // CONSTRUCTOR
    // ================================================================

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
        residual(std::numeric_limits<double>::infinity()),
        iteration(0)
    {
        dU.assign(
            static_cast<std::size_t>(
                (mesh.getNx() + 1) * mesh.getNy()),
            0.0);

        dV.assign(
            static_cast<std::size_t>(
                mesh.getNx() * (mesh.getNy() + 1)),
            0.0);
    }


    // ================================================================
    // SETTINGS
    // ================================================================

    void StaggeredSIMPLE::setPressureRelaxation(double value)
    {
        if (!std::isfinite(value) ||
            value <= 0.0 ||
            value > 1.0)
        {
            throw std::runtime_error(
                "StaggeredSIMPLE: pressure relaxation must be between 0 and 1");
        }

        relaxationPressure = value;
    }


    void StaggeredSIMPLE::setVelocityRelaxation(double value)
    {
        if (!std::isfinite(value) ||
            value <= 0.0 ||
            value > 1.0)
        {
            throw std::runtime_error(
                "StaggeredSIMPLE: velocity relaxation must be between 0 and 1");
        }

        relaxationVelocity = value;
    }


    void StaggeredSIMPLE::setConvergenceTolerance(double value)
    {
        if (!std::isfinite(value) || value <= 0.0)
        {
            throw std::runtime_error(
                "StaggeredSIMPLE: convergence tolerance must be positive");
        }

        convergenceTolerance = value;
    }


    void StaggeredSIMPLE::setMaxIterations(std::size_t value)
    {
        if (value == 0)
        {
            throw std::runtime_error(
                "StaggeredSIMPLE: maximum iterations must be greater than zero");
        }

        maxIterations = value;
    }


    void StaggeredSIMPLE::setDensity(double value)
    {
        if (!std::isfinite(value) || value <= 0.0)
        {
            throw std::runtime_error(
                "StaggeredSIMPLE: density must be positive");
        }

        rho = value;
    }


    void StaggeredSIMPLE::setViscosity(double value)
    {
        if (!std::isfinite(value) || value <= 0.0)
        {
            throw std::runtime_error(
                "StaggeredSIMPLE: viscosity must be positive");
        }

        mu = value;
    }


    double StaggeredSIMPLE::getPressureRelaxation() const
    {
        return relaxationPressure;
    }


    double StaggeredSIMPLE::getVelocityRelaxation() const
    {
        return relaxationVelocity;
    }


    double StaggeredSIMPLE::getConvergenceTolerance() const
    {
        return convergenceTolerance;
    }


    double StaggeredSIMPLE::getDensity() const
    {
        return rho;
    }


    double StaggeredSIMPLE::getViscosity() const
    {
        return mu;
    }


    std::size_t StaggeredSIMPLE::getMaxIterations() const
    {
        return maxIterations;
    }


    std::size_t StaggeredSIMPLE::getIteration() const
    {
        return iteration;
    }


    double StaggeredSIMPLE::getResidual() const
    {
        return residual;
    }


    // ================================================================
    // BOUNDARY CONDITIONS
    // ================================================================

    void StaggeredSIMPLE::applyBoundaryConditions()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();

        // West inlet U
        for (int j = 0; j < ny; ++j)
        {
            if (westBC.hasU())
            {
                fields.u(0, j) = westBC.getU();
            }
        }

        // South / North V
        for (int i = 0; i < nx; ++i)
        {
            if (southBC.hasV())
            {
                fields.v(i, 0) = southBC.getV();
            }

            if (northBC.hasV())
            {
                fields.v(i, ny) = northBC.getV();
            }
        }

        // South / North U
        for (int i = 0; i <= nx; ++i)
        {
            if (southBC.hasU() && i < nx)
            {
                fields.u(i, 0) = southBC.getU();
            }

            if (northBC.hasU() && i < nx)
            {
                fields.u(i, ny - 1) = northBC.getU();
            }
        }

        // West / East V
        for (int j = 0; j <= ny; ++j)
        {
            if (westBC.hasV() && j < ny)
            {
                fields.v(0, j) = westBC.getV();
            }

            if (eastBC.hasV() && j < ny)
            {
                fields.v(nx - 1, j) = eastBC.getV();
            }
        }

        // Outlet pressure
        if (eastBC.hasPressure())
        {
            for (int j = 0; j < ny; ++j)
            {
                fields.p(nx - 1, j) =
                    eastBC.getPressure();
            }
        }
    }


    // ================================================================
    // U MOMENTUM
    // ================================================================

    void StaggeredSIMPLE::assembleUMomentum()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();

        const std::size_t nU =
            static_cast<std::size_t>(
                (nx + 1) * ny);

        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> values;

        uRHS = Vector(nU, 0.0);

        const double dx = mesh.getDx();
        const double dy = mesh.getDy();

        const double Ae = mesh.eastWestFaceArea();
        const double An = mesh.northSouthFaceArea();

        rows.reserve(5 * nU);
        cols.reserve(5 * nU);
        values.reserve(5 * nU);

        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i <= nx; ++i)
            {
                const std::size_t P =
                    static_cast<std::size_t>(
                        fields.uIndex(i, j));

                const bool fixed =
                    (i == 0 && westBC.hasU()) ||
                    (j == 0 && southBC.hasU()) ||
                    (j == ny - 1 && northBC.hasU());

                if (fixed)
                {
                    rows.push_back(P);
                    cols.push_back(P);
                    values.push_back(1.0);

                    double prescribed =
                        fields.u(i, j);

                    if (i == 0 && westBC.hasU())
                    {
                        prescribed = westBC.getU();
                    }

                    if (j == 0 && southBC.hasU())
                    {
                        prescribed = southBC.getU();
                    }

                    if (j == ny - 1 && northBC.hasU())
                    {
                        prescribed = northBC.getU();
                    }

                    uRHS[P] = prescribed;
                    dU[P] = 0.0;

                    continue;
                }

                const double aE =
                    (i < nx)
                    ? mu * Ae / dx
                    : 0.0;

                const double aW =
                    (i > 0)
                    ? mu * Ae / dx
                    : 0.0;

                const double aN =
                    (j < ny - 1)
                    ? mu * An / dy
                    : 0.0;

                const double aS =
                    (j > 0)
                    ? mu * An / dy
                    : 0.0;

                const double aPNoRelax =
                    aE + aW + aN + aS + SMALL;

                const double aP =
                    safeDiagonal(
                        aPNoRelax /
                        relaxationVelocity);

                dU[P] =
                    Ae / aP;

                const double pW =
                    (i > 0)
                    ? fields.p(i - 1, j)
                    : fields.p(0, j);

                const double pE =
                    (i < nx)
                    ? fields.p(i, j)
                    : eastBC.hasPressure()
                    ? eastBC.getPressure()
                    : fields.p(nx - 1, j);

                const double source =
                    (pW - pE) * Ae +
                    ((1.0 - relaxationVelocity) /
                        relaxationVelocity) *
                    aPNoRelax *
                    fields.u(i, j);

                rows.push_back(P);
                cols.push_back(P);
                values.push_back(aP);

                if (i < nx)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            fields.uIndex(i + 1, j)));

                    values.push_back(-aE);
                }

                if (i > 0)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            fields.uIndex(i - 1, j)));

                    values.push_back(-aW);
                }

                if (j < ny - 1)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            fields.uIndex(i, j + 1)));

                    values.push_back(-aN);
                }

                if (j > 0)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            fields.uIndex(i, j - 1)));

                    values.push_back(-aS);
                }

                uRHS[P] = source;
            }
        }

        uMatrix =
            SparseMatrix(nU, nU);

        uMatrix.setFromTriplets(
            rows,
            cols,
            values);
    }


    void StaggeredSIMPLE::solveUMomentum()
    {
        assembleUMomentum();

        const std::size_t nU =
            static_cast<std::size_t>(
                (mesh.getNx() + 1) *
                mesh.getNy());

        Vector solution(nU, 0.0);

        for (std::size_t k = 0;
            k < nU;
            ++k)
        {
            solution[k] =
                fields.uData()[k];
        }

        BiCGSTAB solver(
            1.0e-8,
            10000);

        BiCGSTABResult result =
            solver.solve(
                uMatrix,
                uRHS,
                solution,
                false);

        if (!result.converged)
        {
            throw std::runtime_error(
                "StaggeredSIMPLE: U momentum solver failed to converge");
        }

        for (std::size_t k = 0;
            k < nU;
            ++k)
        {
            fields.uData()[k] =
                solution[k];
        }
    }


    // ================================================================
    // V MOMENTUM
    // ================================================================

    void StaggeredSIMPLE::assembleVMomentum()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();

        const std::size_t nV =
            static_cast<std::size_t>(
                nx * (ny + 1));

        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> values;

        vRHS = Vector(nV, 0.0);

        const double dx = mesh.getDx();
        const double dy = mesh.getDy();

        const double Ae = mesh.eastWestFaceArea();
        const double An = mesh.northSouthFaceArea();

        rows.reserve(5 * nV);
        cols.reserve(5 * nV);
        values.reserve(5 * nV);

        for (int j = 0; j <= ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t P =
                    static_cast<std::size_t>(
                        fields.vIndex(i, j));

                const bool fixed =
                    (j == 0 && southBC.hasV()) ||
                    (j == ny && northBC.hasV()) ||
                    (i == 0 && westBC.hasV()) ||
                    (i == nx - 1 && eastBC.hasV());

                if (fixed)
                {
                    rows.push_back(P);
                    cols.push_back(P);
                    values.push_back(1.0);

                    double prescribed =
                        fields.v(i, j);

                    if (j == 0 && southBC.hasV())
                    {
                        prescribed = southBC.getV();
                    }

                    if (j == ny && northBC.hasV())
                    {
                        prescribed = northBC.getV();
                    }

                    if (i == 0 && westBC.hasV())
                    {
                        prescribed = westBC.getV();
                    }

                    if (i == nx - 1 && eastBC.hasV())
                    {
                        prescribed = eastBC.getV();
                    }

                    vRHS[P] = prescribed;
                    dV[P] = 0.0;

                    continue;
                }

                const double aE =
                    (i < nx - 1)
                    ? mu * Ae / dx
                    : 0.0;

                const double aW =
                    (i > 0)
                    ? mu * Ae / dx
                    : 0.0;

                const double aN =
                    (j < ny)
                    ? mu * An / dy
                    : 0.0;

                const double aS =
                    (j > 0)
                    ? mu * An / dy
                    : 0.0;

                const double aPNoRelax =
                    aE + aW + aN + aS + SMALL;

                const double aP =
                    safeDiagonal(
                        aPNoRelax /
                        relaxationVelocity);

                dV[P] =
                    An / aP;

                const double pS =
                    (j > 0)
                    ? fields.p(i, j - 1)
                    : fields.p(i, 0);

                const double pN =
                    (j < ny)
                    ? fields.p(i, j)
                    : fields.p(i, ny - 1);

                const double source =
                    (pS - pN) * An +
                    ((1.0 - relaxationVelocity) /
                        relaxationVelocity) *
                    aPNoRelax *
                    fields.v(i, j);

                rows.push_back(P);
                cols.push_back(P);
                values.push_back(aP);

                if (i < nx - 1)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            fields.vIndex(i + 1, j)));

                    values.push_back(-aE);
                }

                if (i > 0)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            fields.vIndex(i - 1, j)));

                    values.push_back(-aW);
                }

                if (j < ny)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            fields.vIndex(i, j + 1)));

                    values.push_back(-aN);
                }

                if (j > 0)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            fields.vIndex(i, j - 1)));

                    values.push_back(-aS);
                }

                vRHS[P] = source;
            }
        }

        vMatrix =
            SparseMatrix(nV, nV);

        vMatrix.setFromTriplets(
            rows,
            cols,
            values);
    }


    void StaggeredSIMPLE::solveVMomentum()
    {
        assembleVMomentum();

        const std::size_t nV =
            static_cast<std::size_t>(
                mesh.getNx() *
                (mesh.getNy() + 1));

        Vector solution(nV, 0.0);

        for (std::size_t k = 0;
            k < nV;
            ++k)
        {
            solution[k] =
                fields.vData()[k];
        }

        BiCGSTAB solver(
            1.0e-8,
            10000);

        BiCGSTABResult result =
            solver.solve(
                vMatrix,
                vRHS,
                solution,
                false);

        if (!result.converged)
        {
            throw std::runtime_error(
                "StaggeredSIMPLE: V momentum solver failed to converge");
        }

        for (std::size_t k = 0;
            k < nV;
            ++k)
        {
            fields.vData()[k] =
                solution[k];
        }
    }


    // ================================================================
    // PRESSURE CORRECTION
    // ================================================================

    void StaggeredSIMPLE::assemblePressureCorrection()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();

        const std::size_t nP =
            static_cast<std::size_t>(
                nx * ny);

        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> values;

        pressureRHS =
            Vector(nP, 0.0);

        const double dx = mesh.getDx();
        const double dy = mesh.getDy();

        const double Ae =
            mesh.eastWestFaceArea();

        const double An =
            mesh.northSouthFaceArea();

        rows.reserve(5 * nP);
        cols.reserve(5 * nP);
        values.reserve(5 * nP);

        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t P =
                    static_cast<std::size_t>(
                        fields.pIndex(i, j));

                const bool isReference =
                    (i == nx - 1 &&
                        j == ny - 1);

                if (isReference)
                {
                    rows.push_back(P);
                    cols.push_back(P);
                    values.push_back(1.0);

                    pressureRHS[P] = 0.0;

                    continue;
                }

                double aE = 0.0;
                double aW = 0.0;
                double aN = 0.0;
                double aS = 0.0;

                if (i < nx - 1)
                {
                    aE =
                        rho * Ae *
                        dU[
                            static_cast<std::size_t>(
                                fields.uIndex(i + 1, j))
                        ] / dx;
                }
                else if (eastBC.hasPressure())
                {
                    aE =
                        rho * Ae *
                        dU[
                            static_cast<std::size_t>(
                                fields.uIndex(nx, j))
                        ] / dx;
                }

                if (i > 0)
                {
                    aW =
                        rho * Ae *
                        dU[
                            static_cast<std::size_t>(
                                fields.uIndex(i, j))
                        ] / dx;
                }

                if (j < ny - 1)
                {
                    aN =
                        rho * An *
                        dV[
                            static_cast<std::size_t>(
                                fields.vIndex(i, j + 1))
                        ] / dy;
                }

                if (j > 0)
                {
                    aS =
                        rho * An *
                        dV[
                            static_cast<std::size_t>(
                                fields.vIndex(i, j))
                        ] / dy;
                }

                const double aP =
                    safeDiagonal(
                        aE + aW + aN + aS);

                const double continuity =
                    rho * Ae *
                    fields.u(i + 1, j)
                    -
                    rho * Ae *
                    fields.u(i, j)
                    +
                    rho * An *
                    fields.v(i, j + 1)
                    -
                    rho * An *
                    fields.v(i, j);

                rows.push_back(P);
                cols.push_back(P);
                values.push_back(aP);

                if (i < nx - 1 &&
                    aE > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            fields.pIndex(
                                i + 1,
                                j)));

                    values.push_back(-aE);
                }

                if (i > 0 &&
                    aW > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            fields.pIndex(
                                i - 1,
                                j)));

                    values.push_back(-aW);
                }

                if (j < ny - 1 &&
                    aN > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            fields.pIndex(
                                i,
                                j + 1)));

                    values.push_back(-aN);
                }

                if (j > 0 &&
                    aS > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            fields.pIndex(
                                i,
                                j - 1)));

                    values.push_back(-aS);
                }

                pressureRHS[P] =
                    -continuity;
            }
        }

        pressureMatrix =
            SparseMatrix(nP, nP);

        pressureMatrix.setFromTriplets(
            rows,
            cols,
            values);
    }


    void StaggeredSIMPLE::solvePressureCorrection()
    {
        assemblePressureCorrection();

        const std::size_t nP =
            static_cast<std::size_t>(
                mesh.getNx() *
                mesh.getNy());

        pressureCorrection =
            Vector(nP, 0.0);

        BiCGSTAB solver(
            1.0e-8,
            10000);

        BiCGSTABResult result =
            solver.solve(
                pressureMatrix,
                pressureRHS,
                pressureCorrection,
                false);

        if (!result.converged)
        {
            throw std::runtime_error(
                "StaggeredSIMPLE: pressure correction solver failed to converge");
        }
    }


    // ================================================================
    // CORRECTIONS
    // ================================================================

    void StaggeredSIMPLE::correctPressure()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();

        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                fields.p(i, j) +=
                    relaxationPressure *
                    pressureCorrection[
                        static_cast<std::size_t>(
                            fields.pIndex(i, j))];
            }
        }

        if (eastBC.hasPressure())
        {
            for (int j = 0; j < ny; ++j)
            {
                fields.p(nx - 1, j) =
                    eastBC.getPressure();
            }
        }
    }


    void StaggeredSIMPLE::correctVelocities()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();

        // Interior U faces
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 1; i < nx; ++i)
            {
                const std::size_t idx =
                    static_cast<std::size_t>(
                        fields.uIndex(i, j));

                const double pW =
                    pressureCorrection[
                        static_cast<std::size_t>(
                            fields.pIndex(i - 1, j))];

                const double pE =
                    pressureCorrection[
                        static_cast<std::size_t>(
                            fields.pIndex(i, j))];

                fields.u(i, j) +=
                    dU[idx] *
                    (pW - pE);
            }
        }

        // Outlet U faces
        if (eastBC.hasPressure())
        {
            for (int j = 0; j < ny; ++j)
            {
                const std::size_t idx =
                    static_cast<std::size_t>(
                        fields.uIndex(nx, j));

                const double pW =
                    pressureCorrection[
                        static_cast<std::size_t>(
                            fields.pIndex(nx - 1, j))];

                const double pE = 0.0;

                fields.u(nx, j) +=
                    dU[idx] *
                    (pW - pE);
            }
        }

        // Interior V faces
        for (int j = 1; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t idx =
                    static_cast<std::size_t>(
                        fields.vIndex(i, j));

                const double pS =
                    pressureCorrection[
                        static_cast<std::size_t>(
                            fields.pIndex(i, j - 1))];

                const double pN =
                    pressureCorrection[
                        static_cast<std::size_t>(
                            fields.pIndex(i, j))];

                fields.v(i, j) +=
                    dV[idx] *
                    (pS - pN);
            }
        }

        applyBoundaryConditions();
    }


    // ================================================================
    // RESIDUAL
    // ================================================================

    double StaggeredSIMPLE::calculateResidual()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();

        const double Ae =
            mesh.eastWestFaceArea();

        const double An =
            mesh.northSouthFaceArea();

        double totalResidual = 0.0;

        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                const double continuity =
                    rho * Ae *
                    fields.u(i + 1, j)
                    -
                    rho * Ae *
                    fields.u(i, j)
                    +
                    rho * An *
                    fields.v(i, j + 1)
                    -
                    rho * An *
                    fields.v(i, j);

                totalResidual +=
                    std::abs(continuity);
            }
        }

        return totalResidual;
    }


    bool StaggeredSIMPLE::checkConvergence()
    {
        residual =
            calculateResidual();

        return residual <
            convergenceTolerance;
    }


    // ================================================================
    // ONE SIMPLE ITERATION
    // ================================================================

    void StaggeredSIMPLE::step()
    {
        if (finished())
        {
            return;
        }

        ++iteration;

        applyBoundaryConditions();

        // 1. Solve U momentum
        solveUMomentum();

        // 2. Solve V momentum
        solveVMomentum();

        // 3. Solve pressure correction
        solvePressureCorrection();

        // 4. Correct pressure
        correctPressure();

        // 5. Correct velocities
        correctVelocities();

        // 6. Calculate continuity residual
        residual =
            calculateResidual();

        std::cout
            << "SIMPLE iteration "
            << iteration
            << " | Residual = "
            << residual
            << '\n';
    }


    // ================================================================
    // FINISHED?
    // ================================================================

    bool StaggeredSIMPLE::finished() const
    {
        return
            residual < convergenceTolerance ||
            iteration >= maxIterations;
    }


    // ================================================================
    // COMPLETE SOLVE
    // ================================================================

    void StaggeredSIMPLE::solve()
    {
        iteration = 0;

        residual =
            std::numeric_limits<double>::infinity();

        while (!finished())
        {
            step();
        }

        if (residual <
            convergenceTolerance)
        {
            std::cout
                << "\nStaggered SIMPLE converged after "
                << iteration
                << " iterations.\n";
        }
        else
        {
            std::cout
                << "\nStaggered SIMPLE reached maximum iterations.\n"
                << "Final residual = "
                << residual
                << '\n';
        }
    }

}