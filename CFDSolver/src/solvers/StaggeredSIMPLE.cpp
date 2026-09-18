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
            if (!std::isfinite(value) ||
                std::abs(value) < SMALL)
            {
                return SMALL;
            }

            return value;
        }

        double positive(double value)
        {
            return std::max(value, 0.0);
        }
    }


    // ================================================================
    // DIAGONAL EXTRACTION
    // ================================================================

    static double getDiagonal(
        const SparseMatrix& matrix,
        std::size_t P)
    {
        if (P >= matrix.rows())
            return 0.0;

        const auto& rp = matrix.rowPtr();
        const auto& ci = matrix.colIndices();
        const auto& vv = matrix.values();

        if (rp.size() < matrix.rows() + 1)
            return 0.0;

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
                return vv[k];
        }

        return 0.0;
    }


    // ================================================================
    // RELATIVE ALGEBRAIC RESIDUAL
    // ================================================================

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

        for (std::size_t row = 0;
            row < matrix.rows();
            ++row)
        {
            double ax = 0.0;

            for (std::size_t k = rowPtr[row];
                k < rowPtr[row + 1];
                ++k)
            {
                ax += values[k] * solution[colIdx[k]];
            }

            l1Residual += std::abs(ax - rhs[row]);
            l1Rhs += std::abs(rhs[row]);
        }

        return l1Residual / std::max(l1Rhs, 1.0);
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
        residual(
            std::numeric_limits<double>::infinity()),
        continuityResidual_(
            std::numeric_limits<double>::infinity()),
        uMomentumResidual_(
            std::numeric_limits<double>::infinity()),
        vMomentumResidual_(
            std::numeric_limits<double>::infinity())
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();

        dU.assign(
            static_cast<std::size_t>((nx + 1) * ny),
            0.0);

        dV.assign(
            static_cast<std::size_t>(nx * (ny + 1)),
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
        if (!std::isfinite(value) ||
            value <= 0.0)
        {
            throw std::runtime_error(
                "StaggeredSIMPLE: continuity tolerance must be positive");
        }

        convergenceTolerance = value;
    }


    void StaggeredSIMPLE::setMomentumConvergenceTolerance(double value)
    {
        if (!std::isfinite(value) ||
            value <= 0.0)
        {
            throw std::runtime_error(
                "StaggeredSIMPLE: momentum tolerance must be positive");
        }

        momentumTolerance_ = value;
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
        if (!std::isfinite(value) ||
            value <= 0.0)
        {
            throw std::runtime_error(
                "StaggeredSIMPLE: density must be positive");
        }

        rho = value;
    }


    void StaggeredSIMPLE::setViscosity(double value)
    {
        if (!std::isfinite(value) ||
            value <= 0.0)
        {
            throw std::runtime_error(
                "StaggeredSIMPLE: viscosity must be positive");
        }

        mu = value;
    }


    // ================================================================
    // GETTERS
    // ================================================================

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


    double StaggeredSIMPLE::getMomentumConvergenceTolerance() const
    {
        return momentumTolerance_;
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
        return continuityResidual_;
    }


    double StaggeredSIMPLE::getContinuityResidual() const
    {
        return continuityResidual_;
    }


    double StaggeredSIMPLE::getUMomentumResidual() const
    {
        return uMomentumResidual_;
    }


    double StaggeredSIMPLE::getVMomentumResidual() const
    {
        return vMomentumResidual_;
    }


    // ================================================================
    // BOUNDARY CONDITIONS
    // ================================================================

    void StaggeredSIMPLE::applyBoundaryConditions()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();

        // ------------------------------------------------------------
        // WEST / EAST U boundaries
        //
        // U is normal to the west/east walls.
        //
        // This is important for the cavity:
        //
        // west U = 0
        // east U = 0
        //
        // Without explicitly enforcing the east U boundary, the
        // cavity can behave as though the right wall is partially open.
        // ------------------------------------------------------------

        for (int j = 0; j < ny; ++j)
        {
            if (westBC.hasU())
            {
                fields.u(0, j) =
                    westBC.getU();
            }

            if (eastBC.hasU())
            {
                fields.u(nx, j) =
                    eastBC.getU();
            }
        }

        // ------------------------------------------------------------
        // SOUTH / NORTH V boundaries
        //
        // V is normal to the south/north walls.
        // ------------------------------------------------------------

        for (int i = 0; i < nx; ++i)
        {
            if (southBC.hasV())
            {
                fields.v(i, 0) =
                    southBC.getV();
            }

            if (northBC.hasV())
            {
                fields.v(i, ny) =
                    northBC.getV();
            }
        }

        /*
         * Tangential wall velocities are imposed through the
         * half-cell diffusion terms in the momentum equations.
         */
    }


    // ================================================================
    // U MOMENTUM ASSEMBLY
    // ================================================================

    void StaggeredSIMPLE::assembleUMomentum()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();

        const std::size_t nU =
            static_cast<std::size_t>((nx + 1) * ny);

        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> values;

        uRHS = Vector(nU, 0.0);

        const double dx = mesh.getDx();
        const double dy = mesh.getDy();

        const double Ae =
            mesh.eastWestFaceArea();

        const double An =
            mesh.northSouthFaceArea();

        const double De =
            mu * Ae / dx;

        const double Dw =
            mu * Ae / dx;

        const double Dn =
            mu * An / dy;

        const double Ds =
            mu * An / dy;

        rows.reserve(7 * nU);
        cols.reserve(7 * nU);
        values.reserve(7 * nU);

        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i <= nx; ++i)
            {
                const std::size_t P =
                    static_cast<std::size_t>(
                        fields.uIndex(i, j));

                // ----------------------------------------------------
                // WEST NORMAL VELOCITY BOUNDARY
                // ----------------------------------------------------

                const bool westFixed =
                    (i == 0 && westBC.hasU());

                // ----------------------------------------------------
                // EAST NORMAL VELOCITY BOUNDARY
                //
                // This is the critical cavity correction.
                // A stationary or moving wall has prescribed U.
                // An outlet does not.
                // ----------------------------------------------------

                const bool eastFixed =
                    (i == nx && eastBC.hasU());

                if (westFixed || eastFixed)
                {
                    rows.push_back(P);
                    cols.push_back(P);
                    values.push_back(1.0);

                    uRHS[P] =
                        westFixed
                        ? westBC.getU()
                        : eastBC.getU();

                    dU[P] = 0.0;

                    continue;
                }

                // ----------------------------------------------------
                // DIFFUSIVE COEFFICIENTS
                // ----------------------------------------------------

                const double aE_diff =
                    (i < nx) ? De : 0.0;

                const double aW_diff =
                    (i > 0) ? Dw : 0.0;

                const double aN_diff =
                    (j < ny - 1)
                    ? Dn
                    : (northBC.hasU()
                        ? 2.0 * Dn
                        : 0.0);

                const double aS_diff =
                    (j > 0)
                    ? Ds
                    : (southBC.hasU()
                        ? 2.0 * Ds
                        : 0.0);

                // ----------------------------------------------------
                // CONVECTIVE MASS FLUXES
                //
                // First-order upwind discretisation.
                // ----------------------------------------------------

                double Fe = 0.0;
                double Fw = 0.0;
                double Fn = 0.0;
                double Fs = 0.0;

                // East/west fluxes
                if (i < nx)
                {
                    Fe =
                        rho * Ae *
                        0.5 *
                        (
                            fields.u(i, j) +
                            fields.u(i + 1, j)
                            );
                }

                if (i > 0)
                {
                    Fw =
                        rho * Ae *
                        0.5 *
                        (
                            fields.u(i, j) +
                            fields.u(i - 1, j)
                            );
                }

                // North/south fluxes.
                //
                // U control volumes lie around a vertical velocity
                // face, so interpolate V to the north/south U-CV face.

                if (j < ny - 1)
                {
                    double vNorth = 0.0;

                    if (i == 0)
                    {
                        vNorth =
                            fields.v(0, j + 1);
                    }
                    else if (i == nx)
                    {
                        vNorth =
                            fields.v(nx - 1, j + 1);
                    }
                    else
                    {
                        vNorth =
                            0.5 *
                            (
                                fields.v(i - 1, j + 1) +
                                fields.v(i, j + 1)
                                );
                    }

                    Fn =
                        rho * An * vNorth;
                }

                if (j > 0)
                {
                    double vSouth = 0.0;

                    if (i == 0)
                    {
                        vSouth =
                            fields.v(0, j);
                    }
                    else if (i == nx)
                    {
                        vSouth =
                            fields.v(nx - 1, j);
                    }
                    else
                    {
                        vSouth =
                            0.5 *
                            (
                                fields.v(i - 1, j) +
                                fields.v(i, j)
                                );
                    }

                    Fs =
                        rho * An * vSouth;
                }

                // ----------------------------------------------------
                // UPWIND COEFFICIENTS
                // ----------------------------------------------------

                const double aE =
                    aE_diff + positive(-Fe);

                const double aW =
                    aW_diff + positive(Fw);

                const double aN =
                    aN_diff + positive(-Fn);

                const double aS =
                    aS_diff + positive(Fs);

                const double aPNoRelax =
                    aE +
                    aW +
                    aN +
                    aS +
                    (Fe - Fw + Fn - Fs > 0.0
                        ? Fe - Fw + Fn - Fs
                        : 0.0) +
                    SMALL;

                const double aP =
                    safeDiagonal(
                        aPNoRelax /
                        relaxationVelocity);

                dU[P] =
                    Ae / aP;

                // ----------------------------------------------------
                // PRESSURE SOURCE
                // ----------------------------------------------------

                const double pW =
                    (i > 0)
                    ? fields.p(i - 1, j)
                    : fields.p(0, j);

                const double pE =
                    (i < nx)
                    ? fields.p(i, j)
                    : (eastBC.hasPressure()
                        ? eastBC.getPressure()
                        : fields.p(nx - 1, j));

                double source =
                    (pW - pE) * Ae;

                // ----------------------------------------------------
                // UNDER-RELAXATION SOURCE
                // ----------------------------------------------------

                source +=
                    ((1.0 - relaxationVelocity) /
                        relaxationVelocity) *
                    aPNoRelax *
                    fields.u(i, j);

                // ----------------------------------------------------
                // TANGENTIAL WALL VELOCITY
                // ----------------------------------------------------

                if (j == 0 &&
                    southBC.hasU())
                {
                    source +=
                        2.0 *
                        Dn *
                        southBC.getU();
                }

                if (j == ny - 1 &&
                    northBC.hasU())
                {
                    source +=
                        2.0 *
                        Dn *
                        northBC.getU();
                }

                // ----------------------------------------------------
                // MATRIX
                // ----------------------------------------------------

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


    // ================================================================
    // U MOMENTUM SOLVE
    // ================================================================

    void StaggeredSIMPLE::solveUMomentum()
    {
        assembleUMomentum();

        const std::size_t nU =
            static_cast<std::size_t>(
                (mesh.getNx() + 1) *
                mesh.getNy());

        Vector solution(nU, 0.0);

        for (std::size_t k = 0; k < nU; ++k)
        {
            solution[k] =
                fields.uData()[k];
        }

        BiCGSTAB solver(
            1.0e-8,
            10000);

        const BiCGSTABResult result =
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

        for (std::size_t k = 0; k < nU; ++k)
        {
            fields.uData()[k] =
                solution[k];
        }
    }


    // ================================================================
    // V MOMENTUM ASSEMBLY
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

        const double Ae =
            mesh.eastWestFaceArea();

        const double An =
            mesh.northSouthFaceArea();

        const double De =
            mu * An / dx;

        const double Dw =
            mu * An / dx;

        const double Dn =
            mu * Ae / dy;

        const double Ds =
            mu * Ae / dy;

        rows.reserve(7 * nV);
        cols.reserve(7 * nV);
        values.reserve(7 * nV);

        for (int j = 0; j <= ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t P =
                    static_cast<std::size_t>(
                        fields.vIndex(i, j));

                // ----------------------------------------------------
                // SOUTH / NORTH NORMAL VELOCITY BOUNDARIES
                // ----------------------------------------------------

                const bool southFixed =
                    (j == 0 && southBC.hasV());

                const bool northFixed =
                    (j == ny && northBC.hasV());

                if (southFixed || northFixed)
                {
                    rows.push_back(P);
                    cols.push_back(P);
                    values.push_back(1.0);

                    vRHS[P] =
                        southFixed
                        ? southBC.getV()
                        : northBC.getV();

                    dV[P] = 0.0;

                    continue;
                }

                // ----------------------------------------------------
                // DIFFUSION
                // ----------------------------------------------------

                const double aE_diff =
                    (i < nx - 1)
                    ? De
                    : (eastBC.hasV()
                        ? 2.0 * De
                        : 0.0);

                const double aW_diff =
                    (i > 0)
                    ? Dw
                    : (westBC.hasV()
                        ? 2.0 * Dw
                        : 0.0);

                const double aN_diff =
                    (j < ny)
                    ? Dn
                    : 0.0;

                const double aS_diff =
                    (j > 0)
                    ? Ds
                    : 0.0;

                // ----------------------------------------------------
                // CONVECTIVE MASS FLUXES
                // ----------------------------------------------------

                double Fe = 0.0;
                double Fw = 0.0;
                double Fn = 0.0;
                double Fs = 0.0;

                // East/west U fluxes
                if (i < nx - 1)
                {
                    double uEast =
                        0.5 *
                        (
                            fields.u(i + 1, std::max(0, j - 1)) +
                            fields.u(i + 1, std::min(ny - 1, j))
                            );

                    Fe =
                        rho * An * uEast;
                }
                else
                {
                    // At an east boundary the normal U velocity
                    // lives directly on the east face.
                    double uEast =
                        0.5 *
                        (
                            fields.u(nx, std::max(0, j - 1)) +
                            fields.u(nx, std::min(ny - 1, j))
                            );

                    Fe =
                        rho * An * uEast;
                }

                if (i > 0)
                {
                    double uWest =
                        0.5 *
                        (
                            fields.u(i, std::max(0, j - 1)) +
                            fields.u(i, std::min(ny - 1, j))
                            );

                    Fw =
                        rho * An * uWest;
                }
                else
                {
                    double uWest =
                        0.5 *
                        (
                            fields.u(0, std::max(0, j - 1)) +
                            fields.u(0, std::min(ny - 1, j))
                            );

                    Fw =
                        rho * An * uWest;
                }

                // North/south V fluxes
                if (j < ny)
                {
                    Fn =
                        rho * Ae *
                        fields.v(i, j);
                }

                if (j > 0)
                {
                    Fs =
                        rho * Ae *
                        fields.v(i, j);
                }

                // ----------------------------------------------------
                // UPWIND COEFFICIENTS
                // ----------------------------------------------------

                const double aE =
                    aE_diff + positive(-Fe);

                const double aW =
                    aW_diff + positive(Fw);

                const double aN =
                    aN_diff + positive(-Fn);

                const double aS =
                    aS_diff + positive(Fs);

                const double netFlux =
                    Fe - Fw + Fn - Fs;

                const double aPNoRelax =
                    aE +
                    aW +
                    aN +
                    aS +
                    std::max(netFlux, 0.0) +
                    SMALL;

                const double aP =
                    safeDiagonal(
                        aPNoRelax /
                        relaxationVelocity);

                dV[P] =
                    An / aP;

                // ----------------------------------------------------
                // PRESSURE SOURCE
                // ----------------------------------------------------

                const double pS =
                    (j > 0)
                    ? fields.p(i, j - 1)
                    : fields.p(i, 0);

                const double pN =
                    (j < ny)
                    ? fields.p(i, j)
                    : fields.p(i, ny - 1);

                double source =
                    (pS - pN) * An;

                // ----------------------------------------------------
                // UNDER-RELAXATION
                // ----------------------------------------------------

                source +=
                    ((1.0 - relaxationVelocity) /
                        relaxationVelocity) *
                    aPNoRelax *
                    fields.v(i, j);

                // ----------------------------------------------------
                // TANGENTIAL WALL VELOCITIES
                // ----------------------------------------------------

                if (i == 0 &&
                    westBC.hasV())
                {
                    source +=
                        2.0 *
                        Dw *
                        westBC.getV();
                }

                if (i == nx - 1 &&
                    eastBC.hasV())
                {
                    source +=
                        2.0 *
                        De *
                        eastBC.getV();
                }

                // ----------------------------------------------------
                // MATRIX
                // ----------------------------------------------------

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


    // ================================================================
    // V MOMENTUM SOLVE
    // ================================================================

    void StaggeredSIMPLE::solveVMomentum()
    {
        assembleVMomentum();

        const std::size_t nV =
            static_cast<std::size_t>(
                mesh.getNx() *
                (mesh.getNy() + 1));

        Vector solution(nV, 0.0);

        for (std::size_t k = 0; k < nV; ++k)
        {
            solution[k] =
                fields.vData()[k];
        }

        BiCGSTAB solver(
            1.0e-8,
            10000);

        const BiCGSTABResult result =
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

        for (std::size_t k = 0; k < nV; ++k)
        {
            fields.vData()[k] =
                solution[k];
        }
    }


    // ================================================================
    // PRESSURE CORRECTION ASSEMBLY
    // ================================================================

    void StaggeredSIMPLE::assemblePressureCorrection()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();

        const std::size_t nP =
            static_cast<std::size_t>(nx * ny);

        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> values;

        pressureRHS =
            Vector(nP, 0.0);

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

                // ----------------------------------------------------
                // Pressure reference for closed domains
                // ----------------------------------------------------

                const bool isReference =
                    !eastBC.hasPressure() &&
                    i == nx - 1 &&
                    j == ny - 1;

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

                // ----------------------------------------------------
                // EAST
                // ----------------------------------------------------

                if (i < nx - 1)
                {
                    const std::size_t idx =
                        static_cast<std::size_t>(
                            fields.uIndex(i + 1, j));

                    aE =
                        rho * Ae * dU[idx];
                }
                else if (eastBC.hasPressure())
                {
                    const std::size_t idx =
                        static_cast<std::size_t>(
                            fields.uIndex(nx, j));

                    aE =
                        rho * Ae * dU[idx];
                }

                // ----------------------------------------------------
                // WEST
                // ----------------------------------------------------

                if (i > 0)
                {
                    const std::size_t idx =
                        static_cast<std::size_t>(
                            fields.uIndex(i, j));

                    aW =
                        rho * Ae * dU[idx];
                }

                // ----------------------------------------------------
                // NORTH
                // ----------------------------------------------------

                if (j < ny - 1)
                {
                    const std::size_t idx =
                        static_cast<std::size_t>(
                            fields.vIndex(i, j + 1));

                    aN =
                        rho * An * dV[idx];
                }

                // ----------------------------------------------------
                // SOUTH
                // ----------------------------------------------------

                if (j > 0)
                {
                    const std::size_t idx =
                        static_cast<std::size_t>(
                            fields.vIndex(i, j));

                    aS =
                        rho * An * dV[idx];
                }

                const double aP =
                    safeDiagonal(
                        aE +
                        aW +
                        aN +
                        aS);

                // ----------------------------------------------------
                // CONTINUITY
                // ----------------------------------------------------

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

                // ----------------------------------------------------
                // MATRIX
                // ----------------------------------------------------

                rows.push_back(P);
                cols.push_back(P);
                values.push_back(aP);

                if (i < nx - 1 &&
                    aE > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            fields.pIndex(i + 1, j)));

                    values.push_back(-aE);
                }

                if (i > 0 &&
                    aW > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            fields.pIndex(i - 1, j)));

                    values.push_back(-aW);
                }

                if (j < ny - 1 &&
                    aN > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            fields.pIndex(i, j + 1)));

                    values.push_back(-aN);
                }

                if (j > 0 &&
                    aS > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            fields.pIndex(i, j - 1)));

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


    // ================================================================
    // PRESSURE CORRECTION SOLVE
    // ================================================================

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

        const BiCGSTABResult result =
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
    // PRESSURE CORRECTION
    // ================================================================

    void StaggeredSIMPLE::correctPressure()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();

        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t idx =
                    static_cast<std::size_t>(
                        fields.pIndex(i, j));

                fields.p(i, j) +=
                    relaxationPressure *
                    pressureCorrection[idx];
            }
        }
    }


    // ================================================================
    // VELOCITY CORRECTION
    // ================================================================

    void StaggeredSIMPLE::correctVelocities()
    {
        const int nx = mesh.getNx();
        const int ny = mesh.getNy();

        // ------------------------------------------------------------
        // Internal U faces
        // ------------------------------------------------------------

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

        // ------------------------------------------------------------
        // EAST PRESSURE OUTLET
        // ------------------------------------------------------------

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

                fields.u(nx, j) +=
                    dU[idx] *
                    pW;
            }
        }

        // ------------------------------------------------------------
        // Internal V faces
        // ------------------------------------------------------------

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

        // ------------------------------------------------------------
        // Re-apply all explicitly prescribed normal velocities.
        //
        // This is particularly important for the east cavity wall.
        // ------------------------------------------------------------

        applyBoundaryConditions();
    }


    // ================================================================
    // MOMENTUM RESIDUALS
    // ================================================================

    void StaggeredSIMPLE::updateMomentumResiduals()
    {
        // ------------------------------------------------------------
        // U
        // ------------------------------------------------------------

        assembleUMomentum();

        Vector u(
            static_cast<std::size_t>(
                (mesh.getNx() + 1) *
                mesh.getNy()),
            0.0);

        for (std::size_t k = 0;
            k < u.size();
            ++k)
        {
            u[k] =
                fields.uData()[k];
        }

        uMomentumResidual_ =
            relativeAlgebraicResidual(
                uMatrix,
                u,
                uRHS);

        // ------------------------------------------------------------
        // V
        // ------------------------------------------------------------

        assembleVMomentum();

        Vector v(
            static_cast<std::size_t>(
                mesh.getNx() *
                (mesh.getNy() + 1)),
            0.0);

        for (std::size_t k = 0;
            k < v.size();
            ++k)
        {
            v[k] =
                fields.vData()[k];
        }

        vMomentumResidual_ =
            relativeAlgebraicResidual(
                vMatrix,
                v,
                vRHS);
    }


    // ================================================================
    // CONTINUITY RESIDUAL
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


    // ================================================================
    // CONVERGENCE CHECK
    // ================================================================

    bool StaggeredSIMPLE::checkConvergence()
    {
        continuityResidual_ =
            calculateResidual();

        residual =
            continuityResidual_;

        return
            continuityResidual_ <
            convergenceTolerance
            &&
            uMomentumResidual_ <
            momentumTolerance_
            &&
            vMomentumResidual_ <
            momentumTolerance_;
    }


    // ================================================================
    // ONE SIMPLE ITERATION
    // ================================================================

    void StaggeredSIMPLE::step()
    {
        if (finished_)
            return;

        if (iteration >= maxIterations)
        {
            finished_ = true;
            return;
        }

        ++iteration;

        // ------------------------------------------------------------
        // Apply boundary conditions before assembling momentum
        // equations.
        // ------------------------------------------------------------

        applyBoundaryConditions();

        // ------------------------------------------------------------
        // SIMPLE
        // ------------------------------------------------------------

        solveUMomentum();

        solveVMomentum();

        solvePressureCorrection();

        correctPressure();

        correctVelocities();

        updateMomentumResiduals();

        converged_ =
            checkConvergence();

        // ------------------------------------------------------------
        // Console output
        // ------------------------------------------------------------

        std::cout
            << "Staggered SIMPLE iteration "
            << iteration
            << " | continuity = "
            << continuityResidual_
            << " | U momentum = "
            << uMomentumResidual_
            << " | V momentum = "
            << vMomentumResidual_
            << "\n";

        // ------------------------------------------------------------
        // Finished
        // ------------------------------------------------------------

        if (converged_)
        {
            finished_ = true;

            std::cout
                << "\nStaggered SIMPLE converged after "
                << iteration
                << " iterations.\n";
        }
        else if (iteration >= maxIterations)
        {
            finished_ = true;

            std::cout
                << "\nStaggered SIMPLE reached maximum iterations.\n"
                << "Final continuity residual = "
                << continuityResidual_
                << "\n"
                << "Final U-momentum residual = "
                << uMomentumResidual_
                << "\n"
                << "Final V-momentum residual = "
                << vMomentumResidual_
                << "\n";
        }
    }


    // ================================================================
    // FINISHED
    // ================================================================

    bool StaggeredSIMPLE::finished() const
    {
        return finished_;
    }


    // ================================================================
    // CONVERGED
    // ================================================================

    bool StaggeredSIMPLE::converged() const
    {
        return converged_;
    }


    // ================================================================
    // FULL SOLVE
    // ================================================================

    void StaggeredSIMPLE::solve()
    {
        iteration = 0;

        converged_ = false;
        finished_ = false;

        residual =
            std::numeric_limits<double>::infinity();

        continuityResidual_ =
            std::numeric_limits<double>::infinity();

        uMomentumResidual_ =
            std::numeric_limits<double>::infinity();

        vMomentumResidual_ =
            std::numeric_limits<double>::infinity();

        while (!finished_)
        {
            step();
        }
    }
}