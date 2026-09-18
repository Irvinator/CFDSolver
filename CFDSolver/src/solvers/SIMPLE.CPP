#include "solvers/SIMPLE.h"

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

        bool isFixedU(
            const Mesh::Cell& cell,
            const BoundaryCondition& northBC,
            const BoundaryCondition& southBC,
            const BoundaryCondition& westBC,
            const BoundaryCondition& eastBC)
        {
            if (cell.north == -1 && northBC.hasU() && northBC.getType() == BoundaryType::Wall) return true;
            if (cell.south == -1 && southBC.hasU() && southBC.getType() == BoundaryType::Wall) return true;
            if (cell.west == -1 && westBC.hasU() &&
                (westBC.getType() == BoundaryType::Inlet || westBC.getType() == BoundaryType::Wall)) return true;
            if (cell.east == -1 && eastBC.hasU() && eastBC.getType() == BoundaryType::Wall) return true;
            return false;
        }

        bool isFixedV(
            const Mesh::Cell& cell,
            const BoundaryCondition& northBC,
            const BoundaryCondition& southBC,
            const BoundaryCondition& westBC,
            const BoundaryCondition& eastBC)
        {
            if (cell.north == -1 && northBC.hasV() && northBC.getType() == BoundaryType::Wall) return true;
            if (cell.south == -1 && southBC.hasV() && southBC.getType() == BoundaryType::Wall) return true;
            if (cell.west == -1 && westBC.hasV() &&
                (westBC.getType() == BoundaryType::Inlet || westBC.getType() == BoundaryType::Wall)) return true;
            if (cell.east == -1 && eastBC.hasV() && eastBC.getType() == BoundaryType::Wall) return true;
            return false;
        }

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

    SIMPLE::SIMPLE(
        Mesh& mesh,
        Fields& fields,
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
        initializeFluxStorage();
    }

    void SIMPLE::initializeFluxStorage()
    {
        const std::size_t n = static_cast<std::size_t>(mesh.numberOfCells());
        fluxEast.assign(n, 0.0);
        fluxWest.assign(n, 0.0);
        fluxNorth.assign(n, 0.0);
        fluxSouth.assign(n, 0.0);
        invDiagU.assign(n, 0.0);
        invDiagV.assign(n, 0.0);
    }

    void SIMPLE::setPressureRelaxation(double value)
    {
        if (!std::isfinite(value) || value <= 0.0 || value > 1.0)
        {
            throw std::runtime_error("SIMPLE: pressure relaxation must be between 0 and 1");
        }
        relaxationPressure = value;
    }

    void SIMPLE::setVelocityRelaxation(double value)
    {
        if (!std::isfinite(value) || value <= 0.0 || value > 1.0)
        {
            throw std::runtime_error("SIMPLE: velocity relaxation must be between 0 and 1");
        }
        relaxationVelocity = value;
    }

    void SIMPLE::setConvergenceTolerance(double value)
    {
        if (!std::isfinite(value) || value <= 0.0)
        {
            throw std::runtime_error("SIMPLE: convergence tolerance must be positive");
        }
        convergenceTolerance = value;
    }

    void SIMPLE::setMaxIterations(std::size_t value)
    {
        if (value == 0)
        {
            throw std::runtime_error("SIMPLE: maximum iterations must be greater than zero");
        }
        maxIterations = value;
    }

    void SIMPLE::setDensity(double value)
    {
        if (!std::isfinite(value) || value <= 0.0)
        {
            throw std::runtime_error("SIMPLE: density must be positive");
        }
        rho = value;
    }

    void SIMPLE::setViscosity(double value)
    {
        if (!std::isfinite(value) || value <= 0.0)
        {
            throw std::runtime_error("SIMPLE: viscosity must be positive");
        }
        mu = value;
    }

    double SIMPLE::getPressureRelaxation() const { return relaxationPressure; }
    double SIMPLE::getVelocityRelaxation() const { return relaxationVelocity; }
    double SIMPLE::getConvergenceTolerance() const { return convergenceTolerance; }
    double SIMPLE::getDensity() const { return rho; }
    double SIMPLE::getViscosity() const { return mu; }
    std::size_t SIMPLE::getMaxIterations() const { return maxIterations; }
    std::size_t SIMPLE::getIteration() const { return iteration; }
    double SIMPLE::getResidual() const { return residual; }

    void SIMPLE::applyBoundaryConditions()
    {
        const std::size_t n = static_cast<std::size_t>(mesh.numberOfCells());

        for (std::size_t P = 0; P < n; ++P)
        {
            const auto& cell = mesh.getCell(static_cast<int>(P));

            if (cell.westBoundary)
            {
                if (westBC.hasU()) fields.velocity.getx()[P] = westBC.getU();
                if (westBC.hasV()) fields.velocity.gety()[P] = westBC.getV();
                if (westBC.hasPressure()) fields.pressure[P] = westBC.getPressure();
            }
            if (cell.eastBoundary)
            {
                if (eastBC.hasU()) fields.velocity.getx()[P] = eastBC.getU();
                if (eastBC.hasV()) fields.velocity.gety()[P] = eastBC.getV();
                if (eastBC.hasPressure()) fields.pressure[P] = eastBC.getPressure();
            }
            if (cell.northBoundary)
            {
                if (northBC.hasU()) fields.velocity.getx()[P] = northBC.getU();
                if (northBC.hasV()) fields.velocity.gety()[P] = northBC.getV();
                if (northBC.hasPressure()) fields.pressure[P] = northBC.getPressure();
            }
            if (cell.southBoundary)
            {
                if (southBC.hasU()) fields.velocity.getx()[P] = southBC.getU();
                if (southBC.hasV()) fields.velocity.gety()[P] = southBC.getV();
                if (southBC.hasPressure()) fields.pressure[P] = southBC.getPressure();
            }
        }
    }

    void SIMPLE::assembleUMomentum()
    {
        const std::size_t n = static_cast<std::size_t>(mesh.numberOfCells());
        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> values;
        uRHS = Vector(n, 0.0);

        rows.reserve(5 * n);
        cols.reserve(5 * n);
        values.reserve(5 * n);

        const double dx = mesh.getDx();
        const double dy = mesh.getDy();
        const double Ae = mesh.eastWestFaceArea();
        const double An = mesh.northSouthFaceArea();

        for (int j = 0; j < mesh.getNy(); ++j)
        {
            for (int i = 0; i < mesh.getNx(); ++i)
            {
                const std::size_t P = static_cast<std::size_t>(mesh.cellIndex(i, j));
                const auto& cell = mesh.getCell(static_cast<int>(P));

                if (isFixedU(cell, northBC, southBC, westBC, eastBC))
                {
                    rows.push_back(P);
                    cols.push_back(P);
                    values.push_back(1.0);

                    double prescribed = fields.velocity.getx()[P];
                    if (cell.west == -1 && westBC.hasU()) prescribed = westBC.getU();
                    if (cell.east == -1 && eastBC.hasU()) prescribed = eastBC.getU();
                    if (cell.north == -1 && northBC.hasU()) prescribed = northBC.getU();
                    if (cell.south == -1 && southBC.hasU()) prescribed = southBC.getU();
                    uRHS[P] = prescribed;
                    invDiagU[P] = 0.0;
                    continue;
                }

                const double Fe = fluxEast[P];
                const double Fw = fluxWest[P];
                const double Fn = fluxNorth[P];
                const double Fs = fluxSouth[P];

                double aE = mu * Ae / dx + std::max(-Fe, 0.0);
                double aW = mu * Ae / dx + std::max(Fw, 0.0);
                double aN = mu * An / dy + std::max(-Fn, 0.0);
                double aS = mu * An / dy + std::max(Fs, 0.0);

                if (cell.east == -1) aE = 0.0;
                if (cell.west == -1) aW = 0.0;
                if (cell.north == -1) aN = 0.0;
                if (cell.south == -1) aS = 0.0;

                const double aPNoRelax = aE + aW + aN + aS + (Fe - Fw + Fn - Fs);
                const double aP = safeDiagonal(aPNoRelax / relaxationVelocity);
                invDiagU[P] = 1.0 / aP;

                double source = 0.0;
                const double pW = (cell.west != -1) ? fields.pressure[static_cast<std::size_t>(cell.west)] : fields.pressure[P];
                const double pE = (cell.east != -1) ? fields.pressure[static_cast<std::size_t>(cell.east)] : (eastBC.hasPressure() ? eastBC.getPressure() : fields.pressure[P]);
                source += (pW - pE) * Ae;
                source += ((1.0 - relaxationVelocity) / relaxationVelocity) * aPNoRelax * fields.velocity.getx()[P];

                rows.push_back(P); cols.push_back(P); values.push_back(aP);
                if (cell.east != -1 && aE > 0.0) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(cell.east)); values.push_back(-aE); }
                if (cell.west != -1 && aW > 0.0) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(cell.west)); values.push_back(-aW); }
                if (cell.north != -1 && aN > 0.0) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(cell.north)); values.push_back(-aN); }
                if (cell.south != -1 && aS > 0.0) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(cell.south)); values.push_back(-aS); }
                uRHS[P] = source;
            }
        }

        uMatrix = SparseMatrix(n, n);
        uMatrix.setFromTriplets(rows, cols, values);
    }

    void SIMPLE::solveUMomentum()
    {
        assembleUMomentum();
        const std::size_t n = static_cast<std::size_t>(mesh.numberOfCells());
        Vector solution(n, 0.0);
        for (std::size_t P = 0; P < n; ++P) solution[P] = fields.velocity.getx()[P];

        BiCGSTAB solver(1.0e-8, 10000);
        const auto result = solver.solve(uMatrix, uRHS, solution, false);
        if (!result.converged)
        {
            throw std::runtime_error("SIMPLE: U momentum solver failed to converge");
        }

        for (std::size_t P = 0; P < n; ++P)
        {
            if (!isFixedU(mesh.getCell(static_cast<int>(P)), northBC, southBC, westBC, eastBC))
            {
                fields.velocity.getx()[P] = solution[P];
            }
        }
        applyBoundaryConditions();
    }

    void SIMPLE::assembleVMomentum()
    {
        const std::size_t n = static_cast<std::size_t>(mesh.numberOfCells());
        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> values;
        vRHS = Vector(n, 0.0);

        rows.reserve(5 * n);
        cols.reserve(5 * n);
        values.reserve(5 * n);

        const double dx = mesh.getDx();
        const double dy = mesh.getDy();
        const double Ae = mesh.eastWestFaceArea();
        const double An = mesh.northSouthFaceArea();

        for (int j = 0; j < mesh.getNy(); ++j)
        {
            for (int i = 0; i < mesh.getNx(); ++i)
            {
                const std::size_t P = static_cast<std::size_t>(mesh.cellIndex(i, j));
                const auto& cell = mesh.getCell(static_cast<int>(P));

                if (isFixedV(cell, northBC, southBC, westBC, eastBC))
                {
                    rows.push_back(P);
                    cols.push_back(P);
                    values.push_back(1.0);

                    double prescribed = fields.velocity.gety()[P];
                    if (cell.west == -1 && westBC.hasV()) prescribed = westBC.getV();
                    if (cell.east == -1 && eastBC.hasV()) prescribed = eastBC.getV();
                    if (cell.north == -1 && northBC.hasV()) prescribed = northBC.getV();
                    if (cell.south == -1 && southBC.hasV()) prescribed = southBC.getV();
                    vRHS[P] = prescribed;
                    invDiagV[P] = 0.0;
                    continue;
                }

                const double Fe = fluxEast[P];
                const double Fw = fluxWest[P];
                const double Fn = fluxNorth[P];
                const double Fs = fluxSouth[P];

                double aE = mu * Ae / dx + std::max(-Fe, 0.0);
                double aW = mu * Ae / dx + std::max(Fw, 0.0);
                double aN = mu * An / dy + std::max(-Fn, 0.0);
                double aS = mu * An / dy + std::max(Fs, 0.0);

                if (cell.east == -1) aE = 0.0;
                if (cell.west == -1) aW = 0.0;
                if (cell.north == -1) aN = 0.0;
                if (cell.south == -1) aS = 0.0;

                const double aPNoRelax = aE + aW + aN + aS + (Fe - Fw + Fn - Fs);
                const double aP = safeDiagonal(aPNoRelax / relaxationVelocity);
                invDiagV[P] = 1.0 / aP;

                double source = 0.0;
                const double pS = (cell.south != -1) ? fields.pressure[static_cast<std::size_t>(cell.south)] : fields.pressure[P];
                const double pN = (cell.north != -1) ? fields.pressure[static_cast<std::size_t>(cell.north)] : fields.pressure[P];
                source += (pS - pN) * An;
                source += ((1.0 - relaxationVelocity) / relaxationVelocity) * aPNoRelax * fields.velocity.gety()[P];

                rows.push_back(P); cols.push_back(P); values.push_back(aP);
                if (cell.east != -1 && aE > 0.0) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(cell.east)); values.push_back(-aE); }
                if (cell.west != -1 && aW > 0.0) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(cell.west)); values.push_back(-aW); }
                if (cell.north != -1 && aN > 0.0) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(cell.north)); values.push_back(-aN); }
                if (cell.south != -1 && aS > 0.0) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(cell.south)); values.push_back(-aS); }
                vRHS[P] = source;
            }
        }

        vMatrix = SparseMatrix(n, n);
        vMatrix.setFromTriplets(rows, cols, values);
    }

    void SIMPLE::solveVMomentum()
    {
        assembleVMomentum();
        const std::size_t n = static_cast<std::size_t>(mesh.numberOfCells());
        Vector solution(n, 0.0);
        for (std::size_t P = 0; P < n; ++P) solution[P] = fields.velocity.gety()[P];

        BiCGSTAB solver(1.0e-8, 10000);
        const auto result = solver.solve(vMatrix, vRHS, solution, false);
        if (!result.converged)
        {
            throw std::runtime_error("SIMPLE: V momentum solver failed to converge");
        }

        for (std::size_t P = 0; P < n; ++P)
        {
            if (!isFixedV(mesh.getCell(static_cast<int>(P)), northBC, southBC, westBC, eastBC))
            {
                fields.velocity.gety()[P] = solution[P];
            }
        }
        applyBoundaryConditions();
    }

    void SIMPLE::buildRhieChowFaceFluxes()
    {
        const std::size_t n = static_cast<std::size_t>(mesh.numberOfCells());
        (void)n;

        const double dx = mesh.getDx();
        const double dy = mesh.getDy();
        const double Ae = mesh.eastWestFaceArea();
        const double An = mesh.northSouthFaceArea();

        fluxEast.assign(fluxEast.size(), 0.0);
        fluxWest.assign(fluxWest.size(), 0.0);
        fluxNorth.assign(fluxNorth.size(), 0.0);
        fluxSouth.assign(fluxSouth.size(), 0.0);

        auto gradPx = [&](std::size_t P) -> double
            {
                const auto& cell = mesh.getCell(static_cast<int>(P));
                if (cell.west != -1 && cell.east != -1)
                {
                    return (fields.pressure[static_cast<std::size_t>(cell.east)] -
                        fields.pressure[static_cast<std::size_t>(cell.west)]) / (2.0 * dx);
                }
                if (cell.east != -1)
                {
                    return (fields.pressure[static_cast<std::size_t>(cell.east)] - fields.pressure[P]) / dx;
                }
                if (cell.west != -1)
                {
                    return (fields.pressure[P] - fields.pressure[static_cast<std::size_t>(cell.west)]) / dx;
                }
                return 0.0;
            };

        auto gradPy = [&](std::size_t P) -> double
            {
                const auto& cell = mesh.getCell(static_cast<int>(P));
                if (cell.south != -1 && cell.north != -1)
                {
                    return (fields.pressure[static_cast<std::size_t>(cell.north)] -
                        fields.pressure[static_cast<std::size_t>(cell.south)]) / (2.0 * dy);
                }
                if (cell.north != -1)
                {
                    return (fields.pressure[static_cast<std::size_t>(cell.north)] - fields.pressure[P]) / dy;
                }
                if (cell.south != -1)
                {
                    return (fields.pressure[P] - fields.pressure[static_cast<std::size_t>(cell.south)]) / dy;
                }
                return 0.0;
            };

        for (int j = 0; j < mesh.getNy(); ++j)
        {
            for (int i = 0; i < mesh.getNx(); ++i)
            {
                const std::size_t P = static_cast<std::size_t>(mesh.cellIndex(i, j));
                const auto& cell = mesh.getCell(static_cast<int>(P));

                if (cell.east != -1)
                {
                    const std::size_t E = static_cast<std::size_t>(cell.east);
                    const double uBar = 0.5 * (fields.velocity.getx()[P] + fields.velocity.getx()[E]);
                    const double dFace = 0.5 * (Ae * invDiagU[P] + Ae * invDiagU[E]);
                    const double dpdxFace = (fields.pressure[E] - fields.pressure[P]) / dx;
                    const double dpdxInterp = 0.5 * (gradPx(P) + gradPx(E));
                    fluxEast[P] = rho * (uBar - dFace * (dpdxFace - dpdxInterp)) * Ae;
                }
                else
                {
                    if (eastBC.getType() == BoundaryType::Wall)
                    {
                        fluxEast[P] = 0.0;
                    }
                    else if (eastBC.getType() == BoundaryType::Inlet && eastBC.hasU())
                    {
                        fluxEast[P] = rho * eastBC.getU() * Ae;
                    }
                    else
                    {
                        fluxEast[P] = rho * fields.velocity.getx()[P] * Ae;
                    }
                }

                if (cell.west != -1)
                {
                    const std::size_t W = static_cast<std::size_t>(cell.west);
                    const double uBar = 0.5 * (fields.velocity.getx()[P] + fields.velocity.getx()[W]);
                    const double dFace = 0.5 * (Ae * invDiagU[P] + Ae * invDiagU[W]);
                    const double dpdxFace = (fields.pressure[P] - fields.pressure[W]) / dx;
                    const double dpdxInterp = 0.5 * (gradPx(P) + gradPx(W));
                    fluxWest[P] = rho * (uBar - dFace * (dpdxFace - dpdxInterp)) * Ae;
                }
                else
                {
                    if (westBC.getType() == BoundaryType::Wall)
                    {
                        fluxWest[P] = 0.0;
                    }
                    else if (westBC.getType() == BoundaryType::Inlet && westBC.hasU())
                    {
                        fluxWest[P] = rho * westBC.getU() * Ae;
                    }
                    else
                    {
                        fluxWest[P] = rho * fields.velocity.getx()[P] * Ae;
                    }
                }

                if (cell.north != -1)
                {
                    const std::size_t N = static_cast<std::size_t>(cell.north);
                    const double vBar = 0.5 * (fields.velocity.gety()[P] + fields.velocity.gety()[N]);
                    const double dFace = 0.5 * (An * invDiagV[P] + An * invDiagV[N]);
                    const double dpdyFace = (fields.pressure[N] - fields.pressure[P]) / dy;
                    const double dpdyInterp = 0.5 * (gradPy(P) + gradPy(N));
                    fluxNorth[P] = rho * (vBar - dFace * (dpdyFace - dpdyInterp)) * An;
                }
                else
                {
                    if (northBC.getType() == BoundaryType::Wall)
                    {
                        fluxNorth[P] = 0.0;
                    }
                    else if (northBC.getType() == BoundaryType::Inlet && northBC.hasV())
                    {
                        fluxNorth[P] = rho * northBC.getV() * An;
                    }
                    else
                    {
                        fluxNorth[P] = rho * fields.velocity.gety()[P] * An;
                    }
                }

                if (cell.south != -1)
                {
                    const std::size_t S = static_cast<std::size_t>(cell.south);
                    const double vBar = 0.5 * (fields.velocity.gety()[P] + fields.velocity.gety()[S]);
                    const double dFace = 0.5 * (An * invDiagV[P] + An * invDiagV[S]);
                    const double dpdyFace = (fields.pressure[P] - fields.pressure[S]) / dy;
                    const double dpdyInterp = 0.5 * (gradPy(P) + gradPy(S));
                    fluxSouth[P] = rho * (vBar - dFace * (dpdyFace - dpdyInterp)) * An;
                }
                else
                {
                    if (southBC.getType() == BoundaryType::Wall)
                    {
                        fluxSouth[P] = 0.0;
                    }
                    else if (southBC.getType() == BoundaryType::Inlet && southBC.hasV())
                    {
                        fluxSouth[P] = rho * southBC.getV() * An;
                    }
                    else
                    {
                        fluxSouth[P] = rho * fields.velocity.gety()[P] * An;
                    }
                }
            }
        }
    }

    void SIMPLE::assemblePressureCorrection()
    {
        const std::size_t n = static_cast<std::size_t>(mesh.numberOfCells());
        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> values;
        pressureRHS = Vector(n, 0.0);

        rows.reserve(5 * n);
        cols.reserve(5 * n);
        values.reserve(5 * n);

        const double dx = mesh.getDx();
        const double dy = mesh.getDy();
        const double Ae = mesh.eastWestFaceArea();
        const double An = mesh.northSouthFaceArea();

        for (int j = 0; j < mesh.getNy(); ++j)
        {
            for (int i = 0; i < mesh.getNx(); ++i)
            {
                const std::size_t P = static_cast<std::size_t>(mesh.cellIndex(i, j));
                const auto& cell = mesh.getCell(static_cast<int>(P));

                const bool isReferenceCell = (i == mesh.getNx() - 1 && j == mesh.getNy() - 1);
                if (isReferenceCell)
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

                if (cell.east != -1)
                {
                    const std::size_t E = static_cast<std::size_t>(cell.east);
                    aE = rho * Ae * Ae * 0.5 * (invDiagU[P] + invDiagU[E]) / dx;
                }
                else if (eastBC.getType() == BoundaryType::Outlet)
                {
                    aE = rho * Ae * Ae * invDiagU[P] / dx;
                }

                if (cell.west != -1)
                {
                    const std::size_t W = static_cast<std::size_t>(cell.west);
                    aW = rho * Ae * Ae * 0.5 * (invDiagU[P] + invDiagU[W]) / dx;
                }

                if (cell.north != -1)
                {
                    const std::size_t N = static_cast<std::size_t>(cell.north);
                    aN = rho * An * An * 0.5 * (invDiagV[P] + invDiagV[N]) / dy;
                }

                if (cell.south != -1)
                {
                    const std::size_t S = static_cast<std::size_t>(cell.south);
                    aS = rho * An * An * 0.5 * (invDiagV[P] + invDiagV[S]) / dy;
                }

                const double aP = aE + aW + aN + aS;
                const double imbalance = fluxEast[P] - fluxWest[P] + fluxNorth[P] - fluxSouth[P];

                if (!std::isfinite(aP) || aP <= SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(P);
                    values.push_back(1.0);
                    pressureRHS[P] = 0.0;
                    continue;
                }

                rows.push_back(P); cols.push_back(P); values.push_back(aP);
                if (cell.east != -1 && aE > SMALL) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(cell.east)); values.push_back(-aE); }
                if (cell.west != -1 && aW > SMALL) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(cell.west)); values.push_back(-aW); }
                if (cell.north != -1 && aN > SMALL) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(cell.north)); values.push_back(-aN); }
                if (cell.south != -1 && aS > SMALL) { rows.push_back(P); cols.push_back(static_cast<std::size_t>(cell.south)); values.push_back(-aS); }
                pressureRHS[P] = -imbalance;
            }
        }

        pressureMatrix = SparseMatrix(n, n);
        pressureMatrix.setFromTriplets(rows, cols, values);
    }

    void SIMPLE::solvePressureCorrection()
    {
        assemblePressureCorrection();
        const std::size_t n = static_cast<std::size_t>(mesh.numberOfCells());
        pressureCorrection = Vector(n, 0.0);

        BiCGSTAB solver(1.0e-8, 10000);
        const auto result = solver.solve(pressureMatrix, pressureRHS, pressureCorrection, false);
        if (!result.converged)
        {
            throw std::runtime_error("SIMPLE: pressure correction solver failed to converge");
        }
    }

    void SIMPLE::correctPressure()
    {
        const std::size_t n = static_cast<std::size_t>(mesh.numberOfCells());
        for (std::size_t P = 0; P < n; ++P)
        {
            fields.pressure[P] += relaxationPressure * pressureCorrection[P];
        }

        for (std::size_t P = 0; P < n; ++P)
        {
            const auto& cell = mesh.getCell(static_cast<int>(P));
            if (cell.eastBoundary && eastBC.hasPressure()) fields.pressure[P] = eastBC.getPressure();
            if (cell.westBoundary && westBC.hasPressure()) fields.pressure[P] = westBC.getPressure();
            if (cell.northBoundary && northBC.hasPressure()) fields.pressure[P] = northBC.getPressure();
            if (cell.southBoundary && southBC.hasPressure()) fields.pressure[P] = southBC.getPressure();
        }
    }

    void SIMPLE::correctFaceFluxes()
    {
        const double dx = mesh.getDx();
        const double dy = mesh.getDy();
        const double Ae = mesh.eastWestFaceArea();
        const double An = mesh.northSouthFaceArea();

        for (int j = 0; j < mesh.getNy(); ++j)
        {
            for (int i = 0; i < mesh.getNx(); ++i)
            {
                const std::size_t P = static_cast<std::size_t>(mesh.cellIndex(i, j));
                const auto& cell = mesh.getCell(static_cast<int>(P));

                if (cell.east != -1)
                {
                    const std::size_t E = static_cast<std::size_t>(cell.east);
                    const double dFace = 0.5 * (invDiagU[P] + invDiagU[E]);
                    fluxEast[P] -= rho * Ae * Ae * dFace * (pressureCorrection[E] - pressureCorrection[P]) / dx;
                }
                else if (eastBC.getType() == BoundaryType::Outlet)
                {
                    fluxEast[P] -= rho * Ae * Ae * invDiagU[P] * (0.0 - pressureCorrection[P]) / dx;
                }

                if (cell.west != -1)
                {
                    const std::size_t W = static_cast<std::size_t>(cell.west);
                    const double dFace = 0.5 * (invDiagU[P] + invDiagU[W]);
                    fluxWest[P] -= rho * Ae * Ae * dFace * (pressureCorrection[P] - pressureCorrection[W]) / dx;
                }

                if (cell.north != -1)
                {
                    const std::size_t N = static_cast<std::size_t>(cell.north);
                    const double dFace = 0.5 * (invDiagV[P] + invDiagV[N]);
                    fluxNorth[P] -= rho * An * An * dFace * (pressureCorrection[N] - pressureCorrection[P]) / dy;
                }

                if (cell.south != -1)
                {
                    const std::size_t S = static_cast<std::size_t>(cell.south);
                    const double dFace = 0.5 * (invDiagV[P] + invDiagV[S]);
                    fluxSouth[P] -= rho * An * An * dFace * (pressureCorrection[P] - pressureCorrection[S]) / dy;
                }
            }
        }
    }

    void SIMPLE::correctVelocity()
    {
        const std::size_t n = static_cast<std::size_t>(mesh.numberOfCells());
        const double Ae = mesh.eastWestFaceArea();
        const double An = mesh.northSouthFaceArea();

        for (std::size_t P = 0; P < n; ++P)
        {
            const auto& cell = mesh.getCell(static_cast<int>(P));

            if (!isFixedU(cell, northBC, southBC, westBC, eastBC))
            {
                const double pW = (cell.west != -1) ? pressureCorrection[static_cast<std::size_t>(cell.west)] : pressureCorrection[P];
                const double pE = (cell.east != -1) ? pressureCorrection[static_cast<std::size_t>(cell.east)] : 0.0;
                fields.velocity.getx()[P] += Ae * invDiagU[P] * (pW - pE);
            }

            if (!isFixedV(cell, northBC, southBC, westBC, eastBC))
            {
                const double pS = (cell.south != -1) ? pressureCorrection[static_cast<std::size_t>(cell.south)] : pressureCorrection[P];
                const double pN = (cell.north != -1) ? pressureCorrection[static_cast<std::size_t>(cell.north)] : pressureCorrection[P];
                fields.velocity.gety()[P] += An * invDiagV[P] * (pS - pN);
            }
        }

        applyBoundaryConditions();
    }

    double SIMPLE::calculateResidual()
    {
        const std::size_t n = static_cast<std::size_t>(mesh.numberOfCells());
        double total = 0.0;
        for (std::size_t P = 0; P < n; ++P)
        {
            const double imbalance = fluxEast[P] - fluxWest[P] + fluxNorth[P] - fluxSouth[P];
            if (!std::isfinite(imbalance))
            {
                return std::numeric_limits<double>::infinity();
            }
            total += std::abs(imbalance);
        }
        return total;
    }

    bool SIMPLE::checkConvergence()
    {
        residual = calculateResidual();
        return residual < convergenceTolerance;
    }

    void SIMPLE::solve()
    {
        iteration = 0;
        residual = std::numeric_limits<double>::infinity();

        initializeFluxStorage();
        applyBoundaryConditions();

        for (iteration = 1; iteration <= maxIterations; ++iteration)
        {
            solveUMomentum();
            solveVMomentum();
            buildRhieChowFaceFluxes();
            solvePressureCorrection();
            correctPressure();
            correctFaceFluxes();
            correctVelocity();

            residual = calculateResidual();

            std::cout << "SIMPLE iteration " << iteration
                << " | Residual = " << residual << "\n";

            if (checkConvergence())
            {
                std::cout << "\nSIMPLE converged after " << iteration << " iterations.\n";
                return;
            }
        }

        std::cout << "\nSIMPLE reached maximum iterations.\n"
            << "Final residual = " << residual << "\n";
    }
}
