#include "solvers/SIMPLE.h"

#include "linearAlgebra/BiCGSTAB.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace CFD
{

    namespace
    {
        constexpr double SMALL = 1.0e-14;

        constexpr double MAX_VELOCITY =
            10.0;

        constexpr double MAX_PRESSURE_CORRECTION =
            1.0;


        double getMatrixDiagonal(
            const SparseMatrix& matrix,
            std::size_t P)
        {
            for (std::size_t k = matrix.rowPtr()[P];
                k < matrix.rowPtr()[P + 1];
                ++k)
            {
                if (matrix.colIndices()[k] == P)
                {
                    return matrix.values()[k];
                }
            }

            return 0.0;
        }


        double requirePositiveDiagonal(
            const SparseMatrix& matrix,
            std::size_t P,
            const char* equation)
        {
            const double value =
                getMatrixDiagonal(matrix, P);

            if (!std::isfinite(value) ||
                value <= SMALL)
            {
                throw std::runtime_error(
                    std::string("SIMPLE: invalid ") +
                    equation +
                    " diagonal at cell " +
                    std::to_string(P));
            }

            return value;
        }


        double limitPressureCorrection(
            double value)
        {
            if (!std::isfinite(value))
            {
                throw std::runtime_error(
                    "SIMPLE: non-finite pressure correction");
            }

            return std::max(
                -MAX_PRESSURE_CORRECTION,
                std::min(
                    MAX_PRESSURE_CORRECTION,
                    value));
        }
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
        westBC(westBC)
    {
    }


    void SIMPLE::setPressureRelaxation(double value)
    {
        if (!std::isfinite(value) ||
            value <= 0.0 ||
            value > 1.0)
        {
            throw std::runtime_error(
                "SIMPLE: pressure relaxation must be between 0 and 1");
        }

        relaxationPressure = value;
    }


    void SIMPLE::setVelocityRelaxation(double value)
    {
        if (!std::isfinite(value) ||
            value <= 0.0 ||
            value > 1.0)
        {
            throw std::runtime_error(
                "SIMPLE: velocity relaxation must be between 0 and 1");
        }

        relaxationVelocity = value;
    }


    void SIMPLE::setConvergenceTolerance(double value)
    {
        if (!std::isfinite(value) ||
            value <= 0.0)
        {
            throw std::runtime_error(
                "SIMPLE: convergence tolerance must be positive");
        }

        convergenceTolerance = value;
    }


    void SIMPLE::setMaxIterations(std::size_t value)
    {
        if (value == 0)
        {
            throw std::runtime_error(
                "SIMPLE: maximum iterations must be greater than zero");
        }

        maxIterations = value;
    }


    void SIMPLE::setDensity(double value)
    {
        if (!std::isfinite(value) ||
            value <= 0.0)
        {
            throw std::runtime_error(
                "SIMPLE: density must be positive");
        }

        rho = value;
    }


    double SIMPLE::getPressureRelaxation() const
    {
        return relaxationPressure;
    }


    double SIMPLE::getVelocityRelaxation() const
    {
        return relaxationVelocity;
    }


    double SIMPLE::getConvergenceTolerance() const
    {
        return convergenceTolerance;
    }


    double SIMPLE::getDensity() const
    {
        return rho;
    }


    std::size_t SIMPLE::getMaxIterations() const
    {
        return maxIterations;
    }


    std::size_t SIMPLE::getIteration() const
    {
        return iteration;
    }


    double SIMPLE::getResidual() const
    {
        return residual;
    }


    void SIMPLE::applyBoundaryConditions()
    {
        const std::size_t n =
            mesh.numberOfCells();


        /*
            Velocity boundary conditions are imposed on faces,
            not directly on cell-centred velocity values.

            Pressure boundary conditions are stored at the
            boundary cells because pressure is cell centred.
        */

        for (std::size_t P = 0;
            P < n;
            ++P)
        {
            const auto& cell =
                mesh.getCell(
                    static_cast<int>(P));


            if (cell.westBoundary &&
                westBC.hasPressure())
            {
                fields.pressure[P] =
                    westBC.getPressure();
            }


            if (cell.eastBoundary &&
                eastBC.hasPressure())
            {
                fields.pressure[P] =
                    eastBC.getPressure();
            }


            if (cell.northBoundary &&
                northBC.hasPressure())
            {
                fields.pressure[P] =
                    northBC.getPressure();
            }


            if (cell.southBoundary &&
                southBC.hasPressure())
            {
                fields.pressure[P] =
                    southBC.getPressure();
            }
        }
    }


    void SIMPLE::calculateFaceFluxes()
    {
        const std::size_t n =
            mesh.numberOfCells();

        fluxEast.assign(n, 0.0);
        fluxWest.assign(n, 0.0);
        fluxNorth.assign(n, 0.0);
        fluxSouth.assign(n, 0.0);


        const double Ae =
            mesh.eastWestFaceArea();

        const double An =
            mesh.northSouthFaceArea();


        for (std::size_t j = 0;
            j < mesh.getNy();
            ++j)
        {
            for (std::size_t i = 0;
                i < mesh.getNx();
                ++i)
            {
                const std::size_t P =
                    mesh.cellIndex(i, j);

                const auto& cell =
                    mesh.getCell(
                        static_cast<int>(P));


                /*
                    EAST FACE
                */

                if (cell.east != -1)
                {
                    const std::size_t E =
                        static_cast<std::size_t>(
                            cell.east);

                    const double uFace =
                        0.5 *
                        (
                            fields.velocity.getx()[P] +
                            fields.velocity.getx()[E]
                            );

                    fluxEast[P] =
                        rho *
                        uFace *
                        Ae;
                }
                else if (
                    eastBC.getType() ==
                    BoundaryType::Wall)
                {
                    fluxEast[P] = 0.0;
                }
                else if (
                    eastBC.getType() ==
                    BoundaryType::Inlet)
                {
                    fluxEast[P] =
                        rho *
                        eastBC.getU() *
                        Ae;
                }
                else if (
                    eastBC.getType() ==
                    BoundaryType::Outlet)
                {
                    fluxEast[P] =
                        rho *
                        fields.velocity.getx()[P] *
                        Ae;
                }


                /*
                    WEST FACE
                */

                if (cell.west != -1)
                {
                    const std::size_t W =
                        static_cast<std::size_t>(
                            cell.west);

                    const double uFace =
                        0.5 *
                        (
                            fields.velocity.getx()[P] +
                            fields.velocity.getx()[W]
                            );

                    fluxWest[P] =
                        rho *
                        uFace *
                        Ae;
                }
                else if (
                    westBC.getType() ==
                    BoundaryType::Wall)
                {
                    fluxWest[P] = 0.0;
                }
                else if (
                    westBC.getType() ==
                    BoundaryType::Inlet)
                {
                    fluxWest[P] =
                        rho *
                        westBC.getU() *
                        Ae;
                }
                else if (
                    westBC.getType() ==
                    BoundaryType::Outlet)
                {
                    fluxWest[P] =
                        rho *
                        fields.velocity.getx()[P] *
                        Ae;
                }


                /*
                    NORTH FACE
                */

                if (cell.north != -1)
                {
                    const std::size_t N =
                        static_cast<std::size_t>(
                            cell.north);

                    const double vFace =
                        0.5 *
                        (
                            fields.velocity.gety()[P] +
                            fields.velocity.gety()[N]
                            );

                    fluxNorth[P] =
                        rho *
                        vFace *
                        An;
                }
                else if (
                    northBC.getType() ==
                    BoundaryType::Wall)
                {
                    fluxNorth[P] = 0.0;
                }
                else if (
                    northBC.getType() ==
                    BoundaryType::Inlet)
                {
                    fluxNorth[P] =
                        rho *
                        northBC.getV() *
                        An;
                }
                else if (
                    northBC.getType() ==
                    BoundaryType::Outlet)
                {
                    fluxNorth[P] =
                        rho *
                        fields.velocity.gety()[P] *
                        An;
                }


                /*
                    SOUTH FACE
                */

                if (cell.south != -1)
                {
                    const std::size_t S =
                        static_cast<std::size_t>(
                            cell.south);

                    const double vFace =
                        0.5 *
                        (
                            fields.velocity.gety()[P] +
                            fields.velocity.gety()[S]
                            );

                    fluxSouth[P] =
                        rho *
                        vFace *
                        An;
                }
                else if (
                    southBC.getType() ==
                    BoundaryType::Wall)
                {
                    fluxSouth[P] = 0.0;
                }
                else if (
                    southBC.getType() ==
                    BoundaryType::Inlet)
                {
                    fluxSouth[P] =
                        rho *
                        southBC.getV() *
                        An;
                }
                else if (
                    southBC.getType() ==
                    BoundaryType::Outlet)
                {
                    fluxSouth[P] =
                        rho *
                        fields.velocity.gety()[P] *
                        An;
                }


                if (!std::isfinite(fluxEast[P]) ||
                    !std::isfinite(fluxWest[P]) ||
                    !std::isfinite(fluxNorth[P]) ||
                    !std::isfinite(fluxSouth[P]))
                {
                    throw std::runtime_error(
                        "SIMPLE: non-finite face flux");
                }
            }
        }
    }


    void SIMPLE::assembleUMomentum()
    {
        const std::size_t n =
            mesh.numberOfCells();

        const double dx =
            mesh.getDx();

        const double dy =
            mesh.getDy();

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


        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> values;

        rows.reserve(5 * n);
        cols.reserve(5 * n);
        values.reserve(5 * n);

        uRHS =
            Vector(n, 0.0);


        for (std::size_t j = 0;
            j < mesh.getNy();
            ++j)
        {
            for (std::size_t i = 0;
                i < mesh.getNx();
                ++i)
            {
                const std::size_t P =
                    mesh.cellIndex(i, j);

                const auto& cell =
                    mesh.getCell(
                        static_cast<int>(P));


                const double Fe =
                    fluxEast[P];

                const double Fw =
                    fluxWest[P];

                const double Fn =
                    fluxNorth[P];

                const double Fs =
                    fluxSouth[P];


                double aE =
                    De +
                    std::max(-Fe, 0.0);

                double aW =
                    Dw +
                    std::max(Fw, 0.0);

                double aN =
                    Dn +
                    std::max(-Fn, 0.0);

                double aS =
                    Ds +
                    std::max(Fs, 0.0);


                double source =
                    0.0;


                /*
                    Physical boundary faces.

                    A prescribed velocity is represented by
                    a source contribution.

                    The boundary is NOT turned into a fixed
                    cell-centre velocity.
                */

                if (cell.east == -1)
                {
                    aE =
                        2.0 * De +
                        std::max(-Fe, 0.0);

                    if (eastBC.getType() ==
                        BoundaryType::Wall ||
                        eastBC.getType() ==
                        BoundaryType::Inlet)
                    {
                        source +=
                            aE *
                            eastBC.getU();
                    }
                    else if (
                        eastBC.getType() ==
                        BoundaryType::Outlet)
                    {
                        aE = 0.0;
                    }
                }


                if (cell.west == -1)
                {
                    aW =
                        2.0 * Dw +
                        std::max(Fw, 0.0);

                    if (westBC.getType() ==
                        BoundaryType::Wall ||
                        westBC.getType() ==
                        BoundaryType::Inlet)
                    {
                        source +=
                            aW *
                            westBC.getU();
                    }
                    else if (
                        westBC.getType() ==
                        BoundaryType::Outlet)
                    {
                        aW = 0.0;
                    }
                }


                if (cell.north == -1)
                {
                    aN =
                        2.0 * Dn +
                        std::max(-Fn, 0.0);

                    if (northBC.getType() ==
                        BoundaryType::Wall ||
                        northBC.getType() ==
                        BoundaryType::Inlet)
                    {
                        source +=
                            aN *
                            northBC.getU();
                    }
                    else if (
                        northBC.getType() ==
                        BoundaryType::Outlet)
                    {
                        aN = 0.0;
                    }
                }


                if (cell.south == -1)
                {
                    aS =
                        2.0 * Ds +
                        std::max(Fs, 0.0);

                    if (southBC.getType() ==
                        BoundaryType::Wall ||
                        southBC.getType() ==
                        BoundaryType::Inlet)
                    {
                        source +=
                            aS *
                            southBC.getU();
                    }
                    else if (
                        southBC.getType() ==
                        BoundaryType::Outlet)
                    {
                        aS = 0.0;
                    }
                }


                /*
                    Pressure gradient.

                        (pW - pE) Ae
                */

                if (cell.west != -1)
                {
                    source +=
                        fields.pressure[
                            static_cast<std::size_t>(
                                cell.west)] *
                        Ae;
                }

                if (cell.east != -1)
                {
                    source -=
                        fields.pressure[
                            static_cast<std::size_t>(
                                cell.east)] *
                        Ae;
                }


                /*
                    IMPORTANT:

                    Do not include the current continuity defect
                    in aP.

                    The previous formulation:

                        aP = sum(aNB) + imbalance

                    allowed a large continuity error to produce
                    a negative momentum diagonal.

                    The momentum matrix must remain positive
                    while SIMPLE drives the continuity error down.
                */

                const double aP =
                    aE +
                    aW +
                    aN +
                    aS;


                if (!std::isfinite(aP) ||
                    aP <= SMALL)
                {
                    throw std::runtime_error(
                        "SIMPLE: invalid U-momentum diagonal");
                }


                rows.push_back(P);
                cols.push_back(P);
                values.push_back(aP);


                if (cell.east != -1 &&
                    aE > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            cell.east));
                    values.push_back(-aE);
                }


                if (cell.west != -1 &&
                    aW > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            cell.west));
                    values.push_back(-aW);
                }


                if (cell.north != -1 &&
                    aN > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            cell.north));
                    values.push_back(-aN);
                }


                if (cell.south != -1 &&
                    aS > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            cell.south));
                    values.push_back(-aS);
                }


                uRHS[P] =
                    source;
            }
        }


        uMatrix =
            SparseMatrix(n, n);

        uMatrix.setFromTriplets(
            rows,
            cols,
            values);
    }


    void SIMPLE::solveUMomentum()
    {
        assembleUMomentum();


        const std::size_t n =
            mesh.numberOfCells();


        Vector uOld(n, 0.0);
        Vector uSolution(n, 0.0);


        for (std::size_t P = 0;
            P < n;
            ++P)
        {
            uOld[P] =
                fields.velocity.getx()[P];

            uSolution[P] =
                uOld[P];
        }


        BiCGSTAB solver(
            1.0e-8,
            10000);


        BiCGSTABResult result =
            solver.solve(
                uMatrix,
                uRHS,
                uSolution,
                false);


        if (!result.converged)
        {
            throw std::runtime_error(
                "SIMPLE: U momentum solver failed to converge");
        }


        for (std::size_t P = 0;
            P < n;
            ++P)
        {
            const double uNew =
                uOld[P] +
                relaxationVelocity *
                (
                    uSolution[P] -
                    uOld[P]
                    );


            if (!std::isfinite(uNew) ||
                std::abs(uNew) > MAX_VELOCITY)
            {
                throw std::runtime_error(
                    "SIMPLE: unstable U velocity detected");
            }


            fields.velocity.getx()[P] =
                uNew;
        }
    }


    void SIMPLE::assembleVMomentum()
    {
        const std::size_t n =
            mesh.numberOfCells();

        const double dx =
            mesh.getDx();

        const double dy =
            mesh.getDy();

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


        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> values;

        rows.reserve(5 * n);
        cols.reserve(5 * n);
        values.reserve(5 * n);

        vRHS =
            Vector(n, 0.0);


        for (std::size_t j = 0;
            j < mesh.getNy();
            ++j)
        {
            for (std::size_t i = 0;
                i < mesh.getNx();
                ++i)
            {
                const std::size_t P =
                    mesh.cellIndex(i, j);

                const auto& cell =
                    mesh.getCell(
                        static_cast<int>(P));


                const double Fe =
                    fluxEast[P];

                const double Fw =
                    fluxWest[P];

                const double Fn =
                    fluxNorth[P];

                const double Fs =
                    fluxSouth[P];


                double aE =
                    De +
                    std::max(-Fe, 0.0);

                double aW =
                    Dw +
                    std::max(Fw, 0.0);

                double aN =
                    Dn +
                    std::max(-Fn, 0.0);

                double aS =
                    Ds +
                    std::max(Fs, 0.0);


                double source =
                    0.0;


                if (cell.east == -1)
                {
                    aE =
                        2.0 * De +
                        std::max(-Fe, 0.0);

                    if (eastBC.getType() ==
                        BoundaryType::Wall ||
                        eastBC.getType() ==
                        BoundaryType::Inlet)
                    {
                        source +=
                            aE *
                            eastBC.getV();
                    }
                    else if (
                        eastBC.getType() ==
                        BoundaryType::Outlet)
                    {
                        aE = 0.0;
                    }
                }


                if (cell.west == -1)
                {
                    aW =
                        2.0 * Dw +
                        std::max(Fw, 0.0);

                    if (westBC.getType() ==
                        BoundaryType::Wall ||
                        westBC.getType() ==
                        BoundaryType::Inlet)
                    {
                        source +=
                            aW *
                            westBC.getV();
                    }
                    else if (
                        westBC.getType() ==
                        BoundaryType::Outlet)
                    {
                        aW = 0.0;
                    }
                }


                if (cell.north == -1)
                {
                    aN =
                        2.0 * Dn +
                        std::max(-Fn, 0.0);

                    if (northBC.getType() ==
                        BoundaryType::Wall ||
                        northBC.getType() ==
                        BoundaryType::Inlet)
                    {
                        source +=
                            aN *
                            northBC.getV();
                    }
                    else if (
                        northBC.getType() ==
                        BoundaryType::Outlet)
                    {
                        aN = 0.0;
                    }
                }


                if (cell.south == -1)
                {
                    aS =
                        2.0 * Ds +
                        std::max(Fs, 0.0);

                    if (southBC.getType() ==
                        BoundaryType::Wall ||
                        southBC.getType() ==
                        BoundaryType::Inlet)
                    {
                        source +=
                            aS *
                            southBC.getV();
                    }
                    else if (
                        southBC.getType() ==
                        BoundaryType::Outlet)
                    {
                        aS = 0.0;
                    }
                }


                /*
                    Pressure gradient:

                        (pS - pN) An
                */

                if (cell.south != -1)
                {
                    source +=
                        fields.pressure[
                            static_cast<std::size_t>(
                                cell.south)] *
                        An;
                }

                if (cell.north != -1)
                {
                    source -=
                        fields.pressure[
                            static_cast<std::size_t>(
                                cell.north)] *
                        An;
                }


                /*
                    Keep the momentum diagonal positive.

                    The continuity equation is corrected by the
                    pressure-correction equation rather than by
                    allowing a temporary mass imbalance to make
                    aP negative.
                */

                const double aP =
                    aE +
                    aW +
                    aN +
                    aS;


                if (!std::isfinite(aP) ||
                    aP <= SMALL)
                {
                    throw std::runtime_error(
                        "SIMPLE: invalid V-momentum diagonal");
                }


                rows.push_back(P);
                cols.push_back(P);
                values.push_back(aP);


                if (cell.east != -1 &&
                    aE > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            cell.east));
                    values.push_back(-aE);
                }


                if (cell.west != -1 &&
                    aW > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            cell.west));
                    values.push_back(-aW);
                }


                if (cell.north != -1 &&
                    aN > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            cell.north));
                    values.push_back(-aN);
                }


                if (cell.south != -1 &&
                    aS > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            cell.south));
                    values.push_back(-aS);
                }


                vRHS[P] =
                    source;
            }
        }


        vMatrix =
            SparseMatrix(n, n);

        vMatrix.setFromTriplets(
            rows,
            cols,
            values);
    }


    void SIMPLE::solveVMomentum()
    {
        assembleVMomentum();


        const std::size_t n =
            mesh.numberOfCells();


        Vector vOld(n, 0.0);
        Vector vSolution(n, 0.0);


        for (std::size_t P = 0;
            P < n;
            ++P)
        {
            vOld[P] =
                fields.velocity.gety()[P];

            vSolution[P] =
                vOld[P];
        }


        BiCGSTAB solver(
            1.0e-8,
            10000);


        BiCGSTABResult result =
            solver.solve(
                vMatrix,
                vRHS,
                vSolution,
                false);


        if (!result.converged)
        {
            throw std::runtime_error(
                "SIMPLE: V momentum solver failed to converge");
        }


        for (std::size_t P = 0;
            P < n;
            ++P)
        {
            const double vNew =
                vOld[P] +
                relaxationVelocity *
                (
                    vSolution[P] -
                    vOld[P]
                    );


            if (!std::isfinite(vNew) ||
                std::abs(vNew) > MAX_VELOCITY)
            {
                throw std::runtime_error(
                    "SIMPLE: unstable V velocity detected");
            }


            fields.velocity.gety()[P] =
                vNew;
        }
    }


    void SIMPLE::assemblePressureCorrection()
    {
        const std::size_t n =
            mesh.numberOfCells();

        const double Ae =
            mesh.eastWestFaceArea();

        const double An =
            mesh.northSouthFaceArea();


        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> values;

        rows.reserve(5 * n);
        cols.reserve(5 * n);
        values.reserve(5 * n);


        pressureRHS =
            Vector(n, 0.0);


        /*
            SIMPLE d coefficients.

                dU = alpha_u Ae / aPu
                dV = alpha_v An / aPv
        */

        std::vector<double> dU(n, 0.0);
        std::vector<double> dV(n, 0.0);


        for (std::size_t P = 0;
            P < n;
            ++P)
        {
            const double aPu =
                requirePositiveDiagonal(
                    uMatrix,
                    P,
                    "U-momentum");

            const double aPv =
                requirePositiveDiagonal(
                    vMatrix,
                    P,
                    "V-momentum");


            dU[P] =
                relaxationVelocity *
                Ae /
                aPu;

            dV[P] =
                relaxationVelocity *
                An /
                aPv;


            if (!std::isfinite(dU[P]) ||
                !std::isfinite(dV[P]))
            {
                throw std::runtime_error(
                    "SIMPLE: invalid SIMPLE coefficient");
            }
        }


        /*
            Interior face coefficients.
        */

        auto uFaceCoefficient =
            [&](std::size_t P,
                std::size_t Q)
            -> double
            {
                const double dFace =
                    0.5 *
                    (
                        dU[P] +
                        dU[Q]
                        );

                return rho *
                    Ae *
                    dFace;
            };


        auto vFaceCoefficient =
            [&](std::size_t P,
                std::size_t Q)
            -> double
            {
                const double dFace =
                    0.5 *
                    (
                        dV[P] +
                        dV[Q]
                        );

                return rho *
                    An *
                    dFace;
            };


        for (std::size_t j = 0;
            j < mesh.getNy();
            ++j)
        {
            for (std::size_t i = 0;
                i < mesh.getNx();
                ++i)
            {
                const std::size_t P =
                    mesh.cellIndex(i, j);

                const auto& cell =
                    mesh.getCell(
                        static_cast<int>(P));


                /*
                    East outlet provides the pressure reference.
                */

                if (cell.east == -1 &&
                    eastBC.getType() ==
                    BoundaryType::Outlet &&
                    eastBC.hasPressure())
                {
                    rows.push_back(P);
                    cols.push_back(P);
                    values.push_back(1.0);

                    pressureRHS[P] =
                        0.0;

                    continue;
                }


                double aE =
                    0.0;

                double aW =
                    0.0;

                double aN =
                    0.0;

                double aS =
                    0.0;


                /*
                    EAST
                */

                if (cell.east != -1)
                {
                    const std::size_t E =
                        static_cast<std::size_t>(
                            cell.east);

                    aE =
                        uFaceCoefficient(
                            P,
                            E);
                }
                else if (
                    eastBC.getType() ==
                    BoundaryType::Outlet)
                {
                    aE =
                        rho *
                        Ae *
                        dU[P];
                }


                /*
                    WEST
                */

                if (cell.west != -1)
                {
                    const std::size_t W =
                        static_cast<std::size_t>(
                            cell.west);

                    aW =
                        uFaceCoefficient(
                            P,
                            W);
                }


                /*
                    NORTH
                */

                if (cell.north != -1)
                {
                    const std::size_t N =
                        static_cast<std::size_t>(
                            cell.north);

                    aN =
                        vFaceCoefficient(
                            P,
                            N);
                }
                else if (
                    northBC.getType() ==
                    BoundaryType::Outlet)
                {
                    aN =
                        rho *
                        An *
                        dV[P];
                }


                /*
                    SOUTH
                */

                if (cell.south != -1)
                {
                    const std::size_t S =
                        static_cast<std::size_t>(
                            cell.south);

                    aS =
                        vFaceCoefficient(
                            P,
                            S);
                }
                else if (
                    southBC.getType() ==
                    BoundaryType::Outlet)
                {
                    aS =
                        rho *
                        An *
                        dV[P];
                }


                const double imbalance =
                    fluxEast[P]
                    - fluxWest[P]
                    + fluxNorth[P]
                    - fluxSouth[P];


                const double aP =
                    aE +
                    aW +
                    aN +
                    aS;


                if (!std::isfinite(aP) ||
                    aP <= SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(P);
                    values.push_back(1.0);

                    pressureRHS[P] =
                        0.0;

                    continue;
                }


                /*
                    Pressure correction equation:

                        aP p'P
                        -
                        sum(aNB p'NB)
                        =
                        -continuityError
                */

                pressureRHS[P] =
                    -imbalance;


                rows.push_back(P);
                cols.push_back(P);
                values.push_back(aP);


                if (cell.east != -1 &&
                    aE > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            cell.east));
                    values.push_back(-aE);
                }


                if (cell.west != -1 &&
                    aW > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            cell.west));
                    values.push_back(-aW);
                }


                if (cell.north != -1 &&
                    aN > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            cell.north));
                    values.push_back(-aN);
                }


                if (cell.south != -1 &&
                    aS > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(
                        static_cast<std::size_t>(
                            cell.south));
                    values.push_back(-aS);
                }
            }
        }


        pressureMatrix =
            SparseMatrix(n, n);

        pressureMatrix.setFromTriplets(
            rows,
            cols,
            values);
    }


    void SIMPLE::solvePressureCorrection()
    {
        const std::size_t n =
            mesh.numberOfCells();


        pressureCorrection =
            Vector(n, 0.0);


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
                "SIMPLE: pressure correction solver failed to converge");
        }


        /*
            Limit pressure corrections before they can
            destabilise the velocity field.

            The physical pressure relaxation is still applied
            separately in correctPressure().
        */

        for (std::size_t P = 0;
            P < n;
            ++P)
        {
            pressureCorrection[P] =
                limitPressureCorrection(
                    pressureCorrection[P]);
        }
    }


    void SIMPLE::correctPressure()
    {
        const std::size_t n =
            mesh.numberOfCells();


        for (std::size_t P = 0;
            P < n;
            ++P)
        {
            fields.pressure[P] +=
                relaxationPressure *
                pressureCorrection[P];


            if (!std::isfinite(
                fields.pressure[P]))
            {
                throw std::runtime_error(
                    "SIMPLE: non-finite pressure");
            }
        }


        /*
            Reapply physical pressure boundary conditions.
        */

        for (std::size_t P = 0;
            P < n;
            ++P)
        {
            const auto& cell =
                mesh.getCell(
                    static_cast<int>(P));


            if (cell.westBoundary &&
                westBC.hasPressure())
            {
                fields.pressure[P] =
                    westBC.getPressure();
            }


            if (cell.eastBoundary &&
                eastBC.hasPressure())
            {
                fields.pressure[P] =
                    eastBC.getPressure();
            }


            if (cell.northBoundary &&
                northBC.hasPressure())
            {
                fields.pressure[P] =
                    northBC.getPressure();
            }


            if (cell.southBoundary &&
                southBC.hasPressure())
            {
                fields.pressure[P] =
                    southBC.getPressure();
            }
        }
    }


    void SIMPLE::correctVelocity()
    {
        const std::size_t n =
            mesh.numberOfCells();

        const double Ae =
            mesh.eastWestFaceArea();

        const double An =
            mesh.northSouthFaceArea();


        for (std::size_t P = 0;
            P < n;
            ++P)
        {
            const auto& cell =
                mesh.getCell(
                    static_cast<int>(P));


            const double aPu =
                requirePositiveDiagonal(
                    uMatrix,
                    P,
                    "U-momentum");

            const double aPv =
                requirePositiveDiagonal(
                    vMatrix,
                    P,
                    "V-momentum");


            const double dU =
                relaxationVelocity *
                Ae /
                aPu;

            const double dV =
                relaxationVelocity *
                An /
                aPv;


            /*
                U CORRECTION
            */

            double pPrimeW =
                pressureCorrection[P];

            double pPrimeE =
                pressureCorrection[P];


            /*
                WEST FACE

                Fixed velocity boundary:

                    p'_face = p'_P

                Outlet:

                    p'_face = 0
            */

            if (cell.west != -1)
            {
                pPrimeW =
                    pressureCorrection[
                        static_cast<std::size_t>(
                            cell.west)];
            }
            else if (
                westBC.getType() ==
                BoundaryType::Outlet)
            {
                pPrimeW = 0.0;
            }
            else
            {
                pPrimeW =
                    pressureCorrection[P];
            }


            /*
                EAST FACE
            */

            if (cell.east != -1)
            {
                pPrimeE =
                    pressureCorrection[
                        static_cast<std::size_t>(
                            cell.east)];
            }
            else if (
                eastBC.getType() ==
                BoundaryType::Outlet)
            {
                pPrimeE = 0.0;
            }
            else
            {
                pPrimeE =
                    pressureCorrection[P];
            }


            double uCorrection =
                dU *
                (
                    pPrimeW -
                    pPrimeE
                    );


            /*
                Limit the actual velocity correction.

                This prevents one bad pressure-correction
                iteration from destroying the solution.
            */

            if (!std::isfinite(uCorrection))
            {
                throw std::runtime_error(
                    "SIMPLE: non-finite U correction");
            }


            const double maxCorrection =
                0.25 *
                std::max(
                    1.0,
                    std::abs(
                        fields.velocity.getx()[P]));


            uCorrection =
                std::max(
                    -maxCorrection,
                    std::min(
                        maxCorrection,
                        uCorrection));


            fields.velocity.getx()[P] +=
                uCorrection;


            /*
                V CORRECTION
            */

            double pPrimeS =
                pressureCorrection[P];

            double pPrimeN =
                pressureCorrection[P];


            if (cell.south != -1)
            {
                pPrimeS =
                    pressureCorrection[
                        static_cast<std::size_t>(
                            cell.south)];
            }
            else if (
                southBC.getType() ==
                BoundaryType::Outlet)
            {
                pPrimeS = 0.0;
            }
            else
            {
                pPrimeS =
                    pressureCorrection[P];
            }


            if (cell.north != -1)
            {
                pPrimeN =
                    pressureCorrection[
                        static_cast<std::size_t>(
                            cell.north)];
            }
            else if (
                northBC.getType() ==
                BoundaryType::Outlet)
            {
                pPrimeN = 0.0;
            }
            else
            {
                pPrimeN =
                    pressureCorrection[P];
            }


            double vCorrection =
                dV *
                (
                    pPrimeS -
                    pPrimeN
                    );


            if (!std::isfinite(vCorrection))
            {
                throw std::runtime_error(
                    "SIMPLE: non-finite V correction");
            }


            const double maxVCorrection =
                0.25 *
                std::max(
                    1.0,
                    std::abs(
                        fields.velocity.gety()[P]));


            vCorrection =
                std::max(
                    -maxVCorrection,
                    std::min(
                        maxVCorrection,
                        vCorrection));


            fields.velocity.gety()[P] +=
                vCorrection;


            if (!std::isfinite(
                fields.velocity.getx()[P]) ||
                !std::isfinite(
                    fields.velocity.gety()[P]))
            {
                throw std::runtime_error(
                    "SIMPLE: non-finite velocity after correction");
            }


            if (std::abs(
                fields.velocity.getx()[P])
            > MAX_VELOCITY ||
                std::abs(
                    fields.velocity.gety()[P])
            > MAX_VELOCITY)
            {
                throw std::runtime_error(
                    "SIMPLE: velocity became unstable during correction");
            }
        }
    }


    double SIMPLE::calculateResidual()
    {
        const std::size_t n =
            mesh.numberOfCells();


        double totalResidual =
            0.0;


        for (std::size_t P = 0;
            P < n;
            ++P)
        {
            const double imbalance =
                fluxEast[P]
                - fluxWest[P]
                + fluxNorth[P]
                - fluxSouth[P];


            if (!std::isfinite(imbalance))
            {
                return std::numeric_limits<double>::infinity();
            }


            totalResidual +=
                std::abs(imbalance);
        }


        return totalResidual;
    }


    bool SIMPLE::checkConvergence()
    {
        residual =
            calculateResidual();


        return residual <
            convergenceTolerance;
    }


    void SIMPLE::solve()
    {
        iteration = 0;

        residual =
            std::numeric_limits<double>::infinity();


        for (iteration = 1;
            iteration <= maxIterations;
            ++iteration)
        {
            /*
                1. Physical pressure boundary conditions
            */

            applyBoundaryConditions();


            /*
                2. Current face mass fluxes
            */

            calculateFaceFluxes();


            /*
                3. U momentum
            */

            solveUMomentum();


            /*
                4. V momentum
            */

            solveVMomentum();


            /*
                5. Recalculate fluxes using momentum solution
            */

            calculateFaceFluxes();


            /*
                6. Pressure correction equation
            */

            assemblePressureCorrection();


            /*
                7. Pressure correction solution
            */

            solvePressureCorrection();


            /*
                8. Pressure update
            */

            correctPressure();


            /*
                9. Velocity correction

                The relaxation factor is already included in
                dU and dV.
            */

            correctVelocity();


            /*
                10. Recalculate physical fluxes
            */

            calculateFaceFluxes();


            /*
                11. Continuity residual
            */

            residual =
                calculateResidual();


            std::cout
                << "SIMPLE iteration "
                << iteration
                << " | Residual = "
                << residual
                << "\n";


            /*
                12. Convergence
            */

            if (residual <
                convergenceTolerance)
            {
                std::cout
                    << "\nSIMPLE converged after "
                    << iteration
                    << " iterations.\n";

                return;
            }
        }


        std::cout
            << "\nSIMPLE reached maximum iterations."
            << "\nFinal residual = "
            << residual
            << "\n";
    }

}