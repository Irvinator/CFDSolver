
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
            const BoundaryCondition& westBC)
        {
            if (cell.north == -1 &&
                northBC.getType() == BoundaryType::Wall &&
                northBC.hasU())
            {
                return true;
            }

            if (cell.south == -1 &&
                southBC.getType() == BoundaryType::Wall &&
                southBC.hasU())
            {
                return true;
            }

            if (cell.west == -1 &&
                westBC.getType() == BoundaryType::Inlet &&
                westBC.hasU())
            {
                return true;
            }

            return false;
        }


        bool isFixedV(
            const Mesh::Cell& cell,
            const BoundaryCondition& northBC,
            const BoundaryCondition& southBC,
            const BoundaryCondition& westBC)
        {
            if (cell.north == -1 &&
                northBC.getType() == BoundaryType::Wall &&
                northBC.hasV())
            {
                return true;
            }

            if (cell.south == -1 &&
                southBC.getType() == BoundaryType::Wall &&
                southBC.hasV())
            {
                return true;
            }

            if (cell.west == -1 &&
                westBC.getType() == BoundaryType::Inlet &&
                westBC.hasV())
            {
                return true;
            }

            return false;
        }


        double safeDiagonal(double value)
        {
            if (!std::isfinite(value) ||
                std::abs(value) < SMALL)
            {
                return SMALL;
            }

            return value;
        }
    }


    // ================================================================
    // CONSTRUCTOR
    // ================================================================

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


    // ================================================================
    // SETTERS
    // ================================================================

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


    // ================================================================
    // GETTERS
    // ================================================================

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


    // ================================================================
    // SOLVER STATE GETTERS
    // ================================================================

    std::size_t SIMPLE::getIteration() const
    {
        return iteration;
    }


    double SIMPLE::getResidual() const
    {
        return residual;
    }


    // ================================================================
    // APPLY BOUNDARY CONDITIONS
    // ================================================================

    void SIMPLE::applyBoundaryConditions()
    {
        const std::size_t n =
            mesh.numberOfCells();

        for (std::size_t P = 0;
            P < n;
            ++P)
        {
            const auto& cell =
                mesh.getCell(static_cast<int>(P));


            // --------------------------------------------------------
            // WEST
            // --------------------------------------------------------

            if (cell.westBoundary)
            {
                if (westBC.hasU())
                {
                    fields.velocity.getx()[P] =
                        westBC.getU();
                }

                if (westBC.hasV())
                {
                    fields.velocity.gety()[P] =
                        westBC.getV();
                }

                if (westBC.hasPressure())
                {
                    fields.pressure[P] =
                        westBC.getPressure();
                }
            }


            // --------------------------------------------------------
            // EAST
            // --------------------------------------------------------

            if (cell.eastBoundary)
            {
                if (eastBC.hasU())
                {
                    fields.velocity.getx()[P] =
                        eastBC.getU();
                }

                if (eastBC.hasV())
                {
                    fields.velocity.gety()[P] =
                        eastBC.getV();
                }

                if (eastBC.hasPressure())
                {
                    fields.pressure[P] =
                        eastBC.getPressure();
                }
            }


            // --------------------------------------------------------
            // NORTH
            //
            // Wall takes precedence over inlet/outlet at a corner.
            // --------------------------------------------------------

            if (cell.northBoundary)
            {
                if (northBC.hasU())
                {
                    fields.velocity.getx()[P] =
                        northBC.getU();
                }

                if (northBC.hasV())
                {
                    fields.velocity.gety()[P] =
                        northBC.getV();
                }

                if (northBC.hasPressure())
                {
                    fields.pressure[P] =
                        northBC.getPressure();
                }
            }


            // --------------------------------------------------------
            // SOUTH
            // --------------------------------------------------------

            if (cell.southBoundary)
            {
                if (southBC.hasU())
                {
                    fields.velocity.getx()[P] =
                        southBC.getU();
                }

                if (southBC.hasV())
                {
                    fields.velocity.gety()[P] =
                        southBC.getV();
                }

                if (southBC.hasPressure())
                {
                    fields.pressure[P] =
                        southBC.getPressure();
                }
            }
        }
    }


    // ================================================================
    // CALCULATE FACE FLUXES
    //
    // Sign convention:
    //
    // fluxEast  = positive velocity in +x direction
    // fluxWest  = positive velocity in +x direction
    // fluxNorth = positive velocity in +y direction
    // fluxSouth = positive velocity in +y direction
    //
    // Therefore:
    //
    // continuity residual =
    //
    //     Fe - Fw + Fn - Fs
    //
    // ================================================================

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

        std::cout << "\nINLET FACE FLUXES\n";


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
                    mesh.getCell(static_cast<int>(P));


                // ----------------------------------------------------
                // EAST FACE
                // ----------------------------------------------------

                if (cell.east != -1)
                {
                    const double uP =
                        fields.velocity.getx()[P];

                    const double uE =
                        fields.velocity.getx()[cell.east];

                    const double uFace =
                        0.5 * (uP + uE);

                    fluxEast[P] =
                        rho * uFace * Ae;
                }
                else
                {
                    if (eastBC.getType() == BoundaryType::Wall)
                    {
                        fluxEast[P] = 0.0;
                    }
                    else if (eastBC.getType() == BoundaryType::Inlet)
                    {
                        fluxEast[P] =
                            rho *
                            eastBC.getU() *
                            Ae;
                    }
                    else if (eastBC.getType() == BoundaryType::Outlet)
                    {
                        fluxEast[P] =
                            rho *
                            fields.velocity.getx()[P] *
                            Ae;
                    }
                }


                // ----------------------------------------------------
                // WEST FACE
                // ----------------------------------------------------

                if (cell.west != -1)
                {
                    const double uP =
                        fields.velocity.getx()[P];

                    const double uW =
                        fields.velocity.getx()[cell.west];

                    const double uFace =
                        0.5 * (uP + uW);

                    fluxWest[P] =
                        rho * uFace * Ae;
                }
                else
                {
                    if (westBC.getType() == BoundaryType::Wall)
                    {
                        fluxWest[P] = 0.0;
                    }
                    else if (westBC.getType() == BoundaryType::Inlet)
                    {
                        fluxWest[P] =
                            rho *
                            westBC.getU() *
                            Ae;
                    }
                    else if (westBC.getType() == BoundaryType::Outlet)
                    {
                        fluxWest[P] =
                            rho *
                            fields.velocity.getx()[P] *
                            Ae;
                    }
                }


                // ----------------------------------------------------
                // NORTH FACE
                // ----------------------------------------------------

                if (cell.north != -1)
                {
                    const double vP =
                        fields.velocity.gety()[P];

                    const double vN =
                        fields.velocity.gety()[cell.north];

                    const double vFace =
                        0.5 * (vP + vN);

                    fluxNorth[P] =
                        rho * vFace * An;
                }
                else
                {
                    if (northBC.getType() == BoundaryType::Wall)
                    {
                        fluxNorth[P] = 0.0;
                    }
                    else if (northBC.getType() == BoundaryType::Inlet)
                    {
                        fluxNorth[P] =
                            rho *
                            northBC.getV() *
                            An;
                    }
                    else if (northBC.getType() == BoundaryType::Outlet)
                    {
                        fluxNorth[P] =
                            rho *
                            fields.velocity.gety()[P] *
                            An;
                    }
                }


                // ----------------------------------------------------
                // SOUTH FACE
                // ----------------------------------------------------

                if (cell.south != -1)
                {
                    const double vP =
                        fields.velocity.gety()[P];

                    const double vS =
                        fields.velocity.gety()[cell.south];

                    const double vFace =
                        0.5 * (vP + vS);

                    fluxSouth[P] =
                        rho * vFace * An;
                }
                else
                {
                    if (southBC.getType() == BoundaryType::Wall)
                    {
                        fluxSouth[P] = 0.0;
                    }
                    else if (southBC.getType() == BoundaryType::Inlet)
                    {
                        fluxSouth[P] =
                            rho *
                            southBC.getV() *
                            An;
                    }
                    else if (southBC.getType() == BoundaryType::Outlet)
                    {
                        fluxSouth[P] =
                            rho *
                            fields.velocity.gety()[P] *
                            An;
                    }
                }
            }
        }
    }


    // ================================================================
    // ASSEMBLE U MOMENTUM
    //
    // First-order upwind:
    //
    // aE = De + max(-Fe, 0)
    // aW = Dw + max( Fw, 0)
    // aN = Dn + max(-Fn, 0)
    // aS = Ds + max( Fs, 0)
    //
    // The diagonal contains the net convection contribution:
    //
    // aP =
    //     aE + aW + aN + aS
    //     + Fe - Fw + Fn - Fs
    //
    // This is important because continuity is NOT satisfied during
    // the momentum solve.
    // ================================================================

    void SIMPLE::assembleUMomentum()
    {
        const std::size_t n =
            mesh.numberOfCells();

        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> values;

        uRHS =
            Vector(n, 0.0);

        rows.reserve(5 * n);
        cols.reserve(5 * n);
        values.reserve(5 * n);


        const double dx =
            mesh.getDx();

        const double dy =
            mesh.getDy();

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
                    mesh.getCell(static_cast<int>(P));


                // ----------------------------------------------------
                // Fixed velocity cells
                // ----------------------------------------------------

                if (isFixedU(
                    cell,
                    northBC,
                    southBC,
                    westBC))
                {
                    rows.push_back(P);
                    cols.push_back(P);
                    values.push_back(1.0);

                    double value = 0.0;

                    if (cell.west == -1 &&
                        westBC.getType() == BoundaryType::Inlet &&
                        westBC.hasU())
                    {
                        value = westBC.getU();
                    }

                    if (cell.north == -1 &&
                        northBC.getType() == BoundaryType::Wall &&
                        northBC.hasU())
                    {
                        value = northBC.getU();
                    }

                    if (cell.south == -1 &&
                        southBC.getType() == BoundaryType::Wall &&
                        southBC.hasU())
                    {
                        value = southBC.getU();
                    }

                    uRHS[P] = value;

                    continue;
                }


                // ----------------------------------------------------
                // Face mass fluxes
                // ----------------------------------------------------

                const double Fe = fluxEast[P];
                const double Fw = fluxWest[P];
                const double Fn = fluxNorth[P];
                const double Fs = fluxSouth[P];


                // ----------------------------------------------------
                // Diffusion
                // ----------------------------------------------------

                const double De =
                    mu * Ae / dx;

                const double Dw =
                    mu * Ae / dx;

                const double Dn =
                    mu * An / dy;

                const double Ds =
                    mu * An / dy;


                // ----------------------------------------------------
                // Upwind coefficients
                // ----------------------------------------------------

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


                // ----------------------------------------------------
                // Boundary faces do not couple to another unknown.
                // ----------------------------------------------------

                if (cell.east == -1)
                    aE = 0.0;

                if (cell.west == -1)
                    aW = 0.0;

                if (cell.north == -1)
                    aN = 0.0;

                if (cell.south == -1)
                    aS = 0.0;


                // ----------------------------------------------------
                // Pressure source
                //
                //     (pW - pE) Ae
                // ----------------------------------------------------

                double source = 0.0;

                if (cell.west != -1)
                {
                    source +=
                        fields.pressure[cell.west] * Ae;
                }

                if (cell.east != -1)
                {
                    source -=
                        fields.pressure[cell.east] * Ae;
                }


                // ----------------------------------------------------
                // Momentum diagonal
                //
                // Include net mass flux because continuity is not
                // exactly satisfied while momentum is being solved.
                // ----------------------------------------------------

                double aP =
                    aE +
                    aW +
                    aN +
                    aS +
                    (Fe - Fw + Fn - Fs);


                if (!std::isfinite(aP) ||
                    aP <= SMALL)
                {
                    throw std::runtime_error(
                        "SIMPLE: invalid U-momentum diagonal");
                }


                // ----------------------------------------------------
                // Matrix
                // ----------------------------------------------------

                rows.push_back(P);
                cols.push_back(P);
                values.push_back(aP);


                if (cell.east != -1 &&
                    aE > 0.0)
                {
                    rows.push_back(P);
                    cols.push_back(cell.east);
                    values.push_back(-aE);
                }

                if (cell.west != -1 &&
                    aW > 0.0)
                {
                    rows.push_back(P);
                    cols.push_back(cell.west);
                    values.push_back(-aW);
                }

                if (cell.north != -1 &&
                    aN > 0.0)
                {
                    rows.push_back(P);
                    cols.push_back(cell.north);
                    values.push_back(-aN);
                }

                if (cell.south != -1 &&
                    aS > 0.0)
                {
                    rows.push_back(P);
                    cols.push_back(cell.south);
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


    // ================================================================
    // SOLVE U MOMENTUM
    // ================================================================

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


        // ------------------------------------------------------------
        // Under-relax
        // ------------------------------------------------------------

        for (std::size_t P = 0;
            P < n;
            ++P)
        {
            if (isFixedU(
                mesh.getCell(static_cast<int>(P)),
                northBC,
                southBC,
                westBC))
            {
                continue;
            }

            fields.velocity.getx()[P] =
                uOld[P] +
                relaxationVelocity *
                (uSolution[P] - uOld[P]);
        }


        applyBoundaryConditions();
    }


    // ================================================================
    // ASSEMBLE V MOMENTUM
    // ================================================================

    void SIMPLE::assembleVMomentum()
    {
        const std::size_t n =
            mesh.numberOfCells();

        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> values;

        vRHS =
            Vector(n, 0.0);

        rows.reserve(5 * n);
        cols.reserve(5 * n);
        values.reserve(5 * n);


        const double dx =
            mesh.getDx();

        const double dy =
            mesh.getDy();

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
                    mesh.getCell(static_cast<int>(P));


                // ----------------------------------------------------
                // Fixed velocity cells
                // ----------------------------------------------------

                if (isFixedV(
                    cell,
                    northBC,
                    southBC,
                    westBC))
                {
                    rows.push_back(P);
                    cols.push_back(P);
                    values.push_back(1.0);

                    double value = 0.0;

                    if (cell.west == -1 &&
                        westBC.getType() == BoundaryType::Inlet &&
                        westBC.hasV())
                    {
                        value = westBC.getV();
                    }

                    if (cell.north == -1 &&
                        northBC.getType() == BoundaryType::Wall &&
                        northBC.hasV())
                    {
                        value = northBC.getV();
                    }

                    if (cell.south == -1 &&
                        southBC.getType() == BoundaryType::Wall &&
                        southBC.hasV())
                    {
                        value = southBC.getV();
                    }

                    vRHS[P] = value;

                    continue;
                }


                const double Fe = fluxEast[P];
                const double Fw = fluxWest[P];
                const double Fn = fluxNorth[P];
                const double Fs = fluxSouth[P];


                // ----------------------------------------------------
                // Diffusion
                // ----------------------------------------------------

                const double De =
                    mu * Ae / dx;

                const double Dw =
                    mu * Ae / dx;

                const double Dn =
                    mu * An / dy;

                const double Ds =
                    mu * An / dy;


                // ----------------------------------------------------
                // Upwind coefficients
                // ----------------------------------------------------

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


                if (cell.east == -1)
                    aE = 0.0;

                if (cell.west == -1)
                    aW = 0.0;

                if (cell.north == -1)
                    aN = 0.0;

                if (cell.south == -1)
                    aS = 0.0;


                // ----------------------------------------------------
                // Pressure source
                //
                //     (pS - pN) An
                // ----------------------------------------------------

                double source = 0.0;

                if (cell.south != -1)
                {
                    source +=
                        fields.pressure[cell.south] * An;
                }

                if (cell.north != -1)
                {
                    source -=
                        fields.pressure[cell.north] * An;
                }


                // ----------------------------------------------------
                // Diagonal
                // ----------------------------------------------------

                const double aP =
                    aE +
                    aW +
                    aN +
                    aS +
                    (Fe - Fw + Fn - Fs);


                if (!std::isfinite(aP) ||
                    aP <= SMALL)
                {
                    throw std::runtime_error(
                        "SIMPLE: invalid V-momentum diagonal");
                }


                // ----------------------------------------------------
                // Matrix
                // ----------------------------------------------------

                rows.push_back(P);
                cols.push_back(P);
                values.push_back(aP);


                if (cell.east != -1 &&
                    aE > 0.0)
                {
                    rows.push_back(P);
                    cols.push_back(cell.east);
                    values.push_back(-aE);
                }

                if (cell.west != -1 &&
                    aW > 0.0)
                {
                    rows.push_back(P);
                    cols.push_back(cell.west);
                    values.push_back(-aW);
                }

                if (cell.north != -1 &&
                    aN > 0.0)
                {
                    rows.push_back(P);
                    cols.push_back(cell.north);
                    values.push_back(-aN);
                }

                if (cell.south != -1 &&
                    aS > 0.0)
                {
                    rows.push_back(P);
                    cols.push_back(cell.south);
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


    // ================================================================
    // SOLVE V MOMENTUM
    // ================================================================

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
            if (isFixedV(
                mesh.getCell(static_cast<int>(P)),
                northBC,
                southBC,
                westBC))
            {
                continue;
            }

            fields.velocity.gety()[P] =
                vOld[P] +
                relaxationVelocity *
                (vSolution[P] - vOld[P]);
        }


        applyBoundaryConditions();
    }


    // ================================================================
    // GET MOMENTUM DIAGONAL
    //
    // Extract the actual diagonal coefficient from the matrix.
    // This deliberately does NOT use the identity coefficient of a
    // prescribed velocity cell as a physical SIMPLE coefficient.
    // ================================================================

    static double getDiagonal(
        const SparseMatrix& matrix,
        std::size_t P)
    {
        for (std::size_t k =
            matrix.rowPtr()[P];
            k <
            matrix.rowPtr()[P + 1];
            ++k)
        {
            if (matrix.colIndices()[k] == P)
            {
                return matrix.values()[k];
            }
        }

        return 0.0;
    }


    // ================================================================
    // ASSEMBLE PRESSURE CORRECTION
    //
    // The pressure correction is derived from the velocity response:
    //
    //     u' = d_u (p'_W - p'_E)
    //
    //     v' = d_v (p'_S - p'_N)
    //
    // where
    //
    //     d_u = Ae / aPu
    //     d_v = An / aPv
    //
    // For an internal face, the coefficient is based on the momentum
    // sensitivity of the cells adjacent to that face.
    //
    // For a face touching a fixed-velocity boundary, there is no
    // velocity response from that boundary side, so the pressure
    // correction is not allowed to alter that prescribed velocity.
    //
    // Walls therefore have zero pressure-correction flux.
    //
    // At an outlet, the outlet pressure correction is fixed to zero,
    // while the outlet velocity is allowed to respond.
    // ================================================================

    void SIMPLE::assemblePressureCorrection()
    {
        const std::size_t n =
            mesh.numberOfCells();

        std::vector<std::size_t> rows;
        std::vector<std::size_t> cols;
        std::vector<double> values;

        pressureRHS =
            Vector(n, 0.0);

        rows.reserve(5 * n);
        cols.reserve(5 * n);
        values.reserve(5 * n);


        const double Ae =
            mesh.eastWestFaceArea();

        const double An =
            mesh.northSouthFaceArea();


        // ------------------------------------------------------------
        // Momentum diagonals
        // ------------------------------------------------------------

        std::vector<double> aPu(n, SMALL);
        std::vector<double> aPv(n, SMALL);

        std::vector<bool> fixedU(n, false);
        std::vector<bool> fixedV(n, false);


        for (std::size_t P = 0;
            P < n;
            ++P)
        {
            const auto& cell =
                mesh.getCell(static_cast<int>(P));

            fixedU[P] =
                isFixedU(
                    cell,
                    northBC,
                    southBC,
                    westBC);

            fixedV[P] =
                isFixedV(
                    cell,
                    northBC,
                    southBC,
                    westBC);


            if (!fixedU[P])
            {
                aPu[P] =
                    safeDiagonal(
                        getDiagonal(uMatrix, P));
            }

            if (!fixedV[P])
            {
                aPv[P] =
                    safeDiagonal(
                        getDiagonal(vMatrix, P));
            }
        }


        // ------------------------------------------------------------
        // Face coefficient helpers
        // ------------------------------------------------------------

        auto uFaceCoefficient =
            [&](std::size_t P,
                std::size_t Q) -> double
            {
                if (fixedU[P] &&
                    fixedU[Q])
                {
                    return 0.0;
                }

                if (fixedU[P])
                {
                    return rho * Ae * Ae / aPu[Q];
                }

                if (fixedU[Q])
                {
                    return rho * Ae * Ae / aPu[P];
                }

                // Harmonic-style sensitivity combination.
                //
                // The two momentum equations provide two velocity
                // responses. Combining them through the sum of inverse
                // diagonals gives a symmetric face coefficient.
                const double sensitivity =
                    0.5 *
                    (
                        1.0 / aPu[P] +
                        1.0 / aPu[Q]
                        );

                return rho *
                    Ae *
                    Ae *
                    sensitivity;
            };


        auto vFaceCoefficient =
            [&](std::size_t P,
                std::size_t Q) -> double
            {
                if (fixedV[P] &&
                    fixedV[Q])
                {
                    return 0.0;
                }

                if (fixedV[P])
                {
                    return rho * An * An / aPv[Q];
                }

                if (fixedV[Q])
                {
                    return rho * An * An / aPv[P];
                }

                const double sensitivity =
                    0.5 *
                    (
                        1.0 / aPv[P] +
                        1.0 / aPv[Q]
                        );

                return rho *
                    An *
                    An *
                    sensitivity;
            };


        // ------------------------------------------------------------
        // Assemble one pressure equation per cell.
        // ------------------------------------------------------------

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
                    mesh.getCell(static_cast<int>(P));


                // ----------------------------------------------------
                // Pressure reference at outlet
                // ----------------------------------------------------

                if (cell.east == -1 &&
                    eastBC.getType() == BoundaryType::Outlet &&
                    eastBC.hasPressure())
                {
                    rows.push_back(P);
                    cols.push_back(P);
                    values.push_back(1.0);

                    pressureRHS[P] =
                        0.0;

                    continue;
                }


                // ----------------------------------------------------
                // Pressure-correction coefficients
                // ----------------------------------------------------

                double aE = 0.0;
                double aW = 0.0;
                double aN = 0.0;
                double aS = 0.0;


                // ----------------------------------------------------
                // EAST
                // ----------------------------------------------------

                if (cell.east != -1)
                {
                    const std::size_t E =
                        static_cast<std::size_t>(
                            cell.east);

                    aE =
                        uFaceCoefficient(P, E);
                }
                else if (
                    eastBC.getType() ==
                    BoundaryType::Outlet)
                {
                    // Outlet velocity is allowed to respond to
                    // pressure correction.
                    if (!fixedU[P])
                    {
                        aE =
                            rho *
                            Ae *
                            Ae /
                            aPu[P];
                    }
                }


                // ----------------------------------------------------
                // WEST
                // ----------------------------------------------------

                if (cell.west != -1)
                {
                    const std::size_t W =
                        static_cast<std::size_t>(
                            cell.west);

                    aW =
                        uFaceCoefficient(P, W);
                }
                else
                {
                    // Fixed inlet and wall velocities cannot be
                    // pressure-corrected.
                    aW = 0.0;
                }


                // ----------------------------------------------------
                // NORTH
                // ----------------------------------------------------

                if (cell.north != -1)
                {
                    const std::size_t N =
                        static_cast<std::size_t>(
                            cell.north);

                    aN =
                        vFaceCoefficient(P, N);
                }
                else
                {
                    // Wall velocity is fixed.
                    aN = 0.0;
                }


                // ----------------------------------------------------
                // SOUTH
                // ----------------------------------------------------

                if (cell.south != -1)
                {
                    const std::size_t S =
                        static_cast<std::size_t>(
                            cell.south);

                    aS =
                        vFaceCoefficient(P, S);
                }
                else
                {
                    // Wall velocity is fixed.
                    aS = 0.0;
                }


                // ----------------------------------------------------
                // Current continuity imbalance
                // ----------------------------------------------------

                const double imbalance =
                    fluxEast[P]
                    - fluxWest[P]
                    + fluxNorth[P]
                    - fluxSouth[P];


                // ----------------------------------------------------
                // Pressure correction equation
                //
                //     aP p'P
                //       - aE p'E
                //       - aW p'W
                //       - aN p'N
                //       - aS p'S
                //
                //       = -R_P
                // ----------------------------------------------------

                const double aP =
                    aE +
                    aW +
                    aN +
                    aS;


                if (!std::isfinite(aP) ||
                    aP < SMALL)
                {
                    // A cell with no pressure coupling can only occur
                    // when all its surrounding velocities are fixed.
                    // Such a pressure correction is irrelevant, so
                    // use a reference-style equation.
                    rows.push_back(P);
                    cols.push_back(P);
                    values.push_back(1.0);

                    pressureRHS[P] =
                        0.0;

                    continue;
                }


                pressureRHS[P] =
                    -imbalance;


                rows.push_back(P);
                cols.push_back(P);
                values.push_back(aP);


                if (cell.east != -1 &&
                    aE > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(cell.east);
                    values.push_back(-aE);
                }

                if (cell.west != -1 &&
                    aW > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(cell.west);
                    values.push_back(-aW);
                }

                if (cell.north != -1 &&
                    aN > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(cell.north);
                    values.push_back(-aN);
                }

                if (cell.south != -1 &&
                    aS > SMALL)
                {
                    rows.push_back(P);
                    cols.push_back(cell.south);
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


    // ================================================================
    // SOLVE PRESSURE CORRECTION
    // ================================================================

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
    }


    // ================================================================
    // CORRECT PRESSURE
    // ================================================================

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
        }


        // ------------------------------------------------------------
        // Reapply physical pressure boundary conditions.
        // ------------------------------------------------------------

        for (std::size_t P = 0;
            P < n;
            ++P)
        {
            const auto& cell =
                mesh.getCell(static_cast<int>(P));


            if (cell.eastBoundary &&
                eastBC.hasPressure())
            {
                fields.pressure[P] =
                    eastBC.getPressure();
            }

            if (cell.westBoundary &&
                westBC.hasPressure())
            {
                fields.pressure[P] =
                    westBC.getPressure();
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


    // ================================================================
    // CORRECT VELOCITY
    //
    // The correction is based on the SAME momentum diagonals used
    // when constructing the pressure-correction equation.
    //
    // U:
    //
    //     u' = Ae/aPu * (p'_W - p'_E)
    //
    // V:
    //
    //     v' = An/aPv * (p'_S - p'_N)
    //
    // Fixed velocity cells are never corrected.
    // ================================================================

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
                mesh.getCell(static_cast<int>(P));


            // --------------------------------------------------------
            // U correction
            // --------------------------------------------------------

            if (!isFixedU(
                cell,
                northBC,
                southBC,
                westBC))
            {
                const double aPu =
                    safeDiagonal(
                        getDiagonal(
                            uMatrix,
                            P));


                double pPrimeW =
                    pressureCorrection[P];

                double pPrimeE =
                    pressureCorrection[P];


                bool hasWest =
                    cell.west != -1;

                bool hasEast =
                    cell.east != -1;


                if (hasWest)
                {
                    pPrimeW =
                        pressureCorrection[
                            static_cast<std::size_t>(
                                cell.west)];
                }

                if (hasEast)
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
                    // Outlet pressure correction is zero because the
                    // outlet pressure is the reference.
                    pPrimeE = 0.0;
                }


                // For a west boundary with prescribed inlet velocity,
                // there is no velocity correction through that face.
                //
                // For an interior cell, use the pressure gradient.
                if (hasWest)
                {
                    const double correction =
                        (Ae / aPu) *
                        (pPrimeW - pPrimeE);

                    fields.velocity.getx()[P] +=
                        relaxationVelocity *
                        correction;
                }
                else if (hasEast &&
                    eastBC.getType() ==
                    BoundaryType::Outlet)
                {
                    const double correction =
                        -(Ae / aPu) *
                        pPrimeE;

                    fields.velocity.getx()[P] +=
                        relaxationVelocity *
                        correction;
                }
            }


            // --------------------------------------------------------
            // V correction
            // --------------------------------------------------------

            if (!isFixedV(
                cell,
                northBC,
                southBC,
                westBC))
            {
                const double aPv =
                    safeDiagonal(
                        getDiagonal(
                            vMatrix,
                            P));


                double pPrimeS =
                    pressureCorrection[P];

                double pPrimeN =
                    pressureCorrection[P];


                bool hasSouth =
                    cell.south != -1;

                bool hasNorth =
                    cell.north != -1;


                if (hasSouth)
                {
                    pPrimeS =
                        pressureCorrection[
                            static_cast<std::size_t>(
                                cell.south)];
                }

                if (hasNorth)
                {
                    pPrimeN =
                        pressureCorrection[
                            static_cast<std::size_t>(
                                cell.north)];
                }


                // Interior north/south correction.
                //
                // At a wall, the corresponding boundary velocity is
                // fixed, so there is no velocity correction through
                // the wall.
                if (hasSouth &&
                    hasNorth)
                {
                    const double correction =
                        (An / aPv) *
                        (pPrimeS - pPrimeN);

                    fields.velocity.gety()[P] +=
                        relaxationVelocity *
                        correction;
                }
                else if (hasSouth)
                {
                    // North wall.
                    const double correction =
                        (An / aPv) *
                        pPrimeS;

                    fields.velocity.gety()[P] +=
                        relaxationVelocity *
                        correction;
                }
                else if (hasNorth)
                {
                    // South wall.
                    const double correction =
                        -(An / aPv) *
                        pPrimeN;

                    fields.velocity.gety()[P] +=
                        relaxationVelocity *
                        correction;
                }
            }
        }


        // ------------------------------------------------------------
        // Never allow correction to overwrite prescribed BCs.
        // ------------------------------------------------------------

        applyBoundaryConditions();
    }


    // ================================================================
    // CALCULATE RESIDUAL
    //
    // The residual is the absolute sum of continuity errors over the
    // complete set of control volumes.
    //
    // We do NOT simply discard inlet/wall cells anymore.
    //
    // Instead, the boundary fluxes themselves are defined so that
    // prescribed walls have zero mass flux and prescribed inlets have
    // the correct specified mass flux.
    // ================================================================

    double SIMPLE::calculateResidual()
    {
        const std::size_t n =
            mesh.numberOfCells();

        double totalResidual = 0.0;


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


    // ================================================================
    // CHECK CONVERGENCE
    // ================================================================

    bool SIMPLE::checkConvergence()
    {
        residual =
            calculateResidual();

        return residual <
            convergenceTolerance;
    }


    // ================================================================
    // MAIN SIMPLE SOLVER
    // ================================================================

    void SIMPLE::solve()
    {
        iteration = 0;
        residual = std::numeric_limits<double>::infinity();


        for (iteration = 1;
            iteration <= maxIterations;
            ++iteration)
        {
            // ========================================================
            // 1. Apply physical boundary conditions
            // ========================================================

            applyBoundaryConditions();


            // ========================================================
            // 2. Calculate current face mass fluxes
            // ========================================================

            calculateFaceFluxes();


            // ========================================================
            // 3. Solve U momentum
            // ========================================================

            solveUMomentum();


            // ========================================================
            // 4. Solve V momentum
            // ========================================================

            solveVMomentum();


            // ========================================================
            // 5. Momentum solution changed the velocities.
            //    Therefore the mass fluxes MUST be recalculated.
            // ========================================================

            calculateFaceFluxes();


            // ========================================================
            // 6. Assemble pressure correction
            // ========================================================

            assemblePressureCorrection();


            // ========================================================
            // 7. Solve pressure correction
            // ========================================================

            solvePressureCorrection();


            // ========================================================
            // 8. Correct pressure
            // ========================================================

            correctPressure();


            // ========================================================
            // 9. Correct velocity using the SAME pressure correction
            //    relationship used to derive the pressure equation.
            // ========================================================

            correctVelocity();


            // ========================================================
            // 10. Velocity changed again.
            //     Recalculate mass fluxes.
            // ========================================================

            calculateFaceFluxes();


            // ========================================================
            // 11. Calculate continuity residual
            // ========================================================

            residual =
                calculateResidual();


            std::cout
                << "SIMPLE iteration "
                << iteration
                << " | Residual = "
                << residual
                << "\n";


            // ========================================================
            // 12. Check convergence
            // ========================================================

            if (checkConvergence())
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

