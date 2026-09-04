
#include "solvers/SIMPLE.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace CFD;


// ============================================================
// Print one cell in detail
// ============================================================

void printCellDiagnostics(
    const Mesh& mesh,
    const Fields& fields,
    std::size_t P,
    double rho,
    double mu)
{
    const auto& cell =
        mesh.getCell(static_cast<int>(P));

    const double dx = mesh.getDx();
    const double dy = mesh.getDy();

    const double Ae =
        mesh.eastWestFaceArea();

    const double An =
        mesh.northSouthFaceArea();

    const double uP =
        fields.velocity.getx()[P];

    const double vP =
        fields.velocity.gety()[P];

    const double pP =
        fields.pressure[P];


    std::cout << "\n";
    std::cout << "============================================================\n";
    std::cout << "DETAILED CELL DIAGNOSTICS\n";
    std::cout << "============================================================\n";

    std::cout
        << "Cell P = " << P
        << " (i=" << cell.i
        << ", j=" << cell.j << ")\n";

    std::cout
        << "U_P = " << uP
        << " | V_P = " << vP
        << " | P_P = " << pP
        << "\n";


    // ========================================================
    // NEIGHBOUR VELOCITIES
    // ========================================================

    double uE = uP;
    double uW = uP;
    double vN = vP;
    double vS = vP;

    double pE = pP;
    double pW = pP;
    double pN = pP;
    double pS = pP;


    if (cell.east != -1)
    {
        uE = fields.velocity.getx()[cell.east];
        pE = fields.pressure[cell.east];
    }

    if (cell.west != -1)
    {
        uW = fields.velocity.getx()[cell.west];
        pW = fields.pressure[cell.west];
    }

    if (cell.north != -1)
    {
        vN = fields.velocity.gety()[cell.north];
        pN = fields.pressure[cell.north];
    }

    if (cell.south != -1)
    {
        vS = fields.velocity.gety()[cell.south];
        pS = fields.pressure[cell.south];
    }


    // ========================================================
    // FACE VELOCITIES
    // ========================================================

    const double uEface =
        0.5 * (uP + uE);

    const double uWface =
        0.5 * (uP + uW);

    const double vNface =
        0.5 * (vP + vN);

    const double vSface =
        0.5 * (vP + vS);


    std::cout << "\nFACE VELOCITIES\n";
    std::cout << "------------------------------------------------------------\n";

    std::cout
        << "u_E face = " << uEface
        << "\n";

    std::cout
        << "u_W face = " << uWface
        << "\n";

    std::cout
        << "v_N face = " << vNface
        << "\n";

    std::cout
        << "v_S face = " << vSface
        << "\n";


    // ========================================================
    // MASS FLUXES
    //
    // Positive values here mean velocity in the positive
    // coordinate direction.
    // ========================================================

    double Fe =
        rho * uEface * Ae;

    double Fw =
        rho * uWface * Ae;

    double Fn =
        rho * vNface * An;

    double Fs =
        rho * vSface * An;


    std::cout << "\nMASS FLUXES\n";
    std::cout << "------------------------------------------------------------\n";

    std::cout
        << "F_E = " << Fe << "\n";

    std::cout
        << "F_W = " << Fw << "\n";

    std::cout
        << "F_N = " << Fn << "\n";

    std::cout
        << "F_S = " << Fs << "\n";


    // ========================================================
    // CONTINUITY
    //
    // Net outward mass flux:
    //
    // Fe - Fw + Fn - Fs
    // ========================================================

    const double continuity =
        Fe - Fw + Fn - Fs;


    std::cout << "\nCONTINUITY\n";
    std::cout << "------------------------------------------------------------\n";

    std::cout
        << "Fe - Fw + Fn - Fs = "
        << continuity
        << "\n";

    std::cout
        << "Absolute imbalance = "
        << std::abs(continuity)
        << "\n";


    // ========================================================
    // DIFFUSION COEFFICIENTS
    // ========================================================

    const double De =
        mu * Ae / dx;

    const double Dw =
        mu * Ae / dx;

    const double Dn =
        mu * An / dy;

    const double Ds =
        mu * An / dy;


    std::cout << "\nDIFFUSION COEFFICIENTS\n";
    std::cout << "------------------------------------------------------------\n";

    std::cout
        << "D_E = " << De << "\n";

    std::cout
        << "D_W = " << Dw << "\n";

    std::cout
        << "D_N = " << Dn << "\n";

    std::cout
        << "D_S = " << Ds << "\n";


    // ========================================================
    // UPWIND MOMENTUM COEFFICIENTS
    // ========================================================

    const double aE =
        De + std::max(-Fe, 0.0);

    const double aW =
        Dw + std::max(Fw, 0.0);

    const double aN =
        Dn + std::max(-Fn, 0.0);

    const double aS =
        Ds + std::max(Fs, 0.0);


    const double aP =
        aE + aW + aN + aS;


    std::cout << "\nMOMENTUM COEFFICIENTS\n";
    std::cout << "------------------------------------------------------------\n";

    std::cout
        << "a_E = " << aE << "\n";

    std::cout
        << "a_W = " << aW << "\n";

    std::cout
        << "a_N = " << aN << "\n";

    std::cout
        << "a_S = " << aS << "\n";

    std::cout
        << "a_P = " << aP << "\n";


    // ========================================================
    // PRESSURE GRADIENT TERMS
    // ========================================================

    double pressureGradientU = 0.0;
    double pressureGradientV = 0.0;

    if (cell.east != -1 &&
        cell.west != -1)
    {
        pressureGradientU =
            (pW - pE) * An;
    }

    if (cell.north != -1 &&
        cell.south != -1)
    {
        pressureGradientV =
            (pS - pN) * Ae;
    }


    std::cout << "\nPRESSURE TERMS\n";
    std::cout << "------------------------------------------------------------\n";

    std::cout
        << "U pressure source = "
        << pressureGradientU
        << "\n";

    std::cout
        << "V pressure source = "
        << pressureGradientV
        << "\n";


    // ========================================================
    // CONVECTIVE CONTRIBUTIONS
    // ========================================================

    const double convectionU =
        Fe * uEface -
        Fw * uWface +
        Fn * uP -
        Fs * uP;

    const double convectionV =
        Fe * vP -
        Fw * vP +
        Fn * vNface -
        Fs * vSface;


    std::cout << "\nCONVECTIVE TERMS\n";
    std::cout << "------------------------------------------------------------\n";

    std::cout
        << "U convection = "
        << convectionU
        << "\n";

    std::cout
        << "V convection = "
        << convectionV
        << "\n";


    // ========================================================
    // DIFFUSIVE CONTRIBUTIONS
    // ========================================================

    double diffusionU = 0.0;
    double diffusionV = 0.0;

    if (cell.east != -1)
        diffusionU += De * (uE - uP);

    if (cell.west != -1)
        diffusionU += Dw * (uW - uP);

    if (cell.north != -1)
        diffusionU += Dn * (uP - uP);

    if (cell.south != -1)
        diffusionU += Ds * (uP - uP);


    if (cell.north != -1)
        diffusionV += Dn * (vN - vP);

    if (cell.south != -1)
        diffusionV += Ds * (vS - vP);


    std::cout << "\nDIFFUSIVE TERMS\n";
    std::cout << "------------------------------------------------------------\n";

    std::cout
        << "U diffusion = "
        << diffusionU
        << "\n";

    std::cout
        << "V diffusion = "
        << diffusionV
        << "\n";


    // ========================================================
    // CELL DIAGNOSTIC SUMMARY
    // ========================================================

    std::cout << "\nSUMMARY\n";
    std::cout << "------------------------------------------------------------\n";

    std::cout
        << "Continuity imbalance = "
        << std::abs(continuity)
        << "\n";

    std::cout
        << "U convection         = "
        << convectionU
        << "\n";

    std::cout
        << "U diffusion          = "
        << diffusionU
        << "\n";

    std::cout
        << "U pressure           = "
        << pressureGradientU
        << "\n";

    std::cout
        << "V convection         = "
        << convectionV
        << "\n";

    std::cout
        << "V diffusion          = "
        << diffusionV
        << "\n";

    std::cout
        << "V pressure           = "
        << pressureGradientV
        << "\n";
}


// ============================================================
// Global diagnostics
// ============================================================

void printGlobalDiagnostics(
    const Mesh& mesh,
    const Fields& fields,
    double rho,
    double mu)
{
    const std::size_t n =
        mesh.numberOfCells();

    double maxContinuity = 0.0;

    double totalContinuity = 0.0;

    double minU =
        std::numeric_limits<double>::max();

    double maxU =
        std::numeric_limits<double>::lowest();

    double minV =
        std::numeric_limits<double>::max();

    double maxV =
        std::numeric_limits<double>::lowest();

    double minP =
        std::numeric_limits<double>::max();

    double maxP =
        std::numeric_limits<double>::lowest();


    for (std::size_t P = 0;
        P < n;
        ++P)
    {
        const auto& cell =
            mesh.getCell(static_cast<int>(P));

        const double uP =
            fields.velocity.getx()[P];

        const double vP =
            fields.velocity.gety()[P];

        const double pP =
            fields.pressure[P];


        minU = std::min(minU, uP);
        maxU = std::max(maxU, uP);

        minV = std::min(minV, vP);
        maxV = std::max(maxV, vP);

        minP = std::min(minP, pP);
        maxP = std::max(maxP, pP);


        // ----------------------------------------------------
        // Face velocities
        // ----------------------------------------------------

        double uE = uP;
        double uW = uP;
        double vN = vP;
        double vS = vP;


        if (cell.east != -1)
            uE =
            fields.velocity.getx()[cell.east];

        if (cell.west != -1)
            uW =
            fields.velocity.getx()[cell.west];

        if (cell.north != -1)
            vN =
            fields.velocity.gety()[cell.north];

        if (cell.south != -1)
            vS =
            fields.velocity.gety()[cell.south];


        double Fe =
            rho *
            0.5 * (uP + uE) *
            mesh.eastWestFaceArea();

        double Fw =
            rho *
            0.5 * (uP + uW) *
            mesh.eastWestFaceArea();

        double Fn =
            rho *
            0.5 * (vP + vN) *
            mesh.northSouthFaceArea();

        double Fs =
            rho *
            0.5 * (vP + vS) *
            mesh.northSouthFaceArea();


        // ----------------------------------------------------
        // Boundary flux corrections
        // ----------------------------------------------------

        if (cell.west == -1)
        {
            // Inlet
            Fe = Fe;

            Fw =
                rho *
                uP *
                mesh.eastWestFaceArea();
        }

        if (cell.east == -1)
        {
            Fe =
                rho *
                uP *
                mesh.eastWestFaceArea();
        }

        if (cell.north == -1)
        {
            Fn = 0.0;
        }

        if (cell.south == -1)
        {
            Fs = 0.0;
        }


        const double imbalance =
            Fe - Fw + Fn - Fs;

        totalContinuity +=
            std::abs(imbalance);

        maxContinuity =
            std::max(
                maxContinuity,
                std::abs(imbalance));
    }


    // ========================================================
    // Mass flow through inlet
    // ========================================================

    double inletMassFlow = 0.0;

    for (std::size_t j = 0;
        j < mesh.getNy();
        ++j)
    {
        const std::size_t P =
            mesh.cellIndex(0, j);

        const double u =
            fields.velocity.getx()[P];

        inletMassFlow +=
            rho *
            u *
            mesh.eastWestFaceArea();
    }


    // ========================================================
    // Mass flow through outlet
    // ========================================================

    double outletMassFlow = 0.0;

    for (std::size_t j = 0;
        j < mesh.getNy();
        ++j)
    {
        const std::size_t P =
            mesh.cellIndex(
                mesh.getNx() - 1,
                j);

        const double u =
            fields.velocity.getx()[P];

        outletMassFlow +=
            rho *
            u *
            mesh.eastWestFaceArea();
    }


    std::cout << "\n";
    std::cout << "============================================================\n";
    std::cout << "GLOBAL DIAGNOSTICS\n";
    std::cout << "============================================================\n";

    std::cout
        << "Max continuity imbalance = "
        << maxContinuity
        << "\n";

    std::cout
        << "Total continuity error   = "
        << totalContinuity
        << "\n";

    std::cout
        << "Inlet mass flow          = "
        << inletMassFlow
        << "\n";

    std::cout
        << "Outlet mass flow         = "
        << outletMassFlow
        << "\n";

    std::cout
        << "Mass flow difference     = "
        << inletMassFlow - outletMassFlow
        << "\n";

    std::cout
        << "U min                    = "
        << minU
        << "\n";

    std::cout
        << "U max                    = "
        << maxU
        << "\n";

    std::cout
        << "V min                    = "
        << minV
        << "\n";

    std::cout
        << "V max                    = "
        << maxV
        << "\n";

    std::cout
        << "P min                    = "
        << minP
        << "\n";

    std::cout
        << "P max                    = "
        << maxP
        << "\n";
}


// ============================================================
// Print velocity grid
// ============================================================

void printVelocityGrid(
    const Mesh& mesh,
    const Fields& fields)
{
    std::cout << "\n";
    std::cout << "U VELOCITY GRID\n";
    std::cout << "------------------------------------------------------------\n";

    for (int j =
        static_cast<int>(mesh.getNy()) - 1;
        j >= 0;
        --j)
    {
        std::cout
            << "j=" << j << " : ";

        for (std::size_t i = 0;
            i < mesh.getNx();
            ++i)
        {
            const std::size_t P =
                mesh.cellIndex(
                    i,
                    static_cast<std::size_t>(j));

            std::cout
                << std::setw(10)
                << std::fixed
                << std::setprecision(5)
                << fields.velocity.getx()[P];
        }

        std::cout << "\n";
    }


    std::cout << "\n";
    std::cout << "V VELOCITY GRID\n";
    std::cout << "------------------------------------------------------------\n";

    for (int j =
        static_cast<int>(mesh.getNy()) - 1;
        j >= 0;
        --j)
    {
        std::cout
            << "j=" << j << " : ";

        for (std::size_t i = 0;
            i < mesh.getNx();
            ++i)
        {
            const std::size_t P =
                mesh.cellIndex(
                    i,
                    static_cast<std::size_t>(j));

            std::cout
                << std::setw(10)
                << std::fixed
                << std::setprecision(5)
                << fields.velocity.gety()[P];
        }

        std::cout << "\n";
    }
}


// ============================================================
// MAIN
// ============================================================

int main()
{
    try
    {
        std::cout
            << "============================================================\n"
            << "           SIMPLE 5x5 DEBUG PIPE TEST\n"
            << "============================================================\n";


        // ========================================================
        // PHYSICAL PARAMETERS
        // ========================================================

        const double rho = 1.0;
        const double mu = 0.01;

        const double inletU = 1.0;
        const double inletV = 0.0;

        const double outletPressure = 0.0;


        // ========================================================
        // MESH
        // ========================================================

        const std::size_t nx = 5;
        const std::size_t ny = 5;

        const double width = 5.0;
        const double height = 1.0;


        Mesh mesh(
            nx,
            ny,
            width,
            height);


        // ========================================================
        // FIELDS
        // ========================================================

        Fields fields(
            mesh.numberOfCells());


        for (std::size_t P = 0;
            P < mesh.numberOfCells();
            ++P)
        {
            fields.velocity.getx()[P] = 0.0;
            fields.velocity.gety()[P] = 0.0;
            fields.pressure[P] = 0.0;
        }


        // ========================================================
        // BOUNDARY CONDITIONS
        // ========================================================

        BoundaryCondition westBC(
            BoundarySide::west,
            BoundaryType::Inlet);

        westBC.setVelocity(
            inletU,
            inletV);


        BoundaryCondition eastBC(
            BoundarySide::east,
            BoundaryType::Outlet);

        eastBC.setPressure(
            outletPressure);


        BoundaryCondition northBC(
            BoundarySide::north,
            BoundaryType::Wall);

        northBC.setVelocity(
            0.0,
            0.0);


        BoundaryCondition southBC(
            BoundarySide::south,
            BoundaryType::Wall);

        southBC.setVelocity(
            0.0,
            0.0);


        // ========================================================
        // SIMPLE
        // ========================================================

        SIMPLE simple(
            mesh,
            fields,
            northBC,
            southBC,
            eastBC,
            westBC);


        simple.setDensity(rho);

        simple.setVelocityRelaxation(0.7);

        simple.setPressureRelaxation(0.3);

        simple.setConvergenceTolerance(1e-6);


        // ========================================================
        // IMPORTANT:
        //
        // Run ONE SIMPLE iteration at a time.
        //
        // SIMPLE::solve() resets its iteration counter, but the
        // field variables remain where the previous call left them.
        //
        // This lets us inspect the field after every iteration.
        // ========================================================

        const std::size_t debugIterations = 25;


        std::cout << "\n";
        std::cout
            << "============================================================\n"
            << "STARTING ITERATION-BY-ITERATION DEBUGGING\n"
            << "============================================================\n";


        for (std::size_t iteration = 1;
            iteration <= debugIterations;
            ++iteration)
        {
            std::cout << "\n\n";
            std::cout
                << "############################################################\n";

            std::cout
                << "                    ITERATION "
                << iteration
                << "\n";

            std::cout
                << "############################################################\n";


            // ----------------------------------------------------
            // Run exactly one SIMPLE iteration
            // ----------------------------------------------------

            simple.setMaxIterations(1);

            simple.solve();


            // ----------------------------------------------------
            // Global diagnostics
            // ----------------------------------------------------

            printGlobalDiagnostics(
                mesh,
                fields,
                rho,
                mu);


            // ----------------------------------------------------
            // Detailed diagnostics for selected cells
            //
            // Cell 12 is the centre of the 5x5 grid.
            // Cell 6 is near the lower wall.
            // Cell 13 is near the centre/right.
            // ----------------------------------------------------

            std::cout << "\n";
            std::cout
                << "============================================================\n"
                << "SELECTED CELL DEBUGGING\n"
                << "============================================================\n";


            std::cout
                << "\n*** CENTRE CELL 12 ***\n";

            printCellDiagnostics(
                mesh,
                fields,
                12,
                rho,
                mu);


            std::cout
                << "\n*** CELL 6 ***\n";

            printCellDiagnostics(
                mesh,
                fields,
                6,
                rho,
                mu);


            std::cout
                << "\n*** CELL 13 ***\n";

            printCellDiagnostics(
                mesh,
                fields,
                13,
                rho,
                mu);


            // ----------------------------------------------------
            // Velocity grids
            // ----------------------------------------------------

            printVelocityGrid(
                mesh,
                fields);


            // ----------------------------------------------------
            // Stop if NaN / Inf appears
            // ----------------------------------------------------

            bool finite = true;

            for (std::size_t P = 0;
                P < mesh.numberOfCells();
                ++P)
            {
                if (!std::isfinite(
                    fields.velocity.getx()[P]))
                {
                    std::cout
                        << "\nERROR: U became NaN/Inf in cell "
                        << P
                        << "\n";

                    finite = false;
                }


                if (!std::isfinite(
                    fields.velocity.gety()[P]))
                {
                    std::cout
                        << "\nERROR: V became NaN/Inf in cell "
                        << P
                        << "\n";

                    finite = false;
                }


                if (!std::isfinite(
                    fields.pressure[P]))
                {
                    std::cout
                        << "\nERROR: P became NaN/Inf in cell "
                        << P
                        << "\n";

                    finite = false;
                }
            }


            if (!finite)
            {
                std::cout
                    << "\nSTOPPING DEBUGGING BECAUSE "
                    << "NaN/Inf WAS DETECTED.\n";

                return 1;
            }
        }


        // ========================================================
        // FINAL
        // ========================================================

        std::cout << "\n\n";
        std::cout
            << "============================================================\n"
            << "DEBUG TEST COMPLETE\n"
            << "============================================================\n";


        printGlobalDiagnostics(
            mesh,
            fields,
            rho,
            mu);


        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "\n";
        std::cerr
            << "============================================================\n"
            << "TEST FAILED WITH EXCEPTION\n"
            << "============================================================\n";

        std::cerr
            << e.what()
            << "\n";

        return 1;
    }
}

