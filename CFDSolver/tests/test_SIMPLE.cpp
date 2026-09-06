#include "solvers/SIMPLE.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace CFD;

namespace
{
    bool allFieldsFinite(const Mesh& mesh, const Fields& fields)
    {
        for (std::size_t P = 0; P < static_cast<std::size_t>(mesh.numberOfCells()); ++P)
        {
            if (!std::isfinite(fields.velocity.getx()[P]) ||
                !std::isfinite(fields.velocity.gety()[P]) ||
                !std::isfinite(fields.pressure[P]))
            {
                return false;
            }
        }

        return true;
    }

    void printFieldRanges(const Mesh& mesh, const Fields& fields)
    {
        double minU = std::numeric_limits<double>::max();
        double maxU = std::numeric_limits<double>::lowest();
        double minV = std::numeric_limits<double>::max();
        double maxV = std::numeric_limits<double>::lowest();
        double minP = std::numeric_limits<double>::max();
        double maxP = std::numeric_limits<double>::lowest();

        for (std::size_t P = 0; P < static_cast<std::size_t>(mesh.numberOfCells()); ++P)
        {
            minU = std::min(minU, static_cast<double>(fields.velocity.getx()[P]));
            maxU = std::max(maxU, static_cast<double>(fields.velocity.getx()[P]));
            minV = std::min(minV, static_cast<double>(fields.velocity.gety()[P]));
            maxV = std::max(maxV, static_cast<double>(fields.velocity.gety()[P]));
            minP = std::min(minP, static_cast<double>(fields.pressure[P]));
            maxP = std::max(maxP, static_cast<double>(fields.pressure[P]));
        }

        std::cout << "Field ranges: "
            << "u=[" << minU << ", " << maxU << "]  "
            << "v=[" << minV << ", " << maxV << "]  "
            << "p=[" << minP << ", " << maxP << "]\n";
    }

    double inletMassFlow(const Mesh& mesh, const Fields& fields, double rho, double prescribedU)
    {
        double massFlow = 0.0;
        const double area = mesh.eastWestFaceArea();

        for (int j = 0; j < mesh.getNy(); ++j)
        {
            const int P = mesh.cellIndex(0, j);
            const auto& cell = mesh.getCell(P);

            if (cell.westBoundary)
            {
                (void)fields;
                massFlow += rho * prescribedU * area;
            }
        }

        return massFlow;
    }

    double outletMassFlow(const Mesh& mesh, const Fields& fields, double rho)
    {
        double massFlow = 0.0;
        const double area = mesh.eastWestFaceArea();

        for (int j = 0; j < mesh.getNy(); ++j)
        {
            const int P = mesh.cellIndex(mesh.getNx() - 1, j);
            const auto& cell = mesh.getCell(P);

            if (cell.eastBoundary)
            {
                massFlow += rho * fields.velocity.getx()[P] * area;
            }
        }

        return massFlow;
    }

    void printVelocityGrid(const Mesh& mesh, const Fields& fields)
    {
        std::cout << "\nU velocity grid\n";
        for (int j = mesh.getNy() - 1; j >= 0; --j)
        {
            std::cout << "j=" << j << " : ";
            for (int i = 0; i < mesh.getNx(); ++i)
            {
                const int P = mesh.cellIndex(i, j);
                std::cout << std::setw(11)
                    << std::fixed
                    << std::setprecision(5)
                    << fields.velocity.getx()[P];
            }
            std::cout << "\n";
        }

        std::cout << "\nV velocity grid\n";
        for (int j = mesh.getNy() - 1; j >= 0; --j)
        {
            std::cout << "j=" << j << " : ";
            for (int i = 0; i < mesh.getNx(); ++i)
            {
                const int P = mesh.cellIndex(i, j);
                std::cout << std::setw(11)
                    << std::fixed
                    << std::setprecision(5)
                    << fields.velocity.gety()[P];
            }
            std::cout << "\n";
        }
    }
}

int main()
{
    try
    {
        std::cout << "============================================================\n";
        std::cout << "        SIMPLE Structured-Pipe Regression / Debug Test      \n";
        std::cout << "============================================================\n";

        const double rho = 1.0;
        const double mu = 0.01;
        const double inletU = 1.0;
        const double inletV = 0.0;
        const double outletPressure = 0.0;

        const int nx = 20;
        const int ny = 10;
        const double width = 5.0;
        const double height = 1.0;

        Mesh mesh(nx, ny, width, height);
        Fields fields(static_cast<std::size_t>(mesh.numberOfCells()));
        fields.initialise(0.0, 0.0, 0.0);

        BoundaryCondition westBC(BoundarySide::west, BoundaryType::Inlet);
        westBC.setVelocity(inletU, inletV);

        BoundaryCondition eastBC(BoundarySide::east, BoundaryType::Outlet);
        eastBC.setPressure(outletPressure);

        BoundaryCondition northBC(BoundarySide::north, BoundaryType::Wall);
        northBC.setVelocity(0.0, 0.0);

        BoundaryCondition southBC(BoundarySide::south, BoundaryType::Wall);
        southBC.setVelocity(0.0, 0.0);

        SIMPLE simple(mesh, fields, northBC, southBC, eastBC, westBC);
        simple.setDensity(rho);
        simple.setViscosity(mu);
        simple.setVelocityRelaxation(0.2);
        simple.setPressureRelaxation(0.1);
        simple.setConvergenceTolerance(1.0e-5);
        simple.setMaxIterations(1);

        const std::size_t maxOuterSteps = 400;
        const std::size_t printEvery = 10;

        std::vector<double> residualHistory;
        residualHistory.reserve(maxOuterSteps);

        std::cout << "\nStarting one-SIMPLE-step-at-a-time debug run...\n\n";

        for (std::size_t step = 1; step <= maxOuterSteps; ++step)
        {
            simple.solve();

            const double residual = simple.getResidual();
            residualHistory.push_back(residual);

            if (!std::isfinite(residual))
            {
                throw std::runtime_error("Residual became NaN/Inf");
            }

            if (!allFieldsFinite(mesh, fields))
            {
                throw std::runtime_error("Field variable became NaN/Inf");
            }

            const double mIn = inletMassFlow(mesh, fields, rho, inletU);
            const double mOut = outletMassFlow(mesh, fields, rho);
            const double massMismatch = std::abs(mIn - mOut);

            if (step == 1 || step % printEvery == 0 || residual < simple.getConvergenceTolerance())
            {
                std::cout << "Step " << std::setw(4) << step
                    << " | residual = " << std::scientific << residual
                    << " | inlet = " << mIn
                    << " | outlet = " << mOut
                    << " | |m_in-m_out| = " << massMismatch
                    << "\n";

                printFieldRanges(mesh, fields);
            }

            if (residual < simple.getConvergenceTolerance())
            {
                std::cout << "\nConverged after " << step << " outer SIMPLE steps.\n";
                printVelocityGrid(mesh, fields);
                return 0;
            }
        }

        std::cout << "\nSolver did not reach tolerance within " << maxOuterSteps << " outer SIMPLE steps.\n";

        if (!residualHistory.empty())
        {
            std::cout << "Initial residual = " << residualHistory.front() << "\n";
            std::cout << "Final residual   = " << residualHistory.back() << "\n";
        }

        printVelocityGrid(mesh, fields);
        return 1;
    }
    catch (const std::exception& e)
    {
        std::cerr << "\n============================================================\n";
        std::cerr << "TEST FAILED WITH EXCEPTION\n";
        std::cerr << "============================================================\n";
        std::cerr << e.what() << "\n";
        return 1;
    }
}
