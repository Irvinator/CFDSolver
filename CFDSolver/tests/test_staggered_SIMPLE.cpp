#include "solvers/StaggeredSIMPLE.h"
#include "fields/StaggeredFields.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace CFD;

namespace
{
    bool allFieldsFinite(const StaggeredFields& fields)
    {
        for (double value : fields.pressureData())
        {
            if (!std::isfinite(value)) return false;
        }
        for (double value : fields.uData())
        {
            if (!std::isfinite(value)) return false;
        }
        for (double value : fields.vData())
        {
            if (!std::isfinite(value)) return false;
        }
        return true;
    }

    double inletMassFlow(const Mesh& mesh, const StaggeredFields& fields, double rho)
    {
        double massFlow = 0.0;
        const double area = mesh.eastWestFaceArea();
        for (int j = 0; j < mesh.getNy(); ++j)
        {
            massFlow += rho * fields.u(0, j) * area;
        }
        return massFlow;
    }

    double outletMassFlow(const Mesh& mesh, const StaggeredFields& fields, double rho)
    {
        double massFlow = 0.0;
        const double area = mesh.eastWestFaceArea();
        for (int j = 0; j < mesh.getNy(); ++j)
        {
            massFlow += rho * fields.u(mesh.getNx(), j) * area;
        }
        return massFlow;
    }

    void printFaceVelocityGrids(const Mesh& mesh, const StaggeredFields& fields)
    {
        std::cout << "\nU velocity on vertical faces\n";
        std::cout << "------------------------------------------------------------\n";

        for (int j = mesh.getNy() - 1; j >= 0; --j)
        {
            std::cout << "j=" << j << " : ";
            for (int i = 0; i <= mesh.getNx(); ++i)
            {
                std::cout << std::setw(11)
                    << std::fixed
                    << std::setprecision(5)
                    << fields.u(i, j);
            }
            std::cout << "\n";
        }

        std::cout << "\nV velocity on horizontal faces\n";
        std::cout << "------------------------------------------------------------\n";

        for (int j = mesh.getNy(); j >= 0; --j)
        {
            std::cout << "j=" << j << " : ";
            for (int i = 0; i < mesh.getNx(); ++i)
            {
                std::cout << std::setw(11)
                    << std::fixed
                    << std::setprecision(5)
                    << fields.v(i, j);
            }
            std::cout << "\n";
        }
    }

    void printCellCenteredVelocityTable(const Mesh& mesh, const StaggeredFields& fields)
    {
        std::cout << "\nCell-centred U velocity table (compare with Fluent)\n";
        std::cout << "============================================================\n";
        std::cout << "y\\x";
        for (int i = 0; i < mesh.getNx(); ++i)
        {
            const double x = (static_cast<double>(i) + 0.5) * mesh.getDx();
            std::cout << std::setw(12) << std::fixed << std::setprecision(3) << x;
        }
        std::cout << "\n";

        for (int j = mesh.getNy() - 1; j >= 0; --j)
        {
            const double y = (static_cast<double>(j) + 0.5) * mesh.getDy();
            std::cout << std::setw(4) << std::fixed << std::setprecision(3) << y;
            for (int i = 0; i < mesh.getNx(); ++i)
            {
                const double uCell = 0.5 * (fields.u(i, j) + fields.u(i + 1, j));
                std::cout << std::setw(12) << std::fixed << std::setprecision(6) << uCell;
            }
            std::cout << "\n";
        }

        std::cout << "\nCell-centred V velocity table (compare with Fluent)\n";
        std::cout << "============================================================\n";
        std::cout << "y\\x";
        for (int i = 0; i < mesh.getNx(); ++i)
        {
            const double x = (static_cast<double>(i) + 0.5) * mesh.getDx();
            std::cout << std::setw(12) << std::fixed << std::setprecision(3) << x;
        }
        std::cout << "\n";

        for (int j = mesh.getNy() - 1; j >= 0; --j)
        {
            const double y = (static_cast<double>(j) + 0.5) * mesh.getDy();
            std::cout << std::setw(4) << std::fixed << std::setprecision(3) << y;
            for (int i = 0; i < mesh.getNx(); ++i)
            {
                const double vCell = 0.5 * (fields.v(i, j) + fields.v(i, j + 1));
                std::cout << std::setw(12) << std::fixed << std::setprecision(6) << vCell;
            }
            std::cout << "\n";
        }
    }

    void writeCellCenteredVelocityCSV(const Mesh& mesh, const StaggeredFields& fields, const std::string& filename)
    {
        std::ofstream file(filename);
        if (!file)
        {
            throw std::runtime_error("Failed to open output file: " + filename);
        }

        file << "x,y,u,v,p\n";

        for (int j = 0; j < mesh.getNy(); ++j)
        {
            for (int i = 0; i < mesh.getNx(); ++i)
            {
                const double x = (static_cast<double>(i) + 0.5) * mesh.getDx();
                const double y = (static_cast<double>(j) + 0.5) * mesh.getDy();
                const double uCell = 0.5 * (fields.u(i, j) + fields.u(i + 1, j));
                const double vCell = 0.5 * (fields.v(i, j) + fields.v(i, j + 1));
                const double pCell = fields.p(i, j);

                file << x << ','
                    << y << ','
                    << uCell << ','
                    << vCell << ','
                    << pCell << '\n';
            }
        }
    }
}

int main()
{
    try
    {
        std::cout << "============================================================\n";
        std::cout << "      STAGGERED SIMPLE Structured-Pipe Regression Test      \n";
        std::cout << "============================================================\n";

        const double rho = 1.0;
        const double mu = 0.01;
        const double inletU = 1.0;
        const double inletV = 0.0;
        const double outletPressure = 0.0;

        const int nx = 8;
        const int ny = 8;
        const double width = 5.0;
        const double height = 1.0;

        Mesh mesh(nx, ny, width, height);
        StaggeredFields fields(nx, ny);
        fields.initialise(0.0, 0.0, 0.0);

        BoundaryCondition westBC(BoundarySide::west, BoundaryType::Inlet);
        westBC.setVelocity(inletU, inletV);

        BoundaryCondition eastBC(BoundarySide::east, BoundaryType::Outlet);
        eastBC.setPressure(outletPressure);

        BoundaryCondition northBC(BoundarySide::north, BoundaryType::Wall);
        northBC.setVelocity(0.0, 0.0);

        BoundaryCondition southBC(BoundarySide::south, BoundaryType::Wall);
        southBC.setVelocity(0.0, 0.0);

        StaggeredSIMPLE simple(mesh, fields, northBC, southBC, eastBC, westBC);
        simple.setDensity(rho);
        simple.setViscosity(mu);
        simple.setVelocityRelaxation(0.5);
        simple.setPressureRelaxation(0.3);
        simple.setConvergenceTolerance(1.0e-5);
        simple.setMaxIterations(1);

        const std::size_t maxOuterSteps = 400;

        for (std::size_t step = 1; step <= maxOuterSteps; ++step)
        {
            simple.solve();

            const double residual = simple.getResidual();
            if (!std::isfinite(residual))
            {
                throw std::runtime_error("Residual became NaN/Inf");
            }
            if (!allFieldsFinite(fields))
            {
                throw std::runtime_error("Field variable became NaN/Inf");
            }

            const double mIn = inletMassFlow(mesh, fields, rho);
            const double mOut = outletMassFlow(mesh, fields, rho);

            if (step == 1 || step % 10 == 0 || residual < simple.getConvergenceTolerance())
            {
                std::cout << "Step " << std::setw(4) << step
                    << " | residual = " << std::scientific << residual
                    << " | inlet = " << mIn
                    << " | outlet = " << mOut
                    << " | |m_in-m_out| = " << std::abs(mIn - mOut)
                    << "\n";
            }

            if (residual < simple.getConvergenceTolerance())
            {
                std::cout << "\nConverged after " << step << " outer SIMPLE steps.\n";
                printFaceVelocityGrids(mesh, fields);
                printCellCenteredVelocityTable(mesh, fields);
                writeCellCenteredVelocityCSV(mesh, fields, "staggered_velocity_output.csv");
                std::cout << "\nWrote cell-centred velocity data to staggered_velocity_output.csv\n";
                return 0;
            }
        }

        std::cout << "\nSolver did not reach tolerance within " << maxOuterSteps << " outer SIMPLE steps.\n";
        printFaceVelocityGrids(mesh, fields);
        printCellCenteredVelocityTable(mesh, fields);
        writeCellCenteredVelocityCSV(mesh, fields, "staggered_velocity_output.csv");
        std::cout << "\nWrote cell-centred velocity data to staggered_velocity_output.csv\n";
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
