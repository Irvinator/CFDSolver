#include "solvers/StaggeredSIMPLE.h"
#include "fields/StaggeredFields.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace CFD;

namespace
{
    bool allFieldsFinite(const StaggeredFields& fields)
    {
        for (double value : fields.pressureData()) if (!std::isfinite(value)) return false;
        for (double value : fields.uData()) if (!std::isfinite(value)) return false;
        for (double value : fields.vData()) if (!std::isfinite(value)) return false;
        return true;
    }

    double inletMassFlow(const Mesh& mesh, const StaggeredFields& fields, double rho)
    {
        double massFlow = 0.0;
        const double area = mesh.eastWestFaceArea();
        for (int j = 0; j < mesh.getNy(); ++j) massFlow += rho * fields.u(0, j) * area;
        return massFlow;
    }

    double outletMassFlow(const Mesh& mesh, const StaggeredFields& fields, double rho)
    {
        double massFlow = 0.0;
        const double area = mesh.eastWestFaceArea();
        for (int j = 0; j < mesh.getNy(); ++j) massFlow += rho * fields.u(mesh.getNx(), j) * area;
        return massFlow;
    }
}

int main()
{
    try
    {
        constexpr double rho = 1.0;
        constexpr double mu = 0.01;
        constexpr double inletU = 1.0;
        constexpr double outletPressure = 0.0;
        constexpr int nx = 8;
        constexpr int ny = 8;
        constexpr double width = 5.0;
        constexpr double height = 1.0;
        constexpr double continuityTolerance = 1.0e-5;
        constexpr double momentumTolerance = 1.0e-6;
        constexpr std::size_t maxOuterSteps = 400;

        Mesh mesh(nx, ny, width, height);
        StaggeredFields fields(nx, ny);
        fields.initialise(0.0, 0.0, 0.0);

        BoundaryCondition westBC(BoundarySide::west, BoundaryType::Inlet);
        westBC.setVelocity(inletU, 0.0);

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
        simple.setConvergenceTolerance(continuityTolerance);
        simple.setMomentumConvergenceTolerance(momentumTolerance);

        // One SIMPLE sweep per call lets this test print one complete residual
        // triplet per outer step while retaining the current field state.
        simple.setMaxIterations(1);

        std::cout << "Step | continuity | U momentum | V momentum | inlet | outlet | imbalance\n";
        for (std::size_t step = 1; step <= maxOuterSteps; ++step)
        {
            simple.solve();

            const double continuity = simple.getContinuityResidual();
            const double uMomentum = simple.getUMomentumResidual();
            const double vMomentum = simple.getVMomentumResidual();
            if (!std::isfinite(continuity) || !std::isfinite(uMomentum) ||
                !std::isfinite(vMomentum) || !allFieldsFinite(fields))
            {
                throw std::runtime_error("Solver produced a non-finite residual or field value");
            }

            const double mIn = inletMassFlow(mesh, fields, rho);
            const double mOut = outletMassFlow(mesh, fields, rho);
            const bool converged =
                continuity < continuityTolerance &&
                uMomentum < momentumTolerance &&
                vMomentum < momentumTolerance;

            if (step == 1 || step % 10 == 0 || converged)
            {
                std::cout << std::setw(4) << step << " | "
                    << std::scientific << std::setprecision(4)
                    << continuity << " | " << uMomentum << " | " << vMomentum << " | "
                    << mIn << " | " << mOut << " | " << std::abs(mIn - mOut) << '\n';
            }

            if (converged)
            {
                std::cout << "\nConverged after " << step << " outer SIMPLE steps.\n";
                return 0;
            }
        }

        std::cerr << "\nSolver did not satisfy all three convergence criteria within "
            << maxOuterSteps << " outer steps.\n";
        return 1;
    }
    catch (const std::exception& exception)
    {
        std::cerr << "Test failed: " << exception.what() << '\n';
        return 1;
    }
}