#include "solvers/HeatSolver2D.h"
#include "mesh/mesh2D.h"

#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

    void expect(bool condition, const std::string& message)
    {
        if (!condition) {
            throw std::runtime_error("Test failed: " + message);
        }
    }

    void expectNear(double actual, double expected, double tol, const std::string& message)
    {
        if (std::abs(actual - expected) > tol) {
            throw std::runtime_error(
                "Test failed: " + message +
                " | expected=" + std::to_string(expected) +
                " actual=" + std::to_string(actual));
        }
    }

    double analyticalLinearX(double x, double width, double T_west, double T_east)
    {
        return T_west * (1.0 - x / width) + T_east * (x / width);
    }

    void testConstructorAndAccessors()
    {
        Mesh mesh(8, 6, 2.0, 1.0);
        CFD::HeatSolver2D solver(mesh, 1.0e-3, 0.1, 1.0);

        expectNear(solver.time(), 0.0, 1e-14, "Initial time should be zero");
        expect(solver.steps() == 0, "Initial step count should be zero");
        expect(!solver.finished(), "Solver should not be finished initially");
        expect(solver.T().size() == mesh.numberOfCells(), "Temperature field size mismatch");
    }

    void testUniformInitialCondition()
    {
        Mesh mesh(4, 3, 2.0, 1.5);
        CFD::HeatSolver2D solver(mesh, 1.0e-3, 0.1, 0.1);

        solver.setIC(2.5);

        for (int i = 0; i < mesh.numberOfCells(); ++i) {
            expectNear(solver.T()[i], 2.5, 1e-14, "Uniform IC not applied correctly");
        }
    }

    void testCustomInitialCondition()
    {
        Mesh mesh(5, 4, 2.0, 1.0);
        CFD::HeatSolver2D solver(mesh, 1.0e-3, 0.1, 0.1);

        solver.setIC([](double x, double y) {
            return 2.0 * x + 3.0 * y;
            });

        for (int i = 0; i < mesh.numberOfCells(); ++i) {
            const Mesh::Cell& c = mesh.getCell(i);
            const double expected = 2.0 * c.x + 3.0 * c.y;
            expectNear(solver.T()[i], expected, 1e-14, "Custom IC not applied correctly");
        }
    }

    void testSteadyDirichletCaseRunsReasonably()
    {
        Mesh mesh(30, 20, 3.0, 2.0);
        CFD::HeatSolver2D solver(mesh, 1.0e-2, 1.0, 2000.0);

        CFD::HeatSolver2D::BoundaryConditions bcs;
        bcs.T_west = 1.0;
        bcs.T_east = 0.0;
        bcs.T_south = 0.5;
        bcs.T_north = 0.5;

        solver.setBCs(bcs);
        solver.setIC(0.5);
        solver.setOutputFreq(1000000);
        solver.setSteadyTolerance(1e-10);
        solver.enableSteadyStop(true);

        const auto result = solver.run();

        expect(result.reachedSteady,
            "Solver should reach steady state for this Dirichlet case");
        expect(result.totalSteps > 0,
            "Solver should take at least one step");
        expect(result.avgCgIterations > 0.0,
            "Average CG iterations should be positive");

        for (int i = 0; i < mesh.numberOfCells(); ++i) {
            const double T = solver.T()[i];
            expect(T >= -1e-10, "Temperature should not drop below minimum BC significantly");
            expect(T <= 1.0 + 1e-10, "Temperature should not exceed maximum BC significantly");
        }
    }

    void testInvalidParameters()
    {
        bool threw = false;
        try {
            Mesh mesh(4, 4, 1.0, 1.0);
            CFD::HeatSolver2D solver(mesh, -1.0, 0.1, 1.0);
        }
        catch (const std::invalid_argument&) {
            threw = true;
        }
        expect(threw, "Negative alpha should throw");

        threw = false;
        try {
            Mesh mesh(4, 4, 1.0, 1.0);
            CFD::HeatSolver2D solver(mesh, 1.0, 0.0, 1.0);
        }
        catch (const std::invalid_argument&) {
            threw = true;
        }
        expect(threw, "Zero dt should throw");

        threw = false;
        try {
            Mesh mesh(4, 4, 1.0, 1.0);
            CFD::HeatSolver2D solver(mesh, 1.0, 0.1, 0.0);
        }
        catch (const std::invalid_argument&) {
            threw = true;
        }
        expect(threw, "Zero tEnd should throw");
    }

    void testInvalidOutputFrequencyAndTolerance()
    {
        Mesh mesh(4, 4, 1.0, 1.0);
        CFD::HeatSolver2D solver(mesh, 1.0e-3, 0.1, 1.0);

        bool threw = false;
        try {
            solver.setOutputFreq(0);
        }
        catch (const std::invalid_argument&) {
            threw = true;
        }
        expect(threw, "Output frequency of zero should throw");

        threw = false;
        try {
            solver.setSteadyTolerance(0.0);
        }
        catch (const std::invalid_argument&) {
            threw = true;
        }
        expect(threw, "Non-positive steady tolerance should throw");
    }

} // namespace

int main()
{
    try {
        testConstructorAndAccessors();
        testUniformInitialCondition();
        testCustomInitialCondition();
        testSteadyDirichletCaseRunsReasonably();
        testInvalidParameters();
        testInvalidOutputFrequencyAndTolerance();

        std::cout << "All HeatSolver2D tests passed.\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
