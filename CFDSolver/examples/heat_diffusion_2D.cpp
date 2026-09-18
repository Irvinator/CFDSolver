#include "solvers/HeatSolver2D.h"
#include "mesh/mesh2D.h"

#include <exception>
#include <iostream>
#include <string>

int main()
{
    try {
        Mesh mesh(30, 20, 3.0, 2.0);

        CFD::HeatSolver2D solver(mesh, 1.0e-2, 1.0, 2000.0);

        CFD::HeatSolver2D::BoundaryConditions bcs;
        bcs.T_west = 1.0;
        bcs.T_east = 0.0;
        bcs.T_south = 0.5;
        bcs.T_north = 0.5;

        solver.setBCs(bcs);
        solver.setIC(0.5);
        solver.setOutputFreq(50);
        solver.setSteadyTolerance(1e-10);
        solver.enableSteadyStop(true);

        const auto result = solver.run("output/heat2D_example");

        std::cout << "\n=== Example Summary ================================\n";
        std::cout << "Runtime [s]        : " << result.runtime << "\n";
        std::cout << "Steps              : " << result.totalSteps << "\n";
        std::cout << "Reached steady     : " << (result.reachedSteady ? "yes" : "no") << "\n";
        std::cout << "Final dTinf        : " << result.finalDeltaInf << "\n";
        std::cout << "Average CG iters   : " << result.avgCgIterations << "\n";
        std::cout << "Max error          : " << result.maxError << "\n";
        std::cout << "L2 error           : " << result.l2Error << "\n";
        std::cout << "Validated          : " << (result.validated ? "yes" : "no") << "\n";
        std::cout << "Output directory   : output/heat2D_example\n";
        std::cout << "====================================================\n";

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << '\n';
        return 1;
    }
}
