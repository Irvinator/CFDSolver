#include "BCs/BC.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{
    void expect(bool condition, const std::string& message)
    {
        if (!condition)
        {
            throw std::runtime_error("Test failed: " + message);
        }
    }

    void printHorizontalLine()
    {
        std::cout
            << "+------------------+"
            << "------------------+"
            << "------------------+"
            << "------------------+\n";
    }
}


int main()
{
    try
    {
        // ============================================================
        // 1. CREATE BOUNDARY CONDITIONS
        // ============================================================

        // ------------------------------------------------------------
        // West boundary = inlet
        // Velocity = 2 m/s
        // ------------------------------------------------------------

        CFD::BoundaryCondition west(
            CFD::BoundarySide::west,
            CFD::BoundaryType::Inlet
        );

        west.setVelocity(2.0, 0.0);


        // ------------------------------------------------------------
        // East boundary = outlet
        // Velocity = 5 m/s
        // ------------------------------------------------------------

        CFD::BoundaryCondition east(
            CFD::BoundarySide::east,
            CFD::BoundaryType::Outlet
        );

        east.setVelocity(5.0, 0.0);


        // ------------------------------------------------------------
        // North boundary = wall
        // Velocity = 0 m/s
        // ------------------------------------------------------------

        CFD::BoundaryCondition north(
            CFD::BoundarySide::north,
            CFD::BoundaryType::Wall
        );

        north.setVelocity(0.0, 0.0);


        // ------------------------------------------------------------
        // South boundary = wall
        // Velocity = 0 m/s
        // ------------------------------------------------------------

        CFD::BoundaryCondition south(
            CFD::BoundarySide::south,
            CFD::BoundaryType::Wall
        );

        south.setVelocity(0.0, 0.0);


        // ============================================================
        // 2. TEST BOUNDARY SIDES
        // ============================================================

        expect(
            west.getSide() == CFD::BoundarySide::west,
            "West boundary side is incorrect"
        );

        expect(
            east.getSide() == CFD::BoundarySide::east,
            "East boundary side is incorrect"
        );

        expect(
            north.getSide() == CFD::BoundarySide::north,
            "North boundary side is incorrect"
        );

        expect(
            south.getSide() == CFD::BoundarySide::south,
            "South boundary side is incorrect"
        );

        std::cout << "Boundary side tests passed.\n";


        // ============================================================
        // 3. TEST BOUNDARY TYPES
        // ============================================================

        expect(
            west.getType() == CFD::BoundaryType::Inlet,
            "West boundary should be an Inlet"
        );

        expect(
            east.getType() == CFD::BoundaryType::Outlet,
            "East boundary should be an Outlet"
        );

        expect(
            north.getType() == CFD::BoundaryType::Wall,
            "North boundary should be a Wall"
        );

        expect(
            south.getType() == CFD::BoundaryType::Wall,
            "South boundary should be a Wall"
        );

        std::cout << "Boundary type tests passed.\n";


        // ============================================================
        // 4. TEST VELOCITIES
        // ============================================================

        // West inlet
        expect(
            west.getU() == 2.0,
            "West inlet U velocity should be 2.0 m/s"
        );

        expect(
            west.getV() == 0.0,
            "West inlet V velocity should be 0.0 m/s"
        );


        // East outlet
        expect(
            east.getU() == 5.0,
            "East outlet U velocity should be 5.0 m/s"
        );

        expect(
            east.getV() == 0.0,
            "East outlet V velocity should be 0.0 m/s"
        );


        // North wall
        expect(
            north.getU() == 0.0,
            "North wall U velocity should be 0.0 m/s"
        );

        expect(
            north.getV() == 0.0,
            "North wall V velocity should be 0.0 m/s"
        );


        // South wall
        expect(
            south.getU() == 0.0,
            "South wall U velocity should be 0.0 m/s"
        );

        expect(
            south.getV() == 0.0,
            "South wall V velocity should be 0.0 m/s"
        );

        std::cout << "Velocity tests passed.\n";


        // ============================================================
        // 5. TEST SPECIFIED FLAGS
        // ============================================================

        expect(
            west.hasU(),
            "West U should be specified"
        );

        expect(
            west.hasV(),
            "West V should be specified"
        );


        expect(
            east.hasU(),
            "East U should be specified"
        );

        expect(
            east.hasV(),
            "East V should be specified"
        );


        expect(
            north.hasU(),
            "North U should be specified"
        );

        expect(
            north.hasV(),
            "North V should be specified"
        );


        expect(
            south.hasU(),
            "South U should be specified"
        );

        expect(
            south.hasV(),
            "South V should be specified"
        );


        // Pressure has not been specified
        expect(
            !west.hasPressure(),
            "West pressure should not be specified"
        );

        expect(
            !east.hasPressure(),
            "East pressure should not be specified"
        );

        expect(
            !north.hasPressure(),
            "North pressure should not be specified"
        );

        expect(
            !south.hasPressure(),
            "South pressure should not be specified"
        );

        std::cout << "Specified flag tests passed.\n";


        // ============================================================
        // 6. TEST NAME FUNCTIONS
        // ============================================================

        expect(
            west.getSideName() == "West",
            "West side name is incorrect"
        );

        expect(
            east.getSideName() == "East",
            "East side name is incorrect"
        );

        expect(
            north.getSideName() == "North",
            "North side name is incorrect"
        );

        expect(
            south.getSideName() == "South",
            "South side name is incorrect"
        );


        expect(
            west.getTypeName() == "Inlet",
            "West type name is incorrect"
        );

        expect(
            east.getTypeName() == "Outlet",
            "East type name is incorrect"
        );

        expect(
            north.getTypeName() == "Wall",
            "North type name is incorrect"
        );

        expect(
            south.getTypeName() == "Wall",
            "South type name is incorrect"
        );

        std::cout << "Name tests passed.\n";


        // ============================================================
        // 7. CREATE 4 x 4 DOMAIN
        // ============================================================

        const int nx = 4;
        const int ny = 4;

        const double width = 4.0;
        const double height = 4.0;

        const double dx = width / nx;
        const double dy = height / ny;


        // ============================================================
        // 8. DISPLAY DOMAIN INFORMATION
        // ============================================================

        std::cout << "\n";
        std::cout << "============================================================\n";
        std::cout << "                  4 x 4 CFD DOMAIN\n";
        std::cout << "============================================================\n\n";

        std::cout << "Domain width  = "
            << width
            << " m\n";

        std::cout << "Domain height = "
            << height
            << " m\n";

        std::cout << "Number of cells = "
            << nx
            << " x "
            << ny
            << "\n";

        std::cout << "Cell size = "
            << dx
            << " x "
            << dy
            << " m\n\n";


        // ============================================================
        // 9. DISPLAY BOUNDARY CONDITIONS
        // ============================================================

        std::cout << "Boundary conditions:\n\n";

        std::cout << "West  : "
            << west.getTypeName()
            << " | U = "
            << west.getU()
            << " m/s"
            << " | V = "
            << west.getV()
            << " m/s\n";

        std::cout << "East  : "
            << east.getTypeName()
            << " | U = "
            << east.getU()
            << " m/s"
            << " | V = "
            << east.getV()
            << " m/s\n";

        std::cout << "North : "
            << north.getTypeName()
            << " | U = "
            << north.getU()
            << " m/s"
            << " | V = "
            << north.getV()
            << " m/s\n";

        std::cout << "South : "
            << south.getTypeName()
            << " | U = "
            << south.getU()
            << " m/s"
            << " | V = "
            << south.getV()
            << " m/s\n";


        // ============================================================
        // 10. VISUALISE 4 x 4 MATRIX
        // ============================================================

        std::cout << "\n";
        std::cout << "============================================================\n";
        std::cout << "             BOUNDARY VELOCITY MATRIX\n";
        std::cout << "============================================================\n\n";

        std::cout << "Each cell displays:\n";
        std::cout << "  (i,j) = cell index\n";
        std::cout << "  x,y   = cell position\n";
        std::cout << "  U     = boundary velocity\n";
        std::cout << "  N/A   = interior cell\n\n";


        // Print from north to south
        for (int j = ny - 1; j >= 0; --j)
        {
            printHorizontalLine();


            // --------------------------------------------------------
            // Cell indices
            // --------------------------------------------------------

            for (int i = 0; i < nx; ++i)
            {
                std::ostringstream cellIndex;

                cellIndex
                    << "("
                    << i
                    << ","
                    << j
                    << ")";

                std::cout
                    << "| "
                    << std::left
                    << std::setw(16)
                    << cellIndex.str();
            }

            std::cout << "|\n";


            // --------------------------------------------------------
            // Cell positions
            // --------------------------------------------------------

            for (int i = 0; i < nx; ++i)
            {
                double x = (i + 0.5) * dx;
                double y = (j + 0.5) * dy;

                std::ostringstream position;

                position
                    << "x="
                    << std::fixed
                    << std::setprecision(1)
                    << x
                    << ", y="
                    << y;

                std::cout
                    << "| "
                    << std::left
                    << std::setw(16)
                    << position.str();
            }

            std::cout << "|\n";


            // --------------------------------------------------------
            // Boundary velocity
            // --------------------------------------------------------

            for (int i = 0; i < nx; ++i)
            {
                std::string velocity;

                // North boundary
                if (j == ny - 1)
                {
                    velocity = "U = 0.0 m/s";
                }

                // South boundary
                else if (j == 0)
                {
                    velocity = "U = 0.0 m/s";
                }

                // West boundary
                else if (i == 0)
                {
                    velocity = "U = 2.0 m/s";
                }

                // East boundary
                else if (i == nx - 1)
                {
                    velocity = "U = 5.0 m/s";
                }

                // Interior
                else
                {
                    velocity = "U = N/A";
                }

                std::cout
                    << "| "
                    << std::left
                    << std::setw(16)
                    << velocity;
            }

            std::cout << "|\n";
        }

        printHorizontalLine();


        // ============================================================
        // 11. EXPECTED MATRIX
        // ============================================================

        std::cout << "\n";
        std::cout << "Expected U-velocity boundary pattern:\n\n";

        std::cout
            << "        "
            << std::setw(12) << "i=0"
            << std::setw(12) << "i=1"
            << std::setw(12) << "i=2"
            << std::setw(12) << "i=3"
            << "\n";

        std::cout
            << "j=3    "
            << std::setw(12) << "0.0"
            << std::setw(12) << "0.0"
            << std::setw(12) << "0.0"
            << std::setw(12) << "0.0"
            << "\n";

        std::cout
            << "j=2    "
            << std::setw(12) << "2.0"
            << std::setw(12) << "N/A"
            << std::setw(12) << "N/A"
            << std::setw(12) << "5.0"
            << "\n";

        std::cout
            << "j=1    "
            << std::setw(12) << "2.0"
            << std::setw(12) << "N/A"
            << std::setw(12) << "N/A"
            << std::setw(12) << "5.0"
            << "\n";

        std::cout
            << "j=0    "
            << std::setw(12) << "0.0"
            << std::setw(12) << "0.0"
            << std::setw(12) << "0.0"
            << std::setw(12) << "0.0"
            << "\n";


        // ============================================================
        // 12. TEST 4 x 4 BOUNDARY LOGIC
        // ============================================================

        // Check North
        for (int i = 0; i < nx; ++i)
        {
            expect(
                north.getU() == 0.0,
                "North boundary velocity is incorrect"
            );
        }


        // Check South
        for (int i = 0; i < nx; ++i)
        {
            expect(
                south.getU() == 0.0,
                "South boundary velocity is incorrect"
            );
        }


        // Check West
        for (int j = 1; j < ny - 1; ++j)
        {
            expect(
                west.getU() == 2.0,
                "West inlet velocity is incorrect"
            );
        }


        // Check East
        for (int j = 1; j < ny - 1; ++j)
        {
            expect(
                east.getU() == 5.0,
                "East outlet velocity is incorrect"
            );
        }

        std::cout << "\n4 x 4 boundary layout tests passed.\n";


        // ============================================================
        // 13. FINAL RESULT
        // ============================================================

        std::cout << "\n";
        std::cout << "============================================================\n";
        std::cout << "          ALL BOUNDARY CONDITION TESTS PASSED\n";
        std::cout << "============================================================\n";

        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "\n";
        std::cerr << "============================================================\n";
        std::cerr << "TEST FAILED\n";
        std::cerr << error.what() << "\n";
        std::cerr << "============================================================\n";

        return 1;
    }
}