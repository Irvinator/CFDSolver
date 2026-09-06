// 

#include "fields/Fields.h"
#include "mesh/mesh2D.h"

#include <iostream>
#include <stdexcept>

int main()
{
    try
    {
        // Create a 20 x 20 mesh
        Mesh mesh(20, 20, 1.0, 1.0);

        // Create all fields using the mesh
        CFD::Fields fields(mesh);

        std::cout << "Number of cells: "
            << fields.size() << "\n";

        // Check pressure
        fields.pressure[0] = 101325.0;

        // Check velocity
        fields.velocity.getx()[0] = 1.0;
        fields.velocity.gety()[0] = 0.5;

        std::cout << "Pressure at cell 0: "
            << fields.pressure[0] << "\n";

        std::cout << "u velocity at cell 0: "
            << fields.velocity.getx()[0] << "\n";

        std::cout << "v velocity at cell 0: "
            << fields.velocity.gety()[0] << "\n";

        // Check pressure correction
        fields.pressureCorrection[0] = 10.0;

        std::cout << "Pressure correction at cell 0: "
            << fields.pressureCorrection[0] << "\n";

        // Check initialisation
        
        std::cout << "\nAfter initialise():\n";

		for (std::size_t j = 0; j < mesh.getNx(); ++j)
		{
			for (std::size_t i = 0; i < mesh.getNy(); ++i)
			{
                fields.initialise(500.0, i, j);

				std::cout << mesh.cellIndex(i, j) << ":\n";

                std::cout << "Pressure: "
                    << fields.pressure[mesh.cellIndex(i, j)] << "\n";

                std::cout << "u velocity: "
                    << fields.velocity.getx()[mesh.cellIndex(i, j)] << "\n";

                std::cout << "v velocity: "
                    << fields.velocity.gety()[mesh.cellIndex(i, j)] << "\n";

                std::cout << "\nFields test passed.\n";
				
			}
		}
       
    }
    catch (const std::exception& e)
    {
        std::cerr << "Fields test failed: "
            << e.what() << "\n";

        return 1;
    }

    return 0;
}