#include "mesh/mesh2D.h"

#include <sstream>
#include <cmath>
#include <exception>
#include <iomanip>
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

    void printMeshSummary(const Mesh& mesh)
    {
        std::cout << "\n=== Mesh Summary ======================================\n";
        std::cout << "Nx                 : " << mesh.getNx() << "\n";
        std::cout << "Ny                 : " << mesh.getNy() << "\n";
        std::cout << "Width              : " << mesh.getWidth() << "\n";
        std::cout << "Height             : " << mesh.getHeight() << "\n";
        std::cout << "dx                 : " << mesh.getDx() << "\n";
        std::cout << "dy                 : " << mesh.getDy() << "\n";
        std::cout << "Cell volume        : " << mesh.cellVolume() << "\n";
        std::cout << "EW face area       : " << mesh.eastWestFaceArea() << "\n";
        std::cout << "NS face area       : " << mesh.northSouthFaceArea() << "\n";
        std::cout << "Number of cells    : " << mesh.numberOfCells() << "\n";
    }

    void printCellTable(const Mesh& mesh)
    {
        std::cout << "\n=== Cell Connectivity Table ===========================\n";
        std::cout << std::left
            << std::setw(6) << "idx"
            << std::setw(8) << "(i,j)"
            << std::setw(12) << "center"
            << std::setw(8) << "W"
            << std::setw(8) << "E"
            << std::setw(8) << "S"
            << std::setw(8) << "N"
            << std::setw(18) << "boundaries"
            << "\n";
        std::cout << std::string(76, '-') << "\n";

        for (int idx = 0; idx < mesh.numberOfCells(); ++idx) {
            const Mesh::Cell& c = mesh.getCell(idx);

            std::string boundaryFlags;
            boundaryFlags += c.westBoundary ? "W" : "-";
            boundaryFlags += c.eastBoundary ? "E" : "-";
            boundaryFlags += c.southBoundary ? "S" : "-";
            boundaryFlags += c.northBoundary ? "N" : "-";

            std::ostringstream center;
            center << '(' << std::fixed << std::setprecision(1) << c.x << ',' << c.y << ')';

            std::ostringstream logical;
            logical << '(' << c.i << ',' << c.j << ')';

            std::cout << std::left
                << std::setw(6) << idx
                << std::setw(8) << logical.str()
                << std::setw(12) << center.str()
                << std::setw(8) << c.west
                << std::setw(8) << c.east
                << std::setw(8) << c.south
                << std::setw(8) << c.north
                << std::setw(18) << boundaryFlags
                << "\n";
        }
    }

    void printLogicalLayout(const Mesh& mesh)
    {
        std::cout << "\n=== Logical Layout (top row printed first) ============\n";
        for (int j = mesh.getNy() - 1; j >= 0; --j) {
            for (int i = 0; i < mesh.getNx(); ++i) {
                const int idx = mesh.cellIndex(i, j);
                std::cout << '[' << std::setw(2) << idx << ':' << i << ',' << j << "] ";
            }
            std::cout << "\n";
        }
    }

    void testBasicProperties()
    {
        Mesh mesh(3, 2, 6.0, 4.0);

        expect(mesh.getNx() == 3, "Nx should be 3");
        expect(mesh.getNy() == 2, "Ny should be 2");
        expectNear(mesh.getWidth(), 6.0, 1e-12, "Width should be 6.0");
        expectNear(mesh.getHeight(), 4.0, 1e-12, "Height should be 4.0");
        expectNear(mesh.getDx(), 2.0, 1e-12, "dx should be 2.0");
        expectNear(mesh.getDy(), 2.0, 1e-12, "dy should be 2.0");
        expectNear(mesh.cellVolume(), 4.0, 1e-12, "Cell volume should be 4.0");
        expectNear(mesh.eastWestFaceArea(), 2.0, 1e-12, "East/west face area should be 2.0");
        expectNear(mesh.northSouthFaceArea(), 2.0, 1e-12, "North/south face area should be 2.0");
        expect(mesh.numberOfCells() == 6, "Number of cells should be 6");
    }

    void testCellIndexing()
    {
        Mesh mesh(3, 2, 6.0, 4.0);

        expect(mesh.cellIndex(0, 0) == 0, "cellIndex(0,0) should be 0");
        expect(mesh.cellIndex(1, 0) == 1, "cellIndex(1,0) should be 1");
        expect(mesh.cellIndex(2, 0) == 2, "cellIndex(2,0) should be 2");
        expect(mesh.cellIndex(0, 1) == 3, "cellIndex(0,1) should be 3");
        expect(mesh.cellIndex(1, 1) == 4, "cellIndex(1,1) should be 4");
        expect(mesh.cellIndex(2, 1) == 5, "cellIndex(2,1) should be 5");
    }

    void testCellLogicalIndicesAndCenters()
    {
        Mesh mesh(3, 2, 6.0, 4.0);

        const Mesh::Cell& c00 = mesh.getCell(mesh.cellIndex(0, 0));
        const Mesh::Cell& c10 = mesh.getCell(mesh.cellIndex(1, 0));
        const Mesh::Cell& c21 = mesh.getCell(mesh.cellIndex(2, 1));

        expect(c00.i == 0 && c00.j == 0, "Cell (0,0) logical indices incorrect");
        expect(c10.i == 1 && c10.j == 0, "Cell (1,0) logical indices incorrect");
        expect(c21.i == 2 && c21.j == 1, "Cell (2,1) logical indices incorrect");

        expectNear(c00.x, 1.0, 1e-12, "Cell (0,0) x center");
        expectNear(c00.y, 1.0, 1e-12, "Cell (0,0) y center");

        expectNear(c10.x, 3.0, 1e-12, "Cell (1,0) x center");
        expectNear(c10.y, 1.0, 1e-12, "Cell (1,0) y center");

        expectNear(c21.x, 5.0, 1e-12, "Cell (2,1) x center");
        expectNear(c21.y, 3.0, 1e-12, "Cell (2,1) y center");
    }

    void testInteriorConnectivity()
    {
        Mesh mesh(3, 3, 3.0, 3.0);
        const Mesh::Cell& c = mesh.getCell(mesh.cellIndex(1, 1));

        expect(c.west == mesh.cellIndex(0, 1), "Interior west neighbor incorrect");
        expect(c.east == mesh.cellIndex(2, 1), "Interior east neighbor incorrect");
        expect(c.south == mesh.cellIndex(1, 0), "Interior south neighbor incorrect");
        expect(c.north == mesh.cellIndex(1, 2), "Interior north neighbor incorrect");

        expect(!c.westBoundary, "Interior cell should not be west boundary");
        expect(!c.eastBoundary, "Interior cell should not be east boundary");
        expect(!c.southBoundary, "Interior cell should not be south boundary");
        expect(!c.northBoundary, "Interior cell should not be north boundary");
    }

    void testBoundaryConnectivity()
    {
        Mesh mesh(3, 2, 6.0, 4.0);

        const Mesh::Cell& leftBottom = mesh.getCell(mesh.cellIndex(0, 0));
        expect(leftBottom.west == -1, "Left boundary west neighbor should be -1");
        expect(leftBottom.south == -1, "Bottom boundary south neighbor should be -1");
        expect(leftBottom.east == mesh.cellIndex(1, 0), "Left-bottom east neighbor incorrect");
        expect(leftBottom.north == mesh.cellIndex(0, 1), "Left-bottom north neighbor incorrect");
        expect(leftBottom.westBoundary, "Left-bottom should be west boundary");
        expect(leftBottom.southBoundary, "Left-bottom should be south boundary");

        const Mesh::Cell& rightTop = mesh.getCell(mesh.cellIndex(2, 1));
        expect(rightTop.east == -1, "Right boundary east neighbor should be -1");
        expect(rightTop.north == -1, "Top boundary north neighbor should be -1");
        expect(rightTop.west == mesh.cellIndex(1, 1), "Right-top west neighbor incorrect");
        expect(rightTop.south == mesh.cellIndex(2, 0), "Right-top south neighbor incorrect");
        expect(rightTop.eastBoundary, "Right-top should be east boundary");
        expect(rightTop.northBoundary, "Right-top should be north boundary");
    }

    void testInvalidConstruction()
    {
        bool threw = false;
        try {
            Mesh mesh(0, 2, 1.0, 1.0);
        }
        catch (const std::invalid_argument&) {
            threw = true;
        }
        expect(threw, "Constructing mesh with nx=0 should throw");

        threw = false;
        try {
            Mesh mesh(2, -1, 1.0, 1.0);
        }
        catch (const std::invalid_argument&) {
            threw = true;
        }
        expect(threw, "Constructing mesh with ny<0 should throw");

        threw = false;
        try {
            Mesh mesh(2, 2, 0.0, 1.0);
        }
        catch (const std::invalid_argument&) {
            threw = true;
        }
        expect(threw, "Constructing mesh with zero width should throw");
    }

    void testInvalidAccess()
    {
        Mesh mesh(3, 2, 6.0, 4.0);

        bool threw = false;
        try {
            (void)mesh.getCell(-1);
        }
        catch (const std::out_of_range&) {
            threw = true;
        }
        expect(threw, "getCell(-1) should throw");

        threw = false;
        try {
            (void)mesh.getCell(mesh.numberOfCells());
        }
        catch (const std::out_of_range&) {
            threw = true;
        }
        expect(threw, "getCell(numberOfCells) should throw");

        threw = false;
        try {
            (void)mesh.cellIndex(3, 0);
        }
        catch (const std::out_of_range&) {
            threw = true;
        }
        expect(threw, "cellIndex with out-of-range i should throw");

        threw = false;
        try {
            (void)mesh.cellIndex(0, 2);
        }
        catch (const std::out_of_range&) {
            threw = true;
        }
        expect(threw, "cellIndex with out-of-range j should throw");
    }

    void runVisualDemo()
    {
        const Mesh mesh(4, 3, 8.0, 6.0);
        printMeshSummary(mesh);
        printLogicalLayout(mesh);
        printCellTable(mesh);
    }

} // namespace

int main()
{
    try {
        testBasicProperties();
        testCellIndexing();
        testCellLogicalIndicesAndCenters();
        testInteriorConnectivity();
        testBoundaryConnectivity();
        testInvalidConstruction();
        testInvalidAccess();

        std::cout << "All mesh tests passed.\n";
        runVisualDemo();
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}