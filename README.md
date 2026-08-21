MACK:
CFDSolver/
└── CFDSolver/
    └── src/
        └── mesh/
            └── Mesh2D.hpp    ← HERE
#pragma once
#include <vector>
#include <stdexcept>
#include <iostream>
#include <string>
#include <cmath>

namespace CFD {

/**
 * ═══════════════════════════════════════════
 * Mesh2D — 2D Structured Cartesian Mesh
 * ═══════════════════════════════════════════
 *
 * Divides domain [xMin,xMax] x [yMin,yMax]
 * into nx * ny rectangular cells
 *
 * Cell layout:
 *
 *  j=2  [0,2][1,2][2,2][3,2]
 *  j=1  [0,1][1,1][2,1][3,1]
 *  j=0  [0,0][1,0][2,0][3,0]
 *        i=0  i=1  i=2  i=3
 *
 * Cell (i,j) has 1D index: i + j*nx
 *
 * Cell centre positions:
 *   x(i) = xMin + (i + 0.5) * dx
 *   y(j) = yMin + (j + 0.5) * dy
 *
 * Face positions:
 *   xFace(i) = xMin + i * dx
 *   yFace(j) = yMin + j * dy
 *
 * Neighbours:
 *   Left  (W): cellIndex(i-1, j)
 *   Right (E): cellIndex(i+1, j)
 *   Bottom(S): cellIndex(i,   j-1)
 *   Top   (N): cellIndex(i,   j+1)
 *   Returns -1 if on boundary
 * ═══════════════════════════════════════════
 */
class Mesh2D {
public:

    // ── Constructor ───────────────────────────────────────────
    Mesh2D(double xMin, double xMax,
           double yMin, double yMax,
           int    nx,   int    ny)
        : xMin_(xMin), xMax_(xMax)
        , yMin_(yMin), yMax_(yMax)
        , nx_(nx), ny_(ny)
    {
        // Validate inputs
        if (nx <= 0 || ny <= 0)
            throw std::runtime_error(
                "Mesh2D: nx and ny must "
                "be positive");
        if (xMax <= xMin)
            throw std::runtime_error(
                "Mesh2D: xMax must be > xMin");
        if (yMax <= yMin)
            throw std::runtime_error(
                "Mesh2D: yMax must be > yMin");

        // Compute cell sizes
        dx_ = (xMax - xMin) / nx;
        dy_ = (yMax - yMin) / ny;

        // Build mesh arrays
        buildMesh();
    }

    // ══════════════════════════════════════════════════════════
    // BASIC INFO
    // ══════════════════════════════════════════════════════════

    int    nx()      const { return nx_; }
    int    ny()      const { return ny_; }
    int    nCells()  const { return nx_ * ny_; }
    int    nFacesX() const { return (nx_+1)*ny_; }
    int    nFacesY() const { return nx_*(ny_+1); }
    int    nFaces()  const {
        return nFacesX() + nFacesY();
    }

    double dx()      const { return dx_; }
    double dy()      const { return dy_; }
    double cellVol() const { return dx_ * dy_; }

    double xMin()    const { return xMin_; }
    double xMax()    const { return xMax_; }
    double yMin()    const { return yMin_; }
    double yMax()    const { return yMax_; }

    double Lx()      const { return xMax_-xMin_; }
    double Ly()      const { return yMax_-yMin_; }

    // ══════════════════════════════════════════════════════════
    // CELL INDEXING
    // ══════════════════════════════════════════════════════════

    /**
     * Convert 2D index (i,j) to 1D index
     * This is how your solver accesses fields:
     * T[mesh.cellIndex(i,j)]
     */
    int cellIndex(int i, int j) const {
        return i + j * nx_;
    }

    /**
     * Convert 1D index back to 2D
     * Useful for loops over all cells
     */
    void cellIJ(int idx, int& i, int& j) const {
        i = idx % nx_;
        j = idx / nx_;
    }

    // ══════════════════════════════════════════════════════════
    // CELL POSITIONS
    // ══════════════════════════════════════════════════════════

    // Cell centre x coordinate for column i
    double cellX(int i) const {
        return xMin_ + (i + 0.5) * dx_;
    }

    // Cell centre y coordinate for row j
    double cellY(int j) const {
        return yMin_ + (j + 0.5) * dy_;
    }

    // Face x position
    double xFace(int i) const {
        return xMin_ + i * dx_;
    }

    // Face y position
    double yFace(int j) const {
        return yMin_ + j * dy_;
    }

    // All cell centre x coords (for output)
    const std::vector<double>& cellCentresX()
        const { return cellCentresX_; }

    // All cell centre y coords (for output)
    const std::vector<double>& cellCentresY()
        const { return cellCentresY_; }

    // ══════════════════════════════════════════════════════════
    // NEIGHBOUR QUERIES
    // Returns -1 if neighbour is outside domain
    // ══════════════════════════════════════════════════════════

    /**
     * Left neighbour (West)
     * Cell to the left of (i,j)
     */
    int leftNeighbour(int i, int j) const {
        return (i > 0) ?
            cellIndex(i-1, j) : -1;
    }

    /**
     * Right neighbour (East)
     * Cell to the right of (i,j)
     */
    int rightNeighbour(int i, int j) const {
        return (i < nx_-1) ?
            cellIndex(i+1, j) : -1;
    }

    /**
     * Bottom neighbour (South)
     * Cell below (i,j)
     */
    int bottomNeighbour(int i, int j) const {
        return (j > 0) ?
            cellIndex(i, j-1) : -1;
    }

    /**
     * Top neighbour (North)
     * Cell above (i,j)
     */
    int topNeighbour(int i, int j) const {
        return (j < ny_-1) ?
            cellIndex(i, j+1) : -1;
    }

    // ══════════════════════════════════════════════════════════
    // BOUNDARY CHECKS
    // ══════════════════════════════════════════════════════════

    bool isLeftBoundary  (int i) const {
        return i == 0;
    }
    bool isRightBoundary (int i) const {
        return i == nx_ - 1;
    }
    bool isBottomBoundary(int j) const {
        return j == 0;
    }
    bool isTopBoundary   (int j) const {
        return j == ny_ - 1;
    }

    // Check if cell (i,j) is ANY boundary
    bool isBoundaryCell(int i, int j) const {
        return isLeftBoundary(i)   ||
               isRightBoundary(i)  ||
               isBottomBoundary(j) ||
               isTopBoundary(j);
    }

    // ══════════════════════════════════════════════════════════
    // DISTANCES (for gradient computation)
    // ══════════════════════════════════════════════════════════

    // Distance between cell centres
    // in x direction (uniform = dx)
    double distX(int i1, int i2) const {
        return std::abs(cellX(i1) - cellX(i2));
    }

    // Distance between cell centres
    // in y direction (uniform = dy)
    double distY(int j1, int j2) const {
        return std::abs(cellY(j1) - cellY(j2));
    }

    // ══════════════════════════════════════════════════════════
    // MESH REFINEMENT (for boundary layers)
    // ══════════════════════════════════════════════════════════

    /**
     * Apply hyperbolic tangent stretching
     * in y direction to cluster cells
     * near bottom wall
     *
     * Used for boundary layer resolution
     * in NS solver!
     *
     * beta = 1.0 → uniform (no stretching)
     * beta = 3.0 → strong wall clustering
     */
    void applyYstretching(double beta = 2.0) {
        if (beta <= 0)
            throw std::runtime_error(
                "Mesh2D: beta must be > 0");

        for (int j = 0; j < ny_; ++j) {
            double eta = (j + 0.5) / ny_;
            double yNew = yMin_ +
                (yMax_ - yMin_) *
                std::tanh(beta * eta) /
                std::tanh(beta);

            for (int i = 0; i < nx_; ++i)
                cellCentresY_[
                    cellIndex(i,j)] = yNew;
        }
        stretched_ = true;
        std::cout << "Y-stretching applied"
                  << " (beta=" << beta << ")\n";
    }

    // ══════════════════════════════════════════════════════════
    // PRINT INFO
    // ══════════════════════════════════════════════════════════

    void printInfo() const {
        std::cout
            << "\n── Mesh2D ─────────────────────"
            << "────────\n"
            << "Domain : ["
            << xMin_ << ", " << xMax_
            << "] x ["
            << yMin_ << ", " << yMax_
            << "]\n"
            << "Grid   : "
            << nx_ << " x " << ny_ << "\n"
            << "Cells  : " << nCells() << "\n"
            << "Faces  : " << nFaces() << "\n"
            << "dx     : " << dx_ << "\n"
            << "dy     : " << dy_ << "\n"
            << "Vol    : " << cellVol() << "\n"
            << "Aspect : "
            << dx_/dy_ << "\n"
            << "Stretched: "
            << (stretched_ ? "yes" : "no")
            << "\n"
            << "────────────────────────────────"
            << "────\n";
    }

private:
    // Domain bounds
    double xMin_, xMax_;
    double yMin_, yMax_;

    // Grid dimensions
    int    nx_, ny_;

    // Cell sizes
    double dx_, dy_;

    // Cell centre positions
    std::vector<double> cellCentresX_;
    std::vector<double> cellCentresY_;

    // Flags
    bool stretched_ = false;

    // ── Build mesh arrays ─────────────────────────────────────
    void buildMesh() {
        int N = nCells();
        cellCentresX_.resize(N);
        cellCentresY_.resize(N);

        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                int idx = cellIndex(i, j);
                cellCentresX_[idx] =
                    xMin_ + (i + 0.5) * dx_;
                cellCentresY_[idx] =
                    yMin_ + (j + 0.5) * dy_;
            }
        }
    }
};

} // namespace CFD




Step 3 — Create the Test File
Right-click tests folder
Add → New Item → C++ File
Name it test_mesh2D.cpp
Paste this:

#include "mesh/Mesh2D.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace CFD;

void test_construction() {
    Mesh2D mesh(0.0, 1.0,
                0.0, 1.0,
                4, 3);

    assert(mesh.nx()     == 4);
    assert(mesh.ny()     == 3);
    assert(mesh.nCells() == 12);
    assert(std::abs(mesh.dx() - 0.25) < 1e-10);
    assert(std::abs(mesh.dy() - 1.0/3.0) < 1e-10);
    std::cout << "PASSED: construction\n";
}

void test_cell_index() {
    Mesh2D mesh(0, 1, 0, 1, 4, 3);

    // Cell (0,0) → index 0
    assert(mesh.cellIndex(0, 0) == 0);
    // Cell (1,0) → index 1
    assert(mesh.cellIndex(1, 0) == 1);
    // Cell (0,1) → index 4
    assert(mesh.cellIndex(0, 1) == 4);
    // Cell (3,2) → index 11
    assert(mesh.cellIndex(3, 2) == 11);

    std::cout << "PASSED: cell indexing\n";
}

void test_cell_positions() {
    Mesh2D mesh(0, 1, 0, 1, 4, 4);
    double dx = 0.25;

    // Cell centres at (i+0.5)*dx
    assert(std::abs(
        mesh.cellX(0) - 0.125) < 1e-10);
    assert(std::abs(
        mesh.cellX(1) - 0.375) < 1e-10);
    assert(std::abs(
        mesh.cellX(3) - 0.875) < 1e-10);
    assert(std::abs(
        mesh.cellY(0) - 0.125) < 1e-10);
    assert(std::abs(
        mesh.cellY(2) - 0.625) < 1e-10);

    std::cout << "PASSED: cell positions\n";
}

void test_neighbours() {
    Mesh2D mesh(0, 1, 0, 1, 4, 4);

    // Interior cell (2,2)
    assert(mesh.leftNeighbour  (2,2) ==
           mesh.cellIndex(1,2));
    assert(mesh.rightNeighbour (2,2) ==
           mesh.cellIndex(3,2));
    assert(mesh.bottomNeighbour(2,2) ==
           mesh.cellIndex(2,1));
    assert(mesh.topNeighbour   (2,2) ==
           mesh.cellIndex(2,3));

    // Boundary cells return -1
    assert(mesh.leftNeighbour  (0,2) == -1);
    assert(mesh.rightNeighbour (3,2) == -1);
    assert(mesh.bottomNeighbour(2,0) == -1);
    assert(mesh.topNeighbour   (2,3) == -1);

    std::cout << "PASSED: neighbours\n";
}

void test_boundaries() {
    Mesh2D mesh(0, 1, 0, 1, 4, 4);

    assert( mesh.isLeftBoundary(0));
    assert(!mesh.isLeftBoundary(1));
    assert( mesh.isRightBoundary(3));
    assert(!mesh.isRightBoundary(2));
    assert( mesh.isBottomBoundary(0));
    assert(!mesh.isBottomBoundary(1));
    assert( mesh.isTopBoundary(3));
    assert(!mesh.isTopBoundary(2));

    // Corner cell is boundary
    assert(mesh.isBoundaryCell(0,0));
    assert(mesh.isBoundaryCell(3,3));
    // Interior cell is not
    assert(!mesh.isBoundaryCell(1,1));
    assert(!mesh.isBoundaryCell(2,2));

    std::cout << "PASSED: boundaries\n";
}

void test_non_square() {
    // Non-square domain and grid
    Mesh2D mesh(0, 2, 0, 1, 10, 5);

    assert(mesh.nx()     == 10);
    assert(mesh.ny()     == 5);
    assert(mesh.nCells() == 50);
    assert(std::abs(mesh.dx() - 0.2) < 1e-10);
    assert(std::abs(mesh.dy() - 0.2) < 1e-10);

    std::cout << "PASSED: non-square mesh\n";
}

void test_cell_volume() {
    Mesh2D mesh(0, 1, 0, 1, 5, 5);
    double expected = 0.2 * 0.2;
    assert(std::abs(
        mesh.cellVol() - expected) < 1e-10);
    std::cout << "PASSED: cell volume\n";
}

void test_large_mesh() {
    // Test a realistic solver mesh
    Mesh2D mesh(0, 1, 0, 1, 100, 100);
    assert(mesh.nCells() == 10000);

    // Check corner cells
    assert(mesh.leftNeighbour(0, 50)   == -1);
    assert(mesh.rightNeighbour(99, 50) == -1);
    assert(mesh.bottomNeighbour(50, 0) == -1);
    assert(mesh.topNeighbour(50, 99)   == -1);

    // Check interior cell
    int P = mesh.cellIndex(50, 50);
    int E = mesh.rightNeighbour(50, 50);
    int W = mesh.leftNeighbour(50, 50);
    assert(E == mesh.cellIndex(51, 50));
    assert(W == mesh.cellIndex(49, 50));

    mesh.printInfo();
    std::cout << "PASSED: large mesh\n";
}

int main() {
    std::cout << "\n== Mesh2D Tests ==\n";

    test_construction();
    test_cell_index();
    test_cell_positions();
    test_neighbours();
    test_boundaries();
    test_non_square();
    test_cell_volume();
    test_large_mesh();

    std::cout << "\n✅ All Mesh2D tests passed!\n";
    return 0;
}



Step 4 — Update CMakeLists.txt
Add the test to your outer CMakeLists.txt:

cmake
# Add this alongside test_vector
add_executable(test_mesh2D
    CFDSolver/tests/test_mesh2D.cpp)
target_link_libraries(test_mesh2D
    PRIVATE CFDCore)
add_test(NAME Mesh2DTest
    COMMAND test_mesh2D)
    
