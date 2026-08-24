#pragma once

#include <vector>

class Mesh {
public:
    struct Cell {
        int i = -1;
        int j = -1;

        double x = 0.0;
        double y = 0.0;

        int north = -1;
        int east = -1;
        int south = -1;
        int west = -1;

        bool northBoundary = false;
        bool eastBoundary = false;
        bool southBoundary = false;
        bool westBoundary = false;
    };

    Mesh(int nx, int ny, double width, double height);

    void generate();

    int getNx() const;
    int getNy() const;

    double getWidth() const;
    double getHeight() const;

    double getDx() const;
    double getDy() const;

    double cellVolume() const;
    double eastWestFaceArea() const;
    double northSouthFaceArea() const;

    const Cell& getCell(int index) const;
    Cell& getCell(int index);

    int cellIndex(int i, int j) const;
    int numberOfCells() const;

    const std::vector<Cell>& cells() const;

private:
    int nx_ = 0;
    int ny_ = 0;

    double width_ = 0.0;
    double height_ = 0.0;

    double dx_ = 0.0;
    double dy_ = 0.0;

    std::vector<Cell> cells_;
};