#include "mesh2D.h"

#include <stdexcept>

Mesh::Mesh(int nx, int ny, double width, double height)
    : nx_(nx), ny_(ny), width_(width), height_(height)
{
    if (nx_ <= 0 || ny_ <= 0) {
        throw std::invalid_argument("Mesh: number of cells must be positive");
    }
    if (width_ <= 0.0 || height_ <= 0.0) {
        throw std::invalid_argument("Mesh: domain dimensions must be positive");
    }

    dx_ = width_ / static_cast<double>(nx_);
    dy_ = height_ / static_cast<double>(ny_);

    generate();
}

void Mesh::generate()
{
    cells_.clear();
    cells_.resize(static_cast<std::size_t>(nx_ * ny_));

    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            const int index = cellIndex(i, j);
            Cell& cell = cells_[static_cast<std::size_t>(index)];

            cell.i = i;
            cell.j = j;

            cell.x = (static_cast<double>(i) + 0.5) * dx_;
            cell.y = (static_cast<double>(j) + 0.5) * dy_;

            cell.west = (i > 0) ? cellIndex(i - 1, j) : -1;
            cell.east = (i < nx_ - 1) ? cellIndex(i + 1, j) : -1;
            cell.south = (j > 0) ? cellIndex(i, j - 1) : -1;
            cell.north = (j < ny_ - 1) ? cellIndex(i, j + 1) : -1;

            cell.westBoundary = (i == 0);
            cell.eastBoundary = (i == nx_ - 1);
            cell.southBoundary = (j == 0);
            cell.northBoundary = (j == ny_ - 1);
        }
    }
}

int Mesh::getNx() const
{
    return nx_;
}

int Mesh::getNy() const
{
    return ny_;
}

double Mesh::getWidth() const
{
    return width_;
}

double Mesh::getHeight() const
{
    return height_;
}

double Mesh::getDx() const
{
    return dx_;
}

double Mesh::getDy() const
{
    return dy_;
}

double Mesh::cellVolume() const
{
    return dx_ * dy_;
}

double Mesh::eastWestFaceArea() const
{
    return dy_;
}

double Mesh::northSouthFaceArea() const
{
    return dx_;
}

const Mesh::Cell& Mesh::getCell(int index) const
{
    if (index < 0 || index >= static_cast<int>(cells_.size())) {
        throw std::out_of_range("Mesh::getCell index out of range");
    }
    return cells_[static_cast<std::size_t>(index)];
}

Mesh::Cell& Mesh::getCell(int index)
{
    if (index < 0 || index >= static_cast<int>(cells_.size())) {
        throw std::out_of_range("Mesh::getCell index out of range");
    }
    return cells_[static_cast<std::size_t>(index)];
}

int Mesh::cellIndex(int i, int j) const
{
    if (i < 0 || i >= nx_ || j < 0 || j >= ny_) {
        throw std::out_of_range("Mesh::cellIndex indices out of range");
    }
    return j * nx_ + i;
}

int Mesh::numberOfCells() const
{
    return static_cast<int>(cells_.size());
}

const std::vector<Mesh::Cell>& Mesh::cells() const
{
    return cells_;
}