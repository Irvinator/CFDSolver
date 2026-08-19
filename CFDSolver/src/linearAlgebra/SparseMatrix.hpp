#pragma once
#include "Vector.hpp"
#include <vector>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <string>

namespace CFD
{

	/**
 * Sparse Matrix in Compressed Sparse Row (CSR) format
 *
 * Memory layout:
 *
 *  values_     → stores every non-zero value
 *  colIndices_ → stores the column index of each non-zero
 *  rowPtr_     → stores where each row starts in values_
 *
 * Example — this 4x4 matrix:
 *
 *  | 4  0  0  2 |
 *  | 0  3  0  0 |
 *  | 0  0  7  1 |
 *  | 0  0  0  5 |
 *
 * Is stored as:
 *  values_     = [ 4, 2, 3, 7, 1, 5 ]
 *  colIndices_ = [ 0, 3, 1, 2, 3, 3 ]
 *  rowPtr_     = [ 0, 2, 3, 5, 6 ]
 */
    class SparseMatrix {
    public:

        // ── Constructors ─────────────────────────────────────────────
        SparseMatrix() = default;

        SparseMatrix(std::size_t rows, std::size_t cols)
            : rows_(rows), cols_(cols), rowPtr_(rows + 1, 0) {
        }

        // ── Build from triplets (COO format → CSR format) ─────────────
        /**
         * Assembles the matrix from three arrays:
         *   rowIdx[k], colIdx[k], vals[k]
         * representing the k-th non-zero entry A(row, col) = val
         *
         * This is how FVM assembly works — you loop over cells,
         * compute coefficients, and add triplets
         */
        void setFromTriplets(
            const std::vector<std::size_t>& rowIdx,
            const std::vector<std::size_t>& colIdx,
            const std::vector<double>& vals)
        {
            // Validate inputs
            if (rowIdx.size() != colIdx.size() ||
                rowIdx.size() != vals.size())
                throw std::runtime_error(
                    "SparseMatrix: triplet arrays must have equal size");

            std::size_t nnz = vals.size();
            values_.resize(nnz);
            colIndices_.resize(nnz);
            rowPtr_.assign(rows_ + 1, 0);

            // Step 1 — count how many entries per row
            for (std::size_t k = 0; k < nnz; ++k) {
                if (rowIdx[k] >= rows_ || colIdx[k] >= cols_)
                    throw std::runtime_error(
                        "SparseMatrix: index out of bounds");
                rowPtr_[rowIdx[k] + 1]++;
            }

            // Step 2 — prefix sum to get row pointers
            for (std::size_t i = 1; i <= rows_; ++i)
                rowPtr_[i] += rowPtr_[i - 1];

            // Step 3 — fill values and column indices
            std::vector<std::size_t> pos(rowPtr_.begin(),
                rowPtr_.end());
            for (std::size_t k = 0; k < nnz; ++k) {
                std::size_t dest = pos[rowIdx[k]]++;
                values_[dest] = vals[k];
                colIndices_[dest] = colIdx[k];
            }
        }

        // ── Matrix-Vector multiply: y = A * x ────────────────────────
        /**
         * This is the most called function in your entire solver.
         * The conjugate gradient solver calls this every iteration.
         *
         * For FVM: this computes A*x where A contains your
         * discretised fluxes and x is your field (T, p, u, v)
         */
        Vector multiply(const Vector& x) const {
            if (x.size() != cols_)
                throw std::runtime_error(
                    "SparseMatrix::multiply — size mismatch");

            Vector y(rows_, 0.0);
            for (std::size_t i = 0; i < rows_; ++i)
                for (std::size_t k = rowPtr_[i];
                    k < rowPtr_[i + 1]; ++k)
                    y[i] += values_[k] * x[colIndices_[k]];
            return y;
        }

        // ── operator* overload (cleaner syntax) ──────────────────────
        Vector operator*(const Vector& x) const {
            return multiply(x);
        }

        // ── Accessors ─────────────────────────────────────────────────
        std::size_t rows() const { return rows_; }
        std::size_t cols() const { return cols_; }
        std::size_t nnz()  const { return values_.size(); }

        const std::vector<double>& values()     const { return values_; }
        const std::vector<std::size_t>& colIndices() const { return colIndices_; }
        const std::vector<std::size_t>& rowPtr()     const { return rowPtr_; }

        // ── Diagonal extraction (useful for preconditioners) ──────────
        Vector diagonal() const {
            Vector d(rows_, 0.0);
            for (std::size_t i = 0; i < rows_; ++i)
                for (std::size_t k = rowPtr_[i];
                    k < rowPtr_[i + 1]; ++k)
                    if (colIndices_[k] == i)
                        d[i] = values_[k];
            return d;
        }

        // ── Print info ────────────────────────────────────────────────
        void printInfo() const {
            std::cout << "SparseMatrix ["
                << rows_ << " x " << cols_ << "]"
                << " | NNZ = " << nnz()
                << " | Density = "
                << (100.0 * nnz()) / (rows_ * cols_)
                << "%\n";
        }

        // ── Print full matrix (only use for small matrices!) ──────────
        void printDense() const {
            std::cout << "\nMatrix [" << rows_
                << "x" << cols_ << "]:\n";
            for (std::size_t i = 0; i < rows_; ++i) {
                for (std::size_t j = 0; j < cols_; ++j) {
                    double val = 0.0;
                    for (std::size_t k = rowPtr_[i];
                        k < rowPtr_[i + 1]; ++k)
                        if (colIndices_[k] == j)
                            val = values_[k];
                    std::cout << val << "\t";
                }
                std::cout << "\n";
            }
        }

    private:
        std::size_t rows_ = 0;
        std::size_t cols_ = 0;

        std::vector<double>      values_;      // Non-zero values
        std::vector<std::size_t> colIndices_;  // Column indices
        std::vector<std::size_t> rowPtr_;      // Row pointers
    };
}