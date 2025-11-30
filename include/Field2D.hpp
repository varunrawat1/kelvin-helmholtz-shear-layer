#pragma once

#include <vector>
#include "Grid.hpp"

// Simple 2D scalar field with periodic indexing.
class Field2D {
public:
    explicit Field2D(const Grid& grid, double initialValue = 0.0);

    double& operator()(int i, int j);
    const double& operator()(int i, int j) const;

    const Grid& grid() const { return grid_; }

    void fill(double value);

private:
    Grid grid_;
    std::vector<double> data_;

    int index(int i, int j) const;
};
