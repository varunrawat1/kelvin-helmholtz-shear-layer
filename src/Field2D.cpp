#include "Field2D.hpp"

#include <algorithm>

namespace {
    inline int wrapIndex(int i, int n) {
        if (i >= n) {
            i %= n;
        } else if (i < 0) {
            i = (i % n + n) % n;
        }
        return i;
    }
}

Field2D::Field2D(const Grid& grid, double initialValue)
    : grid_(grid), data_(static_cast<std::size_t>(grid.nx * grid.ny), initialValue)
{
}

int Field2D::index(int i, int j) const {
    const int ii = wrapIndex(i, grid_.nx);
    const int jj = wrapIndex(j, grid_.ny);
    return jj * grid_.nx + ii;
}

double& Field2D::operator()(int i, int j) {
    return data_[index(i, j)];
}

const double& Field2D::operator()(int i, int j) const {
    return data_[index(i, j)];
}

void Field2D::fill(double value) {
    std::fill(data_.begin(), data_.end(), value);
}
