#pragma once

#include "Grid.hpp"
#include "Field2D.hpp"

// Solves ∇²ψ = -ω with periodic BC using SOR iterations.
class PoissonSolver {
public:
    PoissonSolver(const Grid& grid, int maxIterations, double tolerance);

    void solve(const Field2D& omega, Field2D& psi);

private:
    Grid grid_;
    int maxIterations_;
    double tolerance_;
    Field2D buffer_; // no longer used, but harmless to keep
};
