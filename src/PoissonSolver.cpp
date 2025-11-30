#include "PoissonSolver.hpp"

#include <cmath>

PoissonSolver::PoissonSolver(const Grid& grid, int maxIterations, double tolerance)
    : grid_(grid),
      maxIterations_(maxIterations),
      tolerance_(tolerance),
      buffer_(grid, 0.0)
{
}

// In-place SOR for ∇²ψ = -ω, periodic in x and y.
// Discrete equation:
//   (ψ_{i+1,j} + ψ_{i-1,j} + ψ_{i,j+1} + ψ_{i,j-1} - 4 ψ_{i,j}) / dx² = -ω_{i,j}
// Rearranged:
//   ψ_{i,j} = 1/4 [neighbors + dx² ω_{i,j}]
void PoissonSolver::solve(const Field2D& omega, Field2D& psi) {
    const double dx2 = grid_.dx * grid_.dx;
    const double coeff = 0.25;
    const double omegaSOR = 1.7;   // relaxation factor (1 < ω < 2)

    for (int iter = 0; iter < maxIterations_; ++iter) {
        double maxDiff = 0.0;

        for (int j = 0; j < grid_.ny; ++j) {
            for (int i = 0; i < grid_.nx; ++i) {
                const double neighbors =
                    psi(i + 1, j) + psi(i - 1, j) +
                    psi(i, j + 1) + psi(i, j - 1);

                const double newPsi = coeff * (neighbors + dx2 * omega(i, j));
                const double updated = psi(i, j) + omegaSOR * (newPsi - psi(i, j));

                const double diff = std::fabs(updated - psi(i, j));
                if (diff > maxDiff) {
                    maxDiff = diff;
                }

                psi(i, j) = updated;
            }
        }

        if (maxDiff < tolerance_) {
            break;
        }
    }
}
