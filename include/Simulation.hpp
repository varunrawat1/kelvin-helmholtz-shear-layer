#pragma once

#include <string>
#include "Grid.hpp"
#include "Field2D.hpp"
#include "PoissonSolver.hpp"

struct SimulationParameters {
    double viscosity;
    double dt;
    int totalSteps;
    int outputEvery;
};

// 2D incompressible Navier–Stokes in vorticity–streamfunction form:
//   ∇²ψ = -ω
//   ∂ω/∂t + u·∇ω = ν ∇²ω
//   u = ∂ψ/∂y, v = -∂ψ/∂x
class Simulation {
public:
    Simulation(const Grid& grid,
               const SimulationParameters& params,
               const std::string& outputDir);

    void run();

private:
    Grid grid_;
    SimulationParameters params_;
    std::string outputDir_;

    Field2D omega_;
    Field2D omegaTmp_;
    Field2D rhs_;
    Field2D rhsTmp_;
    Field2D psi_;
    Field2D u_;
    Field2D v_;

    PoissonSolver poisson_;

    void initialiseShearLayer();
    void computeVelocity(const Field2D& psi, Field2D& u, Field2D& v);
    void computeRHS(const Field2D& omega,
                    const Field2D& u,
                    const Field2D& v,
                    Field2D& rhs,
                    double viscosity);

    void writeVorticityToFile(int step) const;
    std::string makeFilename(int step) const;

    // Diagnostics
    void writeRunInfo() const;
    void computeDiagnostics(double& maxAbsOmega,
                            double& maxSpeed,
                            double& cflAdv,
                            double& cflDiff) const;
    bool fieldIsFinite(const Field2D& f) const;
};
