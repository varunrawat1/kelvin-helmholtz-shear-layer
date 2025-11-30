#include "Simulation.hpp"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

Simulation::Simulation(const Grid& grid,
                       const SimulationParameters& params,
                       const std::string& outputDir)
    : grid_(grid),
      params_(params),
      outputDir_(outputDir),
      omega_(grid, 0.0),
      omegaTmp_(grid, 0.0),
      rhs_(grid, 0.0),
      rhsTmp_(grid, 0.0),
      psi_(grid, 0.0),
      u_(grid, 0.0),
      v_(grid, 0.0),
      poisson_(grid, /*maxIterations*/ 300, /*tolerance*/ 1e-6)
{
    std::error_code ec;
    std::filesystem::create_directories(outputDir_, ec);
    if (ec) {
        std::cerr << "Warning: could not create output directory '"
                  << outputDir_ << "': " << ec.message() << '\n';
    }

    initialiseShearLayer();
    writeRunInfo();
}

void Simulation::initialiseShearLayer() {
    const double U0 = 1.0;
    const double delta = 0.01 * grid_.Ly;
    const double eps   = 0.08;
    const double twoPiOverLx = 2.0 * M_PI / grid_.Lx;

    for (int j = 0; j < grid_.ny; ++j) {
        const double y = (j + 0.5) * grid_.dy;

        for (int i = 0; i < grid_.nx; ++i) {
            const double x = (i + 0.5) * grid_.dx;

            const double y1 = (y - 0.25 * grid_.Ly) / delta;
            const double y2 = (y - 0.75 * grid_.Ly) / delta;

            const double shear =
                std::tanh(y1) - std::tanh(y2) - 1.0;

            const double gauss1 = std::exp(-y1 * y1);
            const double gauss2 = std::exp(-y2 * y2);

            const double u_val = U0 * shear;
            const double v_val = eps * std::sin(twoPiOverLx * x) * (gauss1 + gauss2);

            u_(i, j) = u_val;
            v_(i, j) = v_val;
        }
    }

    const double inv2dx = 1.0 / (2.0 * grid_.dx);
    const double inv2dy = 1.0 / (2.0 * grid_.dy);

    for (int j = 0; j < grid_.ny; ++j) {
        for (int i = 0; i < grid_.nx; ++i) {
            const double dv_dx = (v_(i + 1, j) - v_(i - 1, j)) * inv2dx;
            const double du_dy = (u_(i, j + 1) - u_(i, j - 1)) * inv2dy;
            omega_(i, j) = dv_dx - du_dy;
        }
    }

    poisson_.solve(omega_, psi_);
}

void Simulation::computeVelocity(const Field2D& psi, Field2D& u, Field2D& v) {
    const double inv2dx = 1.0 / (2.0 * grid_.dx);
    const double inv2dy = 1.0 / (2.0 * grid_.dy);

    for (int j = 0; j < grid_.ny; ++j) {
        for (int i = 0; i < grid_.nx; ++i) {
            const double dpsi_dy = (psi(i, j + 1) - psi(i, j - 1)) * inv2dy;
            const double dpsi_dx = (psi(i + 1, j) - psi(i - 1, j)) * inv2dx;

            u(i, j) = dpsi_dy;      // u = ∂ψ/∂y
            v(i, j) = -dpsi_dx;     // v = -∂ψ/∂x
        }
    }
}

void Simulation::computeRHS(const Field2D& omega,
                            const Field2D& u,
                            const Field2D& v,
                            Field2D& rhs,
                            double viscosity)
{
    const double inv2dx = 1.0 / (2.0 * grid_.dx);
    const double inv2dy = 1.0 / (2.0 * grid_.dy);
    const double invDx2 = 1.0 / (grid_.dx * grid_.dx);

    for (int j = 0; j < grid_.ny; ++j) {
        for (int i = 0; i < grid_.nx; ++i) {
            const double domega_dx =
                (omega(i + 1, j) - omega(i - 1, j)) * inv2dx;
            const double domega_dy =
                (omega(i, j + 1) - omega(i, j - 1)) * inv2dy;

            const double lap =
                (omega(i + 1, j) + omega(i - 1, j) +
                 omega(i, j + 1) + omega(i, j - 1) -
                 4.0 * omega(i, j)) * invDx2;

            const double adv = u(i, j) * domega_dx + v(i, j) * domega_dy;

            rhs(i, j) = -adv + viscosity * lap;
        }
    }
}

std::string Simulation::makeFilename(int step) const {
    std::ostringstream oss;
    oss << outputDir_ << "/vorticity_"
        << std::setw(6) << std::setfill('0') << step << ".dat";
    return oss.str();
}

void Simulation::writeVorticityToFile(int step) const {
    const std::string filename = makeFilename(step);
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Warning: could not open " << filename << " for writing\n";
        return;
    }

    for (int j = 0; j < grid_.ny; ++j) {
        for (int i = 0; i < grid_.nx; ++i) {
            out << omega_(i, j);
            if (i + 1 < grid_.nx) {
                out << ' ';
            }
        }
        out << '\n';
    }
}

void Simulation::writeRunInfo() const {
    const std::string filename = outputDir_ + "/run_info.txt";
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Warning: could not open " << filename
                  << " for writing run info\n";
        return;
    }
    out << "# Kelvin–Helmholtz shear-layer run\n";
    out << "nx " << grid_.nx << "\n";
    out << "ny " << grid_.ny << "\n";
    out << "Lx " << grid_.Lx << "\n";
    out << "Ly " << grid_.Ly << "\n";
    out << "dx " << grid_.dx << "\n";
    out << "dy " << grid_.dy << "\n";
    out << "viscosity " << params_.viscosity << "\n";
    out << "dt " << params_.dt << "\n";
    out << "totalSteps " << params_.totalSteps << "\n";
    out << "outputEvery " << params_.outputEvery << "\n";
}

bool Simulation::fieldIsFinite(const Field2D& f) const {
    for (int j = 0; j < grid_.ny; ++j) {
        for (int i = 0; i < grid_.nx; ++i) {
            const double val = f(i, j);
            if (!std::isfinite(val)) {
                return false;
            }
        }
    }
    return true;
}

void Simulation::computeDiagnostics(double& maxAbsOmega,
                                    double& maxSpeed,
                                    double& cflAdv,
                                    double& cflDiff) const
{
    maxAbsOmega = 0.0;
    maxSpeed    = 0.0;

    for (int j = 0; j < grid_.ny; ++j) {
        for (int i = 0; i < grid_.nx; ++i) {
            const double w = omega_(i, j);
            const double absW = std::fabs(w);
            if (absW > maxAbsOmega) {
                maxAbsOmega = absW;
            }

            const double ux = u_(i, j);
            const double vy = v_(i, j);
            const double speed = std::sqrt(ux * ux + vy * vy);
            if (speed > maxSpeed) {
                maxSpeed = speed;
            }
        }
    }

    const double h = std::min(grid_.dx, grid_.dy);
    cflAdv  = maxSpeed * params_.dt / h;
    cflDiff = params_.viscosity * params_.dt / (h * h);
}

void Simulation::run() {
    using clock = std::chrono::steady_clock;
    const auto tStart = clock::now();

    std::cout << "Running KH simulation on "
              << grid_.nx << "x" << grid_.ny
              << " grid, dt = " << params_.dt
              << ", ν = " << params_.viscosity << '\n';
    std::cout << "Output directory: " << outputDir_ << "\n\n";

    const int N = params_.totalSteps;
    const int widthBar = 40;

    for (int step = 0; step <= N; ++step) {
        const double progress = static_cast<double>(step) / static_cast<double>(N);
        const int filled = static_cast<int>(progress * widthBar);

        auto tNow = clock::now();
        double elapsed = std::chrono::duration<double>(tNow - tStart).count();
        double eta = (progress > 0.0) ? elapsed * (1.0 - progress) / progress : 0.0;

        // Progress bar updated each step (single line).
        std::cout << "\r[";
        for (int k = 0; k < widthBar; ++k) {
            std::cout << (k < filled ? '=' : ' ');
        }
        std::cout << "] "
                  << std::setw(6) << std::fixed << std::setprecision(1)
                  << (progress * 100.0) << "% "
                  << "step " << std::setw(6) << step << "/" << N
                  << "  elapsed " << std::setw(6) << std::setprecision(1) << elapsed << "s"
                  << "  ETA "     << std::setw(6) << std::setprecision(1) << eta << "s"
                  << std::flush;

        // Diagnostics + output every outputEvery steps.
        if (step % params_.outputEvery == 0) {
            poisson_.solve(omega_, psi_);
            computeVelocity(psi_, u_, v_);

            double maxAbsOmega, maxSpeed, cflAdv, cflDiff;
            computeDiagnostics(maxAbsOmega, maxSpeed, cflAdv, cflDiff);

            std::cout << "\n"
                      << "  [diag] step " << step
                      << "  max|ω| = " << std::setprecision(3) << maxAbsOmega
                      << "  max|u| = " << maxSpeed
                      << "  CFL_adv = " << cflAdv
                      << "  CFL_diff = " << cflDiff
                      << '\n';

            writeVorticityToFile(step);

            if (!fieldIsFinite(omega_)) {
                std::cout << "  [diag] Non-finite values detected in ω. Aborting.\n";
                break;
            }
        }

        if (step == N) {
            break;
        }

        // Heun's method (RK2) for dω/dt = F(ω).
        poisson_.solve(omega_, psi_);
        computeVelocity(psi_, u_, v_);
        computeRHS(omega_, u_, v_, rhs_, params_.viscosity);

        for (int j = 0; j < grid_.ny; ++j) {
            for (int i = 0; i < grid_.nx; ++i) {
                omegaTmp_(i, j) = omega_(i, j) + params_.dt * rhs_(i, j);
            }
        }

        poisson_.solve(omegaTmp_, psi_);
        computeVelocity(psi_, u_, v_);
        computeRHS(omegaTmp_, u_, v_, rhsTmp_, params_.viscosity);

        for (int j = 0; j < grid_.ny; ++j) {
            for (int i = 0; i < grid_.nx; ++i) {
                const double rhsAvg = 0.5 * (rhs_(i, j) + rhsTmp_(i, j));
                omega_(i, j) += params_.dt * rhsAvg;
            }
        }
    }

    auto tEnd = clock::now();
    double elapsedTotal = std::chrono::duration<double>(tEnd - tStart).count();
    std::cout << "\nSimulation finished in " << elapsedTotal << " seconds.\n";
}
