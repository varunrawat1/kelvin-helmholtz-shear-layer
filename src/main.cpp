#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "Grid.hpp"
#include "Simulation.hpp"

static std::string makeRunDirectoryName() {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
    localtime_r(&t, &tm);

    std::ostringstream oss;
    oss << "output/run_"
        << (tm.tm_year + 1900) << "-"
        << std::setw(2) << std::setfill('0') << (tm.tm_mon + 1) << "-"
        << std::setw(2) << std::setfill('0') << tm.tm_mday << "_"
        << std::setw(2) << std::setfill('0') << tm.tm_hour
        << std::setw(2) << std::setfill('0') << tm.tm_min
        << std::setw(2) << std::setfill('0') << tm.tm_sec;
    return oss.str();
}

int main() {
    const int    nx = 256;
    const int    ny = 256;
    const double Lx = 1.0;
    const double Ly = 1.0;

    Grid grid(nx, ny, Lx, Ly);

    SimulationParameters params;
    params.viscosity   = 1.0e-6;
    params.dt          = 5.0e-4;
    params.totalSteps  = 8000;
    params.outputEvery = 100;

    const std::string outputDir = makeRunDirectoryName();

    Simulation sim(grid, params, outputDir);
    sim.run();

    return 0;
}
