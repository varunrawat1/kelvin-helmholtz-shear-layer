#pragma once

struct Grid {
    int nx;
    int ny;
    double Lx;
    double Ly;
    double dx;
    double dy;

    Grid(int nx_, int ny_, double Lx_, double Ly_)
        : nx(nx_), ny(ny_), Lx(Lx_), Ly(Ly_),
          dx(Lx_ / static_cast<double>(nx_)),
          dy(Ly_ / static_cast<double>(ny_)) {}
};
