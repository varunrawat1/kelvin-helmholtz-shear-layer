# Kelvin–Helmholtz Shear-Layer Solver (2-D incompressible)

A C++ code for the 2-D incompressible Navier–Stokes equations in vorticity–streamfunction form, targeting the double shear-layer Kelvin–Helmholtz instability on a periodic box.

- Equations: \nabla^2 \psi = -\omega,\ \partial_t \omega + u \partial_x \omega + v \partial_y \omega = \nu \nabla^2 \omega,\ with u = \partial_y \psi,\ v = -\partial_x \psi.
- Numerics: finite differences on a uniform grid, explicit RK2 (Heun) in time, FFT- or SOR-based Poisson solver for the streamfunction, and periodic boundary conditions in both directions.
- Use case: generate swirl gallery animations of KH roll-up and vortex pairing; benchmark advection/elliptic solvers and CFL behaviour.



## Build & run
```bash
mkdir -p build
cd build
cmake ..
cmake --build . -j
cd ..
./build/kh_sim      # writes vorticity snapshots into output/run_.../

## Post-processing
#From the project root:

python3 scripts/plot_latest_run.py
./scripts/make_gif_latest_run.sh

#This produces a GIF of the vorticity field in the latest output/run_* directory.
