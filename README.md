# Kelvin–Helmholtz Shear-Layer Solver (2-D Incompressible)

This repository contains a C++ code for the 2-D incompressible Navier–Stokes equations in vorticity–streamfunction form, targeting the classic double shear-layer Kelvin–Helmholtz instability on a periodic box.

The goal is to produce animations of vortex roll-up and pairing while keeping the numerics and code structure transparent and extensible (parallel Poisson, higher-order advection, GPU, etc.).

---
## 1. Mathematical model

We solve the 2-D incompressible Navier–Stokes equations in vorticity–streamfunction form on a periodic square domain $[0,L_x]\times[0,L_y]$:

- **Streamfunction–vorticity relation (Poisson)**  

  ```math
  \nabla^2 \psi = -\,\omega


* **Vorticity transport**

  ```math
  \frac{\partial \omega}{\partial t}
  + u \frac{\partial \omega}{\partial x}
  + v \frac{\partial \omega}{\partial y}
  = \nu \nabla^2 \omega
  ```

* **Velocity from streamfunction**

  ```math
  u = \frac{\partial \psi}{\partial y}, \qquad
  v = -\,\frac{\partial \psi}{\partial x}
  ```

### Double shear-layer initial condition

The initial velocity field is a **double shear layer**:

* Two counter-flowing $\tanh$ shear layers at $y/L_y \approx 0.25$ and $0.75$.
* A small sinusoidal vertical perturbation localised near each layer to trigger Kelvin–Helmholtz roll-up.

From $(u(x,y,0), v(x,y,0))$ we compute the initial vorticity
$\omega(x,y,0) = \partial_x v - \partial_y u$
and solve the Poisson equation to obtain a consistent initial $\psi(x,y,0)$.



---

## 2. Code structure

- `CMakeLists.txt` – CMake build configuration.
- `include/`
  - `Grid.hpp` – grid dimensions and spacings.
  - `Field2D.hpp` – periodic 2-D scalar field wrapper.
  - `PoissonSolver.hpp` – SOR Poisson solver for \(\nabla^2\psi=-\omega\).
  - `Simulation.hpp` – main time-integration driver (initial condition, RK2 loop, diagnostics, output).
- `src/`
  - `Field2D.cpp` – implementation of the 2-D field and periodic indexing.
  - `PoissonSolver.cpp` – in-place SOR iteration.
  - `Simulation.cpp` – vorticity–streamfunction solver, double shear-layer setup, diagnostics and file output.
  - `main.cpp` – chooses grid/parameters, creates a timestamped run directory, and calls `Simulation::run()`.
- `scripts/`
  - `plot_latest_run.py` – converts vorticity `.dat` snapshots from the latest run into PNG images using Matplotlib.
  - `make_gif_latest_run.sh` – assembles the PNG sequence into an animated GIF using ImageMagick.
- `run.sh` – convenience script to configure, build, and run the simulation in one step.
- `.gitignore` – excludes `build/` and `output/` from version control.

Each simulation run is written into a dedicated directory:

```text
output/run_YYYY-MM-DD_HHMMSS/
  run_info.txt
  vorticity_000000.dat
  vorticity_000100.dat
  ...
  vorticity_008000.dat
  vorticity_*.png        # after plotting
  kh_vorticity.gif       # after GIF generation
````

---

## 3. Building and running

### Prerequisites

* A C++17 compiler (e.g. `g++`),
* CMake ≥ 3.16,
* Python 3 with Matplotlib (`python3-matplotlib`),
* ImageMagick (`convert`) for GIF generation.

On Debian/Ubuntu:

```bash
sudo apt-get update
sudo apt-get install g++ cmake python3 python3-matplotlib imagemagick
```

### Build and run (standard)

From the project root:

```bash
mkdir -p build
cd build
cmake ..
cmake --build . -j
cd ..
./build/kh_sim
```

### Build and run (one-liner convenience)

Alternatively, use the provided script:

```bash
./run.sh
```

which configures CMake (if needed), builds the code, and runs `kh_sim`.

---

## 4. Post-processing: plots and “swirl gallery” GIF

After a simulation has finished:

1. **Generate PNG frames from the latest run**

   ```bash
   python3 scripts/plot_latest_run.py
   ```

   This finds the newest `output/run_*` directory and writes `vorticity_*.png` files in it.

2. **Assemble PNGs into a GIF**

   ```bash
   ./scripts/make_gif_latest_run.sh
   ```

   This uses ImageMagick’s `convert` to create `kh_vorticity.gif` in the same run directory.

---

