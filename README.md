# COVOF 2D Fortran 2018 — High-Performance Multi-Phase CFD Solver

[![Standard](https://img.shields.io/badge/Fortran-2018-blue.svg)](https://en.wikipedia.org/wiki/Fortran)
[![MPI](https://img.shields.io/badge/MPI-OpenMPI%20%2F%20MPICH-green.svg)](https://www.open-mpi.org/)
[![License](https://img.shields.io/badge/License-Proprietary-red.svg)]()
[![Dual-Code](https://img.shields.io/badge/Dual--Code%20Parity-C%2B%2B17-success.svg)]()

**COVOF 2D (Fortran 2018)** is a modern, modular, massively parallel two-dimensional Computational Fluid Dynamics (CFD) code written in pure Fortran 2018 with OpenMPI point-to-point communications. Dedicated to direct numerical simulation (DNS) of multi-phase flows with extreme density and viscosity ratios, it couples high-order interface tracking (VOF) with both monolithic and fractional-step Navier-Stokes formulations, Krylov linear solvers, and automated VTK parallel output.

COVOF 2D Fortran 2018 is developed in strict mathematical and algorithmic parity with its sister implementation **COVOF 2D C++17**.

---

## Key Features and Capabilities

### 1. Navier-Stokes Formulation & Velocity-Pressure Coupling
- **Governing Equations**: Incompressible and low-Mach multi-phase Navier-Stokes with surface tension:
  $$\nabla \cdot \mathbf{u} = 0$$
  $$\rho(C)\left(\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u}\right) = -\nabla p + \rho(C)\mathbf{g} + \nabla \cdot \left[\mu(C)\left(\nabla \mathbf{u} + \nabla \mathbf{u}^T\right)\right] + \mathbf{f}_\sigma$$
- **Velocity-Pressure Solvers**:
  - `MONOLITHIC` (`FC`): All-coupled saddle-point velocity-pressure formulation solved simultaneously.
  - `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (`AAL`): Spatially localized divergence penalization $r_{i,j}(t) = K_{\mathrm{aal}}\,\mu_2\,|\nabla C|$ with Uzawa projection for robust convergence under stiff viscosity contrasts ($\mu_2/\mu_1 > 10^5$).
  - `STANDARD_PROJECTION` (`SP`): Incremental fractional-step projection method with variable-density Poisson solve.
- **Temporal Discretization**: Second-order backward differentiation (Gear / BDF2) and first-order Euler, with Adams-Bashforth convective linearization.
- **Convective Formulations**: Centered second-order, upwind, and TVD schemes.

### 2. High-Fidelity Interface Tracking (VOF)
- **Geometric Advection**:
  - `WEYMOUTH_YUE`: Unsplit Eulerian-Lagrangian flux integration preserving volume conservation in solenoidal fields.
  - `CIAM`: Direction-split Classical Interface Advection Method.
- **Geometric Reconstruction**:
  - `PLIC`: Piecewise Linear Interface Calculation with analytic normal evaluation.
  - `ELVIRA`: Efficient Least-Squares Volume-of-Fluid Interface Reconstruction.
  - `LVIRA`: Least-Squares Volume-of-Fluid Interface Reconstruction.
  - `SLIC`: Simple Line Interface Calculation.
- **Algebraic TVD VOF**: High-resolution TVD limiters (Superbee, Van Leer, MUSCL, HRIC).
- **Smoothed VOF (sVOF)**: PDE-based normalized indicator smoothing for continuous fluid property transitions across thin interfaces.
- **Surface Tension**: Continuum Surface Force (CSF) model with Height-Function (HF) curvature calculation.

### 3. Linear Algebra & Preconditioning
- **Sparse Matrix Representation**: Compressed Sparse Row (CSR) storage format.
- **Iterative Solvers**: BiCGSTAB, BiCGSTAB(2), GMRES, and Conjugate Gradient (CG).
- **Preconditioners**:
  - `BLOCK_TRIANGULAR`: Block Schur-complement preconditioner for monolithic systems.
  - `ILU0`: Incomplete LU factorization with zero fill-in.
  - `SSOR`: Symmetric Successive Over-Relaxation.
  - `JACOBI`: Diagonal preconditioner.

### 4. Parallelization & I/O
- **Parallel Paradigm**: 2D Cartesian grid domain decomposition via OpenMPI / MPICH.
- **Ghost-Cell Synchronization**: Non-blocking asynchronous communications (`MPI_Irecv` / `MPI_Isend`).
- **Data Export**: Serial and parallel VTK XML Rectilinear Grid format (`.vtr` / `.pvtr`) and CSV history monitoring.

---

## Directory Structure

```text
COVOF_2D-FORTRAN/
├── CMakeLists.txt                      # CMake build configuration
├── .gitignore                          # Git ignore rules
└── src/
    ├── module_precision.f90            # Precision definitions (wp = real64)
    ├── module_constants.f90            # Physical and mathematical constants
    ├── module_eos.f90                  # Multi-phase equations of state and properties
    ├── module_json.f90                 # High-performance JSON parser
    ├── module_config.f90               # Configuration and parameter data structures
    ├── module_mesh.f90                 # 2D Cartesian staggered grid and metrics
    ├── module_fields.f90               # 2D scalar and vector staggered field arrays
    ├── module_mpi_comm.f90             # MPI topology, domain splitting, and halo exchanges
    ├── module_interface_reconstruction.f90 # PLIC, ELVIRA, LVIRA, SLIC reconstruction
    ├── module_interface_vof.f90        # Weymouth-Yue, CIAM, and algebraic TVD advection
    ├── module_navier_stokes.f90        # Navier-Stokes assembly and projection solvers
    ├── module_linear_solvers.f90       # CSR matrix operations, BiCGSTAB2, GMRES, preconds
    ├── module_vtk_writer.f90           # Parallel XML VTK (.vtr / .pvtr) export
    ├── module_simulation.f90           # Main simulation loop and physical time stepping
    └── covof2d_main.f90                # Program entry point
```

---

## Build and Installation

### Prerequisites
- **Fortran Compiler**: Modern Fortran 2018 compliant compiler (`gfortran` $\ge 9.0$ or Intel `ifort`/`ifx`).
- **CMake**: Version $\ge 3.18$.
- **MPI**: OpenMPI $\ge 4.0$ or MPICH $\ge 3.3$.
- **BLAS / LAPACK**: Apple Accelerate Framework (macOS) or OpenBLAS / Intel MKL (Linux).

### Compilation Steps

```bash
# Clone the repository
git clone https://github.com/EL-I-Thiam/COVOF_2D-FORTRAN.git
cd COVOF_2D-FORTRAN

# Create build directory
mkdir -p build && cd build

# Configure and compile
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
```

The compiled binary is generated at `bin/covof2d_fortran`.

---

## Running a Simulation

### Serial Execution
```bash
./bin/covof2d_fortran path/to/case.json
```

### Parallel Execution (OpenMPI)
```bash
mpirun -np 4 ./bin/covof2d_fortran path/to/case.json
```

---

## Configuration File Structure (`case.json`)

```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 160,
    "ny": 160,
    "lx": 1.0,
    "ly": 1.0
  },
  "physics": {
    "enable_two_phase": true,
    "fluid1": { "rho": 1.1768, "mu": 0.00185 },
    "fluid2": { "rho": 1800.0, "mu": 500.0 },
    "gx": 0.0,
    "gy": -9.81,
    "sigma": 0.0,
    "enable_surface_tension": false,
    "convective_scheme": "CENTERED",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_svof_smoothing": true,
    "svof_filter_width": 1.0
  },
  "time": {
    "dt": 0.0005,
    "t_final": 2.0,
    "time_scheme": "GEAR",
    "adaptive_dt": false,
    "cfl_max": 0.3
  },
  "initial_condition": {
    "type": "VISCOUS_JET",
    "xc": 0.5,
    "yc": 1.0,
    "r0": 0.04,
    "u0": 0.0,
    "v0": -0.6
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP" },
    "right":  { "type": "FREE_SLIP" },
    "bottom": { "type": "FREE_SLIP" },
    "top": {
      "type": "DIRICHLET",
      "profile": "JET_INLET",
      "jet_center": 0.5,
      "jet_width": 0.08,
      "u": 0.0,
      "v": -0.6
    }
  },
  "solver": {
    "type": "MONOLITHIC",
    "r_al": 0.0,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-06,
    "max_iter": 1000,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "SUPERBEE",
    "cfl_vof": 0.3
  },
  "output": {
    "output_dir": "output/jet_buckling",
    "case_name": "jet_buckling",
    "vtk_interval": 100,
    "log_interval": 20,
    "write_csv_history": true
  }
}
```

---

## Validation Benchmark Suite

COVOF 2D Fortran 2018 is verified against analytical solutions and published DNS data:

1. **Viscous Jet Buckling Instability**: Periodic left/right folding and air entrapment under low Reynolds numbers ($Re = 0.1728 < 0.52$, $H/D > 3\pi$) with exact flow-rate compensation (`JET_INLET`).
2. **Falling Dense Cylinder**: Free fall of a heavy cylinder ($\rho_2/\rho_1 = 1000$) through ambient gas with exact terminal settling rate.
3. **Liquid Sheet Deformation**: Incompressible liquid sheet dynamics under strong density contrasts ($\rho_2/\rho_1 \in \{100, 1000\}$).
4. **Vortex Shear Advection**: PLIC and TVD interface reconstruction under single-vortex shear field and time reversal.
5. **Static Laplace Droplet**: Verification of Laplace-Young pressure jump $\Delta p = \sigma / R$ and damping of parasitic currents.
6. **Poiseuille Flow Suite**: Steady, transient (Szymanski), and variable-viscosity channel flows with $L_2$ error norm $< 0.01\,\%$.
7. **Lid-Driven Cavity**: Recirculating cavity flows up to $Re = 1000$.

---

## Dual-Code Parity Certification

COVOF 2D Fortran 2018 achieves exact mathematical parity with **COVOF 2D C++17**:
- Velocity field max difference: $\max |u^{\mathrm{Fort}} - u^{\mathrm{CPP}}| < 10^{-6}\,\mathrm{m/s}$
- Conserved mass relative error: $|V^{\mathrm{Fort}}(t) - V^{\mathrm{CPP}}(t)| / V_0 < 0.005\,\%$
- Divergence constraint: $\|\nabla\cdot\mathbf{u}\|_\infty \le 10^{-7}\,\mathrm{s}^{-1}$ on matching grids

---

## Authors & Maintenance

- **Author**: Elhadji Ibrahima Thiam
- **Platform**: COVOF 2D Multi-Phase CFD Platform
