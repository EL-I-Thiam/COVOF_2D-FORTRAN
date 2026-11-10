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
  - `MONOLITHIC` (`FC` / `FULLY_COUPLED`): All-coupled saddle-point velocity-pressure formulation solved simultaneously (incompressible and low-Mach weakly compressible $\chi_T(C)\,\mathrm{d}_t p + \nabla\cdot\mathbf{u} = 0$).
  - `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (`AAL` / `C-AAL`): Spatially and temporally adaptive divergence & acoustic penalization $r_{\mathrm{eff},i,j}(t)$ combining inertial (Gear/BDF2), convective, viscous, capillary, Darcy–Brinkman porous drag ($\mu\Delta V/K$), and low-Mach compressibility scales with Uzawa pressure relaxation under extreme viscosity contrasts ($\mu_2/\mu_1 > 10^9$).
  - `STANDARD_PROJECTION` (`SP`): Incremental fractional-step projection method with variable-density Poisson solve.
- **Temporal Discretization**: Second-order backward differentiation (`GEAR` / BDF2) and first-order `EULER`, with full implicit `NEWTON` convective Jacobian linearization, `PICARD`, and `ADAMS_BASHFORTH` schemes.
- **Darcy–Brinkman Volume Penalization**: Implicit porous and solid immersed obstacle penalization $-\frac{\mu(C)}{K}\mathbf{u}$ ($K \in [10^{-14}, 10^{-3}]\;\mathrm{m}^2$).
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
  - `BLOCK_TRIANGULAR` (`BT_SBP` / `BT-SBP`): Block-Triangular Schur–Brinkman–Compressible preconditioner incorporating Augmented Lagrangian $(r_E+r_W)\Delta y/\Delta x$, Darcy–Brinkman resistance $(\mu/K)\Delta V$, Newton convective Jacobian diagonal, and low-Mach acoustic compressibility $-\frac{\chi_T}{\Delta t}\Delta V$.
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

COVOF 2D Fortran 2018 is verified and validated against an exhaustive suite of **24 numerical benchmarks**, covering the full spectrum of computational fluid dynamics: Poisson/Laplace elliptic operators, transient laminar flows, high-contrast multiphase interfaces, low-Mach compressible gas dynamics, Brinkman immersed-boundary penalization, coupled Boussinesq thermal convection, Darcy-Brinkman porous media, and severe hydrodynamic instabilities.

Each benchmark is structured within a standardized verification protocol (`VALIDATION_TOTALE/validation_cas_<k>_<nom>/`) with dedicated input cases, output logs, automated Python/ParaView post-processing scripts, and compiled validation reports.

### Overview Matrix of the 24 Validation Cases

| N° | Directory Name | Physical Mode | Core Operators, Methods & Solvers Validated | Analytical / Reference Benchmark |
| :-: | :--- | :--- | :--- | :--- |
| **1** | `validation_cas_1_laplacian` | `LAPLACIAN` | 2D Poisson operator ($\Delta\Psi = 1$), 5-point FVM stencil, Dirichlet BCs, solvers (`DIRECT`, `CG`, `BICGSTAB2`, `GMRES`, `MINRES`), preconditioners (`NONE`, `JACOBI`, `SSOR`, `ILU0`, `BLOCK_TRIANGULAR`). | 2D double Fourier series expansion on odd eigenmodes ($p_{L_2} = 2.00$). |
| **2** | `validation_cas_2_neumann` | `LAPLACIAN` | 2D Laplace operator ($\Delta\Psi = 0$) with inhomogeneous Neumann BCs, singular matrix gauge regularization, zero-mean pinning. | Harmonic polynomial $\Psi_{\mathrm{ex}} = x^2 - y^2$, exact truncation cancellation ($\|\tau\| \equiv 0$, machine precision $\sim 10^{-14}$). |
| **3** | `validation_cas_3_poiseuille` | `NAVIER_STOKES` | Steady plane Poiseuille flow, 3-point quadratic wall ghost cells, linear pressure drop BC, formulations (`DIRECT`, `AL`, `AAL`, `MONOLITHIC`, `SP`). | Exact parabolic profile $u(y) = \frac{G}{2\mu}y(H-y)$, machine-level error, incompressibility $\|\nabla\cdot\mathbf{u}\|_\infty < 10^{-12}$. |
| **4** | `validation_cas_4_poiseuille_instat` | `NAVIER_STOKES` | Transient start-up Poiseuille flow from rest, cross-convergence study $\Delta t \times \Delta y$, time stepping comparison (`EULER` 1st-order vs `GEAR`/BDF2 2nd-order). | Szymanski analytical Fourier series solution ($p_{\mathrm{time}}=1.0$ for Euler, $p_{\mathrm{time}}=2.0$ for Gear). |
| **5** | `validation_cas_5_poiseuille_visc_variable` | `NAVIER_STOKES` | Variable-viscosity channel flow $\mu(z) = \mu_0 e^{a z^2 / e^2}$, full viscous stress tensor $\nabla\cdot(2\mu\mathbf{D}(\mathbf{u}))$, harmonic viscosity averaging. | Closed-form integral velocity profile, 2nd-order spatial convergence under steep viscosity gradients. |
| **6** | `validation_cas_6_laplace_droplet` | `NAVIER_STOKES` + CSF | Static droplet in capillary equilibrium ($\Delta P = \sigma / R_0$), balanced-force CSF formulation, parabolic curvature smoothing, parasitic current damping. | Laplace-Young pressure jump $\Delta P = \sigma/R_0$, spurious velocity damping $|\mathbf{u}_{\max}| < 10^{-6}\,\mathrm{m/s}$, $\Delta V/V_0 < 10^{-5}$. |
| **7** | `validation_cas_7_vortex_shear` | `VORTEX_ADVECTION` | Circular disk advection in single-vortex shear flow ($T=8\,\mathrm{s}$), geometric unsplit `WEYMOUTH_YUE`, operator-split `CIAM`, algebraic `VOF_TVD`, reconstructions (`PLIC`, `ELVIRA`). | Reversible theoretical shape recovery at $t=T$, strict mass conservation $|\Delta V|/V_0 < 10^{-12}$, boundedness $C \in [0, 1]$. |
| **8** | `validation_cas_8_dense_cylinder` | `NAVIER_STOKES` Multi-Phase | Free fall of a dense cylinder in air under gravity ($r_\rho \approx 85\,000, r_\mu \approx 5.4\times 10^9$), fine meshes $100\times 200$ and $200\times 400$, Gear/BDF2 + AB2 + SVOF smoothing (`FC`, `AL`, `AAL`, `SP`). | Aerodynamic drag balance and terminal velocity trajectory under extreme near-infinite phase contrast. |
| **9** | `validation_cas_9_50x100_newton_dual_code` | `NAVIER_STOKES` Multi-Phase | Dense cylinder on $50\times 100$ grid: Adaptive Augmented Lagrangian (`AAL`) vs Fully Coupled (`FC`), convective linearization (Newton vs AB2 vs Picard), Euler vs Gear/BDF2, SVOF. | Dual-code cross-benchmark: bit-wise parity $\|\mathbf{u}^{\mathrm{Fort}} - \mathbf{u}^{\mathrm{CPP}}\| < 10^{-11}$, Uzawa convergence speedup. |
| **10** | `validation_cas_10_liquid_sheet` | `NAVIER_STOKES` Multi-Phase | Periodic shearing of an inviscid liquid sheet under large density contrasts ($r_\rho \in \{100, 1000\}$), energy-conserving centered convective scheme (`CENTERED`). | Exact conservation of total kinetic energy $E_k(t)/E_k(0) \approx 1.0$, zero numerical dissipation, interface filament stability. |
| **11** | `validation_cas_11_piston_compressible` | `NAVIER_STOKES` Compressible | Isothermal air compression by a water piston in a closed cavity, low-Mach single-field model with acoustic pressure correction $\frac{\chi_T}{\Delta t}(p^{n+1}-p^n)$ (`MONOLITHIC`, `DIRECT`, `BICGSTAB2`). | Analytical thermodynamic compression trajectory $\rho_{\mathrm{air}}(t) = \rho_0 / (1 - V_0 t / H_0)$, total trapped gas mass conservation. |
| **12** | `validation_cas_12_oscillation_bulle_capillaire` | `NAVIER_STOKES` + CSF | Damped capillary quadrupole oscillations ($n=2$) of an elliptical gas bubble in liquid metal ($r_\rho \approx 5950$), balanced-force surface tension, CSF curvature smoothing. | Lamb & Prosperetti dispersion relation ($f_{\mathrm{th}} = 5.70632\,\mathrm{Hz}$, error $< 0.5\%$), exponential viscous damping rate verification. |
| **13** | `validation_cas_13_compression_bulle` | `NAVIER_STOKES` Compressible | Isotropic centripetal compression of a gas bubble by 4-wall normal inflow ($V_0 = 2.5\times 10^{-3}\,\mathrm{m/s}$), low-Mach single-field formulation with distributed acoustic source. | Analytical bubble radius and density evolution $\rho_{\mathrm{air}}(t) = \rho_0 / (1 - q t / (\pi R_0^2))$, collapse time $t_{\mathrm{collapse}} = 2.8274\,\mathrm{s}$. |
| **14** | `validation_cas_14_marche_descendante` | `NAVIER_STOKES` + Brinkman | Isothermal laminar flow over a backward-facing step ($1:2$ expansion, $Re \in [100, 1000]$), Brinkman volume penalization ($K = 10^{-14}\,\mathrm{m^2}$), inflow `POISEUILLE_STEP`. | Armaly et al. & Gartling benchmark data for primary reattachment $X_r/H$ and upper-wall secondary eddy $(X_s/H, X_{rs}/H)$. |
| **15** | `validation_cas_15_obstacle_fin` | `NAVIER_STOKES` + Brinkman | Confined channel flow past a thin vertical rib obstacle on bottom wall ($Re = 82.5$), Brinkman penalization of singular razor-thin solid geometries. | Reference recirculation eddy lengths: upstream corner $x_1/S \approx 4.99$, rib roof $x_2/S \approx 4.01$, downstream $x_3/S \approx 9.19$. |
| **16** | `validation_cas_16_sillage_cylindre` | `NAVIER_STOKES` + Brinkman | Steady symmetric wake flow past an immersed circular cylinder (`CYLINDER`, $Re = 20$), smooth Brinkman solid volume mask. | Coutanceau & Bouard / Tritton benchmark data: wake length $L_r/D \approx 0.87$, separation angle $\theta_s \approx 126.5^\circ$, drag $C_d \approx 2.05$. |
| **17** | `validation_cas_17_cavite_entrainee` | `NAVIER_STOKES` | 2D lid-driven square cavity flow at $Re = 1000$, staggered MAC scheme, convective discretizations (`UPWIND`, `CENTERED`, `TVD_SUPERBEE`), Uzawa and block-triangular preconditioning. | Ghia, Ghia & Shin (1982) benchmark: primary vortex $\psi_{\max} = 0.1189$, corner eddies $\psi_{\min}$, centerline velocity profiles $u(0.5, y)$ and $v(x, 0.5)$. |
| **18** | `validation_cas_18_cavite_chauffage_differentiel` | `NAVIER_STOKES` + Energy | Boussinesq laminar natural convection in a differentially heated cavity ($Pr = 0.71$, $Ra \in [10^3, 10^6]$), staggered scalar transport, coupled thermal buoyancy. | De Vahl Davis (1983) benchmark: average Nusselt numbers $\overline{Nu} \in \{1.118, 2.245, 4.522, 8.825\}$ with 2nd-order spatial accuracy. |
| **19** | `validation_cas_19_convection_poreux_darcy_brinkman` | `NAVIER_STOKES` + Porous | Natural convection in a Darcy-Brinkman saturated porous cavity ($K = 10^{-8}\,\mathrm{m^2}$, $Da = 10^{-6}$, $Ra_f = 10^7$, $Ra^* = 10$, $Pr = 1$), coupled thermal-porous equations. | Analytical boundary layer asymptotic limit $\overline{Nu}_\infty = 1.072157$, 2nd-order convergence of heat transfer rate. |
| **20** | `validation_cas_20_rupture_barrage` | `NAVIER_STOKES` Multi-Phase | 2D dam-break water column collapse over dry bed, sharp density contrast ($\rho_w/\rho_a \approx 1000$), VOF interface tracking with hydrostatic pressure recovery. | Ritter / Saint-Venant analytical rarefaction wave, surge front velocity $x_{\mathrm{front}}(t)$, water height profiles at $t = 0.24\,\mathrm{s}$ and $0.40\,\mathrm{s}$. |
| **21** | `validation_cas_21_ascension_bulle_2d` | `NAVIER_STOKES` + CSF | Buoyancy-driven 2D bubble rise under surface tension in viscous liquid, Bond numbers $Bo = 10$ (ellipsoidal) and $Bo = 125$ (skirted dimpled cap). | Hysing et al. (2009) quantitative CFD benchmark: peak rise velocity $V_{\max} \approx 0.2417\,\mathrm{m/s}$, centroid $y_c(t), circularity tracking. |
| **22** | `validation_cas_22_sedimentation_faxen` | `NAVIER_STOKES` Multi-Phase AAL | Confined sedimentation of a rigid cylinder at low Reynolds ($k_x = a/l = 0.125$), high viscosity penalization ($\mu_s/\mu_f = 10^6$), Adaptive Augmented Lagrangian (`AAL`). | Faxén's analytical wall-correction law $U_{c,\mathrm{th}} = -0.0506\,\mathrm{m/s}$ (error $< 0.8\%$), incompressibility within inclusion $\|\nabla\cdot\mathbf{u}\| < 10^{-10}$. |
| **23** | `validation_cas_23_instabilite_interface_gravitaire` | `NAVIER_STOKES` + CSF | Multi-phase Rayleigh-Taylor instability under gravity and interfacial surface tension, modal perturbation, mushroom roll-up vortex dynamics. | Linear stability dispersion relation $\omega_{\mathrm{th}}^2 = [(\rho_1-\rho_2)gk - \sigma k^3]/(\rho_1+\rho_2)$, capillary cutoff wavenumber $k_c$, exponential growth. |
| **24** | `validation_cas_24_viscous_jet_buckling` | `NAVIER_STOKES` Multi-Phase | 2D viscous liquid jet buckling instability under severe fluid contrast ($r_\rho \approx 1530$, $r_\mu \sim 10^5$), `JET_INLET` boundary with exact flow-rate compensation. | Cruickshank & Munson criterion ($Re \le 0.52$, $H/D > 3\pi$), regime transitions across $Re \in \{0.1728, 0.60, 3.60\}$, divergence control $\|\nabla\cdot\mathbf{u}\|_\infty < 10^{-8}$. |

### Detailed Thematic Breakdown

#### 1. Elliptic Operators & Single-Phase Poiseuille Benchmarks
- **Case 1 (`validation_cas_1_laplacian`)**: 2D Poisson operator with homogeneous Dirichlet conditions ($\Delta\Psi = 1$). Compares direct solver against `CG`, `BICGSTAB2`, `GMRES`, and `MINRES`. Achieves theoretical 2nd-order spatial convergence ($p \approx 2.00$) against the analytical double Fourier series.
- **Case 2 (`validation_cas_2_neumann`)**: 2D Laplace operator with inhomogeneous Neumann conditions on all boundaries ($\Delta\Psi = 0$). Validates singular Laplacian gauge pinning and zero-mean projection, achieving exact machine-precision cancellation ($\sim 10^{-14}$) on the harmonic polynomial $\Psi = x^2 - y^2$.
- **Case 3 (`validation_cas_3_poiseuille`)**: Steady plane laminar Poiseuille channel flow ($Re = 0.125$). Validates 3-point quadratic wall ghost cells and velocity-pressure coupling formulations (`DIRECT`, `MONOLITHIC`, `AL`, `AAL`, `SP`), matching the exact parabolic profile down to machine error.
- **Case 4 (`validation_cas_4_poiseuille_instat`)**: Transient start-up Poiseuille flow from rest. Rigorous cross-convergence analysis ($\Delta t \times \Delta y$) validating 1st-order convergence for Implicit Euler and 2nd-order for Gear/BDF2 against the analytical Szymanski Fourier series.
- **Case 5 (`validation_cas_5_poiseuille_visc_variable`)**: Plane Poiseuille flow with exponential viscosity $\mu(z) = \mu_0 e^{a z^2 / e^2}$. Validates the full variable-viscosity tensor discretization and harmonic face averaging, achieving 2nd-order convergence.

#### 2. Interfacial Kinematics, Capillarity & High-Contrast Multiphase Flows
- **Case 6 (`validation_cas_6_laplace_droplet`)**: Static droplet capillary equilibrium. Verifies the balanced-force CSF formulation and parabolic curvature smoothing, eliminating spurious currents ($|\mathbf{u}_{\max}| < 10^{-6}\,\mathrm{m/s}$) and recovering the exact Laplace-Young pressure jump $\Delta P = \sigma / R_0$.
- **Case 7 (`validation_cas_7_vortex_shear`)**: Advection of a circular fluid disk in a reversible single-vortex shear flow ($T=8\,\mathrm{s}$). Benchmarks geometric unsplit `WEYMOUTH_YUE`, split `CIAM`, and algebraic `VOF_TVD` (`SUPERBEE`, `VAN_LEER`, `HRIC`) with `PLIC`/`ELVIRA`, ensuring exact volume conservation ($|\Delta V|/V_0 < 10^{-12}$) and boundedness $C \in [0, 1]$.
- **Case 8 (`validation_cas_8_dense_cylinder`)**: Free fall of a heavy fluid cylinder through air under gravity ($r_\rho \approx 85\,000, r_\mu \approx 5.4\times 10^9$). Tests mesh convergence ($100\times 200$, $200\times 400$), Gear/BDF2, Adams-Bashforth, and SVOF Helmholtz filtering across `FC`, `AL`, `AAL`, and `SP` formulations.
- **Case 9 (`validation_cas_9_50x100_newton_dual_code`)**: Dense cylinder benchmark on $50\times 100$ grid. In-depth comparative assessment of Adaptive Augmented Lagrangian (`AAL`) vs Fully Coupled (`FC`), convective linearization (Newton vs Adams-Bashforth vs Picard), and time integration (Euler vs Gear/BDF2), verifying machine-level dual-code parity.
- **Case 10 (`validation_cas_10_liquid_sheet`)**: Periodic shearing and kinetic energy conservation of an inviscid liquid sheet under extreme density ratios ($r_\rho \in \{100, 1000\}$). Proves total kinetic energy preservation ($E_k(t)/E_k(0) \approx 1.0$) using the energy-conserving centered convective scheme.

#### 3. Low-Mach Compressible Single-Field Dynamics
- **Case 11 (`validation_cas_11_piston_compressible`)**: Isothermal compression of trapped air by an incoming water piston in a closed cavity. Validates the low-Mach single-field compressible formulation and acoustic pressure correction $\frac{\chi_T}{\Delta t}(p^{n+1}-p^n)$, matching the analytical density evolution $\rho_{\mathrm{air}}(t) = \rho_0 / (1 - V_0 t / H_0)$.
- **Case 12 (`validation_cas_12_oscillation_bulle_capillaire`)**: Damped capillary quadrupole oscillations ($n=2$) of an elliptical gas bubble in liquid metal ($r_\rho \approx 5950$). Matches Lamb & Prosperetti's dispersion frequency ($f_{\mathrm{th}} = 5.70632\,\mathrm{Hz}$) within $0.5\%$ and confirms analytical viscous damping rates.
- **Case 13 (`validation_cas_13_compression_bulle`)**: Isotropic centripetal compression of a gas bubble by four-wall uniform liquid inflow. Accurately predicts the analytical collapse trajectory $\rho_{\mathrm{air}}(t) = \rho_0 / (1 - q t / (\pi R_0^2))$ and collapse time $t_{\mathrm{collapse}} = 2.8274\,\mathrm{s}$ without acoustic CFL limitations.

#### 4. Immersed Boundaries, Obstacles & Brinkman Penalization
- **Case 14 (`validation_cas_14_marche_descendante`)**: Laminar backward-facing step flow ($Re \in [100, 1000]$). Solid step modeled via Brinkman porous penalization ($K = 10^{-14}\,\mathrm{m^2}$) with `POISEUILLE_STEP` inflow profile. Accurately reproduces Armaly et al. & Gartling benchmark reattachment lengths $X_r(Re)$ and upper-wall secondary recirculation zones.
- **Case 15 (`validation_cas_15_obstacle_fin`)**: Confined channel flow past a thin vertical rib obstacle ($Re = 82.5$). Successfully captures all three distinct separation and recirculation bubbles (upstream corner $x_1/S \approx 4.99$, rib roof $x_2/S \approx 4.01$, and downstream eddy $x_3/S \approx 9.19$) and velocity profiles $u(y/H)$.
- **Case 16 (`validation_cas_16_sillage_cylindre`)**: Steady laminar cross-flow past an immersed circular cylinder at $Re = 20$ using a circular Brinkman mask (`CYLINDER`). Accurately recovers wake length ($L_r/D \approx 0.87$), separation angle ($\theta_s \approx 126.5^\circ$), and total drag coefficient ($C_d \approx 2.05$).

#### 5. Classical Recirculating Flows, Thermal Convection & Porous Media
- **Case 17 (`validation_cas_17_cavite_entrainee`)**: 2D lid-driven cavity flow at $Re = 1000$. Matches the benchmark centerline profiles $u(0.5, y)$ and $v(x, 0.5)$ of Ghia, Ghia & Shin (1982), primary vortex strength $\psi_{\max} = 0.1189366$, and bottom corner secondary vortices.
- **Case 18 (`validation_cas_18_cavite_chauffage_differentiel`)**: Boussinesq laminar natural convection in a differentially heated cavity ($Pr = 0.71$, $Ra \in [10^3, 10^6]$). Validates the coupled energy equation and thermal buoyancy source term, yielding 2nd-order spatial convergence of mean Nusselt numbers $\overline{Nu}$ against De Vahl Davis (1983).
- **Case 19 (`validation_cas_19_convection_poreux_darcy_brinkman`)**: Natural convection in a Darcy-Brinkman saturated porous cavity ($K = 10^{-8}\,\mathrm{m^2}, Da = 10^{-6}, Ra_f = 10^7$). Confirms 2nd-order convergence towards the analytical asymptotic limit $\overline{Nu}_\infty = 1.072157$.

#### 6. Hydrodynamic Instabilities, Capillary Dynamics & Free-Surface Flows
- **Case 20 (`validation_cas_20_rupture_barrage`)**: 2D dam-break water column collapse over a dry bed. Validates VOF sharp interface advection under gravity, matching Ritter / Saint-Venant analytical rarefaction wave dynamics, front surge velocity, and free-surface profiles at $t = 0.24\,\mathrm{s}$ and $0.40\,\mathrm{s}$.
- **Case 21 (`validation_cas_21_ascension_bulle_2d`)**: Buoyancy-driven 2D bubble rise in a viscous liquid column ($Bo = 10$ ellipsoidal and $Bo = 125$ skirted dimpled cap). Accurately tracks terminal velocity $V_{\max} \approx 0.2417\,\mathrm{m/s}$, centroid position $y_c(t)$, and interface circularity against the quantitative benchmark of Hysing et al. (2009).
- **Case 22 (`validation_cas_22_sedimentation_faxen`)**: Confined low-Reynolds sedimentation of a solid cylinder ($k_x = 0.125$) modeled via viscosity penalization ($\mu_s/\mu_f = 10^6$) and Adaptive Augmented Lagrangian (`AAL`). Matches Faxén's analytical wall-correction terminal velocity $U_{c,\mathrm{th}} = -0.0506\,\mathrm{m/s}$ within $0.8\%$.
- **Case 23 (`validation_cas_23_instabilite_interface_gravitaire`)**: Multi-phase Rayleigh-Taylor instability under gravity and surface tension. Verifies the exponential growth rate $\omega_{\mathrm{th}}$ in the linear regime, the capillary cutoff wavenumber $k_c$, and subsequent non-linear mushroom roll-up and spike formation.
- **Case 24 (`validation_cas_24_viscous_jet_buckling`)**: Viscous liquid jet buckling instability under severe fluid contrast ($r_\rho \approx 1530$, $r_\mu \sim 10^5$). Validates the Cruickshank & Munson buckling criterion ($Re \le 0.52, H/D > 3\pi$), reproduces the physical transition from periodic folding ($Re = 0.1728$) to symmetric spreading ($Re = 0.60$) and inertial stagnation ($Re = 3.60$), and demonstrates velocity-pressure coupling robustness (`FC`, `AAL`, `SP`) with exact divergence control $\|\nabla\cdot\mathbf{u}\|_\infty < 10^{-8}\,\mathrm{s}^{-1}$.

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
