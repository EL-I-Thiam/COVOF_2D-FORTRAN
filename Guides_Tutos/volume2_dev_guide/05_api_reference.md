# Référence Exhaustive de l'API Interne

Ce chapitre documente de manière exhaustive les interfaces logicielles internes réellement définies dans les **15 fichiers d'en-tête C++17** (`COVOF_2D-C++/include/covof/*.hpp`) et les **16 fichiers sources Fortran 2018** (`COVOF_2D-FORTRAN/src/*.f90`).

---

## 5.1. Référence Complète de l'API C++17 (Espace de Nommage `covof`)

### 5.1.1. `Types.hpp` — Énumérations et Structures de Données de Configuration

- **Énumérations :**
  - `enum class BCType { DIRICHLET, NEUMANN, FREE_SLIP, PERIODIC };`
  - `enum class InitialConditionType { SINGLE_PHASE, DROPLET, BUBBLE, CYLINDER, DAM_BREAK, SHEAR_LAYER, LIQUID_SHEET, VISCOUS_JET, RAYLEIGH_TAYLOR, VORTEX_SHEAR };`
- **Structures :**
  - `BoundaryCondition` : `type`, `u_val`, `v_val`, `profile`, `u_mean`, `y_step`, `jet_width`, `jet_center`, `T_type`, `T_val`.
  - `Obstacle` : `type` (`"BOX"`/`"CYLINDER"`), `x0`, `y0`, `x1`, `y1`, `xc`, `yc`, `radius`, `permeability`.
  - `DomainBC` : `left`, `right`, `bottom`, `top` (`BoundaryCondition`).
  - `FluidProperties` : `rho`, `mu`, `Cp`, `k_thermal`, `beta`, `eos`, `chi_T`, `gamma`, `r_gas`, `cv`, `p_inf`, `b_nasg`, `q_nasg`.
  - `InitialConditionConfig` : `type`, `u`, `v`, `center_x`, `center_y`, `radius`, `axis_a`, `axis_b`, `box_x0`, `box_y0`, `box_x1`, `box_y1`, `pert_amp`, `pert_k`, `thickness`, `jet_width`, `p_init`.
  - `SolverOptions` : `type`, `linear_solver`, `tol`, `max_iter`, `r_al`, `adaptative_al`, `k_aal`, `uzawa_tol`, `l0_ref`, `t0_ref`, `u0_ref`, `p0_ref`, `uzawa_iter`, `use_precond`, `precond_type`.
  - `OutputConfig` : `output_dir`, `case_name`, `vtk_interval`, `log_interval`, `write_csv_history`.
  - `VofConfig` : `reconstruction_method`, `advection_method`, `tvd_limiter`, `cfl_vof`, `use_markers`, `markers_per_cell`, `marker_advection_scheme`, `marker_curvature`, `marker_normals_plic`, `marker_max_per_cell`, `write_markers_vtk`.

### 5.1.2. `Config.hpp` — Classe `Config`

Stocke tous les paramètres de simulation (`nx`, `ny`, `lx`, `ly`, `mode`, `test_case`, drapeaux physiques, `Ra`, `Pr`, `T_ref`, `beta`, modèles de viscosité, `source_term`, `std::vector<Obstacle> obstacles`, `fluid1`, `fluid2`, `sigma`, `gx`, `gy`, paramètres temporels, `bc`, `init_cond`, `solver`, `output`, `vof`).
- `bool load_from_json(const std::string& filepath);`
- `void print() const;`

### 5.1.3. `EOS.hpp` — Classe `EOS` (*Header-Only*)

- `enum class Type { INCOMPRESSIBLE, LINEAR, IDEAL_GAS, NASG };`
- `static Type parse_type(const std::string& name);`
- `static double compute_chi_T(Type type, const FluidProperties& props, double p, double T, double rho);`
- `static double compute_density(Type type, const FluidProperties& props, double p, double T, double p_ref = 1.01325e5);`
- `static double compute_temperature(Type type, const FluidProperties& props, double e, double v_sq);`
- `static double compute_sound_speed(Type type, const FluidProperties& props, double p, double rho, double T);`

### 5.1.4. `Mesh.hpp` — Classe `Mesh2D`

Encapsule la géométrie et la décomposition MPI 2D :
- **Membres publics :** `nx_global`, `ny_global`, `lx`, `ly`, `dx`, `dy`, `comm`, `rank`, `num_procs`, `dims[2]`, `coords[2]`, `rank_left`, `rank_right`, `rank_bottom`, `rank_top`, `is_boundary_left`, `is_boundary_right`, `is_boundary_bottom`, `is_boundary_top`, `nx_local`, `ny_local`, `i_offset`, `j_offset`, `nx_ghost` (`nx_local + 2`), `ny_ghost` (`ny_local + 2`).
- **Méthodes :**
  - `void init(MPI_Comm parent_comm, int nx_g, int ny_g, double length_x, double length_y, bool periodic_x = false, bool periodic_y = false);`
  - `inline int cell_idx(int i, int j) const;` (`j * nx_ghost + i`)
  - `inline int u_idx(int i, int j) const;` (`j * (nx_ghost + 1) + i`)
  - `inline int v_idx(int i, int j) const;` (`j * nx_ghost + i`)
  - `inline int node_idx(int i, int j) const;` (`j * (nx_ghost + 1) + i`)
  - `inline double xc(int i) const;`, `yc(int j) const;`, `xu(int i) const;`, `xv(int i) const;`, `yv(int j) const;`, `yu(int j) const;`

### 5.1.5. `Field.hpp` — Classe `FlowFields` (*Header-Only*)

Conteneur SoA (`std::vector<double>` / `std::vector<uint8_t>`) de l'ensemble des champs discrets :
- `u`, `v`, `p`, `C`, `T`, `rho`, `rho_u`, `rho_v`, `mu`, `mu_node`, `kappa`, `norm_x`, `norm_y`, `F_sigma_x`, `F_sigma_y`, `div_u`, `r_aal`, `chi_T`, `p_old`, `T_old`, `u_old`, `v_old`, `u_nm1`, `v_nm1`, `rho_old`, `rho_u_old`, `rho_v_old`, `rho_u_nm1`, `rho_v_nm1`, `C_s`, `has_nm1` (`bool`), `energy`, `rho_e`, `rhou_star`, `rhov_star`, `rhoe_star`, `C_tilde`, `marker_nx`, `marker_ny`, `marker_kappa`, `has_marker_normal`, `has_marker_kappa`.
- `void allocate(const Mesh2D& mesh);`

### 5.1.6. `MpiComm.hpp` — Classe `MpiComm` (Méthodes Statiques)

- `static void sync_cell_scalar(const Mesh2D& mesh, std::vector<double>& field);`
- `static void sync_u(const Mesh2D& mesh, std::vector<double>& u);`
- `static void sync_v(const Mesh2D& mesh, std::vector<double>& v);`
- `static void apply_boundary_conditions(const Mesh2D& mesh, const DomainBC& bc, FlowFields& fields, bool homogeneous = false);`
- `static void apply_thermal_boundary_conditions(const Mesh2D& mesh, const DomainBC& bc, FlowFields& fields);`
- `static double global_dot_product(const Mesh2D& mesh, const std::vector<double>& a, const std::vector<double>& b, int size);`
- `static double global_norm2(const Mesh2D& mesh, const std::vector<double>& v, int size);`
- `static double global_max(const Mesh2D& mesh, double local_val);`
- `static double global_sum(const Mesh2D& mesh, double local_val);`

### 5.1.7. `CSRMatrix.hpp` — Classe `CSRMatrix`

- **Membres :** `nrows`, `ncols`, `row_ptr`, `col_idx`, `values`, `diag`.
- **Méthodes :** `init(int r, int c)`, `set_from_triplets(...)`, `matvec(const std::vector<double>& x, std::vector<double>& y) const`, `solve_jacobi(const std::vector<double>& b, std::vector<double>& x, int num_iters = 5) const`.

### 5.1.8. `InterfaceReconstruction.hpp` — Classe `InterfaceReconstruction`

- **Types :** `ReconstructionMethod` (`SLIC`, `PLIC`, `LVIRA`, `ELVIRA`), `AdvectionMethod` (`VOF_TVD`, `CIAM`, `WEYMOUTH_YUE`), `TvdLimiter` (`SUPERBEE`, `VAN_LEER`, `MUSCL`, `HRIC`), `InterfaceSegment` (`nx`, `ny`, `alpha`).
- **Méthodes statiques :** `parse_reconstruction_method`, `parse_advection_method`, `parse_tvd_limiter`, `cut_box_area`, `compute_alpha`, `clip_segment_in_cell`, `reconstruct_slic`, `reconstruct_plic`, `reconstruct_lvira`, `reconstruct_elvira`, `reconstruct`, `evaluate_limiter`.

### 5.1.9. `InterfaceMarkers.hpp` — Structure `Marker` & Classe `InterfaceMarkers`

- **`struct Marker` :** `x`, `y`, `nx`, `ny`, `kappa`, `cell_i`, `cell_j`, `active`.
- **`class InterfaceMarkers` :** membre `std::vector<Marker> cloud;`, méthodes publiques :
  - `void init_from_vof(const Mesh2D& mesh, const VofConfig& vof_cfg, const InitialConditionConfig& ic, const DomainBC& bc, FlowFields& fields);`
  - `void interpolate_velocity_cvi(const Mesh2D& mesh, const FlowFields& fields, double xp, double yp, double& up, double& vp) const;`
  - `void interpolate_velocity_bilinear(const Mesh2D& mesh, const FlowFields& fields, double xp, double yp, double& up, double& vp) const;`
  - `void advect_rk2(const Mesh2D& mesh, const VofConfig& vof_cfg, const DomainBC& bc, const FlowFields& fields, double dt);` *(+ surcharge sans `bc`)*
  - `void migrate_mpi(const Mesh2D& mesh, const VofConfig& vof_cfg, const DomainBC& bc);`
  - `void prune_and_reseed(const Mesh2D& mesh, const VofConfig& vof_cfg, const DomainBC& bc, FlowFields& fields);`
  - `void compute_curvature_from_markers(const Mesh2D& mesh, const VofConfig& vof_cfg, const DomainBC& bc, FlowFields& fields);`
  - `void write_vtk(const Mesh2D& mesh, const std::string& output_dir, const std::string& case_name, int step) const;`

### 5.1.10. `InterfaceTracker.hpp` — Classe `InterfaceTracker`

- **Membres :** `fluid1`, `fluid2`, `sigma`, `vof_cfg`, `interface_segments`, `markers` (`InterfaceMarkers`), `step_count`, `bc`, indicateurs SVOF/diffusion/compressible.
- **Méthodes :** `init(const Config& cfg)`, `smooth_svof_helmholtz`, `smooth_diffusion`, `set_initial_condition`, `advect`, `advect_tvd`, `advect_ciam`, `advect_weymouth_yue`, `update_physical_properties`, `compute_curvature_and_CSF`, `apply_marker_normals_to_segments`.

### 5.1.11. `NavierStokesSolver.hpp` — Classe `NavierStokesSolver`

- **Membres :** `mesh`, `bc`, paramètres physiques/solveurs, `obstacles`, `n_owned_u`, `n_owned_v`, `n_owned_p`, `n_owned_total`, matrices locales `Fuu`, `Fvv`, `SchurLaplace`, vecteurs `diag_u`, `diag_v`, `diag_p`, tampons `tf_scratch`, `div_cell_scratch`.
- **Méthodes principales :** `init`, `apply_system_matrix`, `assemble_rhs`, `assemble_preconditioner_blocks`, `pack_vector`, `unpack_vector`, `compute_adaptative_r_al`, `update_pressure_al`, `solve_standard_projection`, `solve_energy`, `advect_ssp_rk3_momentum`, `solve_total_energy`, `compute_mean_nusselt_hot_wall`, prédicats `is_*_in_obstacle` / `get_*_obstacle_perm`, indexeurs `local_{u,v,p}_dof`.

### 5.1.12. `Preconditioner.hpp` & `LinearSolver.hpp` — Hiérarchies Polymorphes

- **`PreconditionerBase` :** `init(const NavierStokesSolver*)`, `update(const FlowFields&, double dt)`, `apply(const std::vector<double>& z, std::vector<double>& r) const`, `get_name()`.
  Dérivées : `JacobiPreconditioner`, `SSORPreconditioner`, `ILU0Preconditioner`, `BlockTriangularPreconditioner`, `CustomPreconditioner` (avec `set_custom_operator`). Fabrique : `PreconditionerFactory::create`.
- **`LinearSolverBase` :** `init(const SolverOptions&)`, `solve(const NavierStokesSolver&, const PreconditionerBase&, const FlowFields&, double dt, const std::vector<double>& b, std::vector<double>& x) const`, `get_name()`.
  Dérivées : `BiCGSTAB2Solver`, `BiCGSTABSolver`, `CGSolver`, `GMRESSolver`, `MINRESSolver`, `DirectSolver` (`A_fact`, `ipiv_cached`, `reset_cache`), `MumpsSolverBridge`, `PetscSolverBridge`, `HypreSolverBridge`. Fabrique : `SolverFactory::create`.

### 5.1.13. `VTKWriter.hpp` & `Simulation.hpp`

- **`VTKWriter` :** `init(mesh, out_dir, name)`, `write_step(mesh, fields, step, time)`, `write_diagnostics(mesh, fields, step, time, dt)`, `finalize()`.
- **`Simulation` :** possède `config`, `mesh`, `fields`, `tracker`, `ns_solver`, `precond` (`unique_ptr`), `linear_solver` (`unique_ptr`), `vtk_writer`. Méthodes : `init`, `run`, `run_laplacian`, `run_laplacian_cas3_2`, `export_poiseuille_results`, `export_poiseuille_unsteady_results`, `run_vortex_advection`, `compute_time_step`.

---

## 5.2. Référence Complète de l'API Fortran 2018 (16 Fichiers `.f90`)

1. **`module_precision.f90` :** Exporte `dp = selected_real_kind(15, 307)`, `sp`, `i4`, `ip`, `i8`.
2. **`module_constants.f90` :** Exporte `PI`, `EPS_SMALL`, `PERMEABILITY_SOLID_DEFAULT`, `BC_*` (`DIRICHLET`, `NEUMANN`, `PERIODIC`, `FREESLIP`), `THERMAL_*`, `SOLVER_*`, `LINEAR_*`, `PRECOND_*`, `CONV_*`, `VISC_*`, `IC_*`, `OBS_*`, `EOS_*`, `TIME_SCHEME_*`, `CONV_LIN_*`.
3. **`module_eos.f90` :** Fonctions pures `eos_parse_type`, `eos_compute_chi_T`, `eos_compute_density`, `eos_compute_temperature`, `eos_compute_sound_speed`.
4. **`module_json.f90` :** `json_read_file(filename, buffer, ierr)`, `json_get_string`, `json_get_real`, `json_get_int`, `json_get_bool`, `json_get_obstacle`.
5. **`module_config.f90` :** Types `bc_t`, `config_t` ; procédure `load_config_from_json(filename, cfg, ierr)`.
6. **`module_mesh.f90` :** Type `mesh_2d_t` (`nx_global`, `ny_global`, `lx`, `ly`, `dx`, `dy`, `rank`, `size`, `comm`, `cart_comm`, `px`, `py`, `coord_x`, `coord_y`, `neighbor_{west,east,south,north}`, `is_boundary_*`, `nx_local`, `ny_local`, `i_start_global`, `j_start_global`) ; `init_mesh`, `destroy_mesh`, fonctions pures `x_face`, `y_face`, `x_center`, `y_center`, `is_in_obstacle`.
7. **`module_fields.f90` :** Type `flow_fields_t` avec tableaux 2D `(0:nx+1, 0:ny+1)` (`u`, `v`, `p`, `div`, `r_aal`, `C`, `T`, `rho`, `rho_u`, `rho_v`, `mu`, `mu_node`, `fx_sigma`, `fy_sigma`, `curvature`, `solid_mask`, `chi_T`, `*_old`, `*_nm1`, `C_s`, `energy`, `rho_e`, `rhou_star`, `rhov_star`, `rhoe_star`, `C_tilde`, `marker_nx`, `marker_ny`, `marker_kappa`, `has_marker_normal`, `has_marker_kappa`) ; méthodes liées `allocate => fields_allocate_mesh`, `free => fields_free` ; routines `allocate_fields`, `deallocate_fields`, `reset_fields`.
8. **`module_mpi_comm.f90` :** `exchange_halos_scalar` (`sync_cell_scalar`), `exchange_halos_velocity` (`sync_u`, `sync_v`), `apply_velocity_bc` (`apply_boundary_conditions`), `apply_thermal_bc` (`apply_thermal_boundary_conditions`), `apply_vof_bc`, `global_max`, `global_min`, `global_sum`, `global_norm2`, `global_dot_product`.
9. **`module_interface_reconstruction.f90` :** Type `interface_segment_t` (`nx`, `ny`, `alpha`) ; `cut_box_area`, `compute_alpha`, `clip_segment_in_cell`, `reconstruct_slic`, `reconstruct_plic`, `reconstruct_elvira`, `reconstruct_lvira`, `reconstruct_interface`, `evaluate_limiter`.
10. **`module_markers.f90` :** Types `marker_t`, `marker_cloud_t` ; `init_markers_from_vof(cloud, mesh, cfg, fields)`, `interpolate_velocity_cvi(mesh, u, v, xp, yp, scheme, up, vp)` (`pure`), `advect_markers_rk2(cloud, mesh, cfg, fields, dt)`, `migrate_markers_mpi(cloud, mesh, cfg)`, `prune_and_reseed_markers(cloud, mesh, cfg, fields)`, `compute_curvature_from_markers(cloud, mesh, cfg, fields)`, `write_markers_vtk(cloud, mesh, output_dir, case_name, step)`, `free_markers(cloud)`.
11. **`module_interface_vof.f90` :** Type `interface_vof_t` (méthodes liées `init`, `set_initial_condition`, `advect`), `global_marker_cloud`, `set_initial_conditions`, `advect_vof`, `update_physical_properties`, `compute_curvature_and_csf`, `smooth_svof_helmholtz`, `write_vof_markers_vtk`.
12. **`module_navier_stokes.f90` :** Type `navier_stokes_t` (procédures liées `init`, `pack_vector`, `unpack_vector`, `apply_system_matrix`, `assemble_rhs`, `compute_adaptative_r_al`, `update_pressure_al`, `assemble_preconditioner_blocks`, `solve_energy`, `advect_ssp_rk3_momentum`, `solve_total_energy`, `compute_mean_nusselt_hot_wall`), fonctions pures `is_u_face_in_obstacle`, `is_v_face_in_obstacle`.
13. **`module_linear_solvers.f90` :** Interface abstraite `custom_precond_proc` et pointeur `custom_precond_hook` ; `solve_linear_system`, `solve_bicgstab2`, `solve_bicgstab`, `solve_cg`, `solve_gmres`, `solve_minres`, `solve_direct`, `apply_preconditioner`.
14. **`module_vtk_writer.f90` :** Type `vtk_writer_t` (`init => vtk_init`, `write_step => vtk_write_step`, `write_diagnostics => vtk_write_diagnostics`, `finalize => vtk_finalize`).
15. **`module_simulation.f90` :** Type `simulation_t` (`init => sim_init`, `run => sim_run`, `compute_time_step => sim_compute_time_step`, `run_laplacian`, `run_laplacian_cas3_2`, `export_poiseuille_results`, `export_poiseuille_unsteady_results`, `run_vortex_advection`, `solve_standard_projection`).
16. **`covof2d_main.f90` :** Programme principal `covof2d_main` (CLI, MPI Init/Finalize/Abort, lancement `sim%init` et `sim%run`).
