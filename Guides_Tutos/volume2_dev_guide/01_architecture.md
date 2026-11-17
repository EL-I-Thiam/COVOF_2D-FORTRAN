# Architecture Logicielle Duale C++17 / Fortran 2018

## Vision d'Ensemble et Paradigmes de Conception

Le projet **COVOF 2D** repose sur une architecture logicielle **gémellaire symétrique** (*dual-mirror architecture*) entre deux implémentations parallèles indépendantes et fonctionnellement isomorphes :

1. **La branche C++17 (`COVOF_2D-C++`) :** Constituée de **15 fichiers d'en-tête (`.hpp`)** dans `include/covof/` et **13 unités de compilation (`.cpp`)** dans `src/`, produisant le binaire `bin/covof2d_cpp` (et le lien `bin/covof2d`). Elle adopte une conception orientée objet moderne dans l'espace de nommage `covof` :
   - Encapsulation stricte des entités (`Config`, `Mesh2D`, `FlowFields`, `EOS`, `CSRMatrix`, `InterfaceReconstruction`, `InterfaceMarkers`, `InterfaceTracker`, `NavierStokesSolver`, `VTKWriter`, `Simulation`).
   - Polymorphisme dynamique via les hiérarchies `LinearSolverBase` (`SolverFactory`) et `PreconditionerBase` (`PreconditionerFactory`) gérées par `std::unique_ptr`.
   - Gestion déterministe de la mémoire par RAII (`std::vector<double>`, `std::ofstream`) sans allocation manuelle `new`/`delete`.
2. **La branche Fortran 2018 (`COVOF_2D-FORTRAN`) :** Constituée de **16 fichiers sources (`.f90`)** dans `src/` (15 modules + le programme principal `covof2d_main.f90`), produisant le binaire `bin/covof2d_fortran`. Elle applique un style modulaire strict (`implicit none`, encapsulation `private` par défaut avec liste explicite `public :: ...`), des types dérivés à procédures liées (`type :: simulation_t`, `navier_stokes_t`, `interface_vof_t`, `vtk_writer_t`, `flow_fields_t`), et un stockage 2D contigu `(0:nx+1, 0:ny+1)` optimisé pour l'ordonnancement *column-major*.

Les deux branches fonctionnent exclusivement en **parallélisme distribué MPI pur** (OpenMPI / MPI-3, sans OpenMP), lisent exactement le même format JSON et génèrent des sorties identiques (`.vtr`, `.pvtr`, `_markers_*.vtk`, `_history.csv` et fichiers CSV analytiques).

---

## Cartographie Intégrale de l'Arborescence des Dépôts

### 1. Inventaire Complet de la Branche C++17 (`COVOF_2D-C++`)

#### Fichiers d'En-tête (`include/covof/` — 15 fichiers `.hpp`)

| Fichier `.hpp` | Classes / Structures / Enums Principaux | Responsabilité Logicielle |
| :--- | :--- | :--- |
| `Types.hpp` | `BCType`, `BoundaryCondition`, `Obstacle`, `DomainBC`, `FluidProperties`, `InitialConditionType`, `InitialConditionConfig`, `SolverOptions`, `OutputConfig`, `VofConfig` | Structures de paramètres typées et énumérations partagées. |
| `Config.hpp` | `Config` | Conteneur global des paramètres JSON, méthodes `load_from_json()` et `print()`. |
| `EOS.hpp` | `EOS`, `EOS::Type` | Équations d'état (`INCOMPRESSIBLE`, `LINEAR`, `IDEAL_GAS`, `NASG`) : calcul pur de $\chi_T, \rho(p,T), T(e,\|\mathbf{u}\|^2), c_s$. |
| `Mesh.hpp` | `Mesh2D` | Topologie cartésienne MPI (`dims`, `coords`, rangs cardinaux), métriques MAC et accesseurs `cell_idx`, `u_idx`, `v_idx`, `node_idx`, `xc`, `yc`, `xu`, `yv`. |
| `Field.hpp` | `FlowFields` | Conteneur SoA (*Structure of Arrays*) de tous les champs MAC (`u`, `v`, `p`, `C`, `C_s`, `C_tilde`, `T`, `rho`, `mu`, `mu_node`, `kappa`, `F_sigma_*`, `chi_T`, `r_aal`, `energy`, `rho_e`, `marker_*`, `has_marker_*`). |
| `MpiComm.hpp` | `MpiComm` | Primitives statiques MPI : échange de halos (`sync_cell_scalar`, `sync_u`, `sync_v`), CL hydrodynamiques/thermiques et réductions (`global_dot_product`, `global_norm2`, `global_max`, `global_sum`). |
| `CSRMatrix.hpp` | `CSRMatrix` | Matrice creuse locale CSR (`row_ptr`, `col_idx`, `values`, `diag`) et méthodes `set_from_triplets`, `matvec`, `solve_jacobi`. |
| `InterfaceReconstruction.hpp` | `ReconstructionMethod`, `AdvectionMethod`, `TvdLimiter`, `InterfaceSegment`, `InterfaceReconstruction` | Découpe géométrique (`cut_box_area`, `compute_alpha`, `clip_segment_in_cell`), reconstructions (`SLIC`, `PLIC`, `LVIRA`, `ELVIRA`) et limiteurs TVD. |
| `InterfaceMarkers.hpp` | `Marker`, `InterfaceMarkers` | Nuage de marqueurs lagrangiens non connectés : `init_from_vof`, `interpolate_velocity_cvi`, `interpolate_velocity_bilinear`, `advect_rk2`, `migrate_mpi`, `prune_and_reseed`, `compute_curvature_from_markers`, `write_vtk`. |
| `InterfaceTracker.hpp` | `InterfaceTracker` | Orchestrateur diphasique VOF + Marqueurs : `set_initial_condition`, `advect` (`advect_tvd`, `advect_ciam`, `advect_weymouth_yue`), `smooth_svof_helmholtz`, `smooth_diffusion`, `update_physical_properties`, `compute_curvature_and_CSF`, `apply_marker_normals_to_segments`. |
| `NavierStokesSolver.hpp` | `NavierStokesSolver` | Opérateur *matrix-free* `apply_system_matrix`, `assemble_rhs`, `compute_adaptative_r_al`, `update_pressure_al`, `solve_standard_projection`, `assemble_preconditioner_blocks`, `solve_energy`, `advect_ssp_rk3_momentum`, `solve_total_energy`, Brinkman. |
| `Preconditioner.hpp` | `PreconditionerBase`, `JacobiPreconditioner`, `SSORPreconditioner`, `ILU0Preconditioner`, `BlockTriangularPreconditioner`, `CustomPreconditioner`, `PreconditionerFactory` | Hiérarchie polymorphe des préconditionneurs parallèles (incluant BT-SBP à complément de Schur). |
| `LinearSolver.hpp` | `LinearSolverBase`, `BiCGSTAB2Solver`, `BiCGSTABSolver`, `CGSolver`, `GMRESSolver`, `MINRESSolver`, `DirectSolver`, `MumpsSolverBridge`, `PetscSolverBridge`, `HypreSolverBridge`, `SolverFactory` | Solveurs de Krylov *matrix-free*, solveur direct LAPACK `dgetrf`/`dgetrs` et passerelles d'extension. |
| `VTKWriter.hpp` | `VTKWriter` | Export parallèle `.vtr`/`.pvtr` et écriture du fichier `<case_name>_history.csv`. |
| `Simulation.hpp` | `Simulation` | Chef d'orchestre global : `init`, `run`, `compute_time_step`, `run_laplacian`, `run_laplacian_cas3_2`, `run_vortex_advection`, `export_poiseuille_*`. |

#### Fichiers Sources (`src/` — 13 fichiers `.cpp`)

| Fichier `.cpp` | En-tête Associé | Rôle Algorithmique Principal |
| :--- | :--- | :--- |
| `Config.cpp` | `Config.hpp` | Analyseur JSON natif tolérant `//`, extraction de toutes les sections et `Config::print()`. |
| `Mesh.cpp` | `Mesh.hpp` | Création de la topologie cartésienne `MPI_Cart_create`, calcul de `nx_local`, `ny_local`, `i_offset`, `j_offset` et détection des bords physiques. |
| `MpiComm.cpp` | `MpiComm.hpp` | Échanges non bloquants en 2 passes (X puis Y) des halos scalaires/MAC, application des CL et `MPI_Allreduce`. |
| `CSRMatrix.cpp` | `CSRMatrix.hpp` | Construction CSR triée depuis triplets, SpMV local et solveur Jacobi. |
| `Preconditioner.cpp` | `Preconditioner.hpp` | Mise à jour des diagonales/blocs et application (notamment 8 balayages SSOR sur $S$ + 4 sur $F_u, F_v$ dans `BlockTriangularPreconditioner`). |
| `LinearSolver.cpp` | `LinearSolver.hpp` | Boucles de Krylov distribuées et factorisation LU directe via `dgetrf_`/`dgetrs_`. |
| `InterfaceReconstruction.cpp` | `InterfaceReconstruction.hpp` | Algorithme Sutherland-Hodgman + Shoelace, bisection `compute_alpha`, segment clipping et schémas SLIC/PLIC/ELVIRA/LVIRA. |
| `InterfaceMarkers.cpp` | `InterfaceMarkers.hpp` | CVI Variante B exacte, advection RK2-Heun avec routage MPI 4 passes, migration 2 passes, élagage par sous-octants, régression parabolique $(\tau, \eta)$ avec halo $2.5 h$, export VTK PolyData. |
| `InterfaceTracker.cpp` | `InterfaceTracker.hpp` | Conditions initiales (avec sous-échantillonnage $16\times 16$ si marqueurs/ellipse), lissage Helmholtz SVOF, lissage diffusif, Weymouth-Yue / CIAM / TVD, couplage marqueurs-CSF. |
| `NavierStokesSolver.cpp` | `NavierStokesSolver.hpp` | Pochoirs MAC complets (BDF2, Picard/AB2/Newton, viscosité nodale $\mu_{\text{node}}$, AL/AAL/C-AAL, Monolithique, Brinkman, Projection, SSP-RK3, Énergie). |
| `VTKWriter.cpp` | `VTKWriter.hpp` | Écriture autonome `.vtr`, assemblage `.pvtr` sur rang 0 et calcul des 9 métriques de `_history.csv`. |
| `Simulation.cpp` | `Simulation.hpp` | Boucle temporelle principale, pas adaptatif CFL+capillaire, modes analytiques Laplacien, Poiseuille et Vortex. |
| `main.cpp` | `Simulation.hpp` | Analyse CLI (`--help`), initialisation/clôture MPI et capture des exceptions C++. |

*(Note : `EOS.hpp`, `Field.hpp` et `Types.hpp` sont *header-only* et ne nécessitent pas de `.cpp` dédié.)*

### 2. Inventaire Complet de la Branche Fortran 2018 (`COVOF_2D-FORTRAN/src/` — 16 fichiers `.f90`)

| Fichier `.f90` (Ordre CMake) | Nom du Module / Programme | Équivalent C++ Direct | Responsabilité Principale |
| :--- | :--- | :--- | :--- |
| `module_precision.f90` | `module_precision` | Types natifs (`double`, `int`) | Définition de `dp = selected_real_kind(15, 307)`, `sp`, `i4`, `ip`, `i8`. |
| `module_constants.f90` | `module_constants` | `Types.hpp` (constantes/enums) | Constantes `PI`, `EPS_SMALL`, énumérations `BC_*`, `SOLVER_*`, `LINEAR_*`, `PRECOND_*`, `CONV_*`, `IC_*`, `OBS_*`, `EOS_*`, `TIME_SCHEME_*`. |
| `module_eos.f90` | `module_eos` | `EOS.hpp` | Fonctions pures `eos_parse_type`, `eos_compute_chi_T`, `eos_compute_density`, `eos_compute_temperature`, `eos_compute_sound_speed`. |
| `module_json.f90` | `module_json` | Helpers de `Config.cpp` | Extraction JSON par chemin pointé (`json_read_file`, `json_get_*`, `json_get_obstacle`). |
| `module_config.f90` | `module_config` | `Config.hpp` + `Config.cpp` | Types `bc_t`, `config_t` et routine `load_config_from_json`. |
| `module_mesh.f90` | `module_mesh` | `Mesh.hpp` + `Mesh.cpp` | Type `mesh_2d_t`, `init_mesh`, `destroy_mesh`, `x_face`, `y_face`, `x_center`, `y_center`, `is_in_obstacle`. |
| `module_fields.f90` | `module_fields` | `Field.hpp` | Type `flow_fields_t` avec tableaux 2D allouables `(0:nx+1, 0:ny+1)` et `allocate_fields`/`reset_fields`. |
| `module_mpi_comm.f90` | `module_mpi_comm` | `MpiComm.hpp` + `MpiComm.cpp` | Échanges `exchange_halos_scalar` / `sync_*`, `apply_velocity_bc`, `apply_thermal_bc`, `apply_vof_bc`, réductions globales. |
| `module_interface_reconstruction.f90` | `module_interface_reconstruction` | `InterfaceReconstruction.*` | Type `interface_segment_t`, `cut_box_area`, `compute_alpha`, `clip_segment_in_cell`, SLIC/PLIC/ELVIRA/LVIRA, limiteurs TVD. |
| `module_markers.f90` | `module_markers` | `InterfaceMarkers.*` | Types `marker_t`, `marker_cloud_t`, `interpolate_velocity_cvi`, `advect_markers_rk2`, `migrate_markers_mpi`, `prune_and_reseed_markers`, `compute_curvature_from_markers`, `write_markers_vtk`. |
| `module_interface_vof.f90` | `module_interface_vof` | `InterfaceTracker.*` | Type `interface_vof_t`, `global_marker_cloud`, `set_initial_conditions`, `advect_vof`, `smooth_svof_helmholtz`, `update_physical_properties`, `compute_curvature_and_csf`. |
| `module_navier_stokes.f90` | `module_navier_stokes` | `NavierStokesSolver.*` | Type `navier_stokes_t` (stencils `schur_a*`, `u_a*`, `v_a*` intégrés), `ns_apply_system_matrix`, `ns_assemble_rhs`, `ns_compute_adaptative_r_al`, `ns_update_pressure_al`, `ns_solve_energy`, `ns_advect_ssp_rk3_momentum`, `ns_solve_total_energy`. |
| `module_linear_solvers.f90` | `module_linear_solvers` | `LinearSolver.*` + `Preconditioner.*` | `apply_preconditioner` (`JACOBI`, `SSOR`, `ILU0`, `BLOCK_TRIANGULAR`, `CUSTOM`), `solve_bicgstab2`, `solve_bicgstab`, `solve_cg`, `solve_gmres`, `solve_minres`, `solve_direct`. |
| `module_vtk_writer.f90` | `module_vtk_writer` | `VTKWriter.*` | Type `vtk_writer_t`, écriture `.vtr`/`.pvtr` et `_history.csv` (11 colonnes incluant `avg_rho_gas`, `avg_p_gas`). |
| `module_simulation.f90` | `module_simulation` | `Simulation.*` | Type `simulation_t`, `sim_init`, `sim_run`, `sim_compute_time_step`, projection standard, modes Laplacien, Poiseuille et Vortex. |
| `covof2d_main.f90` | `program covof2d_main` | `main.cpp` | Point d'entrée CLI Fortran 2018 / MPI. |

---

## Graphe de Dépendances et Cycle d'Exécution Complet

### Graphe de Dépendances Modulaire (Identique C++ et Fortran)

```text
           +-----------------------------------------+
           |  Types / Constants / Precision / EOS    |
           +-----------------------------------------+
                                |
           +--------------------+--------------------+
           |                                         |
           v                                         v
+---------------------+                   +---------------------+
|   Mesh2D / Fields   |                   | Config / JSON Parser|
+---------------------+                   +---------------------+
           |                                         |
           +--------------------+--------------------+
                                |
                                v
                   +-------------------------+
                   |        MpiComm          |
                   +-------------------------+
                                |
           +--------------------+--------------------+
           |                                         |
           v                                         v
+---------------------+                   +---------------------+
|InterfaceReconstruct.|                   |  NavierStokesSolver |
+---------------------+                   |     (+ CSRMatrix)   |
           |                              +---------------------+
           v                                         |
+---------------------+                              v
|  InterfaceMarkers   |                   +---------------------+
+---------------------+                   | Preconditioner &    |
           |                              |    LinearSolver     |
           v                              +---------------------+
+---------------------+                              |
|  InterfaceTracker   |                              |
+---------------------+                              |
           |                                         |
           +--------------------+--------------------+
                                |
                                v
                   +-------------------------+
                   |  Simulation (+ VTKWriter|
                   +-------------------------+
                                |
                                v
                   +-------------------------+
                   | main.cpp / covof2d_main |
                   +-------------------------+
```

### Cycle d'Exécution d'un Pas de Temps (`Simulation::run` / `sim_run`)

Lors de chaque itération $t^n \to t^{n+1} = t^n + \Delta t$ en mode `NAVIER_STOKES` :

```text
+------------------------------------------------------------------------+
| 0. Calcul du Pas de Temps & Sauvegarde Temporelle                      |
| - compute_time_step() : min(config.dt, dt_cfl, dt_sigma)               |
| - Décalage états : (u_nm1, v_nm1) <- (u_old, v_old) <- (u, v, p, T...) |
+------------------------------------------------------------------------+
                                    |
                                    v
+------------------------------------------------------------------------+
| 1. Transport Inertiel Conservatif (si enable_ssp_rk3_momentum)         |
| - ns_solver.advect_ssp_rk3_momentum(fields, dt)                        |
+------------------------------------------------------------------------+
                                    |
                                    v
+------------------------------------------------------------------------+
| 2. Capture d'Interface VOF + Marqueurs CVI (si enable_two_phase)       |
|   Inside InterfaceTracker::advect() / advect_vof() :                   |
|   a) Advection C : WEYMOUTH_YUE (ou CIAM / VOF_TVD) avec normales PLIC |
|      optionnellement guidées par marker_nx/marker_ny                   |
|   b) sync_cell_scalar(C)                                               |
|   c) Si enable_svof_smoothing : smooth_svof_helmholtz() -> C_s         |
|   d) update_physical_properties() -> rho, rho_u, rho_v, mu, mu_node,   |
|      chi_T (via EOS si compressible)                                   |
|   e) Si enable_diffusion_smoothing : smooth_diffusion() -> C_tilde     |
|   f) Si use_markers :                                                  |
|      - markers.advect_rk2 (CVI Variante B + Heun + routage MPI 4 passes|
|      - markers.migrate_mpi (2 passes X puis Y)                         |
|      - markers.prune_and_reseed (sous-octants + reseeding si Np == 0)  |
|      - markers.compute_curvature_from_markers (repère (tau,eta), 2.5h) |
|   g) compute_curvature_and_CSF() -> F_sigma_x, F_sigma_y               |
|   h) apply_boundary_conditions() (et apply_vof_bc)                     |
+------------------------------------------------------------------------+
                                    |
                                    v
+------------------------------------------------------------------------+
| 3. Pénalité AAL/C-AAL & Assemblage du Préconditionneur                 |
| - Si AL/AAL/C-AAL : compute_adaptative_r_al(fields, dt)                |
| - assemble_preconditioner_blocks(fields, dt) + precond->update()       |
+------------------------------------------------------------------------+
                                    |
                                    v
+------------------------------------------------------------------------+
| 4. Couplage Vitesse-Pression                                           |
| - Si STANDARD_PROJECTION : solve_standard_projection()                 |
| - Sinon boucle d'Uzawa m = 0 .. uzawa_max-1 (uzawa_max=1 si MONOLITHIC)|
|     assemble_rhs -> pack_vector -> linear_solver->solve ->             |
|     unpack_vector -> [si AL/AAL/C-AAL : update_pressure_al & test tol] |
+------------------------------------------------------------------------+
                                    |
                                    v
+------------------------------------------------------------------------+
| 5. Propriétés Barotropes, Énergie & Diagnostics                        |
| - Si compressible && two_phase : update_physical_properties() (rho(p)) |
| - Si compressible && enable_energy : solve_total_energy()              |
|   Sinon si enable_energy : solve_energy() (Boussinesq)                 |
| - Calcul global_max_u, global_max_div, global_du_max, nu_mean          |
| - Écriture _history.csv, .vtr/.pvtr et _markers_*.vtk selon intervalles|
+------------------------------------------------------------------------+
```
