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
# Parité Duale C++17 / Fortran 2018 et Spécification d'Interopérabilité

## 2.1. Modèle Architectural Réel : Deux Implémentations Miroirs Autonomes

Avant d'aborder les conventions d'interface binaire (`iso_c_binding` / `extern "C"`), il est essentiel d'établir clairement l'architecture actuelle de `PROJET_COVOF-2D` :

- **`COVOF_2D-C++` (`bin/covof2d_cpp`)** et **`COVOF_2D-FORTRAN` (`bin/covof2d_fortran`)** sont **deux exécutables séparés, compilés indépendamment** via leurs propres `CMakeLists.txt` et liés chacun à leur bibliothèque MPI (`MPI::MPI_CXX` et `MPI::MPI_Fortran`) et à LAPACK (`-framework Accelerate` sous macOS).
- Il n'existe donc pas, dans l'état actuel du dépôt, de fichier passerelle `InteropBindings.hpp` ou `module_interop_bindings.f90` reliant les deux langages au sein d'un binaire unique : la gémellité est **algorithmique, mathématique et structurelle**.
- En revanche, les deux codes partagent exactement la même organisation mémoire, les mêmes conventions d'indiçage de grille MAC, la même arithmétique double précision IEEE 754, le même format d'entrée JSON et les mêmes formats VTK/CSV, ce qui permet de comparer leurs résultats bit-à-bit ou de coupler ultérieurement un noyau Fortran à un pilote C++ via l'ABI standard C.

---

## 2.2. Correspondance Exacte des Types et Disposition Mémoire

### 2.2.1. Correspondance des Types Fondamentaux

| Concept | C++17 (`COVOF_2D-C++`) | Taille | Fortran 2018 (`COVOF_2D-FORTRAN`) | Équivalent `iso_c_binding` |
| :--- | :--- | :---: | :--- | :--- |
| Réel double précision | `double` | 64 bits | `real(dp)` (`selected_real_kind(15, 307)`) | `real(c_double)` |
| Entier standard | `int` | 32 bits | `integer` / `integer(i4)` | `integer(c_int)` |
| Entier 64 bits | `int64_t` / `size_t` | 64 bits | `integer(i8)` | `integer(c_int64_t)` |
| Masque booléen compact | `uint8_t` (dans `FlowFields::has_marker_*`) | 8 bits | `logical` (dans `flow_fields_t`) | `logical(c_bool)` |
| Communicateur MPI | `MPI_Comm` | Opaque / Handle | `integer` (`mpif.h`) | `MPI_Comm_c2f` / `MPI_Comm_f2c` |

### 2.2.2. Parité Stricte de l'Ordonnancement Mémoire (Row-Major C++ vs Column-Major Fortran)

Le point le plus critique pour la performance vectorielle et la parité entre C++17 et Fortran 2018 réside dans l'alignement du **premier indice spatial $i$ (direction $X$) comme dimension contiguë en mémoire** :

1. **En C++17 (`Mesh.hpp`, stockage 1D aplati en *row-major*) :**
   Pour toute grille avec 1 couche de halo (`nx_ghost = nx_local + 2`, `ny_ghost = ny_local + 2`), les fonctions d'indexation définies dans `Mesh2D` sont :
   ```cpp
   inline int cell_idx(int i, int j) const { return j * nx_ghost + i; }
   inline int u_idx(int i, int j)    const { return j * (nx_ghost + 1) + i; }
   inline int v_idx(int i, int j)    const { return j * nx_ghost + i; }
   inline int node_idx(int i, int j) const { return j * (nx_ghost + 1) + i; }
   ```
   En parcourant `j` en boucle externe et `i` en boucle interne, `cell_idx(i, j)` et `cell_idx(i + 1, j)` sont **adjacents en mémoire** ($\Delta\text{adresse} = 8$ octets).

2. **En Fortran 2018 (`module_fields.f90`, tableaux 2D *column-major*) :**
   Tous les champs de `flow_fields_t` sont alloués selon `(0:nx+1, 0:ny+1)` :
   ```fortran
   allocate(fields%C(0:nx+1, 0:ny+1))
   allocate(fields%u(0:nx+1, 0:ny+1))
   ```
   En Fortran, le premier indice `i` varie de façon contiguë en mémoire. Ainsi, la boucle canonique dans tout le code Fortran place systématiquement `do j = 1, ny` à l'extérieur et `do i = 1, nx` à l'intérieur.

| Opération | Indexation C++ (`COVOF_2D-C++`) | Indexation Fortran (`COVOF_2D-FORTRAN`) | Contiguïté |
| :--- | :--- | :--- | :---: |
| Cellule $(i, j)$ | `fields.p[mesh.cell_idx(i, j)]` | `fields%p(i, j)` | Stride-1 selon $i$ |
| Face Ouest $u$ de $(i,j)$ | `fields.u[mesh.u_idx(i, j)]` | `fields%u(i, j)` | Stride-1 selon $i$ |
| Face Sud $v$ de $(i,j)$ | `fields.v[mesh.v_idx(i, j)]` | `fields%v(i, j)` | Stride-1 selon $i$ |
| Sommet Sud-Ouest $(i,j)$ | `fields.mu_node[mesh.node_idx(i, j)]` | `fields%mu_node(i, j)` | Stride-1 selon $i$ |

Grâce à cette convention symétrique, un buffer C++ `fields.p.data()` de taille `(nx_local + 2) * (ny_local + 2)` possède exactement la même représentation binaire linéaire que le tableau Fortran `fields%p(0:nx+1, 0:ny+1)` (à l'exception de `u` et `mu_node` qui ont `nx_ghost + 1` colonnes en C++).

---

## 2.3. Table Exhaustive de Correspondance des Procédures Jumelles

Le tableau suivant recense les procédures critiques implémentées en miroir dans les deux dépôts :

| Fonctionnalité / Opérateur | Méthode C++17 (`COVOF_2D-C++`) | Procédure Fortran 2018 (`COVOF_2D-FORTRAN`) |
| :--- | :--- | :--- |
| **Chargement JSON** | `Config::load_from_json` (`Config.cpp`) | `module_config:load_config_from_json` |
| **Initialisation Maillage & MPI** | `Mesh2D::init` (`Mesh.cpp`) | `module_mesh:init_mesh` |
| **Allocation des Champs MAC** | `FlowFields::allocate` (`Field.hpp`) | `module_fields:allocate_fields` |
| **Échange Halos Scalaires** | `MpiComm::sync_cell_scalar` (`MpiComm.cpp`) | `module_mpi_comm:exchange_halos_scalar` (`sync_cell_scalar`) |
| **Échange Halos Vitesses** | `MpiComm::sync_u`, `MpiComm::sync_v` | `module_mpi_comm:exchange_halos_velocity` (`sync_u`, `sync_v`) |
| **CL Hydrodynamiques & Thermiques** | `MpiComm::apply_boundary_conditions`, `apply_thermal_boundary_conditions` | `module_mpi_comm:apply_velocity_bc`, `apply_thermal_bc`, `apply_vof_bc` |
| **Réductions Globales MPI** | `MpiComm::global_dot_product`, `global_norm2`, `global_max`, `global_sum` | `module_mpi_comm:global_dot_product`, `global_norm2`, `global_max`, `global_sum` |
| **Équations d'État (EOS)** | `EOS::compute_chi_T`, `compute_density`, `compute_temperature`, `compute_sound_speed` | `module_eos:eos_compute_chi_T`, `eos_compute_density`, `eos_compute_temperature`, `eos_compute_sound_speed` |
| **Découpe Géométrique VOF** | `InterfaceReconstruction::cut_box_area`, `compute_alpha`, `clip_segment_in_cell` | `module_interface_reconstruction:cut_box_area`, `compute_alpha`, `clip_segment_in_cell` |
| **Reconstructions SLIC/PLIC/ELVIRA/LVIRA** | `InterfaceReconstruction::reconstruct_*` | `module_interface_reconstruction:reconstruct_*` |
| **Interpolation CVI Variante B** | `InterfaceMarkers::interpolate_velocity_cvi` | `module_markers:interpolate_velocity_cvi` |
| **Advection Marqueurs RK2-Heun** | `InterfaceMarkers::advect_rk2` | `module_markers:advect_markers_rk2` |
| **Migration MPI Marqueurs (2 passes)** | `InterfaceMarkers::migrate_mpi` | `module_markers:migrate_markers_mpi` |
| **Élagage & Ré-ensemencement** | `InterfaceMarkers::prune_and_reseed` | `module_markers:prune_and_reseed_markers` |
| **Courbure Parabolique $(\tau, \eta)$** | `InterfaceMarkers::compute_curvature_from_markers` | `module_markers:compute_curvature_from_markers` |
| **Export VTK Marqueurs** | `InterfaceMarkers::write_vtk` | `module_markers:write_markers_vtk` |
| **Lissage SVOF Helmholtz** | `InterfaceTracker::smooth_svof_helmholtz` | `module_interface_vof:smooth_svof_helmholtz` |
| **Advection VOF Weymouth-Yue / CIAM / TVD** | `InterfaceTracker::advect_weymouth_yue`, `advect_ciam`, `advect_tvd` | `module_interface_vof:advect_vof_weymouth_yue`, `advect_vof_ciam`, `advect_vof_tvd` |
| **Propriétés & Force Capillaire CSF** | `InterfaceTracker::update_physical_properties`, `compute_curvature_and_CSF` | `module_interface_vof:update_physical_properties`, `compute_curvature_and_csf` |
| **Opérateur $A \mathbf{x}$ Navier-Stokes** | `NavierStokesSolver::apply_system_matrix` | `module_navier_stokes:ns_apply_system_matrix` |
| **Assemblage Second Membre $\mathbf{b}$** | `NavierStokesSolver::assemble_rhs` | `module_navier_stokes:ns_assemble_rhs` |
| **Pénalité Adaptative AAL / C-AAL** | `NavierStokesSolver::compute_adaptative_r_al`, `update_pressure_al` | `module_navier_stokes:ns_compute_adaptative_r_al`, `ns_update_pressure_al` |
| **Transport Inertiel SSP-RK3 & Énergie** | `NavierStokesSolver::advect_ssp_rk3_momentum`, `solve_energy`, `solve_total_energy` | `module_navier_stokes:ns_advect_ssp_rk3_momentum`, `ns_solve_energy`, `ns_solve_total_energy` |
| **Blocs de Préconditionnement** | `NavierStokesSolver::assemble_preconditioner_blocks` + `Preconditioner::apply` | `module_navier_stokes:ns_assemble_preconditioner_blocks` + `module_linear_solvers:apply_preconditioner` |
| **Solveurs de Krylov & LU LAPACK** | `BiCGSTAB2Solver`, `BiCGSTABSolver`, `CGSolver`, `GMRESSolver`, `MINRESSolver`, `DirectSolver` | `module_linear_solvers:solve_bicgstab2`, `solve_bicgstab`, `solve_cg`, `solve_gmres`, `solve_minres`, `solve_direct` |
| **Sorties `.vtr`/`.pvtr` et `_history.csv`** | `VTKWriter::init`, `write_step`, `write_diagnostics` | `module_vtk_writer:vtk_init`, `vtk_write_step`, `vtk_write_diagnostics` |

---

## 2.4. Spécification d'Interopérabilité Binaire (`iso_c_binding` et `extern "C"`)

Bien que les deux codes soient compilés séparément aujourd'hui, le solveur s'appuie déjà sur la convention d'appel C pour les appels externes à LAPACK (`dgetrf_` / `dgetrs_`) et fournit une base idéale pour un futur couplage binaire direct (par exemple pour appeler un noyau vectorisé Fortran depuis un pilote C++ ou intégrer des bibliothèques externes comme MUMPS/PETSc/HYPRE).

### 2.4.1. Appel Direct des Routines LAPACK depuis C++ (`LinearSolver.cpp`)

Dans `COVOF_2D-C++/src/LinearSolver.cpp`, le solveur direct `DirectSolver` interface directement les routines Fortran 77/90 de LAPACK via une déclaration `extern "C"` respectant l'ajout du tiret bas (*symbol underscoring*) standard :

```cpp
extern "C" {
    void dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
    void dgetrs_(const char* trans, int* n, int* nrhs, const double* a, int* lda,
                 const int* ipiv, double* b, int* ldb, int* info);
}
```

Côté Fortran (`module_linear_solvers.f90`), `solve_direct` appelle directement `dgetrf` et `dgetrs` via l'interface implicite LAPACK liée par CMake (`-framework Accelerate` ou `-llapack`). Dans les deux cas, la matrice dense `A_local` de dimension $N \times N$ est remplie en ordre **column-major** (`col * n + row` en C++, `A_local(row, col)` en Fortran) pour être immédiatement compatible avec `dgetrf_`.

### 2.4.2. Patron Standard `iso_c_binding` pour Passerelle C++ $\leftrightarrow$ Fortran

Pour exposer sans copie mémoire un noyau Fortran 2018 à C++17 (par exemple un opérateur sur champs scalaires centrés `(0:nx+1, 0:ny+1)`), la norme Fortran 2018 impose l'attribut `bind(C, name="...")` et le module intrinsèque `iso_c_binding` :

```fortran
! Patron de liaison binaire sans copie (Zero-Copy ABI)
subroutine covof_f_kernel_cvi(u_ptr, v_ptr, nx, ny, dx, dy, xp, yp, up, vp) &
           bind(C, name="covof_f_kernel_cvi")
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr),    intent(in), value :: u_ptr, v_ptr
  integer(c_int), intent(in), value :: nx, ny
  real(c_double), intent(in), value :: dx, dy, xp, yp
  real(c_double), intent(out)       :: up, vp

  real(c_double), pointer :: u(:, :), v(:, :)
  ! Vue Fortran 2D sur le buffer contigu C++ sans réallocation
  call c_f_pointer(u_ptr, u, [nx + 2, ny + 2])
  call c_f_pointer(v_ptr, v, [nx + 2, ny + 2])
  ! ... calcul direct sur u(1:nx+2, 1:ny+2) ...
end subroutine covof_f_kernel_cvi
```

### 2.4.3. Conversion de Communicateurs MPI (`MPI_Comm_c2f` / `MPI_Comm_f2c`) et Gestion des Erreurs

Dans tout couplage hybride MPI entre C++ et Fortran :
1. **Conversion de communicateurs :** C++ manipule un `MPI_Comm` opaque (pointeur sous OpenMPI, entier sous MPICH) tandis que Fortran manipule un `integer` 32 bits. La conversion portable passe obligatoirement par :
   ```cpp
   MPI_Fint f_comm = MPI_Comm_c2f(mesh.comm); // C++ -> Fortran
   ```
   et côté Fortran `comm_c = MPI_Comm_f2c(f_comm)` si l'on appelle une routine C depuis Fortran.
2. **Propagation des codes d'erreur :** Aucune exception C++ (`throw`) ne doit traverser une frontière `extern "C"`, et aucun `stop` Fortran ne doit contourner `MPI_Finalize()`. Tout appel inter-langage doit retourner un entier `integer(c_int), intent(out) :: ierr` vérifié côté appelant avant de déclencher une exception C++ ou un `MPI_Abort` coordonné.
3. **Synchronisation des flux d'E/S :** Toujours vider les tampons (`std::cout << std::flush;` en C++ et `flush(6)` en Fortran) avant de franchir une frontière de langage afin d'éviter l'entrelacement des journaux MPI sur stdout.
# Structures de Données, Mémoire et Parallélisme HPC

## 3.1. Représentation Interne du Maillage MAC et des Champs

### 3.1.1. Dimensions Locales et Couche Unique de Cellules Fantômes

Dans les deux implémentations (**C++17** et **Fortran 2018**), chaque processus MPI gère un sous-domaine rectangulaire de dimensions intérieures `nx_local` $\times$ `ny_local`, entouré exactement d'**une unique couche périphérique de cellules fantômes** (*halo* d'épaisseur $1$ de chaque côté).

- **Côté C++ (`Mesh.hpp`) :**
  `nx_ghost = nx_local + 2`, `ny_ghost = ny_local + 2`. Les mailles intérieures occupent $i \in [1, \text{nx\_local}]$, $j \in [1, \text{ny\_local}]$ ; les cellules fantômes occupent les bords $i = 0$, $i = \text{nx\_local}+1$, $j = 0$ et $j = \text{ny\_local}+1$.
  Les tailles exactes allouées dans `FlowFields::allocate(const Mesh2D& mesh)` sont :
  - Champs centrés ($p, C, C_s, \widetilde{C}, T, \rho, \mu, \kappa, \chi_T, r_{\text{aal}}, \dots$) : `n_cells = nx_ghost * ny_ghost`
  - Faces verticales $u$ ($\rho_u, F_{\sigma, x}, (\rho u)^*$) : `n_u = (nx_ghost + 1) * ny_ghost` (soit `(nx_local + 3) * (ny_local + 2)`)
  - Faces horizontales $v$ ($\rho_v, F_{\sigma, y}, (\rho v)^*$) : `n_v = nx_ghost * (ny_ghost + 1)`
  - Sommets / nœuds ($\mu_{\text{node}}$) : `n_nodes = (nx_ghost + 1) * (ny_ghost + 1)`

- **Côté Fortran (`module_fields.f90`) :**
  Tous les tableaux 2D de `flow_fields_t` sont alloués uniformément avec les bornes explicites `(0:nx+1, 0:ny+1)` (`allocate_fields(fields, mesh%nx_local, mesh%ny_local)`), garantissant la compatibilité d'indexation directe sur les faces $(i, j)$ et les centres $(i, j)$.

:::note
L'épaisseur de halo de largeur 1 (`0..nx+1, 0..ny+1`) suffit à l'ensemble des opérateurs du code : pochoir visqueux à 9 points (incluant les nœuds $\mu_{\text{node}}$), reconstruction PLIC/ELVIRA $3\times 3$, filtre de Helmholtz SVOF, advection Weymouth-Yue et **interpolation CVI Variante B** (qui n'accède qu'à $j_u \in \{j-1, j\} \subset [0, N_y]$ et $i_v \in \{i-1, i\} \subset [0, N_x]$).
:::

### 3.1.2. Architecture SoA (*Structure of Arrays*) des Champs d'Écoulement

COVOF 2D adopte strictement le modèle **SoA** (*Structure of Arrays*) afin de maximiser le taux de succès des lignes de cache L1/L2 et la vectorisation SIMD lors des parcours de grille : chaque variable physique est stockée dans son propre tableau contigu.

Dans `FlowFields` (`Field.hpp`) et `flow_fields_t` (`module_fields.f90`), les champs sont regroupés en cinq familles :
1. **Champs primaires :** `u`, `v`, `p`, `C`, `T`.
2. **Propriétés matérielles effectives :** `rho`, `rho_u`, `rho_v`, `mu`, `mu_node`.
3. **Géométrie d'interface & CSF :** `kappa` (`curvature` en Fortran), `norm_x`/`norm_y` (C++), `F_sigma_x`/`F_sigma_y` (`fx_sigma`/`fy_sigma` en Fortran), `C_s` (SVOF), `C_tilde` (diffusion capillaire).
4. **Historique temporel, All-Mach & AAL :** `div_u` (`div`), `r_aal`, `chi_T`, `p_old`, `T_old`, `u_old`, `v_old`, `u_nm1`, `v_nm1`, `rho_old`, `rho_u_old`, `rho_v_old`, `rho_u_nm1`, `rho_v_nm1`, `has_nm1`, `energy`, `rho_e`, `rhou_star`, `rhov_star`, `rhoe_star`.
5. **Champs géométriques issus des marqueurs CVI :** `marker_nx`, `marker_ny`, `marker_kappa`, `has_marker_normal`, `has_marker_kappa` (et `solid_mask` en Fortran).

### 3.1.3. Politique *Zero-Allocation* dans la Boucle Temporelle

Pour éviter les allocations dynamiques coûteuses dans la boucle de Krylov :
- C++ préalloue `state_vector`, `rhs`, `u_prev` dans `Simulation::run()`, ainsi que les tampons temporaires `mutable FlowFields tf_scratch` et `mutable std::vector<double> div_cell_scratch` dans `NavierStokesSolver`, réutilisés sans réallocation à chaque appel de `apply_system_matrix()`.
- Fortran préalloue `state_vector`, `rhs`, `u_prev` dans `sim_run`, ainsi que les stencils `diag_u`, `diag_v`, `diag_p`, `schur_a*`, `u_a*`, `v_a*` dans `ns_init`.

---

## 3.2. Structures de Données du Nuage de Marqueurs CVI Non Connectés

Contrairement aux méthodes *Front-Tracking* classiques qui maintiennent des chaînes d'arêtes ou des triangles surfaciques (nécessitant une chirurgie topologique complexe lors des coalescences ou ruptures), COVOF 2D stocke les marqueurs dans un **nuage plat non connecté** :

- **Côté C++ (`InterfaceMarkers.hpp`) :**
  ```cpp
  struct Marker {
      double x = 0.0, y = 0.0;       // Coordonnées globales [m]
      double nx = 0.0, ny = 0.0;     // Normale unitaire (+grad C / |grad C|)
      double kappa = 0.0;            // Courbure locale [1/m]
      int cell_i = 1, cell_j = 1;    // Cellule hôte locale (1..nx_local, 1..ny_local)
      bool active = true;
  };
  // Conteneur dynamique dans InterfaceMarkers : std::vector<Marker> cloud;
  ```
- **Côté Fortran (`module_markers.f90`) :**
  ```fortran
  type :: marker_t
     real(dp) :: x = 0.0_dp, y = 0.0_dp
     real(dp) :: nx = 0.0_dp, ny = 0.0_dp
     real(dp) :: kappa = 0.0_dp
     integer  :: i_cell = 1, j_cell = 1
     logical  :: active = .true.
  end type marker_t

  type :: marker_cloud_t
     type(marker_t), allocatable :: list(:)
     integer :: count = 0
     integer :: capacity = 0
     integer, allocatable :: count_in_cell(:, :)
  end type marker_cloud_t
  ```
  Le redimensionnement dynamique de `cloud%list` en Fortran utilise la sous-routine `ensure_capacity` avec doublement géométrique (`max(min_cap, 2 * cloud%capacity)`) et `move_alloc`, garantissant un coût amorti $\mathcal{O}(1)$ par insertion.

---

## 3.3. Parallélisme Distribué MPI Pur (Cartésien 2D)

COVOF 2D repose **exclusivement sur MPI pur** (OpenMPI / MPICH / Intel MPI) sans parallélisme partagé OpenMP.

### 3.3.1. Décomposition Cartésienne et Échanges de Halos en 2 Passes

Lors de l'initialisation (`Mesh2D::init` / `init_mesh`), `MPI_Dims_create` et `MPI_Cart_create` construisent le communicateur cartésien 2D avec périodicité éventuelle `periods = {periodic_x, periodic_y}`.

Les fonctions `MpiComm::sync_cell_scalar`, `sync_u`, `sync_v` (`exchange_halos_scalar`, `exchange_halos_velocity` en Fortran) procèdent toujours en **deux passes séquentielles** avec requêtes non bloquantes `MPI_Isend`/`MPI_Irecv` et barrière de complétion `MPI_Waitall` :
- **Pour les scalaires centrés et $u$ :** Passe 1 sur l'axe $X$ (`1..ny_local`), puis Passe 2 sur l'axe $Y$ sur toute la largeur incluant les colonnes fantômes (`0..nx_local+1`), ce qui propage automatiquement les **4 cellules de coins diagonaux**.
- **Pour $v$ :** Passe 1 sur l'axe $Y$ (`1..ny_local`), puis Passe 2 sur l'axe $X$ (`0..ny_local+1`) pour assurer la même propagation exacte des 4 coins.
- **Cas séquentiel (`num_procs == 1`) avec périodicité :** Copie directe des bords opposés sans appel MPI.

### 3.3.2. Communications MPI Spécifiques aux Marqueurs CVI

Le suivi lagrangien en environnement distribué utilise trois schémas de communication sans voisinage diagonal explicite :

1. **Routage aller-retour à 4 passes dans `advect_rk2` (`advect_markers_rk2`) :**
   Lorsqu'un marqueur situé près d'une frontière inter-processus a sa position prédite $\mathbf{x}_p^* = \mathbf{x}_p^n + \Delta t \, \mathbf{u}_{\text{CVI}}(\mathbf{x}_p^n; \mathbf{u}^n)$ hors du sous-domaine local, sa coordonnée et son identifiant de requête `(q_id, dir_x, x*, y*)` sont envoyés via `MPI_Sendrecv` en **Passe X aller** puis **Passe Y aller** vers le rang propriétaire. Ce dernier évalue $\mathbf{u}_{\text{CVI}}(\mathbf{x}_p^*; \mathbf{u}^{n+1})$ sur sa propre grille MAC et renvoie `(u2, v2)` par **Passe Y retour** puis **Passe X retour**.
2. **Migration dimensionnelle en 2 passes (`migrate_mpi` / `migrate_markers_mpi`) :**
   Après mise à jour finale $\mathbf{x}_p^{n+1}$, la **Passe X** transfère (via `MPI_Sendrecv`) les marqueurs sortis à gauche ou à droite (avec décalage périodique $\pm L_x$ si `periodic_x`), puis la **Passe Y** balaie le nuage mis à jour et transfère ceux sortis en bas ou en haut (avec décalage $\pm L_y$ si `periodic_y`). Un marqueur traversant un coin de sous-domaine est ainsi automatiquement relayé en 2 étapes sans jamais nécessiter de communicateur à 8 voisins.
3. **Échange de marqueurs fantômes dans une bande $2.5 h$ (`compute_curvature_from_markers`) :**
   Pour que l'ajustement parabolique $3\times 3$ autour d'une cellule de bordure de sous-domaine donne exactement le même résultat qu'en séquentiel (`np=1` vs `np=4`), les rangs échangent en 2 passes (X puis Y) une copie temporaire de leurs marqueurs situés dans une bande frontalière de largeur $2.5 \max(\Delta x, \Delta y)$.

### 3.3.3. Réductions Globales et Gestion des DDLs Possédés

Dans le système linéaire distribué, chaque rang MPI possède sans duplication :
- `n_owned_u = u_nx * ny_local` avec `u_nx = nx_local + (is_boundary_right ? 1 : 0)` (la face frontière partagée entre deux rangs appartient au rang de droite).
- `n_owned_v = nx_local * v_ny` avec `v_ny = ny_local + (is_boundary_top ? 1 : 0)`.
- `n_owned_p = nx_local * ny_local` (uniquement lorsque `is_monolithic()` est vrai ; sinon `n_owned_p = 0`).

Les produits scalaires et normes de Krylov (`global_dot_product`, `global_norm2`) ainsi que les extrema (`global_max`, `global_sum`) utilisent `MPI_Allreduce` sur `mesh.comm`.

---

## 3.4. Stockage Creux CSR en C++ (`CSRMatrix`) et Blocs Locaux Fortran

Contrairement à l'opérateur complet $A \mathbf{x}$ qui est **100 % *matrix-free*** dans les deux langages, les blocs diagonaux de préconditionnement locaux ($F_{uu}$, $F_{vv}$, $S_p$) utilisent :

- **En C++ (`include/covof/CSRMatrix.hpp` & `src/CSRMatrix.cpp`) :**
  La classe `CSRMatrix` stocke une matrice creuse locale au format *Compressed Sparse Row* via quatre vecteurs contigus :
  - `row_ptr` (`std::vector<int>`, taille `nrows + 1`)
  - `col_idx` (`std::vector<int>`)
  - `values` (`std::vector<double>`)
  - `diag` (`std::vector<double>`, diagonale extraite explicitement pour Jacobi / Gauss-Seidel)
  Elle fournit `set_from_triplets(r, c, rows, cols, vals)` (tri et fusion des doublons), `matvec(x, y)` (produit $y = Ax$ séquentiel par sous-domaine) et `solve_jacobi(b, x, num_iters)`. Dans `NavierStokesSolver::assemble_preconditioner_blocks`, elle assemble `Fuu`, `Fvv` et `SchurLaplace`.
- **En Fortran (`module_navier_stokes.f90`) :**
  Il n'existe pas de module `module_csr_matrix.f90` séparé : `navier_stokes_t` stocke directement les stencils structurés à 5 points sous forme de tableaux 1D/2D contigus (`u_aE, u_aW, u_aN, u_aS, u_aP`, `v_aE, ..., v_aP`, `schur_aE, ..., schur_aP`), optimisés pour les balayages SSOR vectorisés dans `apply_preconditioner`.
# Correspondance Théorie $\leftrightarrow$ Discrétisation $\leftrightarrow$ Code Source

Ce chapitre constitue la **section pivot** du Guide de Développement. Il établit le lien explicite entre les équations continues de la mécanique des fluides, leur discrétisation spatiotemporelle sur grille MAC et leur implémentation exacte dans `COVOF_2D-C++` (C++17) et `COVOF_2D-FORTRAN` (Fortran 2018).

Chaque opérateur est présenté selon le triptyque rigoureux :
1. **Équation discrète exacte** (indicée selon les conventions du code).
2. **Localisation exacte** (fichiers, classes/modules, méthodes/subroutines et plages de lignes).
3. **Analyse de l'implémentation** (boucles critiques, régularisations et invariants).

---

## 4.1. Diffusion Visqueuse et Moyenne Harmonique aux Nœuds ($\mu_{\text{node}}$)

### 4.1.1. Équation Discrète Exacte

Sur la grille MAC, la composante $u_{i,j}$ (face Ouest de la cellule $(i,j)$) intègre deux types de flux visqueux sur le volume de contrôle $\Delta x \Delta y$ :
- Les contraintes normales d'élongation $2\mu \partial_x u$ évaluées aux centres $(i,j)$ et $(i-1,j)$ :
  $$\text{Diff}_{\text{elong}} = \left[ 2 \mu_{i,j} \frac{u_{i+1,j} - u_{i,j}}{\Delta x} - 2 \mu_{i-1,j} \frac{u_{i,j} - u_{i-1,j}}{\Delta x} \right] \Delta y$$
- Les contraintes tangentielles de cisaillement $\mu_{\text{node}}(\partial_y u + \partial_x v)$ évaluées aux nœuds Nord $(i, j+1)$ et Sud $(i, j)$ :
  $$\text{Shear} = \left[ \mu_{\text{node}, i, j+1} \left( \frac{u_{i,j+1} - u_{i,j}}{\Delta y} + \frac{v_{i,j+1} - v_{i-1,j+1}}{\Delta x} \right) - \mu_{\text{node}, i, j} \left( \frac{u_{i,j} - u_{i,j-1}}{\Delta y} + \frac{v_{i,j} - v_{i-1,j}}{\Delta x} \right) \right] \Delta x$$

Pour assurer la continuité des contraintes tangentielles à travers une interface diphasique à fort saut $\mu_1/\mu_2$, $\mu_{\text{node}}$ est calculée par **moyenne harmonique** (ou interpolation SVOF harmonique de $C_{s,\text{node}}$ si `enable_svof_smoothing = true`) :

$$\mu_{\text{node}, i, j} = \frac{4}{\frac{1}{\mu_{i-1, j-1}} + \frac{1}{\mu_{i, j-1}} + \frac{1}{\mu_{i-1, j}} + \frac{1}{\mu_{i, j}}} \quad \text{ou} \quad \frac{\mu_1 \mu_2}{\mu_2 C_{s,\text{node}} + \mu_1 (1 - C_{s,\text{node}})}$$

### 4.1.2. Localisation dans le Code Source

- **Calcul de $\mu_{\text{node}}$ :**
  - C++17 : `COVOF_2D-C++/src/InterfaceTracker.cpp` $\rightarrow$ `InterfaceTracker::update_physical_properties` (lignes 558–582).
  - Fortran 2018 : `COVOF_2D-FORTRAN/src/module_interface_vof.f90` $\rightarrow$ `update_physical_properties` (lignes 479–626).
- **Application dans l'opérateur linéaire :**
  - C++17 : `COVOF_2D-C++/src/NavierStokesSolver.cpp` $\rightarrow$ `NavierStokesSolver::apply_system_matrix` (lignes 334–358 pour $u$, lignes 480–510 pour $v$).
  - Fortran 2018 : `COVOF_2D-FORTRAN/src/module_navier_stokes.f90` $\rightarrow$ `ns_apply_system_matrix` (lignes 380–809).

### 4.1.3. Analyse de l'Implémentation C++ (`NavierStokesSolver.cpp`, lignes 336–358)

```cpp
// (a) Viscosité d'élongation 2*mu d2u/dx2
double mu_cell_E = fields_old.mu[mesh.cell_idx(i, j)];
double mu_cell_W = fields_old.mu[mesh.cell_idx(i - 1, j)];
double diff_elong_E = 2.0 * mu_cell_E * (tf.u[mesh.u_idx(i + 1, j)] - tf.u[mesh.u_idx(i, j)]) / dx * dy;
double diff_elong_W = 2.0 * mu_cell_W * (tf.u[mesh.u_idx(i, j)] - tf.u[mesh.u_idx(i - 1, j)]) / dx * dy;
val -= (diff_elong_E - diff_elong_W);

// (b) Viscosité de cisaillement et rotation mu (d2u/dy2 + d2v/dxdy)
double mu_node_N = fields_old.mu_node[mesh.node_idx(i, j + 1)];
double mu_node_S = fields_old.mu_node[mesh.node_idx(i, j)];
double shear_N = mu_node_N * ((tf.u[mesh.u_idx(i, j + 1)] - tf.u[mesh.u_idx(i, j)]) / dy +
                              (tf.v[mesh.v_idx(i, j + 1)] - tf.v[mesh.v_idx(i - 1, j + 1)]) / dx) * dx;
double shear_S = mu_node_S * ((tf.u[mesh.u_idx(i, j)] - tf.u[mesh.u_idx(i, j - 1)]) / dy +
                              (tf.v[mesh.v_idx(i, j)] - tf.v[mesh.v_idx(i - 1, j)]) / dx) * dx;
val -= (shear_N - shear_S);
```

Notice : L'opérateur complet étant multiplié par le volume cellulaire $\Delta x \Delta y$, les flux surfaciques sont directement multipliés par $\Delta y$ (faces $X$) et $\Delta x$ (faces $Y$), garantissant la conservation symétrique.

---

## 4.2. Advection Convective et Linéarisations (Picard, AB2, Newton)

### 4.2.1. Équation Discrète Exacte

Dans l'opérateur implicite `apply_system_matrix`, la forme semi-implicite de la convection sur la composante $u_{i,j}$ avec vitesse advectante de référence $(\mathbf{u}_{\text{ref}}, \mathbf{v}_{\text{ref}})$ s'écrit avec flux aux faces $F_E = \frac{1}{2}(u_{\text{ref}, i,j} + u_{\text{ref}, i+1,j})\Delta y$, $F_W = \frac{1}{2}(u_{\text{ref}, i-1,j} + u_{\text{ref}, i,j})\Delta y$, $F_N = \bar{v}_{\text{ref}, N}\Delta x$, $F_S = \bar{v}_{\text{ref}, S}\Delta x$ :

$$\mathcal{C}_u(\delta \mathbf{u}) = \frac{\rho_{u, i, j}}{2} \left[ F_E (\delta u_{i+1,j} + \delta u_{i,j}) - F_W (\delta u_{i-1,j} + \delta u_{i,j}) + F_N (\delta u_{i,j+1} + \delta u_{i,j}) - F_S (\delta u_{i,j-1} + \delta u_{i,j}) \right]$$

- **Mode `"PICARD"` :** $\mathbf{u}_{\text{ref}} = \mathbf{u}^n$.
- **Mode `"ADAMS_BASHFORTH"` / `"AB2"` :** lorsque `has_nm1 == true`, $\mathbf{u}_{\text{ref}} = 2 \mathbf{u}^n - \mathbf{u}^{n-1}$ (extrapolation à $t^{n+1}$ depuis $n$ et $n-1$).
- **Mode `"NEWTON"` :** ajoute le jacobien croisé $\rho_u \left( \delta u_{i,j} \frac{u^n_{i+1,j} - u^n_{i-1,j}}{2\Delta x} + \bar{\delta v}_{i,j} \frac{u^n_{i,j+1} - u^n_{i,j-1}}{2\Delta y} \right) \Delta x \Delta y$, et compense le terme quadratique dans `assemble_rhs`.

### 4.2.2. Localisation dans le Code Source

- **C++17 :** `COVOF_2D-C++/src/NavierStokesSolver.cpp` $\rightarrow$ `apply_system_matrix` (lignes 230–332) et `assemble_rhs` (lignes 553–839).
- **Fortran 2018 :** `COVOF_2D-FORTRAN/src/module_navier_stokes.f90` $\rightarrow$ `ns_apply_system_matrix` (lignes 380–809) et `ns_assemble_rhs` (lignes 810–1124).

---

## 4.3. Lagrangien Augmenté Adaptatif (AL, AAL, C-AAL) et Mise à Jour d'Uzawa

### 4.3.1. Équation Discrète Exacte

1. **Pénalité spatio-temporelle locale AAL (Vincent et al. 2004, Éq. 14) :**
   Dans chaque cellule $(i, j)$ :
   $$r_{i,j} = K_{\text{aal}} \max\left( \alpha_{\text{time}} \rho_{i,j} \frac{L_0^2}{t_0}, \; \rho_{i,j} u_0 L_0, \; \rho_{i,j} \frac{L_0^2}{u_0}\|\mathbf{g}\|, \; \frac{p_0 L_0}{u_0}, \; \mu_{i,j}, \; \mu_{i,j} \frac{\Delta x \Delta y}{\mathcal{K}_{i,j}}, \; \frac{\sigma}{u_0} \right)$$
   En mode compressible (`C-AAL` ou `compressible && !is_aal_incompressible`), régularisation acoustique :
   $$r_{\text{eff}, i,j} = \frac{r_{i,j}}{1 + r_{i,j} \frac{\chi_{T, i,j}}{\Delta t}}$$

2. **Terme de divergence pénalisée dans $A \mathbf{x}$ (sur face $u_{i,j}$) :**
   $$-\left( r_{\text{eff}, i,j} D_{i,j} - r_{\text{eff}, i-1,j} D_{i-1,j} \right) \Delta y$$
   où $D_{i,j} = \frac{u_{i+1,j} - u_{i,j}}{\Delta x} + \frac{v_{i,j+1} - v_{i,j}}{\Delta y}$ (en mode `C-AAL` ou `MONOLITHIC` compressible, la ligne de pression inclut également $+\frac{\chi_T}{\Delta t} p^{n+1} \Delta x \Delta y$).

3. **Actualisation duale d'Uzawa :**
   $$p_{i,j}^{m+1} = p_{i,j}^m - r_{\text{eff}, i,j} \left[ D_{i,j}^{*, m+1} + \frac{\chi_{T,i,j}}{\Delta t}(p_{i,j}^m - p_{i,j}^n) \right]$$

### 4.3.2. Localisation dans le Code Source

- **C++17 :** `COVOF_2D-C++/src/NavierStokesSolver.cpp` $\rightarrow$ `compute_adaptative_r_al` (lignes 841–927), `update_pressure_al` (lignes 929–951), `apply_system_matrix` (lignes 361–365).
- **Fortran 2018 :** `COVOF_2D-FORTRAN/src/module_navier_stokes.f90` $\rightarrow$ `ns_compute_adaptative_r_al` (lignes 1125–1241), `ns_update_pressure_al` (lignes 1242–1286).

---

## 4.4. Pénalisation d'Obstacles Immergés (Darcy-Brinkman)

### 4.4.1. Équation Discrète Exacte

Pour toute face $u_{i,j}$ (ou $v_{i,j}$) détectée dans un obstacle (`BOX` ou `CYLINDER`) de perméabilité $\mathcal{K}$ :
- Si $\mathcal{K} \le 10^{-12}$ (solide imperméable) : la ligne de la matrice est remplacée par la condition de Dirichlet homogène $\left(\frac{\alpha_{\text{time}} \rho_u}{\Delta t} \Delta x \Delta y\right) \delta u_{i,j} = 0$.
- Si $\mathcal{K} > 10^{-12}$ (milieu poreux) : ajout implicite du terme de Darcy $+\frac{\mu_u}{\mathcal{K}} \Delta x \Delta y \, \delta u_{i,j}$ à l'opérateur et à la diagonale `diag_u` / `u_aP`.

### 4.4.2. Localisation dans le Code Source

- **C++17 :** `COVOF_2D-C++/include/covof/NavierStokesSolver.hpp` (`is_point_in_obstacle`, `is_u_face_in_obstacle`, `get_u_obstacle_perm`, lignes 101–202) et `src/NavierStokesSolver.cpp` (lignes 278–287 et 376–382).
- **Fortran 2018 :** `COVOF_2D-FORTRAN/src/module_mesh.f90` (`is_in_obstacle`, lignes 203–222) et `module_navier_stokes.f90` (`is_u_face_in_obstacle`, `is_v_face_in_obstacle`, lignes 82–131).

---

## 4.5. Suivi d'Interface VOF-PLIC et Advection de Weymouth-Yue

### 4.5.1. Équation Discrète Exacte

1. **Reconstruction affine par morceaux (`InterfaceSegment` : $n_x x + n_y y + \alpha = 0$ orienté vers $C = 0$) :**
   - `reconstruct_plic` : normale de Youngs $3\times 3$ pondérée $[1, 2, 1]$.
   - `reconstruct_elvira` : sélection du candidat minimisant l'erreur volumique $L_2$ sur le bloc $3\times 3$ parmi 6 fonctions de hauteur horizontales et verticales.
   - `compute_alpha` : bisection exacte sur l'intervalle $[\min_{4\text{ coins}} (-\mathbf{n}\cdot\mathbf{x}), \max_{4\text{ coins}} (-\mathbf{n}\cdot\mathbf{x})]$ utilisant l'aire polygonale exacte `cut_box_area` (clipping de Sutherland-Hodgman + formule de Gauss/Shoelace).

2. **Advection conservative de Weymouth & Yue (2010) avec splitting alterné de Strang :**
   - En début de pas : $C_{d,i,j} = 1$ si $C_{i,j}^n \ge 0.5$, sinon $0$.
   - Passe X : $C_{i,j}^* = C_{i,j}^n - \frac{f_{x, i+1, j} - f_{x, i, j}}{\Delta x \Delta y} + C_{d,i,j} \frac{\Delta t}{\Delta x}(u_{i+1,j} - u_{i,j})$ (sans écrêtage intermédiaire : `clamp_out = false`).
   - Passe Y : $C_{i,j}^{n+1} = C_{i,j}^* - \frac{f_{y, i, j+1} - f_{y, i, j}}{\Delta x \Delta y} + C_{d,i,j} \frac{\Delta t}{\Delta y}(v_{i,j+1} - v_{i,j})$ (avec écrêtage final `[0, 1]`).

### 4.5.2. Localisation dans le Code Source

- **C++17 :**
  - `COVOF_2D-C++/src/InterfaceReconstruction.cpp` : `cut_box_area` (l. 37), `compute_alpha` (l. 106), `clip_segment_in_cell` (l. 147), `reconstruct_plic` (l. 252), `reconstruct_elvira` (l. 321), `reconstruct_lvira` (l. 436).
  - `COVOF_2D-C++/src/InterfaceTracker.cpp` : `advect` (l. 849–894), `advect_weymouth_yue` (l. 1148–1327).
- **Fortran 2018 :**
  - `COVOF_2D-FORTRAN/src/module_interface_reconstruction.f90` : `cut_box_area` (l. 33), `compute_alpha` (l. 105), `clip_segment_in_cell` (l. 152), `reconstruct_elvira` (l. 346).
  - `COVOF_2D-FORTRAN/src/module_interface_vof.f90` : `advect_vof` (l. 792), `advect_vof_weymouth_yue` (l. 1162).

---

## 4.6. Méthode Hybride FC-WY-CVI-Marker (Variante B Exacte sur Boîtes MAC)

### 4.6.1. Formulation Mathématique Discrète

Pour tout point matériel $(x_p, y_p)$ localisé dans la cellule $(i, j)$, la **Variante B exacte** (Wang et al., 2015 ; Pusok et al., 2017) définit deux boîtes décalées distinctes restant strictement dans le halo `(0:nx+1, 0:ny+1)` :
- Boîte horizontale $u$ : $j_u = (y_p < y_j^c) \;?\; j-1 : j$, $\xi_u = \frac{x_p - x_i^f}{\Delta x}$, $\eta_u = \frac{y_p - y_{j_u}^c}{\Delta y}$.
- Boîte verticale $v$ : $i_v = (x_p < x_i^c) \;?\; i-1 : i$, $\xi_v = \frac{x_p - x_{i_v}^c}{\Delta x}$, $\eta_v = \frac{y_p - y_j^f}{\Delta y}$.

Avec $\Delta u_S = u_{i+1, j_u} - u_{i, j_u}$, $\Delta u_N = u_{i+1, j_u+1} - u_{i, j_u+1}$, $\Delta v_W = v_{i_v, j+1} - v_{i_v, j}$ et $\Delta v_E = v_{i_v+1, j+1} - v_{i_v+1, j}$ :

$$u_{\text{CVI}} = u_{\text{bilin}}(\xi_u, \eta_u) + \frac{\Delta x}{2 \Delta y} \xi_u (1 - \xi_u)(\Delta v_E - \Delta v_W)$$
$$v_{\text{CVI}} = v_{\text{bilin}}(\xi_v, \eta_v) + \frac{\Delta y}{2 \Delta x} \eta_v (1 - \eta_v)(\Delta u_N - \Delta u_S)$$

L'advection temporelle utilise le schéma RK2 de Heun (prédicteur $\mathbf{x}_p^*$ avec $\mathbf{u}^n$ puis correcteur avec $\mathbf{u}^{n+1}$ ou le champ courant en mode cinématique), suivie de la migration MPI 2 passes, de l'élagage par sous-octants et de l'ajustement parabolique pondéré dans le repère tourné $(\tau, \eta)$ :

$$\min_{a,b,c} \sum_{p=1}^{M} w_p \left(a \tau_p^2 + b \tau_p + c - \eta_p\right)^2 \quad \text{sous } \mathbf{n}_p \cdot \mathbf{n}^{\text{ref}} \ge 0.5 \implies \kappa_{i,j}^{\text{marker}} = \frac{2a}{(1 + b^2)^{3/2}}$$

### 4.6.2. Localisation dans le Code Source

- **C++17 :** `COVOF_2D-C++/src/InterfaceMarkers.cpp` :
  - `init_from_vof` (l. 16–109)
  - `interpolate_velocity_cvi` (l. 123–185)
  - `advect_rk2` (l. 239–495, incluant routage MPI aller-retour à 4 passes)
  - `migrate_mpi` (l. 497–679, migration X puis Y)
  - `prune_and_reseed` (l. 681–935)
  - `compute_curvature_from_markers` (l. 937–1384, incluant halo $2.5 h$ et lissage tangentiel)
  - `write_vtk` (l. 1390–1459)
- **Fortran 2018 :** `COVOF_2D-FORTRAN/src/module_markers.f90` :
  - `init_markers_from_vof` (l. 110), `interpolate_velocity_cvi` (l. 232), `advect_markers_rk2` (l. 327), `migrate_markers_mpi` (l. 763), `prune_and_reseed_markers` (l. 995), `compute_curvature_from_markers` (l. 1277), `write_markers_vtk` (l. 1785).
- **Couplage dans VOF et CSF :**
  - C++17 : `InterfaceTracker::advect` (`InterfaceTracker.cpp`, l. 887–893), `compute_curvature_and_CSF` (l. 658–733), `apply_marker_normals_to_segments` (l. 736–847).
  - Fortran 2018 : `module_interface_vof.f90` $\rightarrow$ `advect_vof` (l. 792+), `compute_curvature_and_csf` (l. 627+).

### 4.6.3. Analyse de l'Implémentation CVI (`InterfaceMarkers.cpp`, lignes 145–185)

```cpp
int j_u = (yp < mesh.yc(j)) ? (j - 1) : j;
int i_v = (xp < mesh.xc(i)) ? (i - 1) : i;

double xi_u  = std::clamp((xp - mesh.xu(i))   / dx, 0.0, 1.0);
double eta_u = std::clamp((yp - mesh.yc(j_u)) / dy, 0.0, 1.0);
double xi_v  = std::clamp((xp - mesh.xc(i_v)) / dx, 0.0, 1.0);
double eta_v = std::clamp((yp - mesh.yv(j))   / dy, 0.0, 1.0);
// ... lecture des 4 sommets MAC u_{i..i+1, j_u..j_u+1} et v_{i_v..i_v+1, j..j+1} ...
double corr_u = 0.5 * (dx / dy) * xi_u  * (1.0 - xi_u)  * (dv_e - dv_w);
double corr_v = 0.5 * (dy / dx) * eta_v * (1.0 - eta_v) * (du_n - du_s);
up = u_bilin + corr_u;
vp = v_bilin + corr_v;
```

Notice sur les conventions de signes : dans `InterfaceReconstruction`, les segments PLIC $n_x^{\text{PLIC}} x + n_y^{\text{PLIC}} y + \alpha = 0$ ont leur normale orientée vers $C = 0$ ($-\nabla C / \|\nabla C\|$), tandis que les marqueurs et le modèle CSF orientent $\mathbf{n}_p$ vers $C = 1$ ($+\nabla C / \|\nabla C\|$). À chaque ensemencement et lors de l'injection `apply_marker_normals_to_segments`, le signe est rigoureusement inversé (`mk.nx = -seg.nx`, `snx = -fields.marker_nx[idx]`).

---

## 4.7. Lissage SVOF Helmholtz et Préconditionneur Bloc Triangulaire BT-SBP

### 4.7.1. Lissage SVOF Helmholtz (`smooth_svof_helmholtz`)

- **Équation discrète :** $(1 + a_E + a_W + a_N + a_S) C_{s,i,j} - \sum_{\text{nb}} a_{\text{nb}} C_{s,\text{nb}} = C_{i,j}$ avec $a_{E,W} = k \Delta y / \Delta x$, $a_{N,S} = k \Delta x / \Delta y$ (annulés sur les bords physiques de Neumann homogène) et 25 balayages Gauss-Seidel avant/arrière synchronisés en MPI.
- **Localisation :** `InterfaceTracker.cpp` (l. 288–380) et `module_interface_vof.f90` (l. 387–478).

### 4.7.2. Préconditionneur Triangulaire par Blocs BT-SBP

- **Équation discrète :** Pour le système monolithique $\begin{pmatrix} F_u & 0 & G_x \\ 0 & F_v & G_y \\ D_x & D_y & -C_p \end{pmatrix}$, `BlockTriangularPreconditioner` résout d'abord le complément de Schur $S \mathbf{r}_p = \mathbf{z}_p$ (8 balayages de Gauss-Seidel symétrique sur les coefficients $a_W, a_E = \Delta y^2 / u_{aP}$, $a_S, a_N = \Delta x^2 / v_{aP}$, $a_P = \sum a_{\text{nb}} + \frac{\chi_T}{\Delta t}\Delta x\Delta y + 10^{-12}\Delta x\Delta y$, avec retrait du mode constant en incompressible), puis calcule $\text{rhs}_u = \mathbf{z}_u - G_x \mathbf{r}_p$, $\text{rhs}_v = \mathbf{z}_v - G_y \mathbf{r}_p$ et résout $F_u \mathbf{r}_u = \text{rhs}_u$ et $F_v \mathbf{r}_v = \text{rhs}_v$ par 4 balayages SSOR.
- **Localisation :**
  - C++17 : `COVOF_2D-C++/src/Preconditioner.cpp` $\rightarrow$ `BlockTriangularPreconditioner::update` (l. 210–382) et `apply` (l. 384–579).
  - Fortran 2018 : `COVOF_2D-FORTRAN/src/module_navier_stokes.f90` $\rightarrow$ `ns_assemble_preconditioner_blocks` (l. 1287–1585) et `module_linear_solvers.f90` $\rightarrow$ `apply_preconditioner` (l. 47–489).
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
# Vérification, Validation (V&V), Parité Dual-Code et Assurance Qualité

Le développement d'un solveur bi-fluide à suivi d'interface (FC-WY-CVI-Marker + Navier-Stokes + Pénalisation de Brinkman) exige une discipline de **Vérification et Validation (V&V)** stricte et reproductible. Dans l'architecture **COVOF 2D**, la vérification repose sur deux piliers complémentaires : l'exécution automatisée du protocole à trois paliers `verify_cvi_markers.py` (Noyau mathématique, Physique capillaire/cinématique, Parité MPI), et les campagnes multi-cas orchestrées sous `VALIDATION_TOTALE/`. L'architecture étant strictement **dual-code** (`COVOF_2D-C++` en C++17 et `COVOF_2D-FORTRAN` en Fortran 2018), toute évolution d'algorithme doit préserver la parité numérique entre les deux implémentations.

---

## 6.1 Stratégie V&V en Trois Niveaux

L'architecture de vérification de COVOF 2D sépare la vérification mathématique pure (solveurs elliptiques, divergence discrète des marqueurs) de la validation physique multi-fluides et de la reproductibilité parallèle MPI :

| Niveau | Objet vérifié | Outils & Cibles réelles | Critère de succès |
| :--- | :--- | :--- | :--- |
| **Niveau 1 : Noyaux & Analytique** | Divergence CVI analytique & DF, opérateurs elliptiques, écoulements exacts | `verify_cvi_markers.py` (Palier 1), modes internes `"LAPLACIAN"`, `"CAS_3_7"`, `"CAS_3_8"`, `"VORTEX_ADVECTION"` | $\|\nabla\cdot\mathbf{u}\|_{\text{ana}} < 10^{-13}$, $\|\nabla\cdot\mathbf{u}\|_{\text{DF}} < 10^{-10}$, convergence $L^2 / L^\infty$ d'ordre 2 |
| **Niveau 2 : Parité Dual-Code & MPI** | Régression capillaire (Young-Laplace), advection non-linéaire (Rider-Kothe), décomposition MPI | `verify_cvi_markers.py` (Paliers 2 & 3), `covof2d_cpp` vs `covof2d_fortran` (`np=1` vs `mpirun -np 4`) | $\Delta\kappa_{\text{MPI}} < 10^{-5}$, invariance exacte `n_markers(np=1) == n_markers(np=4)` |
| **Niveau 3 : Campagnes V&V** | 14 cas canoniques + cas avancés (immersion Brinkman, diphasique compressible, jet visqueux) | `VALIDATION_TOTALE/Run_Valide-CASES/run_all_validations.sh`, `Run_CASES/`, `Run_FC-WY-CVI-Marker/` | Accord avec solutions analytiques / benchmarks (Faxén, Rayleigh-Plesset, Hinze, Popinet) |

> **Note sur l'infrastructure de test du dépôt :** COVOF 2D n'utilise pas de framework de test externe (ni Catch2, ni GoogleTest, ni pFUnit) et ne dispose pas de pipeline GitHub Actions sur le dépôt. Un binaire utilitaire compilé `COVOF_2D-C++/bin/test_solver` existe en complément, mais la suite officielle de non-régression pour le couplage marqueurs-VOF et le solveur dual-code repose sur les scripts Python/Shell et les modes de validation intégrés à `Simulation` / `module_simulation`.

---

## 6.2 Niveau 1 : Vérification Mathématique du Noyau CVI et Modes Analytiques Internes

### 6.2.1 Palier 1 de `verify_cvi_markers.py` : Invariance sans divergence sur boîtes MAC

Le script `Guides_Tutos/Run_FC-WY-CVI-Marker/verify_cvi_markers.py` implémente dans `test_palier_1_cvi_divergence()` une batterie stochastique de **5 000 blocs MAC $3 \times 3$ aléatoires à divergence nulle** (avec $\Delta x = 0{,}03125$ et $\Delta y = 0{,}025$) afin de comparer mathématiquement trois schémas d'interpolation de vitesse aux marqueurs :

1. **Interpolation bilinéaire standard (`BILINEAR`) :** ne satisfait pas $\nabla \cdot \mathbf{u} = 0$ en tout point continu d'une cellule.
2. **Interpolation CVI Variante A :** interpolation quadratique à divergence analytique nulle dans la boîte centrale, mais susceptible de présenter des sauts de composante normale aux bords.
3. **Interpolation CVI Variante B (Wang et al. 2015 / Pusok et al. 2017) :** interpolation utilisée par COVOF 2D (`marker_advection_scheme = "CVI_RK2"`), combinant continuité $C^0$ des composantes normales aux faces et exactitude sans divergence.

Pour chaque tirage aléatoire et chaque coordonnée locale $(\tilde{x}, \tilde{y})$ dans le domaine d'interpolation, le test vérifie simultanément :
$$\left|\frac{\partial u_{\text{CVI}}}{\partial x} + \frac{\partial v_{\text{CVI}}}{\partial y}\right|_{\text{analytique}} < 10^{-13}, \qquad \left|\frac{u_{\text{CVI}}(x+\delta,y)-u_{\text{CVI}}(x-\delta,y)}{2\delta} + \frac{v_{\text{CVI}}(x,y+\delta)-v_{\text{CVI}}(x,y-\delta)}{2\delta}\right| < 10^{-10}$$
où $\delta$ désigne le pas de différenciation centrée utilisé dans `test_palier_1_cvi_divergence()`. L'échec sur un seul des 5 000 blocs interrompt la suite dès le Palier 1.

### 6.2.2 Modes analytiques intégrés à `Simulation` et `module_simulation`

Plusieurs configurations d'étude (`simulation_mode` ou `case_name` dédiés dans le JSON) déclenchent une comparaison directe contre une solution exacte et exportent automatiquement des fichiers CSV métrologiques en fin de calcul :

* **Mode `"LAPLACIAN"` (`CAS_3_1` Dirichlet et `CAS_3_2_NEUMANN` Neumann) :** résout $-\Delta p = f$ sur maillage cartésien avec les solveurs et préconditionneurs du code (`PCG` / `BICGSTAB`, préconditionneurs `DIAG`, `ILU0`, `SSOR`, `AMG`), puis génère :
  * `laplacian_error.csv` : normes d'erreur discrète $L^2$ et $L^\infty$ par rapport à la solution manufacturée.
  * `solver_history.csv` : historique itératif du résidu relatif $\|r_k\|_2 / \|r_0\|_2$ à chaque itération Krylov.
* **Poiseuille stationnaire (`"CAS_3_7"`) et instationnaire (`"CAS_3_8"`) :** compare le profil de vitesse axiale $u(y)$ et le gradient de pression $\partial p / \partial x$ aux solutions analytiques de Poiseuille (régime établi parabolique et développement transitoire par série de Fourier). Exporte :
  * `poiseuille_error.csv` et `poiseuille_unsteady_error.csv` ;
  * `profile_u_y.csv`, `profile_p_x.csv`, `profile_x_mid.csv`, `profile_y_mid.csv` ;
  * `solution_field_2d.csv` et `solution_field.vtk`.
* **Advection tourbillonnaire (`"VORTEX_ADVECTION"`) :** évalue l'erreur géométrique $E_1$ et la conservation de masse $\Delta M / M_0$ de l'interface dans un champ de vitesse cinématique prescrit, en écrivant `vortex_error.csv`.

```bash
# Exécution directe d'un test analytique intégré en C++ et en Fortran
./COVOF_2D-C++/bin/covof2d_cpp VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_1_laplacian/input.json
./COVOF_2D-FORTRAN/bin/covof2d_fortran VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_1_laplacian/input.json
```

---

## 6.3 Niveau 2 : Parité Dual-Code et Invariance Parallèle MPI (`verify_cvi_markers.py`)

Les Paliers 2 et 3 de `verify_cvi_markers.py` pilotent simultanément les exécutables `COVOF_2D-FORTRAN/bin/covof2d_fortran` et `COVOF_2D-C++/bin/covof2d_cpp` sur des configurations générées dynamiquement.

### 6.3.1 Palier 2 : Goutte statique de Young-Laplace ($R = 0{,}25$, $\sigma = 1{,}0$)

Sur un domaine unitaire $[0,1]^2$ discrétisé en $32 \times 32$ puis $64 \times 64$, une goutte circulaire de rayon $R = 0{,}25$ (courbure théorique $\kappa_{\text{ex}} = 1/R = 4{,}0$, saut de pression théorique $\Delta P_{\text{ex}} = \sigma \kappa_{\text{ex}} = 4{,}0$) est simulée selon quatre variantes algorithmiques pour les deux langages :

1. **CSF brut** (`vof.use_markers = false`, lissage désactivé) ;
2. **CSF + SVOF** (`enable_svof_smoothing = true`, largeur de filtre $k = 0{,}5$) ;
3. **CVI + Marqueurs** (`vof.use_markers = true`, `marker_advection_scheme = "CVI_RK2"`, `marker_curvature = true`) ;
4. **CVI + Marqueurs + SVOF** (couplage complet).

Le script parse les sorties console et les fichiers `<case_name>_history.csv` (9 colonnes en C++, 11 colonnes en Fortran) pour vérifier quantitativement la réduction de l'erreur de courbure `k_l2_err` et des courants parasites `u_max` par rapport au CSF brut.

### 6.3.2 Palier 3 : Vortex de Rider-Kothe et invariance MPI (`np=1` vs `mpirun -np 4`)

Le Palier 3 teste d'une part l'advection non-linéaire réversible du vortex de cisaillement de Rider-Kothe ($64 \times 64$, $T = 1{,}0\text{ s}$) et d'autre part la cohérence parallèle stricte sur une grille $2 \times 2$ de sous-domaines MPI (`np=1` contre `mpirun -np 4`) sur Young-Laplace $32 \times 32$. Les assertions d'acceptation suivantes sont imposées :

```python
assert abs(yl_np1["k_mean"] - yl_np4["k_mean"]) < 1e-5
assert yl_np1["n_markers"] == yl_np4["n_markers"]
```

Cette double condition garantit que :
1. L'échange des marqueurs fantômes dans la bande $2{,}5 h$ (`exchange_ghost_markers`) produit rigoureusement le même ajustement parabolique local $(\tau, \eta)$ aux frontières entre rangs MPI.
2. Le protocole d'élagage/ré-ensemencement local et la migration MPI en 2 passes (`migrate_markers`) ne perdent ni ne dupliquent aucun marqueur lors des traversées d'interfaces inter-processeurs.

---

## 6.4 Niveau 3 : Campagnes `VALIDATION_TOTALE/`

### 6.4.1 Les 14 cas canoniques (`Run_Valide-CASES/run_all_validations.sh`)

Le script `VALIDATION_TOTALE/Run_Valide-CASES/run_all_validations.sh` exécute de façon systématique **14 cas de validation analytiques et canoniques** couvrant l'ensemble de la physique du solveur :

| Dossier dans `Run_Valide-CASES/` | Physique / Benchmark vérifié | Grandeurs mesurées |
| :--- | :--- | :--- |
| `validation_cas_1_laplacian` | Poisson avec conditions de Dirichlet (`CAS_3_1`) | Normes $L^2, L^\infty$ dans `laplacian_error.csv` |
| `validation_cas_2_neumann` | Poisson avec conditions de Neumann (`CAS_3_2_NEUMANN`) | Compatibilité noyau pression moyenne nulle |
| `validation_cas_3_poiseuille` | Écoulement de Poiseuille établi | Profil $u(y)$ (`profile_u_y.csv`), $\partial p/\partial x$ |
| `validation_cas_4_poiseuille_instat` | Poiseuille transitoire impulsionnel | Convergence vers solution série temporelle (`poiseuille_unsteady_error.csv`) |
| `validation_cas_5_poiseuille_visc_variable` | Écoulement de Poiseuille stratifié à viscosité variable | Saut de contrainte visqueuse à l'interface |
| `validation_cas_6_laplace_droplet` | Goutte statique de Young-Laplace | Saut $\Delta P = \sigma / R$ et courants parasites $u_{\max}$ |
| `validation_cas_7_vortex_shear` | Vortex de cisaillement (Rider-Kothe) | Conservation de masse et erreur de forme $E_1$ |
| `validation_cas_8_dense_cylinder` | Écoulement autour d'un cylindre rigide (Brinkman) | $C_D$, $C_L$, longueur de recirculation |
| `validation_cas_9_sedimentation_faxen` | Sédimentation confinée & correction de paroi de Faxén | Vitesse terminale vs formule de Faxén |
| `validation_cas_10_piston_compressible` | Piston compressible 1D (EOS compressible / NASG) | Compression isentropique / onde acoustique |
| `validation_cas_11_compression_bulle` | Compression/collapsus de bulle compressible | Dynamique de Rayleigh-Plesset |
| `validation_cas_12_rupture_barrage` | Rupture de barrage diphasique (*dam break*) | Position du front $x_{\text{front}}(t)$ et impact pariétal |
| `validation_cas_13_oscillation_bulle_capillaire` | Oscillation capillaire d'une bulle elliptique | Fréquence propre $\omega_2$ de Lamb/Rayleigh |
| `validation_cas_14_instabilite_interface_gravitaire` | Instabilité gravitaire de Rayleigh-Taylor | Taux de croissance et morphologie champignon |

### 6.4.2 Cas avancés (`Run_CASES/`) et suite dédiée `Run_FC-WY-CVI-Marker/`

En complément des 14 cas de base, le répertoire `VALIDATION_TOTALE/Run_CASES/` regroupe des configurations industrielles et multi-physiques plus complexes, notamment :
* `validation_cas_9_50x100_newton_dual_code` (comparaison directe C++ vs Fortran) ;
* `validation_cas_10_liquid_sheet` (atomisation/étirement de nappe liquide) ;
* `validation_cas_14_marche_descendante` à `validation_cas_19_convection_poreux_darcy_brinkman` (écoulement sur marche descendante, cavité thermique, convection en milieu poreux Darcy-Brinkman) ;
* `validation_cas_21_ascension_bulle_2d` (cas test de Hysing et al. : bulle ascendante, vitesse d'ascension $V_c$ et circularité) ;
* `validation_cas_24_viscous_jet_buckling` (flambement gravitationnel de jet visqueux haute contrainte).

La suite documentée `Guides_Tutos/Run_FC-WY-CVI-Marker/` contient 8 études ciblées (`01_cas_goutte_statique_young_laplace_cvi_markers.md` à `08_cas_flambement_jet_visqueux_cvi_markers.md`) servant de référence pour valider l'activation des marqueurs lagrangiens face au CSF classique.

---

## 6.5 Protocole de Convergence de Maillage (MMS) et Vérification de Parité VTR

Lors de l'ajout d'un nouvel opérateur discret, le développeur applique la **Méthode des Solutions Manufacturées (MMS)** sur une séquence de maillages raffinés par facteur 2 ($N \in \{32, 64, 128, 256\}$). L'ordre de convergence observé est calculé entre deux maillages consécutifs de pas $2h$ et $h$ par :

$$p_{\text{obs}} = \frac{\ln\left( \|e_{2h}\|_{L^2} / \|e_h\|_{L^2} \right)}{\ln 2}$$

Un schéma centré MAC d'ordre 2 doit exhiber $p_{\text{obs}} \in [1{,}90 ; 2{,}05]$ en norme $L^2$ pour des champs lisses.

Pour valider la **parité stricte C++17 $\leftrightarrow$ Fortran 2018** sur un cas diphasique complet, on compare en post-traitement :
1. Les fichiers `<case_name>_history.csv` sur les 9 colonnes communes (`step`, `time`, `dt`, `max_u`, `max_v`, `max_div`, `mass_fluid1`, `kinetic_energy`, `solver_iters`).
2. Les tableaux CellData (`PhaseIndicator`, `PhaseIndicator_Smoothed`, `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Curvature`, `Velocity`) des fichiers `<case_name>_p<rank>_<step:06d>.vtr` à l'aide de `pyvista` ou `vtk` en Python.

---

## 6.6 Audits Mémoire, Sanitizers CMake et Pré-validation Locale

Avant tout commit ou fusion sur le code principal, la vérification de la robustesse mémoire et de l'absence de comportements indéfinis (dépassements de halo `0..nx+1`, fuites de tampons MPI, divisions par zéro dans `InterfaceMarkers`) s'effectue via des configurations CMake dédiées en mode Debug :

### 6.6.1 Build C++17 avec AddressSanitizer et UndefinedBehaviorSanitizer

```bash
cmake -S COVOF_2D-C++ -B build_asan \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -Wall -Wextra"
cmake --build build_asan -j4

# Exécution sous MPI pur (OpenMPI)
mpirun -np 2 ./build_asan/covof2d_cpp VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_6_laplace_droplet/input.json
```

### 6.6.2 Build Fortran 2018 avec vérifications de bornes et IEEE

```bash
cmake -S COVOF_2D-FORTRAN -B build_fdebug \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_Fortran_FLAGS="-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall -Wextra"
cmake --build build_fdebug -j4

mpirun -np 2 ./build_fdebug/covof2d_fortran VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_6_laplace_droplet/input.json
```

### 6.6.3 Audit Valgrind et checklist avant commit

Pour inspecter d'éventuelles lectures hors-limites dans les routines d'échange de halos ou les allocations de marqueurs :

```bash
mpirun -np 2 valgrind --leak-check=full --track-origins=yes \
  ./build_asan/covof2d_cpp VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_1_laplacian/input.json
```

**Checklist minimale avant validation d'une modification :**
1. Compilation sans avertissement bloquant des deux cibles `covof2d_cpp` (C++17, `-Wall -Wextra -Wno-unused-parameter`) et `covof2d_fortran` (Fortran 2018).
2. Passage intégral de `python3 Guides_Tutos/Run_FC-WY-CVI-Marker/verify_cvi_markers.py` (Paliers 1, 2 et 3).
3. Exécution de `bash VALIDATION_TOTALE/Run_Valide-CASES/run_all_validations.sh` sans divergence ni régression sur les 14 cas canoniques.
4. Vérification de non-fuite mémoire sous `mpirun -np 4` en mode Debug (`-fsanitize=address` / `-fcheck=all`) sur au moins un cas avec `vof.use_markers = true`.
# Standards de Codage, Profilage HPC et Guide d'Extension

Ce chapitre rassemble les conventions d'ingénierie logicielle appliquées dans **COVOF 2D** (`COVOF_2D-C++` et `COVOF_2D-FORTRAN`), les bonnes pratiques de profilage sous parallélisme **MPI pur**, ainsi qu'un tutoriel complet pas à pas montrant comment étendre simultanément les deux codes avec une nouvelle équation de transport (un **scalaire passif** $S(\mathbf{x}, t)$).

---

## 7.1 Standards de Codage C++17 et Fortran 2018

### 7.1.1 Règles C++17 (`COVOF_2D-C++`)

Le code C++17 est compilé avec `-std=c++17 -O3 -Wall -Wextra -Wno-unused-parameter` (et `-framework Accelerate` sous macOS). Tous les fichiers d'en-tête résident dans `include/covof/` (15 fichiers `.hpp`) et leurs implémentations dans `src/` (13 fichiers `.cpp`).

1. **Typage explicite et conventions de nommage :**
   * Enveloppe unique sous `namespace covof { ... }`.
   * Types et classes en `PascalCase` (`Mesh2D`, `FlowFields`, `NavierStokesSolver`, `InterfaceMarkers`, `BoundaryCondition`, `SolverConfig`).
   * Méthodes et fonctions en `snake_case` (`cell_idx`, `u_idx`, `v_idx`, `sync_cell_scalar`, `compute_time_step`, `load_from_json`).
   * Membres de configuration en `snake_case` alignés sur les clés JSON (`enable_two_phase`, `use_markers`, `svof_filter_width`).
2. **Gestion mémoire et contiguïté cache :**
   * Zéro `new`/`delete` nus : utilisation systématique de `std::vector<double>`, `std::vector<int>` et pointeurs intelligents (`std::unique_ptr`).
   * Indexation linéarisée avec **une seule couche de cellules fantômes** (`0 .. nx_local+1`, `0 .. ny_local+1`, soit `nx_ghost = nx_local + 2`) via `mesh.cell_idx(i, j)`, `mesh.u_idx(i, j)` et `mesh.v_idx(i, j)`.
   * Toujours boucler avec **`j` en boucle externe et `i` en boucle interne** afin de parcourir contiguëment les tableaux Row-Major (`i + nx_ghost * j`) et favoriser la vectorisation automatique du compilateur (`-O3`).
3. **Absence d'OpenMP et discipline MPI :**
   * Le solveur repose **exclusivement sur MPI pur** (`MPI::MPI_CXX`). N'introduisez pas de directives `#pragma omp` ni de variables globales mutables dans les noyaux numériques.
   * Tout accès aux données voisines inter-domaines passe par `MpiComm` (`sync_cell_scalar`, `sync_u_face`, `sync_v_face`, `exchange_ghost_markers`, `migrate_markers`).

### 7.1.2 Règles Fortran 2018 (`COVOF_2D-FORTRAN`)

Le code Fortran 2018 est compilé avec `-O3 -march=native` et regroupé en 16 fichiers `.f90` dans `src/`, chacun commençant par `module ...` et :

1. **`implicit none` obligatoire** en tête de chaque module et sous-programme.
2. **Précision portable via `module_precision` :** toute variable flottante doit être déclarée `real(dp)` (et non `real*8` ni `double precision`), où `dp` provient de `use module_precision, only: dp`.
3. **Ordre mémoire Column-Major (inverse du C++) :**
   * Les tableaux 2D de champs sont alloués en `(0:nx+1, 0:ny+1)` (avec `nx = mesh%nx_local`, `ny = mesh%ny_local`).
   * En Fortran, le premier indice (`i`) est contigu en mémoire : écrivez systématiquement **`do j = 1, ny` en externe et `do i = 1, nx` en interne**.
4. **Encapsulation et imports explicites :**
   * Toujours utiliser `use module_xyz, only: ...` pour éviter la pollution d'espace de noms.
   * Les routines MPI doivent systématiquement vérifier l'argument `ierr` ou déléguer à `module_mpi_comm` (`exchange_halos_scalar`, `exchange_halos_u`, `exchange_halos_v`, `mpi_global_max`, `mpi_global_sum`).

---

## 7.2 Profilage et Optimisation HPC sous MPI Pur

Lorsqu'une évolution dégrade les temps de calcul, l'analyse de performance suit un protocole rigoureux avant toute micro-optimisation.

### 7.2.1 Chronométrage granulaire intégré et profilage système

Les deux solveurs affichent le temps de calcul par pas de temps et les compteurs d'itérations Krylov dans `<case_name>_history.csv` (`solver_iters`). Pour localiser précisément les points chauds CPU (construction matricielle CSR, solveur PCG/AMG, advection VOF/marqueurs ou barrières MPI), utilisez :

* **Linux (`perf` + FlameGraph) sur un rang MPI ciblé :**
  ```bash
  mpirun -np 4 perf record -F 99 -g -- ./COVOF_2D-C++/bin/covof2d_cpp cases/benchmark.json
  perf report
  ```
* **macOS (`Instruments` / `sample`) sur exécution locale :**
  ```bash
  mpirun -np 4 ./COVOF_2D-C++/bin/covof2d_cpp cases/benchmark.json &
  sample <PID_RANG_0> 10 -file profile_rank0.txt
  ```
* **Analyse des communications MPI (Score-P / Scalasca / OpenMPI timers) :**
  Permet de détecter un déséquilibre de charge (*load imbalance*) lorsque `markers.enabled = true`, par exemple si une interface déformée concentre la majorité des marqueurs CVI sur un seul rang MPI de la grille cartésienne `px * py`.

### 7.2.2 Pièges classiques de performance dans COVOF 2D

| Symptôme observé | Cause fréquente dans COVOF 2D | Remède recommandé |
| :--- | :--- | :--- |
| Chute de bande passante mémoire sur grilles fines | Inversion d'indices de boucle (`i` externe en Fortran ou `j` interne en C++) | Respecter `j` externe / `i` interne dans les deux langages |
| Explosion du coût Krylov (`solver_iters > 500`) | Changement brutal de $\rho_1/\rho_2$ ou préconditionneur `DIAG` utilisé sur AMG disponible | Basculer `solver.preconditioner` vers `"AMG"` ou activer `enable_svof_smoothing` |
| Surcoût temporel dans `InterfaceMarkers` | Ré-allocations répétées ou nombre excessif de marqueurs par cellule | Ajuster `marker_max_per_cell` (défaut `24`) et `markers_per_cell` (défaut `4`) |
| Ralentissement à fort nombre de rangs `np` | Domaine local `nx_local x ny_local` trop petit (< $16 \times 16$) par rapport au coût des halos MPI | Maintenir au moins $32 \times 32$ cellules internes par cœur MPI |

---

## 7.3 Tutoriel d'Extension : Ajouter le Transport d'un Scalaire Passif $S(\mathbf{x}, t)$

Nous illustrons maintenant la procédure complète pour ajouter dans **COVOF 2D** une équation d'advection-diffusion d'un traceur passif (par exemple une concentration chimique adimensionnée $S \in [0, 1]$) :

$$\frac{\partial S}{\partial t} + \nabla \cdot (\mathbf{u} S) = D_S \nabla^2 S$$

Cette extension touche de manière symétrique les 6 couches architecturales du solveur : Configuration JSON, Structures de Champs, Synchronisation MPI, Noyau Numérique, Boucle Temporelle et Sorties VTK.

---

### Étape 1 : Ajouter les paramètres dans la Configuration (`Config` / `module_config`)

#### C++17 (`include/covof/Config.hpp` et `src/Config.cpp`)

Dans `include/covof/Config.hpp`, ajoutez les champs dans `ModelConfig` (ou `SolverConfig` selon la section choisie ; ici dans `ModelConfig`, bloc JSON `"physics"`) :

```cpp
struct ModelConfig {
    // ... membres existants ...
    bool enable_scalar_transport = false;
    double scalar_diffusivity = 1.0e-5;
};
```

Dans `src/Config.cpp`, à l'intérieur de `Config::load_from_json` lors du parsing du bloc `"physics"` (`model_block`), parsez les clés avec les helpers internes de `Config.cpp` (COVOF 2D utilise son propre parseur JSON léger interne, **pas** `nlohmann::json`) :

```cpp
model.enable_scalar_transport = get_bool(model_block, "enable_scalar_transport", false);
model.scalar_diffusivity = get_double(model_block, "scalar_diffusivity", 1.0e-5);
```

Pensez également à afficher ces deux paramètres dans `Config::print()` (sur le rang 0).

#### Fortran 2018 (`src/module_config.f90`)

Dans `type :: config_t`, ajoutez les champs correspondants et lisez-les dans `load_config_from_json` via `module_json` :

```fortran
type :: config_t
    ! ... champs existants ...
    logical  :: enable_scalar_transport = .false.
    real(dp) :: scalar_diffusivity = 1.0e-5_dp
end type config_t
```

```fortran
! Dans load_config_from_json(cfg, filename)
call json_get_bool(buf, 'physics.enable_scalar_transport', cfg%enable_scalar_transport, .false.)
call json_get_real(buf, 'physics.scalar_diffusivity',    cfg%scalar_diffusivity,    1.0e-5_dp)
```

---

### Étape 2 : Allouer les tableaux de champs avec halo de 1 cellule (`Field` / `module_fields`)

#### C++17 (`include/covof/Field.hpp` et `src/Field.cpp`)

Dans la structure `FlowFields`, déclarez les vecteurs linéarisés `S` et `S_old` :

```cpp
struct FlowFields {
    // ... champs existants (u, v, p, C, T, rho, mu, ...) ...
    std::vector<double> S;
    std::vector<double> S_old;

    void allocate(const Mesh2D& mesh);
    void reset();
};
```

Dans `FlowFields::allocate(const Mesh2D& mesh)` (`src/Field.cpp`), allouez exactement `n_cells = mesh.nx_ghost * mesh.ny_ghost` éléments :

```cpp
const std::size_t n_cells = static_cast<std::size_t>(mesh.nx_ghost) * mesh.ny_ghost;
S.assign(n_cells, 0.0);
S_old.assign(n_cells, 0.0);
```

#### Fortran 2018 (`src/module_fields.f90`)

Dans `type :: flow_fields_t`, ajoutez les champs allouables `(0:nx+1, 0:ny+1)` et mettez à jour `allocate_fields`, `deallocate_fields` et `reset_fields` :

```fortran
type :: flow_fields_t
    ! ... champs existants ...
    real(dp), allocatable :: S(:,:)
    real(dp), allocatable :: S_old(:,:)
end type flow_fields_t
```

```fortran
! Dans allocate_fields(fields, mesh)
allocate(fields%S(0:mesh%nx_local+1, 0:mesh%ny_local+1))
allocate(fields%S_old(0:mesh%nx_local+1, 0:mesh%ny_local+1))
fields%S     = 0.0_dp
fields%S_old = 0.0_dp
```

---

### Étape 3 : Implémenter le noyau d'advection-diffusion (`NavierStokesSolver` / `module_navier_stokes`)

Sur le maillage décalé MAC, les vitesses normales aux faces $(i+1/2, j)$ et $(i, j+1/2)$ sont stockées respectivement dans `u(i, j)` (`u_idx(i, j)`) et `v(i, j)` (`v_idx(i, j)`). En conservant un schéma centré pour la diffusion et *upwind* (ou TVD) sur la vitesse faciale, le flux advectif-diffusif garantit la conservation stricte de $\int_\Omega S \, d\Omega$.

#### C++17 (`include/covof/NavierStokesSolver.hpp` et `src/NavierStokesSolver.cpp`)

Déclarez la méthode membre publique dans `NavierStokesSolver` :

```cpp
void advance_passive_scalar(FlowFields& fields, double dt) const;
```

Implémentez-la dans `src/NavierStokesSolver.cpp` en respectant l'ordre cache `j` externe / `i` interne et **sans directive OpenMP** :

```cpp
void NavierStokesSolver::advance_passive_scalar(FlowFields& fields, double dt) const {
    if (!config_.model.enable_scalar_transport) return;

    const double D = config_.model.scalar_diffusivity;
    const double dx = mesh_.dx;
    const double dy = mesh_.dy;
    const double inv_dx = 1.0 / dx;
    const double inv_dy = 1.0 / dy;

    fields.S_old = fields.S;

    for (int j = 1; j <= mesh_.ny_local; ++j) {
        for (int i = 1; i <= mesh_.nx_local; ++i) {
            const int c  = mesh_.cell_idx(i, j);
            const int cw = mesh_.cell_idx(i - 1, j);
            const int ce = mesh_.cell_idx(i + 1, j);
            const int cs = mesh_.cell_idx(i, j - 1);
            const int cn = mesh_.cell_idx(i, j + 1);

            const double u_e = fields.u[mesh_.u_idx(i, j)];
            const double u_w = fields.u[mesh_.u_idx(i - 1, j)];
            const double v_n = fields.v[mesh_.v_idx(i, j)];
            const double v_s = fields.v[mesh_.v_idx(i, j - 1)];

            // Flux advectifs upwind sur faces MAC
            const double fe = (u_e >= 0.0) ? u_e * fields.S_old[c]  : u_e * fields.S_old[ce];
            const double fw = (u_w >= 0.0) ? u_w * fields.S_old[cw] : u_w * fields.S_old[c];
            const double fn = (v_n >= 0.0) ? v_n * fields.S_old[c]  : v_n * fields.S_old[cn];
            const double fs = (v_s >= 0.0) ? v_s * fields.S_old[cs] : v_s * fields.S_old[c];

            const double adv = (fe - fw) * inv_dx + (fn - fs) * inv_dy;
            const double diff = D * (
                (fields.S_old[ce] - 2.0 * fields.S_old[c] + fields.S_old[cw]) * inv_dx * inv_dx +
                (fields.S_old[cn] - 2.0 * fields.S_old[c] + fields.S_old[cs]) * inv_dy * inv_dy
            );

            fields.S[c] = fields.S_old[c] + dt * (diff - adv);
        }
    }
}
```

#### Fortran 2018 (`src/module_navier_stokes.f90`)

```fortran
subroutine ns_advance_passive_scalar(mesh, cfg, fields, dt)
    type(mesh_t),        intent(in)    :: mesh
    type(config_t),      intent(in)    :: cfg
    type(flow_fields_t), intent(inout) :: fields
    real(dp),            intent(in)    :: dt

    integer  :: i, j, nx, ny
    real(dp) :: D, inv_dx, inv_dy, ue, uw, vn, vs, fe, fw, fn, fs, adv, diff

    if (.not. cfg%enable_scalar_transport) return

    nx = mesh%nx_local; ny = mesh%ny_local
    D = cfg%scalar_diffusivity
    inv_dx = 1.0_dp / mesh%dx; inv_dy = 1.0_dp / mesh%dy

    fields%S_old = fields%S

    do j = 1, ny
        do i = 1, nx
            ue = fields%u(i, j);     uw = fields%u(i - 1, j)
            vn = fields%v(i, j);     vs = fields%v(i, j - 1)

            fe = merge(ue * fields%S_old(i, j),     ue * fields%S_old(i + 1, j), ue >= 0.0_dp)
            fw = merge(uw * fields%S_old(i - 1, j), uw * fields%S_old(i, j),     uw >= 0.0_dp)
            fn = merge(vn * fields%S_old(i, j),     vn * fields%S_old(i, j + 1), vn >= 0.0_dp)
            fs = merge(vs * fields%S_old(i, j - 1), vs * fields%S_old(i, j),     vs >= 0.0_dp)

            adv  = (fe - fw) * inv_dx + (fn - fs) * inv_dy
            diff = D * ( (fields%S_old(i+1,j) - 2.0_dp*fields%S_old(i,j) + fields%S_old(i-1,j)) * inv_dx**2 &
                       + (fields%S_old(i,j+1) - 2.0_dp*fields%S_old(i,j) + fields%S_old(i,j-1)) * inv_dy**2 )

            fields%S(i, j) = fields%S_old(i, j) + dt * (diff - adv)
        end do
    end do
end subroutine ns_advance_passive_scalar
```

> **Note de stabilité :** Si `enable_scalar_transport = true`, ajoutez également la contrainte de Fourier explicite $\Delta t \le \frac{1}{2 D_S (\Delta x^{-2} + \Delta y^{-2})}$ dans `Simulation::compute_time_step()` (`src/Simulation.cpp`) et `sim_compute_dt` (`src/module_simulation.f90`).

---

### Étape 4 : Synchronisation MPI et Conditions aux Limites (`Simulation` / `module_simulation`)

Dans la boucle temporelle `Simulation::run()` (`src/Simulation.cpp`) et `sim_run` (`src/module_simulation.f90`), après l'avancement de Navier-Stokes (et avant/après l'avancement de l'énergie $T$), appelez la synchronisation des halos MPI et l'application des conditions aux limites scalaires :

#### C++17 (`src/Simulation.cpp`)

```cpp
if (config_.model.enable_scalar_transport) {
    ns_solver_->advance_passive_scalar(fields_, dt);
    mpi_comm_->sync_cell_scalar(mesh_, fields_.S);
    // Appliquer les conditions aux limites scalaires (Neumann nul ou Dirichlet) sur les bords globaux
}
```

#### Fortran 2018 (`src/module_simulation.f90`)

```fortran
if (cfg%enable_scalar_transport) then
    call ns_advance_passive_scalar(mesh, cfg, fields, dt)
    call exchange_halos_scalar(fields%S, mesh)
    ! Appliquer les conditions aux limites scalaires sur les frontieres physiques
end if
```

---

### Étape 5 : Exporter le champ vers ParaView (`VTKWriter` / `module_vtk_writer`)

Les fichiers de maillage rectiligne distribué utilisent le format XML VTK `<case_name>_p<rank>_<step:06d>.vtr` agrégé par `<case_name>_<step:06d>.pvtr`.

1. **Dans `src/VTKWriter.cpp` :** à l'intérieur de `VTKWriter::write_step(...)`, ajoutez le `<DataArray type="Float64" Name="PassiveScalar" format="ascii">` dans `<CellData>` de chaque fichier `.vtr`, ainsi que `<PDataArray type="Float64" Name="PassiveScalar"/>` dans `<PCellData>` du fichier `.pvtr` écrit par le rang 0.
2. **Dans `src/module_vtk_writer.f90` :** ajoutez la même balise dans `vtk_write_step` (routine `.vtr` et fichier maître `.pvtr`) afin de garantir que la lecture dans ParaView ou `pyvista` produise exactement la même structure de champs entre les exécutables C++ et Fortran.

---

### Étape 6 : Vérifier analytiquement la nouvelle fonctionnalité

Avant de clôturer l'extension, créez un cas de test dans `VALIDATION_TOTALE/` (avec `u = v = 0`, $S(x,y,0) = \sin(\pi x)\sin(\pi y)$ sur $[0,1]^2$, conditions de Dirichlet nulles aux quatre bords). La solution analytique exacte est :

$$S_{\text{ex}}(x, y, t) = \exp\left(-2 \pi^2 D_S t\right) \sin(\pi x)\sin(\pi y)$$

1. Exécutez les deux exécutables sous `mpirun -np 4` sur une grille $64 \times 64$ et $128 \times 128$.
2. Vérifiez par un script Python (`pyvista`) lisant `<case_name>_<step:06d>.pvtr` que :
   * L'ordre de convergence spatiale vérifie $p_{\text{obs}} \approx 2{,}0$ ;
   * L'écart maximal entre `covof2d_cpp` et `covof2d_fortran` sur le champ `PassiveScalar` reste inférieur à $10^{-11}$ ;
   * Les builds Debug (`-fsanitize=address,undefined` en C++ et `-fcheck=all` en Fortran) ne signalent aucun accès hors-bornes sur `0..nx+1, 0..ny+1`.
