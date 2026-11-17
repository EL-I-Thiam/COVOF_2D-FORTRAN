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
