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
