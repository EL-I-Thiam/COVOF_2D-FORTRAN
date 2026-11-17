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
