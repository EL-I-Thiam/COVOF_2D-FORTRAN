# Introduction Générale et Capacités du Solveur

## Présentation et Philosophie de COVOF 2D

Le code **COVOF 2D** (*Coupled Volume-Of-Fluid 2D*) est un environnement logiciel de recherche et d'ingénierie avancée dédié à la simulation numérique directe (DNS) et à la modélisation à haute fidélité des écoulements de fluides complexes, incompressibles ou faiblement compressibles (*All-Mach*), monophasiques ou diphasiques à interface nette. Issu des développements méthodologiques conduits dans le lignage des codes académiques de référence du laboratoire **I2M (UMR CNRS 5295, Université de Bordeaux)** — en particulier les formulations monolithiques couplées, le Lagrangien augmenté adaptatif (**AAL / C-AAL**), le lissage rhéologique **SVOF** (*Smoothed VOF*) et la méthode hybride **FC-WY-CVI-Marker** — COVOF 2D constitue simultanément un banc d'essai algorithmique rigoureux et un outil de production scientifique.

### Paradigme de Conception : L'Architecture Duale Jumelle C++17 / Fortran 2018

Contrairement aux logiciels CFD traditionnels reposant sur un seul langage ou encapsulant du code hérité via des passerelles opaques, COVOF 2D est structuré selon un principe de **double implémentation miroir intégrale et indépendante** :

1. **La branche `COVOF_2D-C++` (15 en-têtes `.hpp`, 13 sources `.cpp`) :** Écrite en C++ moderne conforme au standard **ISO C++17**, elle privilégie une conception modulaire orientée objet au sein de l'espace de nommage `covof`, avec gestion RAII des tampons (`std::vector<double>`), polymorphisme dynamique pour les solveurs de Krylov (`LinearSolverBase`, `SolverFactory`) et les préconditionneurs (`PreconditionerBase`, `PreconditionerFactory`), et matrice creuse locale CSR (`CSRMatrix`) pour les blocs de préconditionnement. Elle produit l'exécutable `bin/covof2d_cpp` (et le lien symbolique `bin/covof2d`).
2. **La branche `COVOF_2D-FORTRAN` (16 fichiers sources `.f90`) :** Écrite en **ISO Fortran 2018**, elle implémente rigoureusement les mêmes opérateurs discrets au sein de 15 modules fortement typés (`implicit none`, encapsulation `private` par défaut, double précision IEEE 754 `dp = selected_real_kind(15, 307)`) et du programme principal `covof2d_main.f90`. Elle exploite l'ordonnancement mémoire natif *column-major* `(0:nx+1, 0:ny+1)` pour la vectorisation SIMD des boucles internes. Elle produit l'exécutable `bin/covof2d_fortran`.

Cette architecture duale vise la **parité numérique stricte** entre deux implémentations autonomes :
- Les deux exécutables lisent le **même fichier de configuration JSON** (tolérant les commentaires `//` en ligne).
- Les deux branches partagent la **même décomposition de domaine cartésienne 2D en MPI pur** (`MPI_Cart_create`, sans OpenMP), avec une couche unique de cellules fantômes (*halo*) de largeur 1 (`0:nx_local+1`, `0:ny_local+1`).
- Les deux solveurs exportent les **mêmes formats de sortie normalisés** : maillages rectilignes XML VTK parallèles (`<case_name>_p<rank>_<step:06d>.vtr` et fichier maître `<case_name>_<step:06d>.pvtr`), nuage de marqueurs VTK PolyData (`<case_name>_markers_<step:06d>.vtk`) et historique temporel CSV (`<case_name>_history.csv`).

```text
              +---------------------------------------------+
              |        Fichier de Configuration JSON        |
              |     (mon_cas.json, commentaires // admis)   |
              +---------------------------------------------+
                                     |
                     +---------------+---------------+
                     |                               |
                     v                               v
  +-------------------------------------+ +-------------------------------------+
  | Exécutable C++17 (MPI pur)          | | Exécutable Fortran 2018 (MPI pur)   |
  | COVOF_2D-C++/bin/covof2d_cpp        | | COVOF_2D-FORTRAN/bin/covof2d_fortran|
  +-------------------------------------+ +-------------------------------------+
  | - Config, Mesh2D, FlowFields, EOS   | | - module_config, module_mesh,       |
  | - NavierStokesSolver + CSRMatrix    | |   module_fields, module_eos         |
  | - InterfaceTracker +                | | - module_navier_stokes            |
  |   InterfaceReconstruction           | | - module_interface_vof +            |
  | - InterfaceMarkers (CVI-RK2)        | |   module_interface_reconstruction   |
  | - LinearSolver & Preconditioner     | | - module_markers (CVI-RK2)          |
  | - MpiComm & VTKWriter               | | - module_linear_solvers             |
  |                                     | | - module_mpi_comm & module_vtk_writer|
  +-------------------------------------+ +-------------------------------------+
                     |                               |
                     +---------------+---------------+
                                     |
                                     v
              +---------------------------------------------+
              | Sorties Numériques Unifiées                 |
              | - Champs volumiques : <case>_<step>.pvtr    |
              |   + sous-domaines   : <case>_p<r>_<step>.vtr|
              | - Marqueurs CVI     : <case>_markers_*.vtk  |
              | - Métriques         : <case>_history.csv    |
              | - Erreurs analytiques (laplacian/poiseuille/|
              |   vortex_error.csv, profiles CSV)           |
              +---------------------------------------------+
```

### Domaines d'Application Visés

COVOF 2D couvre un large spectre de configurations académiques et industrielles en mécanique des fluides 2D et transferts couplés :

- **Hydrodynamique diphasique à forts contrastes ($\rho_1/\rho_2 \sim 10^3$, $\mu_1/\mu_2 \sim 10^5$) :** Rupture de barrage (*dam break*), chute libre d'un cylindre dense, instabilités de Rayleigh-Taylor et de Kelvin-Helmholtz sur nappe liquide mince, ascension de bulles déformables à jupes fines, flambement de jets visqueux injectés (`JET_INLET`).
- **Dynamique capillaire et microfluidique :** Équilibre de Young-Laplace sur goutte statique, oscillations capillaires amorties de Lamb, réduction drastique des courants parasites et suivi sous-maille des ligaments par la méthode hybride **FC-WY-CVI-Marker**.
- **Écoulements faiblement compressibles et acoustique *All-Mach* :** Propagation d'ondes acoustiques et compression de bulles ou poches de gaz sous loi d'état linéaire, gaz parfait ou *Noble-Abel Stiffened Gas* (**NASG**) avec couplage **C-AAL** et transport inertiel **SSP-RK3**.
- **Convection thermique et approximation de Boussinesq :** Cavités différentiellement chauffées, panaches thermiques et instabilités de Rayleigh-Bénard avec calcul automatique du nombre de Nusselt moyen à la paroi chaude.
- **Écoulements en milieux poreux et obstacles immergés :** Modélisation d'obstacles rectangulaires (`BOX`) ou cylindriques (`CYLINDER`) fixes sans remaillage par pénalisation volumique de Darcy-Brinkman (marche descendante, sillage de cylindre, plaque mince, milieu poreux).

---

## Matrice Complète des Capacités Physiques et Numériques

Le tableau ci-dessous recense l'ensemble des modèles, discrétisations, solveurs et options réellement implémentés et vérifiés dans les sources C++17 et Fortran 2018 de COVOF 2D.

| Catégorie | Fonctionnalité / Modèle | Mots-Clés / Options Exactes | Description Technique et Références |
| :--- | :--- | :--- | :--- |
| **Topologie & Maillage** | Grille de calcul | Grille Cartésienne Structurée Uniforme 2D | Maillage décalé MAC (*Marker-and-Cell*, Harlow & Welch 1965) sur $[0, L_x] \times [0, L_y]$ |
| | Cellules fantômes | 1 couche périphérique (`0..nx+1`, `0..ny+1`) | `nx_ghost = nx_local + 2`, `ny_ghost = ny_local + 2` |
| | Parallélisme distribué | MPI pur 2D (`MPI_Cart_create`) | Partitionnement automatique $(P_x, P_y)$ via `MPI_Dims_create` ; 4 voisins cardinaux (`MPI_Isend`/`MPI_Irecv`) |
| **Modes d'Exécution** | Navier-Stokes complet | `mode = "NAVIER_STOKES"` | Boucle instationnaire couplée (monophasique, diphasique, thermique, compressible) |
| | Vérification elliptique | `mode = "LAPLACIAN"` | Résolution analytique $-\Delta \psi = f$ (Dirichlet `CAS_3_1` et Neumann `CAS_3_2_NEUMANN`) |
| | Cinématique d'interface | `mode = "VORTEX_ADVECTION"` | Tourbillon cisaillé réversible de Rider-Kothe sans résolution dynamique |
| **Régimes & Thermodynamique** | Incompressible | `compressible = false` | Contrainte solénoïdale $\nabla \cdot \mathbf{u} = 0$ |
| | Compressible *All-Mach* | `compressible = true` ou `flow_regime = "ALL_MACH"` | Continuité barotrope $\chi_T \partial_t p + \nabla \cdot \mathbf{u} = 0$ + énergie totale $\rho e$ (El Ouafa et al., *Fluids* 2024) |
| | Lois d'état (`eos`) | `INCOMPRESSIBLE`, `LINEAR`, `IDEAL_GAS`, `NASG` | Compressibilité isotherme $\chi_T$, $\rho(p,T)$ et vitesse du son $c_s$ (`EOS.hpp` / `module_eos.f90`) |
| | Thermique (Boussinesq) | `enable_energy = true` | Équation en température $T$ avec flottabilité $-\rho_0 \beta (T - T_{\text{ref}})\mathbf{g}$ (ou énergie totale si compressible) |
| | Obstacles & Milieux poreux | `obstacles = [{type: "BOX" | "CYLINDER"}]` | Pénalisation de Darcy-Brinkman $-\frac{\mu}{K}\mathbf{u}$ ($K \le 10^{-12}$ : solide rigide adhérent) |
| **Suivi d'Interface VOF & Marqueurs** | Reconstruction géométrique | `PLIC`, `ELVIRA`, `LVIRA`, `SLIC`, `NONE` | Reconstruction affine $n_x x + n_y y + \alpha = 0$ sur pochoir $3\times 3$ |
| | Advection VOF | `WEYMOUTH_YUE`, `CIAM`, `VOF_TVD` | Splitting alterné de Strang ; correction exacte de dilatation $C_{d,i,j}$ (Weymouth & Yue 2010) |
| | Limiteurs TVD/NVD | `SUPERBEE`, `VAN_LEER`, `MUSCL`, `HRIC`, `NONE` | Actifs uniquement avec `advection_method = "VOF_TVD"` |
| | Lissage rhéologique SVOF | `enable_svof_smoothing`, `svof_filter_width` | Filtre de Helmholtz $(I - D\Delta) C_s = C$ ($D = k \Delta x \Delta y$) pour $\rho(C_s), \mu(C_s), \mu_{\text{node}}(C_s)$ |
| | Lissage capillaire diffusif | `enable_diffusion_smoothing`, `diffusion_smoothing_steps` | Diffusion parabolique explicite $\widetilde{C}$ pour la courbure CSF |
| | Marqueurs CVI non connectés | `vof.use_markers`, `CVI_RK2`, `BILINEAR_RK2` | Nuage lagrangien non connecté ; CVI Variante B ($\nabla \cdot \mathbf{u}_{\text{CVI}} \equiv D_{i,j}$), RK2-Heun, migration MPI 2 passes, courbure parabolique $(\tau, \eta)$ |
| | Tension superficielle | `enable_surface_tension`, `sigma` | Force volumique CSF équilibrée (*Balanced-Force*) $\mathbf{f}_\sigma = \sigma \kappa \nabla C$ |
| **Couplage Vitesse-Pression** | Lagrangien Augmenté Standard | `AUGMENTED_LAGRANGIAN` | Pénalité scalaire $r_{\text{al}} \nabla(\nabla \cdot \mathbf{u})$ et boucle d'Uzawa |
| | Lagrangien Adaptatif | `AAL` / `ADAPTATIVE_AUGMENTED_LAGRANGIAN` | Champ local $r(t, M)$ selon les 7 échelles physiques (Vincent et al. 2004, Éq. 14) |
| | Lagrangien Compressible | `C-AAL` | Pénalité acoustique régularisée $r_{\text{eff}} = \frac{r}{1 + r \chi_T / \Delta t}$ (El Ouafa et al. 2024) |
| | Monolithique | `MONOLITHIC` / `FULLY_COUPLED` / `FC` | Système point-selle couplé vitesse-pression sans pénalité |
| | Projection Standard | `STANDARD_PROJECTION` / `PROJECTION` | Schéma fractionnaire de Chorin-Temam (prédiction + Poisson en pression + projection) |
| **Convection & Temps** | Schémas convectifs | `UPWIND`, `CENTERED`, `TVD_SUPERBEE` | Discrétisation spatiale des flux de quantité de mouvement |
| | Linéarisation convective | `PICARD`, `ADAMS_BASHFORTH` (`AB2`), `NEWTON` | Traitement du terme $(\mathbf{u}\cdot\nabla)\mathbf{u}$ dans le système implicite |
| | Schémas temporels | `EULER`, `GEAR` / `BDF2` + `enable_ssp_rk3_momentum` | Euler implicite (ordre 1), BDF2/Gear (ordre 2, $\alpha_{\text{time}}=1.5$), transport inertiel SSP-RK3 |
| **Solveurs & Préconditionneurs** | Solveurs linéaires | `BICGSTAB2`, `BICGSTAB`, `CG`, `GMRES`, `MINRES`, `DIRECT` | Implémentations Krylov *matrix-free* + LU directe LAPACK `dgetrf`/`dgetrs` (passerelles `MUMPS`/`PETSC`/`HYPRE` en C++) |
| | Préconditionneurs | `JACOBI`, `SSOR`, `ILU0`, `BLOCK_TRIANGULAR` (`BT_SBP`), `CUSTOM`, `NONE` | Préconditionnement diagonal enrichi, SSOR, ILU(0) et bloc triangulaire à complément de Schur |

---

## Architecture Fonctionnelle Vue par l'Utilisateur

L'exploitation de COVOF 2D suit un flux de travail unifié, identique pour les exécutables C++17 et Fortran 2018 :

```text
+-------------------------------------------------------------------------+
| 1. PRÉ-TRAITEMENT ET CONFIGURATION                                      |
| - Définition du domaine [0, Lx] x [0, Ly] et de la grille (nx, ny)      |
| - Choix du régime (incompressible, compressible All-Mach, Boussinesq)   |
| - Propriétés fluid1 / fluid2 (rho, mu, Cp, k_thermal, eos, chi_T...)    |
| - Paramétrage VOF (WEYMOUTH_YUE + PLIC/ELVIRA, SVOF, Marqueurs CVI-RK2) |
| - Conditions initiales (initial_condition), CL (boundary_conditions)    |
| - Solveur (MONOLITHIC / AAL / C-AAL / AL, BiCGSTAB2, BT_SBP / JACOBI)   |
+-------------------------------------------------------------------------+
                                     |
                                     v
+-------------------------------------------------------------------------+
| 2. EXÉCUTION DU CALCUL (Séquentiel ou MPI pur)                          |
| - C++17       : mpirun -np P ./COVOF_2D-C++/bin/covof2d_cpp cas.json    |
| - Fortran 2018: mpirun -np P ./COVOF_2D-FORTRAN/bin/covof2d_fortran ... |
| - Partitionnement cartésien MPI 2D, allocation des halos (0..nx+1)      |
| - Boucle temporelle : [SSP-RK3] -> VOF+CVI-Markers -> AL/NS -> Énergie  |
+-------------------------------------------------------------------------+
                                     |
                                     v
+-------------------------------------------------------------------------+
| 3. CONTRÔLE EN TEMPS RÉEL (Monitoring)                                  |
| - Bannière synthétique Config::print() au démarrage (rang 0)            |
| - Tableau console périodique (log_interval) : Pas, Temps, dt, Iter,     |
|   Vitesse Max, Div Max, [Nu Moyen]                                      |
| - Série temporelle CSV : <output_dir>/<case_name>_history.csv           |
+-------------------------------------------------------------------------+
                                     |
                                     v
+-------------------------------------------------------------------------+
| 4. POST-TRAITEMENT ET VISUALISATION PARAVIEW / PYTHON                   |
| - Champs volumiques : <case_name>_<step:06d>.pvtr + *_p<rank>_*.vtr     |
| - Nuage de marqueurs : <case_name>_markers_<step:06d>.vtk (POLYDATA)    |
| - Diagnostics analytiques (modes dédiés) : laplacian_error.csv,         |
|   poiseuille_error.csv, vortex_error.csv, profile_*.csv                 |
+-------------------------------------------------------------------------+
```

### Organisation Réelle du Dépôt `PROJET_COVOF-2D`

À la racine `/Users/elhadjiibrahimathiam/Documents/PROJET_COVOF-2D/`, l'arborescence de travail est organisée comme suit :

- `COVOF_2D-C++/` : Code source C++17 (`include/covof/` : 15 fichiers `.hpp` ; `src/` : 13 fichiers `.cpp`), `CMakeLists.txt`, et binaires compilés dans `bin/` (`covof2d_cpp` et lien `covof2d`).
- `COVOF_2D-FORTRAN/` : Code source Fortran 2018 (`src/` : 16 fichiers `.f90`), `CMakeLists.txt`, et binaire compilé dans `bin/covof2d_fortran`.
- `VALIDATION_TOTALE/` : Dossier maître de vérification et validation regroupant :
  - `Run_Valide-CASES/` : 14 cas de validation analytique et canonique pilotés par `run_all_validations.sh` (`validation_cas_1_laplacian` à `validation_cas_14_instabilite_interface_gravitaire`).
  - `Run_CASES/` : Campagnes avancées (`validation_cas_9_50x100_newton_dual_code`, `validation_cas_10_liquid_sheet`, `validation_cas_14_marche_descendante` à `validation_cas_19_convection_poreux_darcy_brinkman`, `validation_cas_21_ascension_bulle_2d`, `validation_cas_24_viscous_jet_buckling`) ainsi que le sous-dossier `Run_FC-WY-CVI-Marker/` dédié aux 8 campagnes hybrides CVI-Marqueurs et au script `verify_cvi_markers.py`.
  - `Theoritical_Description/` : Descriptions théoriques des 24 cas de validation et rapport consolidé `Validation_Description.pdf`.
- `Guides_Tutos/` : Suite documentaire complète (Volumes I, II et III, script `verify_cvi_markers.py` dans `Run_FC-WY-CVI-Marker/`, synthèses théoriques `SYNTHESE_METHODE_FC_WY_CVI_MARQUEURS.pdf` et `FEUILLE_DE_ROUTE_IMPLEMENTATION_COVOF_CVI_MARQUEURS.pdf`).
- `output/` : Répertoire de sortie par défaut généré lors des simulations si `output.output_dir` vaut `"output"`.
# Installation, Compilation et Déploiement HPC

## Prérequis Logiciels et Matrice de Compatibilité

Le déploiement de **COVOF 2D** sur stations de travail, serveurs de calcul ou supercalculateurs massivement parallèles requiert uniquement des outils standards de l'écosystème scientifique. Les deux branches jumelles (`COVOF_2D-C++` et `COVOF_2D-FORTRAN`) ont été conçues sans dépendance JSON externe (analyseurs JSON natifs intégrés dans `Config.cpp` et `module_json.f90`) et reposent exclusivement sur **MPI pur** (*Message Passing Interface*) et une implémentation LAPACK/BLAS pour le solveur direct LU optionnel.

### Matrice de Compatibilité des Compilateurs

Les deux versions adhèrent strictement aux standards internationaux **ISO C++17** et **ISO Fortran 2018** :

| Famille de Compilateurs | Version Minimale | Compilateur C++ | Compilateur Fortran | Drapeaux CMake / Compilation |
| :--- | :--- | :--- | :--- | :--- |
| **GNU Compiler Collection (GCC)** | $\ge 9.3$ (recommandé $\ge 11.2$) | `g++` | `gfortran` | C++ : `-O3 -Wall -Wextra -Wno-unused-parameter` <br> Fortran : `-O3 -march=native` |
| **LLVM / Clang + Flang/GFortran** | $\ge 14.0$ | `clang++` | `gfortran` / `flang` | `-O3 -march=native` |
| **Intel oneAPI** | $\ge 2022.1$ | `icpx` | `ifx` | `-O3 -xHost` |
| **Apple Clang + Homebrew GCC (macOS)** | Xcode CLT $\ge 14$ | `clang++` | `gfortran` (Homebrew) | Liaison automatique `-framework Accelerate` |

### Bibliothèques Requises et Modèle de Parallélisme

1. **MPI pur (OpenMPI / MPICH / Intel MPI) :**
   - COVOF 2D utilise un **parallélisme distribué MPI pur** (`MPI::MPI_CXX` en C++, `MPI::MPI_Fortran` en Fortran), basé sur une topologie cartésienne 2D (`MPI_Cart_create`).
   - **Aucun thread OpenMP n'est utilisé par le solveur** : les fichiers `CMakeLists.txt` n'activent pas `-fopenmp`, et il n'existe aucune région parallèle `#pragma omp` ou `!$omp parallel` dans les noyaux de calcul. Chaque processus MPI gère un sous-domaine mémoire indépendant.
   - Implémentations validées : **OpenMPI** ($\ge 4.0$, référence principale du projet), **MPICH** ($\ge 3.4$), **Intel MPI** ($\ge 2021$).
2. **Algèbre Linéaire Dense (LAPACK) :**
   - Requise uniquement lorsque le solveur direct exact (`solver.linear_solver = "DIRECT"` ou `"MUMPS"` en repli autonome) invoque la factorisation LU dense `dgetrf_` / `dgetrs_`.
   - Sous **macOS (`APPLE`)** : `CMakeLists.txt` lie automatiquement le framework natif `-framework Accelerate`.
   - Sous **Linux / Clusters HPC** : lier `OpenBLAS` / `LAPACK` (`liblapack-dev`, `libopenblas-dev`) ou Intel MKL selon l'environnement système.
3. **CMake et Python (post-traitement / vérification) :**
   - **CMake $\ge 3.16$** pour `COVOF_2D-C++` et **$\ge 3.18$** pour `COVOF_2D-FORTRAN`.
   - **Python 3** avec `numpy`, `matplotlib`, `pyvista` (optionnel, pour le post-traitement et `verify_cvi_markers.py`).

---

## Compilation Pas à Pas avec CMake

Les deux dépôts disposent chacun de leur fichier `CMakeLists.txt` autonome à leur racine respective.

### Compilation de la Branche C++17 (`COVOF_2D-C++`)

Le fichier `COVOF_2D-C++/CMakeLists.txt` impose `CMAKE_CXX_STANDARD 17` (`REQUIRED ON`), active `CMAKE_EXPORT_COMPILE_COMMANDS ON`, recherche `MPI REQUIRED` et compile les 13 unités `.cpp` vers l'exécutable `COVOF_2D-C++/bin/covof2d_cpp` tout en créant automatiquement un lien symbolique de commodité `COVOF_2D-C++/bin/covof2d` :

```bash
# 1. Se positionner dans le répertoire C++
cd /Users/elhadjiibrahimathiam/Documents/PROJET_COVOF-2D/COVOF_2D-C++

# 2. Créer et entrer dans le dossier de construction hors-source
mkdir -p build && cd build

# 3. Configurer avec CMake en mode Release
cmake -DCMAKE_BUILD_TYPE=Release ..

# 4. Compiler en parallèle sur les cœurs disponibles
make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)
```

À l'issue de l'édition de liens, vérifier la présence des binaires :

```bash
ls -la ../bin/
# Doit afficher :
#   covof2d_cpp
#   covof2d -> .../bin/covof2d_cpp
```

#### Options CMake et Drapeaux Utiles pour C++

| Option CMake | Valeur / Comportement | Rôle |
| :--- | :--- | :--- |
| `CMAKE_BUILD_TYPE` | `Release`, `Debug`, `RelWithDebInfo` | `Release` active `-O3` ; `Debug` permet l'ajout de `-g -O0 -fsanitize=address,undefined` |
| `CMAKE_CXX_STANDARD` | `17` (`REQUIRED ON`) | Standard ISO C++17 requis pour `std::filesystem`, `std::clamp`, etc. |
| `CMAKE_CXX_COMPILER` | ex. `mpicxx` ou `g++` | Surcharge explicite du compilateur C++ si plusieurs suites coexistent |
| `CMAKE_EXPORT_COMPILE_COMMANDS` | `ON` | Génère `build/compile_commands.json` pour les analyseurs LSP (`clangd`) |

### Compilation de la Branche Fortran 2018 (`COVOF_2D-FORTRAN`)

Le fichier `COVOF_2D-FORTRAN/CMakeLists.txt` impose `CMAKE_Fortran_STANDARD 2018` (`REQUIRED ON`), ajoute les drapeaux `-O3 -march=native`, recherche `MPI REQUIRED` et compile les 16 fichiers `.f90` dans l'ordre exact des dépendances modulaires vers l'exécutable `COVOF_2D-FORTRAN/bin/covof2d_fortran` :

```bash
# 1. Se positionner dans le répertoire Fortran
cd /Users/elhadjiibrahimathiam/Documents/PROJET_COVOF-2D/COVOF_2D-FORTRAN

# 2. Créer et entrer dans le dossier de construction
mkdir -p build && cd build

# 3. Configurer avec CMake en mode Release
cmake -DCMAKE_BUILD_TYPE=Release ..

# 4. Compiler en parallèle
make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)
```

:::warning
Sur cluster hétérogène (nœuds de compilation frontaux d'architecture différente des nœuds de calcul), le drapeau `-march=native` de `COVOF_2D-FORTRAN/CMakeLists.txt` doit être remplacé par l'architecture cible des nœuds de calcul (par exemple `-march=skylake-avx512` ou `-xCORE-AVX512`) afin d'éviter une erreur `Illegal Instruction (SIGILL)` à l'exécution.
:::

### Drapeaux d'Optimisation SIMD et Résolution de LAPACK sous Linux

Sur architecture Linux x86_64 pour nœuds Xeon ou EPYC, si `dgetrf_`/`dgetrs_` nécessite une liaison explicite avec OpenBLAS ou LAPACK :

```bash
# Exemple C++ et Fortran avec drapeaux AVX-512 ciblés et LAPACK explicite :
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_FLAGS="-O3 -march=skylake-avx512" \
      -DCMAKE_EXE_LINKER_FLAGS="-llapack -lblas" ..
```

---

## Déploiement sur Supercalculateurs (Clusters HPC)

Dans un environnement distribué piloté par un ordonnanceur (*workload manager*), **COVOF 2D étant un code MPI pur**, chaque tâche Slurm ou PBS correspond à **1 processus MPI sur 1 cœur physique** (`--cpus-per-task=1`). L'utilisation de ressources multi-cœurs par rang n'apporte aucun gain, car le code ne crée pas de threads applicatifs OpenMP (seule l'éventuelle bibliothèque BLAS liée lors de l'appel à `DIRECT` doit être contrainte à un thread si elle est multi-threadée, par exemple `OPENBLAS_NUM_THREADS=1` ou `VECLIB_MAXIMUM_THREADS=1` sous macOS).

### Script de Soumission pour Gestionnaire Slurm

Voici un script de soumission complet pour **Slurm** (`run_covof_slurm.sh`), allouant 2 nœuds de 64 cœurs (soit 128 processus MPI purs) :

```bash
#!/bin/bash
#SBATCH --job-name=COVOF2D_PROD
#SBATCH --output=covof2d_%j.log
#SBATCH --error=covof2d_%j.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --partition=standard
#SBATCH --exclusive

# Chargement de l'environnement logiciel
module purge
module load gcc/12.2.0
module load openmpi/4.1.5
module load openblas/0.3.21

echo "=========================================================="
echo "Démarrage COVOF 2D (MPI pur) sur Slurm"
echo "Job ID                     : $SLURM_JOB_ID"
echo "Nombre de nœuds            : $SLURM_NNODES"
echo "Nombre total de rangs MPI  : $SLURM_NTASKS"
echo "Liste des nœuds            : $SLURM_JOB_NODELIST"
echo "=========================================================="

# Neutralisation d'un éventuel multi-threading implicite des bibliothèques BLAS
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Transport réseau haute vitesse UCX (InfiniBand)
export UCX_TLS=rc,mm,sm
export OMPI_MCA_pml=ucx

WORKDIR=$SLURM_SUBMIT_DIR
CASE_CONFIG="$WORKDIR/case_config.json"
cd "$WORKDIR"

# Lancement parallèle MPI pur (C++17 ou Fortran 2018)
srun --mpi=pmix --cpu-bind=cores ./COVOF_2D-C++/bin/covof2d_cpp "$CASE_CONFIG"

echo "=========================================================="
echo "Fin du calcul COVOF 2D - Code retour : $?"
echo "=========================================================="
```

### Script de Soumission pour Gestionnaire PBS / Torque

Pour les infrastructures s'appuyant sur **PBS Professional** ou **Torque** :

```bash
#!/bin/bash
#PBS -N COVOF2D_PBS
#PBS -l select=2:ncpus=64:mpiprocs=64
#PBS -l walltime=08:00:00
#PBS -q standard
#PBS -j oe
#PBS -o covof2d_pbs.log

cd "$PBS_O_WORKDIR"

module purge
module load intel-oneapi-compilers/2023.1.0
module load intel-oneapi-mpi/2021.9.0
module load intel-oneapi-mkl/2023.1.0

export MKL_NUM_THREADS=1
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=core
export I_MPI_FABRICS=shm:ofi

mpirun -np 128 ./COVOF_2D-FORTRAN/bin/covof2d_fortran case_config.json
```

### Gestion de l'Affinité Mémoire et Équilibre de Charge Cartésien

Lors du partitionnement par `MPI_Cart_create` / `MPI_Dims_create`, il est crucial de choisir un nombre de processus $P = P_x \times P_y$ compatible avec la résolution $(N_x, N_y)$ du maillage :

1. **Verrouillage cœur par cœur (*CPU pinning*) :**
   ```bash
   # OpenMPI : 1 rang MPI épinglé sur 1 cœur physique, pas de migration cache
   mpirun -np 8 --bind-to core --map-by core ./COVOF_2D-C++/bin/covof2d_cpp input.json
   ```
2. **Contraintes sur $P_x$ et $P_y$ :**
   - Chaque sous-domaine local doit posséder au moins $N_{x,\text{local}} \ge 2$ et $N_{y,\text{local}} \ge 2$ mailles (idéalement $\ge 8$ à $16$ mailles par côté pour limiter le ratio communication/calcul des halos et des marqueurs CVI).
   - Éviter un nombre premier $P$ élevé (ex. $P = 17$) sur un domaine carré : `MPI_Dims_create` produirait des tranches 1D très allongées $(17 \times 1)$ augmentant la surface d'échange et le nombre de marqueurs migrant entre rangs. Privilégier des facteurs équilibrés ($P \in \{2, 4, 6, 8, 12, 16, 32, 64\}$).

---

## Conteneurisation et Vérification Post-Installation

### Définition Singularity / Apptainer (`Apptainer.def`)

Sur supercalculateurs où la compilation native est délicate, une image **Apptainer/Singularity** (*rootless*) garantit la portabilité complète :

```dockerfile
Bootstrap: docker
From: ubuntu:22.04

%post
    apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        cmake \
        libopenmpi-dev \
        openmpi-bin \
        libopenblas-dev \
        liblapack-dev \
        python3 \
        python3-numpy \
        ca-certificates
    apt-get clean && rm -rf /var/lib/apt/lists/*

%environment
    export LC_ALL=C
    export OPENBLAS_NUM_THREADS=1

%runscript
    exec "$@"
```

Construction et exécution :

```bash
apptainer build covof2d.sif Apptainer.def
mpirun -np 4 apptainer exec covof2d.sif ./COVOF_2D-C++/bin/covof2d_cpp cas.json
```

### Protocole de Vérification Immédiate Post-Installation

Après compilation des deux branches dans `/Users/elhadjiibrahimathiam/Documents/PROJET_COVOF-2D`, exécuter les vérifications suivantes depuis la racine du projet :

1. **Test d'aide en ligne de commande (`--help`) :**
   ```bash
   ./COVOF_2D-C++/bin/covof2d_cpp --help
   ./COVOF_2D-FORTRAN/bin/covof2d_fortran --help
   ```

2. **Test de vérification analytique rapide (Laplacien 2D en séquentiel et MPI) :**
   ```bash
   ./COVOF_2D-C++/bin/covof2d_cpp VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_1_laplacian/cas_1_laplacian_32x32.json
   mpirun -np 4 ./COVOF_2D-FORTRAN/bin/covof2d_fortran VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_1_laplacian/cas_1_laplacian_32x32.json
   ```
   *(Adapter le nom du fichier `.json` présent dans le dossier du cas de validation selon la grille souhaitée).*

3. **Suite de vérification automatique FC-WY-CVI-Marker (C++ & Fortran, `np=1` et `np=4`) :**
   ```bash
   python3 Guides_Tutos/Run_FC-WY-CVI-Marker/verify_cvi_markers.py
   ```
   Ce script exécute automatiquement les 3 paliers de vérification mathématique et physique (divergence analytique nulle de l'interpolation CVI Variante B, équilibre capillaire de Young-Laplace à $32\times 32$ et $64\times 64$, et parité stricte `np=1` vs `np=4` sur le tourbillon de Rider-Kothe et la goutte statique) pour les deux binaires `covof2d_cpp` et `covof2d_fortran`.

Lorsque la simulation Navier-Stokes se termine normalement, la console affiche :
`>>> SIMULATION COVOF 2D TERMINEE AVEC SUCCES A t = ... s <<<`
suivi du chemin du répertoire de sortie contenant les fichiers `.pvtr`, `.vtr` et `_history.csv`.
# Formulation Mathématique et Modèles Physiques

## Système d'Équations Gouvernantes

Le solveur **COVOF 2D** résout les équations fondamentales de la mécanique des milieux continus régissant les écoulements 2D de fluides visqueux, newtoniens, monophasiques ou diphasiques, incompressibles ou faiblement compressibles (*All-Mach*), avec ou sans transfert thermique et en présence éventuelle d'obstacles solides ou poreux immergés.

### Équations de Navier-Stokes sous Forme Continue (*One-Fluid*)

Dans un domaine rectangulaire ouvert et borné $\Omega = [0, L_x] \times [0, L_y] \subset \mathbb{R}^2$ de frontière $\partial\Omega$ et pour $t \in [0, t_{\text{final}}]$, la formulation *One-Fluid* regroupe les deux phases en un milieu continu unique de masse volumique $\rho(\mathbf{x}, t)$ et de viscosité dynamique $\mu(\mathbf{x}, t)$ :

$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0$$

$$\frac{\partial (\rho \mathbf{u})}{\partial t} + \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u}) = -\nabla p + \nabla \cdot \boldsymbol{\tau} + \rho \mathbf{g} + \mathbf{f}_\sigma + \mathbf{f}_{\text{brinkman}} + \mathbf{f}_{\text{boussinesq}}$$

où :
- $\rho = \rho(\mathbf{x}, t)$ est la masse volumique locale du mélange $[\text{kg}\cdot\text{m}^{-3}]$.
- $\mathbf{u} = (u, v)^T$ est le vecteur vitesse eulérien $[\text{m}\cdot\text{s}^{-1}]$.
- $p = p(\mathbf{x}, t)$ est la pression hydrodynamique ou thermodynamique $[\text{Pa}]$.
- $\boldsymbol{\tau} = \mu \left( \nabla \mathbf{u} + (\nabla \mathbf{u})^T \right)$ est le tenseur des contraintes visqueuses newtoniennes en écoulement incompressible (avec contribution de dilatation $\lambda (\nabla \cdot \mathbf{u}) \mathbf{I}$ négligée en régime bas-Mach isotherme ou prise en compte dans les flux d'élongation).
- $\mathbf{g} = (g_x, g_y)^T$ est le vecteur accélération de la pesanteur $[\text{m}\cdot\text{s}^{-2}]$ (`physics.gx`, `physics.gy`).
- $\mathbf{f}_\sigma$ est la force volumique de tension superficielle interfaciale $[\text{N}\cdot\text{m}^{-3}]$.
- $\mathbf{f}_{\text{brinkman}}$ est la force de traînée de pénalisation de Darcy-Brinkman dans les obstacles ou milieux poreux $[\text{N}\cdot\text{m}^{-3}]$.
- $\mathbf{f}_{\text{boussinesq}}$ est la force de flottabilité thermique sous approximation de Boussinesq $[\text{N}\cdot\text{m}^{-3}]$.

### Formulation Advective Utilisée sur la Grille Décalée MAC

En exploitant la conservation de la masse, la forme advective non conservative résolue par le noyau implicite de COVOF 2D s'écrit, projetée selon $x$ et $y$ :

$$\rho \left( \frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} \right) = -\frac{\partial p}{\partial x} + \frac{\partial}{\partial x}\left( 2\mu \frac{\partial u}{\partial x} \right) + \frac{\partial}{\partial y}\left[ \mu \left( \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \right) \right] + \rho g_x + f_{\sigma, x} + f_{\text{brinkman}, x}$$

$$\rho \left( \frac{\partial v}{\partial t} + u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} \right) = -\frac{\partial p}{\partial y} + \frac{\partial}{\partial x}\left[ \mu \left( \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \right) \right] + \frac{\partial}{\partial y}\left( 2\mu \frac{\partial v}{\partial y} \right) + \rho g_y + f_{\sigma, y} + f_{\text{brinkman}, y} + f_{\text{boussinesq}, y}$$

Lorsque l'option `physics.enable_ssp_rk3_momentum = true` est activée (automatique par défaut en régime `compressible = true`), une étape préliminaire explicite **SSP-RK3** assure le transport inertiel conservatif de $(\rho u, \rho v)$ avant le couplage implicite vitesse-pression.

---

## Modèle Diphasique VOF, Lissage SVOF et Marqueurs Hybrides

### Fonction Indicatrice de Phase (Fraction Volumique $C$)

Dans la formulation *One-Fluid*, les deux fluides immiscibles — **Fluide 1** (`physics.fluid1`, généralement le fluide lourd/porteur ou ambiant) et **Fluide 2** (`physics.fluid2`, généralement la goutte, bulle, cylindre dense ou gaz selon le cas initial) — sont repérés par la fonction indicatrice de couleur $C(\mathbf{x}, t) \in [0, 1]$ (`fields.C` en C++, `fields%C` en Fortran) :

$$C(\mathbf{x}, t) = \begin{cases} 1 & \text{dans la Phase 1 (`fluid1`)} \\ 0 & \text{dans la Phase 2 (`fluid2`)} \\ ]0, 1[ & \text{dans les mailles mixtes coupées par l'interface } \Gamma(t) \end{cases}$$

Sans changement de phase, le transport lagrangien de $C$ obéit à :

$$\frac{\partial C}{\partial t} + \mathbf{u} \cdot \nabla C = 0 \quad \Longleftrightarrow \quad \frac{\partial C}{\partial t} + \nabla \cdot (C \mathbf{u}) = C (\nabla \cdot \mathbf{u})$$

### Régularisation Rhéologique SVOF par Filtre de Helmholtz

Conformément à El Ouafa, Vincent & Le Chenadec (*Fluids* 2021), lorsque le contraste de masse volumique ($\rho_1 / \rho_2 \gg 1$) ou de viscosité ($\mu_1 / \mu_2 \gg 1$) est élevé, COVOF 2D distingue rigoureusement :

1. **La fraction volumique advectée $C$** (`fields.C` / `fields%C`), qui demeure **parfaitement raide** (épaisseur d'une maille) afin de préserver la netteté géométrique de l'interface.
2. **La fonction couleur lissée rhéologique $C_s$** (`fields.C_s` / `fields%C_s`), activée par `physics.enable_svof_smoothing = true` et calculée de manière purement diagnostique après l'advection en résolvant l'équation elliptique de Helmholtz avec conditions aux limites de Neumann homogènes ($\partial_n C_s = 0$ sur $\partial\Omega$) :
   $$C_s - D \nabla^2 C_s = C, \qquad D = k \, \Delta x \Delta y$$
   où $k$ est le paramètre `physics.svof_filter_width` (défaut `1.0`). Ce système diagonalement dominant est résolu par 25 itérations de Gauss-Seidel symétrique parallèle avec bornage strict $C_s \in [0, 1]$.
3. **La fonction couleur filtrée capillaire $\widetilde{C}$** (`fields.C_tilde` / `fields%C_tilde`), activée par `physics.enable_diffusion_smoothing = true` (`diffusion_smoothing_steps`, défaut `3`), obtenue par $N_{\text{diff}}$ pas explicites de diffusion laplacienne ($\Delta \tau = 0.2 \min(\Delta x^2, \Delta y^2)$) pour régulariser le calcul classique de la courbure CSF lorsque les marqueurs ne sont pas activés.

### Lois de Mélange des Propriétés Matérielles sur Grille MAC

En notant $C_{\text{src}} = C_s$ si `enable_svof_smoothing = true` (sinon $C_{\text{src}} = C$) :

1. **Centres des cellules de pression $(i, j)$** — moyenne arithmétique :
   $$\rho_{i,j} = C_{\text{src}, i,j} \rho_1 + (1 - C_{\text{src}, i,j}) \rho_2, \qquad \mu_{i,j} = C_{\text{src}, i,j} \mu_1 + (1 - C_{\text{src}, i,j}) \mu_2$$
2. **Faces décalées de vitesse $(i - 1/2, j)$ et $(i, j - 1/2)$** — moyenne arithmétique des cellules adjacentes :
   $$\rho_{u, i, j} = \frac{1}{2}\left(\rho_{i-1, j} + \rho_{i, j}\right), \qquad \rho_{v, i, j} = \frac{1}{2}\left(\rho_{i, j-1} + \rho_{i, j}\right)$$
3. **Sommets / Nœuds de cisaillement $(i - 1/2, j - 1/2)$** — **moyenne harmonique** (`fields.mu_node` / `fields%mu_node`), indispensable pour assurer la continuité de la contrainte tangentielle à travers les interfaces à grand saut de viscosité :
   - En mode VOF standard (`enable_svof_smoothing = false`) :
     $$\mu_{\text{node}, i, j} = \frac{4}{\frac{1}{\mu_{i-1, j-1}} + \frac{1}{\mu_{i, j-1}} + \frac{1}{\mu_{i-1, j}} + \frac{1}{\mu_{i, j}}}$$
   - En mode SVOF (`enable_svof_smoothing = true`, Éqs. 8–9 *Fluids* 2021) :
     $$\mu_{\text{node}, i, j} = \frac{\mu_1 \mu_2}{\mu_2 \, C_{s,\text{node}} + \mu_1 (1 - C_{s,\text{node}})}, \qquad C_{s,\text{node}} = \frac{1}{4}\left( C_{s,i-1,j-1} + C_{s,i,j-1} + C_{s,i-1,j} + C_{s,i,j} \right)$$

En régime **monophasique** (`enable_two_phase = false`), COVOF 2D propose en outre les modèles de viscosité spatiale `physics.viscosity_model` :
- `"CONSTANT"` : $\mu = \mu_1$.
- `"EXPONENTIAL"` : $\mu(y) = \mu_0 \exp\left( a \left( \frac{y - e}{e} \right)^2 \right)$ (avec `visc_mu0`, `visc_a`, `visc_e`).
- `"ASYMMETRIC_EXPONENTIAL"` (en C++) : $\mu(y) = \mu_0 \exp\left( -a \frac{y - e}{e} \right)$.

### Modèle de Tension Superficielle CSF et Couplage Marqueurs Lagrangiens CVI

La tension superficielle est modélisée par la formulation volumique *Balanced-Force CSF* (Brackbill et al., 1992) discrétisée directement aux faces décalées MAC :

$$F_{\sigma, x, i-1/2, j} = \sigma \, \bar{\kappa}_{i-1/2, j} \, \frac{C_{i, j} - C_{i-1, j}}{\Delta x}, \qquad F_{\sigma, y, i, j-1/2} = \sigma \, \bar{\kappa}_{i, j-1/2} \, \frac{C_{i, j} - C_{i, j-1}}{\Delta y}$$

où $\sigma$ est la tension superficielle (`physics.sigma` $[\text{N}\cdot\text{m}^{-1}]$) et $\kappa$ la courbure algébrique plane $[\text{m}^{-1}]$. Deux modes de calcul de $(\mathbf{n}, \kappa)$ sont disponibles :

1. **Mode VOF eulérien standard :**
   À partir du champ couleur $C_{\text{curv}}$ ($\widetilde{C}$ si `enable_diffusion_smoothing`, sinon $C_s$ si `enable_svof_smoothing`, sinon $C$) :
   $$\mathbf{n}_{i,j} = \frac{\nabla_h C_{\text{curv}, i,j}}{\|\nabla_h C_{\text{curv}, i,j}\|}, \qquad \kappa_{i,j} = -\nabla_h \cdot \mathbf{n}_{i,j}$$
2. **Mode hybride FC-WY-CVI-Marker (`vof.use_markers = true` et `vof.marker_curvature = true`) :**
   Un nuage de marqueurs lagrangiens non connectés $(x_p, y_p, n_{x,p}, n_{y,p}, \kappa_p)$ est advecté à l'ordre 2 par le schéma RK2-Heun à l'aide de l'**Interpolation Conservative de la Vitesse (CVI — Variante B exacte sur boîtes décalées MAC)** garantissant $\nabla \cdot \mathbf{u}_{\text{CVI}} \equiv D_{i,j}$ en tout point. Dans chaque cellule mixte, la normale et la courbure $\kappa_{i,j}^{\text{marker}}$ sont obtenues par régression parabolique pondérée dans le repère local tangent-normal $(\tau, \eta)$ avec filtrage angulaire $\mathbf{n}_p \cdot \mathbf{n}^{\text{ref}} \ge 0.5$ pour protéger les ligaments minces. Cette approche élimine la majeure partie des courants parasites et préserve la courbure avec une convergence d'ordre 2.

---

## Modèle Faiblement Compressible, Régime All-Mach et Équations d'État (EOS)

Pour les écoulements diphasiques compressibles ou acoustiques à bas nombre de Mach (compression d'une bulle sous saut de pression, onde de choc acoustique faible, piston isotherme), COVOF 2D implémente la formulation **All-Mach C-AAL** (El Ouafa et al., *Fluids* 2024).

### Équation de Continuité Barotrope et Compressibilité Isotherme

En introduisant la compressibilité isotherme $\chi_T = \frac{1}{\rho}\left(\frac{\partial \rho}{\partial p}\right)_T = \frac{1}{\rho c_s^2}$ (`fields.chi_T` $[\text{Pa}^{-1}]$), la contrainte de continuité discrète s'écrit :

$$\chi_{T, i,j} \frac{p_{i,j}^{n+1} - p_{i,j}^n}{\Delta t} + (\nabla_h \cdot \mathbf{u}^{n+1})_{i,j} = 0$$

La masse volumique cellulaire est actualisée en fonction de la pression et de la température via l'équation d'état propre à chaque phase (`physics.fluid1.eos`, `physics.fluid2.eos`), puis moyennée par $C_{\text{src}}$ : $\rho = C_{\text{src}} \rho_1(p, T) + (1 - C_{\text{src}}) \rho_2(p, T)$ et $\chi_T = C_{\text{src}} \chi_{T,1} + (1 - C_{\text{src}}) \chi_{T,2}$.

### Catalogue des Quatre Lois d'État Implémentées (`EOS.hpp` / `module_eos.f90`)

| Mot-Clé `eos` (et Alias) | Loi de Masse Volumique $\rho(p, T)$ | Compressibilité Isotherme $\chi_T$ | Célérité du Son $c_s$ | Paramètres JSON Associés |
| :--- | :--- | :--- | :--- | :--- |
| `"INCOMPRESSIBLE"` | $\rho = \rho_0$ | $\chi_T = 0$ | $10^5\,\text{m}\cdot\text{s}^{-1}$ | `rho` |
| `"LINEAR"` / `"ACOUSTIC"` | $\rho = \rho_0 \left[ 1 + \chi_{T,0}(p - p_{\text{ref}}) \right]$ | $\chi_T = \chi_{T,0}$ | $1 / \sqrt{\rho \chi_{T,0}}$ | `rho`, `chi_T` ($p_{\text{ref}} = 1.01325\times 10^5\,\text{Pa}$) |
| `"IDEAL_GAS"` / `"PERFECT_GAS"` | $\rho = \frac{p}{r_{\text{gas}} T}$ | $\chi_T = \frac{1}{\max(1, p)}$ | $\sqrt{\gamma \, r_{\text{gas}} T}$ | `r_gas`, `gamma`, `cv` |
| `"NASG"` / `"STIFFENED_GAS"` | $\rho = \frac{p + \gamma p_\infty}{(\gamma - 1) c_v T + b_{\text{nasg}}(p + \gamma p_\infty)}$ | $\chi_T = \frac{1 - \rho b_{\text{nasg}}}{\rho \max(1, p + p_\infty)}$ | $1500\,\text{m}\cdot\text{s}^{-1}$ (repli liquide) | `gamma`, `cv`, `p_inf`, `b_nasg`, `q_nasg` |

En régime compressible non isotherme (`compressible = true` et `enable_energy = true`), la routine `solve_total_energy` résout l'équation de l'énergie totale spécifique $e = u_{\text{int}} + \frac{1}{2}\|\mathbf{u}\|^2$ (`fields.energy`, `fields.rho_e`), puis reconstruit la température via $T = \frac{\gamma - 1}{r_{\text{gas}}}\left(e - \frac{1}{2}\|\mathbf{u}\|^2\right)$.

---

## Modèle de Milieu Poreux et Pénalisation de Darcy-Brinkman

Pour modéliser des obstacles solides ou des matrices poreuses sans déformer la grille cartésienne MAC, COVOF 2D implémente la pénalisation volumique de Darcy-Brinkman (`obstacles` dans le fichier JSON) :

$$\mathbf{f}_{\text{brinkman}} = -\frac{\mu}{\mathcal{K}(\mathbf{x})} \mathbf{u}$$

Deux comportements distincts sont automatiquement sélectionnés selon la valeur de `permeability` ($\mathcal{K}$) de chaque obstacle (`BOX` ou `CYLINDER`) :
- **Solide rigide imperméable ($\mathcal{K} \le 10^{-12}\,\text{m}^2$, défaut $10^{-14}$ en C++ / $10^{-40}$ en Fortran) :** Les faces de vitesse $u$ et $v$ situées dans l'obstacle reçoivent directement la condition de Dirichlet nulle ($u = 0, v = 0$) dans l'opérateur linéaire et le second membre.
- **Milieu poreux perméable ($\mathcal{K} > 10^{-12}\,\text{m}^2$) :** Le terme de traînée implicite $+\frac{\mu}{\mathcal{K}} \Delta x \Delta y \, u$ est ajouté à la diagonale de la matrice de quantité de mouvement (et dans l'échelle $\tau_{\text{brinkman}} = \mu \frac{\Delta x \Delta y}{\mathcal{K}}$ du mode AAL).

---

## Transfert Thermique et Approximation de Boussinesq

En régime incompressible avec `physics.enable_energy = true` (et `compressible = false`), COVOF 2D résout l'équation de transport-diffusion de la température $T$ (`solve_energy` / `ns_solve_energy`) :

$$\rho C_p \left( \frac{\partial T}{\partial t} + \mathbf{u} \cdot \nabla T \right) = \nabla \cdot (k_{\text{th}} \nabla T)$$

Le couplage avec l'équation de quantité de mouvement verticale s'effectue par le terme de flottabilité de Boussinesq utilisant le nombre de Rayleigh `physics.Ra`, le nombre de Prandtl `physics.Pr`, le coefficient de dilatation `physics.beta` et la température de référence `physics.T_ref` :

$$\mathbf{f}_{\text{boussinesq}} = -\rho \, \beta \, (T - T_{\text{ref}}) \, \mathbf{g}$$

À chaque pas de temps, le solveur évalue et affiche le **nombre de Nusselt moyen sur la paroi chaude ($x = 0$)** via `compute_mean_nusselt_hot_wall`.

---

## Nombres Adimensionnels et Table Générale de Nomenclature

| Nombre Adimensionnel | Définition | Signification Physique |
| :--- | :--- | :--- |
| **Reynolds ($Re$)** | $Re = \frac{\rho U L}{\mu}$ | Rapport forces d'inertie / forces visqueuses |
| **Froude ($Fr$)** | $Fr = \frac{U}{\sqrt{g L}}$ | Rapport forces d'inertie / forces gravitationnelles |
| **Weber ($We$)** | $We = \frac{\rho U^2 L}{\sigma}$ | Rapport forces d'inertie / forces capillaires |
| **Bond / Eötvös ($Bo$)** | $Bo = \frac{\Delta \rho g L^2}{\sigma}$ | Rapport flottabilité gravitaire / tension de surface |
| **Capillaire ($Ca$)** | $Ca = \frac{\mu U}{\sigma}$ | Rapport forces visqueuses / tension de surface |
| **Ohnesorge ($Oh$)** | $Oh = \frac{\mu}{\sqrt{\rho \sigma L}}$ | Amortissement visqueux des ondes capillaires |
| **Prandtl ($Pr$)** | $Pr = \frac{\mu C_p}{k_{\text{th}}}$ | Rapport diffusivité cinématique / diffusivité thermique |
| **Rayleigh ($Ra$)** | $Ra = \frac{g \beta \Delta T L^3}{\nu \alpha}$ | Intensité motrice de la convection naturelle |
| **Darcy ($Da$)** | $Da = \frac{\mathcal{K}}{L^2}$ | Perméabilité adimensionnelle du milieu poreux |

### Table Générale de Nomenclature (Symboles $\leftrightarrow$ Variables C++ & Fortran)

| Symbole | Signification Physique | Unités SI | Identifiant C++ (`COVOF_2D-C++`) | Identifiant Fortran (`COVOF_2D-FORTRAN`) |
| :--- | :--- | :--- | :--- | :--- |
| $\mathbf{u} = (u,v)$ | Vitesse sur faces décalées MAC | $\text{m}\cdot\text{s}^{-1}$ | `fields.u`, `fields.v` | `fields%u`, `fields%v` |
| $p, p^n$ | Pression cellulaire et état précédent | $\text{Pa}$ | `fields.p`, `fields.p_old` | `fields%p`, `fields%p_old` |
| $C$ | Fraction volumique VOF raide | $-$ | `fields.C` | `fields%C` |
| $C_s$ | Fraction volumique lissée SVOF (Helmholtz) | $-$ | `fields.C_s` | `fields%C_s` |
| $\widetilde{C}$ | Fraction volumique diffusée pour CSF | $-$ | `fields.C_tilde` | `fields%C_tilde` |
| $\rho, \rho_u, \rho_v$ | Masse volumique (centres et faces) | $\text{kg}\cdot\text{m}^{-3}$ | `fields.rho`, `fields.rho_u`, `fields.rho_v` | `fields%rho`, `fields%rho_u`, `fields%rho_v` |
| $\mu, \mu_{\text{node}}$ | Viscosité dynamique (centres et nœuds) | $\text{Pa}\cdot\text{s}$ | `fields.mu`, `fields.mu_node` | `fields%mu`, `fields%mu_node` |
| $\mathbf{n}, \kappa$ | Normale unitaire et courbure d'interface | $-, \text{m}^{-1}$ | `fields.norm_x`, `fields.norm_y`, `fields.kappa` | `fields%curvature` (et normales locales) |
| $\mathbf{f}_\sigma$ | Force capillaire CSF aux faces | $\text{N}\cdot\text{m}^{-3}$ | `fields.F_sigma_x`, `fields.F_sigma_y` | `fields%fx_sigma`, `fields%fy_sigma` |
| $(x_p, y_p, \mathbf{n}_p, \kappa_p)$ | Marqueur lagrangien non connecté | $\text{m}, -, \text{m}^{-1}$ | `tracker.markers.cloud[m]` (`Marker`) | `global_marker_cloud%list(m)` (`marker_t`) |
| $\chi_T$ | Compressibilité isotherme | $\text{Pa}^{-1}$ | `fields.chi_T` | `fields%chi_T` |
| $r_{\text{al}}, r(t, M)$ | Pénalité AL constante et champ local AAL | $\text{Pa}\cdot\text{s}$ | `config.solver.r_al`, `fields.r_aal` | `cfg%r_al`, `fields%r_aal` |
| $T, C_p, k_{\text{th}}, \beta$ | Température et paramètres thermiques | $\text{K}, \text{J}/(\text{kg}\cdot\text{K}), \text{W}/(\text{m}\cdot\text{K}), \text{K}^{-1}$ | `fields.T`, `config.fluid1.Cp`, `config.fluid1.k_thermal`, `config.beta` | `fields%T`, `cfg%Cp`, `cfg%k_thermal`, `cfg%beta` |
| $e, \rho e$ | Énergie totale spécifique et volumique | $\text{J}\cdot\text{kg}^{-1}, \text{J}\cdot\text{m}^{-3}$ | `fields.energy`, `fields.rho_e` | `fields%energy`, `fields%rho_e` |
| $\mathcal{K}$ | Perméabilité de Darcy-Brinkman | $\text{m}^2$ | `config.obstacles[k].permeability` | `cfg%obs_perm` |
| $\sigma, \mathbf{g}$ | Tension de surface et gravité | $\text{N}\cdot\text{m}^{-1}, \text{m}\cdot\text{s}^{-2}$ | `config.sigma`, `config.gx`, `config.gy` | `cfg%sigma`, `cfg%gx`, `cfg%gy` |
| $\Delta x, \Delta y, \Delta t$ | Pas d'espace et pas de temps | $\text{m}, \text{s}$ | `mesh.dx`, `mesh.dy`, `dt` | `mesh%dx`, `mesh%dy`, `dt` |
# Méthodes Numériques

## Discrétisation Spatiale sur Grille Cartésienne Décalée (MAC)

Le solveur **COVOF 2D** emploie la discrétisation par volumes finis sur grille cartésienne décalée dite *Marker-and-Cell* (MAC, Harlow & Welch 1965). Ce positionnement entrelacé des inconnues garantit la condition inf-sup (Ladyzhenskaya-Babuška-Brezzi discrète) et élimine intrinsèquement tout mode parasite en damier (*checkerboard*) sur la pression sans recourir à l'interpolation de Rhie-Chow ni à une stabilisation artificielle.

### Localisation des Variables Primordiales

Sur une grille uniforme de pas $\Delta x = L_x / N_x$ et $\Delta y = L_y / N_y$, les grandeurs sont stockées aux positions géométriques suivantes (indices locaux $i \in [1, N_{x,\text{loc}}]$, $j \in [1, N_{y,\text{loc}}]$ entourés d'une couche fantôme $0$ et $N+1$) :

1. **Centres des cellules $(i, j)$ :**
   Pression $p_{i,j}$, fraction volumique VOF $C_{i,j}$, versions régularisées $C_{s,i,j}$ et $\widetilde{C}_{i,j}$, masse volumique $\rho_{i,j}$, viscosité $\mu_{i,j}$, température $T_{i,j}$, énergie totale $e_{i,j}$ et $\rho e_{i,j}$, compressibilité $\chi_{T,i,j}$, paramètre AAL $r_{\text{aal},i,j}$, courbure $\kappa_{i,j}$ et normale $\mathbf{n}_{i,j}$.
2. **Faces verticales $(i - \frac{1}{2}, j)$ (notées `u(i, j)` dans le code, face Ouest de la cellule $(i,j)$) :**
   Vitesse horizontale $u_{i,j}$, masse volumique de face $\rho_{u,i,j}$, force capillaire horizontale $F_{\sigma, x, i, j}$, quantité de mouvement intermédiaire $(\rho u)^*_{i,j}$.
3. **Faces horizontales $(i, j - \frac{1}{2})$ (notées `v(i, j)` dans le code, face Sud de la cellule $(i,j)$) :**
   Vitesse verticale $v_{i,j}$, masse volumique de face $\rho_{v,i,j}$, force capillaire verticale $F_{\sigma, y, i, j}$, flottabilité de Boussinesq et quantité de mouvement intermédiaire $(\rho v)^*_{i,j}$.
4. **Sommets / Nœuds $(i - \frac{1}{2}, j - \frac{1}{2})$ (notés `mu_node(i, j)` dans le code) :**
   Viscosité harmonique nodale $\mu_{\text{node}, i, j}$ utilisée pour discrétiser les dérivées croisées de cisaillement $\mu(\partial_y u + \partial_x v)$.

```text
                       v(i, j+1)  [face Nord de (i,j)]
                 +---------^---------+
                 |                   |
                 |  p(i,j), C(i,j)   |
      u(i, j)  ---->  rho, mu, T   ---->  u(i+1, j)
   [face Ouest]  |  kappa, chi_T     |  [face Est]
                 |                   |
                 +---------^---------+
                       v(i, j)    [face Sud de (i,j)]
     mu_node(i, j) au coin Sud-Ouest (i-1/2, j-1/2)
```

### Dualité Divergence-Gradient sur Grille MAC

Avec cette convention d'indiçage exacte du code (`u(i,j)` à gauche et `u(i+1,j)` à droite de la cellule `(i,j)`), la divergence discrète $(\nabla_h \cdot \mathbf{u})_{i,j}$ et le gradient de pression s'écrivent :

$$(\nabla_h \cdot \mathbf{u})_{i,j} = \frac{u_{i+1,j} - u_{i,j}}{\Delta x} + \frac{v_{i,j+1} - v_{i,j}}{\Delta y}$$

$$\left(\frac{\partial p}{\partial x}\right)_{i,j} = \frac{p_{i,j} - p_{i-1,j}}{\Delta x} \quad (\text{sur face } u_{i,j}), \qquad \left(\frac{\partial p}{\partial y}\right)_{i,j} = \frac{p_{i,j} - p_{i,j-1}}{\Delta y} \quad (\text{sur face } v_{i,j})$$

Cette dualité adjointe ($\mathbf{D} = -\mathbf{G}^T$ à un facteur volumique $\Delta x \Delta y$ près) assure la symétrie parfaite de la matrice de Stokes et l'annulation exacte de la divergence dans les solveurs monolithiques et AAL.

---

## Schémas Convectifs et Linéarisations de l'Inertie

Le traitement du terme non linéaire $(\mathbf{u} \cdot \nabla)\mathbf{u}$ est configuré par deux clés complémentaires dans `physics` (ou au niveau racine) :
- `physics.convective_scheme` : `"UPWIND"`, `"CENTERED"`, `"TVD_SUPERBEE"`.
- `physics.convective_linearization` : `"PICARD"`, `"ADAMS_BASHFORTH"` (alias `"AB2"`), `"NEWTON"`.

### Schémas Spatiaux (`convective_scheme`)

1. **Décentré Amont d'Ordre 1 (`"UPWIND"`) :**
   Inconditionnellement stable et monotone même à grand nombre de Reynolds de maille $Re_\Delta > 2$, au prix d'une diffusion numérique $\mathcal{O}(\Delta x)$.
2. **Centré d'Ordre 2 (`"CENTERED"`) :**
   Évalue les vitesses aux faces par moyenne arithmétique symétrique, atteignant l'ordre $\mathcal{O}(\Delta x^2)$ sans dissipation artificielle, idéal pour la DNS lorsque $Re_\Delta \le 2$.
3. **Haute Résolution TVD (`"TVD_SUPERBEE"`) :**
   Combine la précision d'ordre 2 dans les zones régulières et le limiteur non linéaire Superbee près des gradients raides.

### Linéarisations Temporelles du Terme Convectif (`convective_linearization`)

1. **Linéarisation de Picard (`"PICARD"`, défaut) :**
   La vitesse advectante est figée au pas connu $\mathbf{u}^n$ : $(\mathbf{u}^n \cdot \nabla)\mathbf{u}^{n+1}$. Simple et robuste, précise à l'ordre 1 en temps sur l'advection.
2. **Extrapolation d'Adams-Bashforth d'Ordre 2 (`"ADAMS_BASHFORTH"` / `"AB2"`) :**
   Lorsque l'historique $\mathbf{u}^{n-1}$ est disponible (`fields.has_nm1 == true`), la vitesse advectrice est extrapolée à l'ordre 2 selon l'Équation 10 de *Fluids* 2021 :
   $$\mathbf{u}_{\text{adv}} = 2 \mathbf{u}^n - \mathbf{u}^{n-1} \quad \text{ou} \quad \text{flux extrapolés } \frac{3}{2}\mathcal{C}(\mathbf{u}^n) - \frac{1}{2}\mathcal{C}(\mathbf{u}^{n-1})$$
   préservant la linéarité du système implicite tout en atteignant l'ordre 2 temporel.
3. **Linéarisation Complète de Newton (`"NEWTON"`) :**
   Inclut le terme jacobien croisé $(\delta \mathbf{u} \cdot \nabla)\mathbf{u}^n$ en plus de $(\mathbf{u}^n \cdot \nabla)\mathbf{u}^{n+1}$ dans l'opérateur matriciel `apply_system_matrix` (et le terme quadratique $-(\mathbf{u}^n \cdot \nabla)\mathbf{u}^n$ au second membre `assemble_rhs`), accélérant la convergence sur les écoulements stationnaires ou à grands pas de temps.

---

## Discrétisation Temporelle et Pas de Temps Adaptatif

### Schémas Temporels Implémentés (`time.time_scheme`)

1. **Euler Implicite d'Ordre 1 (`"EULER"`, défaut) :**
   $$\frac{\mathbf{u}^{n+1} - \mathbf{u}^n}{\Delta t}, \qquad \alpha_{\text{time}} = 1.0$$
2. **Schéma de Gear / BDF2 d'Ordre 2 (`"GEAR"` ou `"BDF2"`) :**
   Dès le pas 2 (`fields.has_nm1 == true`), le terme instationnaire utilise la formule rétrograde d'ordre 2 à 3 niveaux temporels :
   $$\frac{3 \mathbf{u}^{n+1} - 4 \mathbf{u}^n + \mathbf{u}^{n-1}}{2 \Delta t}, \qquad \alpha_{\text{time}} = 1.5$$
   (avec amorçage automatique par Euler implicite $\alpha_{\text{time}} = 1.0$ au pas 1).
3. **Transport Inertiel Conservatif `SSP-RK3` (`physics.enable_ssp_rk3_momentum = true`) :**
   Activé automatiquement lorsque `compressible = true` (ou manuellement), un schéma de Runge-Kutta TVD à 3 étages transporte la quantité de mouvement $(\rho u, \rho v)$ de façon conservative avant l'étape de couplage vitesse-pression (étape 1 de *Fluids* 2024).

### Pas de Temps Adaptatif (`time.adaptive_dt = true`)

À chaque itération, `Simulation::compute_time_step()` (`sim_compute_time_step` en Fortran) évalue la vitesse maximale globale $U_{\max} = \max_{i,j} \sqrt{|u_{i,j}|^2 + |v_{i,j}|^2}$ via `MPI_Allreduce(..., MPI_MAX)` et borne le pas de temps par :

$$\Delta t_{\text{CFL}} = \text{CFL}_{\max} \frac{\min(\Delta x, \Delta y)}{\max(U_{\max}, 10^{-6})}$$

$$\Delta t_{\sigma} = \begin{cases} \sqrt{\frac{(\rho_1 + \rho_2) \, \min(\Delta x, \Delta y)^3}{4 \pi \sigma}} & \text{si } \text{enable\_surface\_tension} \text{ et } \sigma > 0 \\ +\infty & \text{sinon} \end{cases}$$

$$\Delta t = \min\left(\Delta t_{\text{init}}, \, \Delta t_{\text{CFL}}, \, \Delta t_{\sigma}\right)$$

où $\Delta t_{\text{init}}$ est `time.dt` et $\text{CFL}_{\max}$ est `time.cfl_max` (le pas est également tronqué sur le dernier pas pour atteindre exactement `time.t_final`).

---

## Suivi d'Interface VOF et Méthode Hybride FC-WY-CVI-Marker

### Reconstruction Géométrique et Advection VOF (`vof.*`)

Le module `InterfaceTracker` (`module_interface_vof` en Fortran) met à disposition deux familles d'advection de la fraction volumique $C$ :

1. **Schémas géométriques (`"WEYMOUTH_YUE"` par défaut, et `"CIAM"`) :**
   - Couplés aux algorithmes de reconstruction affine $n_x x + n_y y + \alpha = 0$ : `"PLIC"` (Youngs), `"ELVIRA"` (6 candidats de fonctions de hauteur, ordre 2 exact sur droites), `"LVIRA"` (moindres carrés non linéaires sur $3\times 3$) ou `"SLIC"`.
   - **Schéma conservatif de Weymouth & Yue (2010)** avec splitting alterné de Strang (X puis Y aux pas pairs, Y puis X aux pas impairs) : fixe l'indicatrice binaire de dilatation $C_{d,i,j} = \mathbf{1}_{C_{i,j}^n \ge 0.5}$ au début du pas et transporte :
     $$C_{i,j}^* = C_{i,j}^n - \frac{\Delta V_{i+1,j}^x - \Delta V_{i,j}^x}{\Delta x \Delta y} + C_{d,i,j} \frac{\Delta t}{\Delta x}(u_{i+1,j}^{n+1} - u_{i,j}^{n+1})$$
     $$C_{i,j}^{n+1} = C_{i,j}^* - \frac{\Delta V_{i,j+1}^y - \Delta V_{i,j}^y}{\Delta x \Delta y} + C_{d,i,j} \frac{\Delta t}{\Delta y}(v_{i,j+1}^{n+1} - v_{i,j}^{n+1})$$
     où $\Delta V$ est calculé exactement par découpe géométrique Sutherland-Hodgman + Shoelace (`cut_box_area`). Sous $\text{CFL}_{\text{vof}} < 0.5$ et $\nabla_h \cdot \mathbf{u}^{n+1} = 0$, la somme des flux télescope et **conserve le volume à la précision du solveur linéaire**.
2. **Schéma algébrique (`"VOF_TVD"` / `"TVD"`) :**
   Sans reconstruction géométrique (`reconstruction_method` forcée à `"NONE"`), piloté par les limiteurs `"SUPERBEE"`, `"VAN_LEER"`, `"MUSCL"` ou `"HRIC"`.

### Méthode Hybride FC-WY-CVI-Marker (Variante B Exacte sur Boîtes MAC Décalées)

Lorsqu'on active `vof.use_markers = true` (ou `markers.enabled = true`), COVOF 2D associe au transport VOF-WY un **nuage de marqueurs lagrangiens non connectés** (`InterfaceMarkers` / `module_markers`) sans arêtes ni maillage surfacique, totalement insensible aux ruptures et coalescences.

#### 1. Interpolation Conservative de la Vitesse (CVI — Variante B exacte)
Selon le théorème de Liouville, l'interpolation bilinéaire standard sur grille MAC produit une divergence continue non nulle d'ordre $\mathcal{O}(h)$ qui crée des amas artificiels de marqueurs. COVOF 2D implémente la **Variante B exacte de Wang et al. (2015) / Pusok et al. (2017)** (`interpolate_velocity_cvi`) en interpolant chaque composante sur **sa propre boîte décalée MAC** :
- Boîte pour $u$ : $\Omega_u(i, j_u) = [x_i^f, x_{i+1}^f] \times [y_{j_u}^c, y_{j_u+1}^c]$ avec $j_u = j - 1$ si $y_p < y_j^c$ sinon $j$, et coordonnées réduites $(\xi_u, \eta_u) \in [0, 1]^2$.
- Boîte pour $v$ : $\Omega_v(i_v, j) = [x_{i_v}^c, x_{i_v+1}^c] \times [y_j^f, y_{j+1}^f]$ avec $i_v = i - 1$ si $x_p < x_i^c$ sinon $i$, et coordonnées réduites $(\xi_v, \eta_v) \in [0, 1]^2$.
- Formule CVI exacte avec corrections quadratiques solénoïdales ($\Delta v_E = v_{i_v+1,j+1} - v_{i_v+1,j}$, $\Delta v_W = v_{i_v,j+1} - v_{i_v,j}$, $\Delta u_N = u_{i+1,j_u+1} - u_{i,j_u+1}$, $\Delta u_S = u_{i+1,j_u} - u_{i,j_u}$) :
  $$u_{\text{CVI}}(x_p, y_p) = u_{\text{bilin}}(\xi_u, \eta_u) + \frac{\Delta x}{2 \Delta y} \, \xi_u (1 - \xi_u) (\Delta v_E - \Delta v_W)$$
  $$v_{\text{CVI}}(x_p, y_p) = v_{\text{bilin}}(\xi_v, \eta_v) + \frac{\Delta y}{2 \Delta x} \, \eta_v (1 - \eta_v) (\Delta u_N - \Delta u_S)$$
En dérivant analytiquement dans chacun des 4 sous-quadrants de $(i, j)$, tous les termes linéaires croisés s'annulent identiquement :
$$\frac{\partial u_{\text{CVI}}}{\partial x} + \frac{\partial v_{\text{CVI}}}{\partial y} \equiv \frac{u_{i+1,j} - u_{i,j}}{\Delta x} + \frac{v_{i,j+1} - v_{i,j}}{\Delta y} = D_{i,j}$$
tout en restant strictement dans la couche unique de halo `(0:nx+1, 0:ny+1)`. *(L'option `marker_advection_scheme = "BILINEAR_RK2"` permet de désactiver la correction quadratique à des fins d'étude comparative).*

#### 2. Advection RK2-Heun, Migration MPI, Élagage et Courbure Parabolique
À l'intérieur de `InterfaceTracker::advect` (`advect_vof`), juste après l'advection VOF et le lissage SVOF/diffusif et juste avant `compute_curvature_and_CSF` :
1. **Advection RK2-Heun (`advect_rk2` / `advect_markers_rk2`) :** Prédit $\mathbf{x}_p^* = \mathbf{x}_p^n + \Delta t \, \mathbf{u}_{\text{CVI}}(\mathbf{x}_p^n; \mathbf{u}^n)$, effectue si nécessaire un routage MPI aller-retour à 4 passes pour évaluer $\mathbf{u}_{\text{CVI}}(\mathbf{x}_p^*; \mathbf{u}^{n+1})$ sur le rang propriétaire, puis corrige $\mathbf{x}_p^{n+1} = \mathbf{x}_p^n + \frac{\Delta t}{2}[\mathbf{u}_1 + \mathbf{u}_2]$.
2. **Migration MPI en 2 passes (`migrate_mpi` / `migrate_markers_mpi`) :** Passe X (Ouest $\leftrightarrow$ Est avec repliement périodique éventuel) puis Passe Y (Sud $\leftrightarrow$ Nord) sans nécessiter de voisins diagonaux.
3. **Élagage et ré-ensemencement sélectif (`prune_and_reseed`) :**
   - Supprime les marqueurs tombés en phase pure ($C < 10^{-6}$ ou $C > 1 - 10^{-6}$ dans le voisinage $3\times 3$) et élague par sous-octants si une cellule dépasse `marker_max_per_cell` (défaut `8`).
   - Ensemence uniquement les cellules interfaciales dépeuplées ($N_p = 0$) avec `markers_per_cell` marqueurs (défaut `4`) sans jamais écraser ni projeter les marqueurs actifs existants.
4. **Courbure et normale dans le repère tourné $(\tau, \eta)$ (`compute_curvature_from_markers`) :**
   - Échange les marqueurs fantômes dans une bande $2.5 h$ entre rangs MPI.
   - Filtre les marqueurs voisins par proximité et cohérence angulaire $\mathbf{n}_p \cdot \mathbf{n}_{i,j}^{\text{ref}} \ge 0.5$ (protégeant les filaments et jupes minces contre le mélange des faces opposées).
   - Ajuste une parabole pondérée $\eta(\tau) = a \tau^2 + b \tau + c$ par moindres carrés $3\times 3$ :
     $$\mathbf{n}_{i,j}^{\text{marker}} = \frac{\mathbf{e}_\eta - b \mathbf{e}_\tau}{\sqrt{1 + b^2}}, \qquad \kappa_{i,j}^{\text{marker}} = \frac{2a}{(1 + b^2)^{3/2}}$$
     avec écrêtage Nyquist $|\kappa| \le 2/\min(\Delta x, \Delta y)$, passe conservative de lissage tangentiel et repli sur la courbure CSF si $M < 4$.
   - Si `marker_normals_plic = true`, injecte $\mathbf{n}^{\text{marker}}$ (sous réserve de cohérence angulaire $\ge 0.5$ avec le gradient VOF) dans les segments PLIC pour l'étape suivante (`apply_marker_normals_to_segments`).

---

## Formulations du Couplage Vitesse-Pression

COVOF 2D propose cinq formulations pilotées par `solver.type` :

### 1. Lagrangien Augmenté Standard (`"AUGMENTED_LAGRANGIAN"`, défaut)
Ajoute le terme de pénalité $-r_{\text{al}} \nabla(\nabla \cdot \mathbf{u}^{k+1})$ à l'équation de quantité de mouvement avec un scalaire uniforme `solver.r_al`, puis actualise la pression par boucle d'Uzawa (`solver.uzawa_iter`) :

$$p^{k+1} = p^k - r_{\text{al}} (\nabla_h \cdot \mathbf{u}^{k+1})$$

avec arrêt anticipé dès que $\|\nabla_h \cdot \mathbf{u}^{k+1}\|_\infty \le \text{uzawa\_tol}$.

### 2. Lagrangien Augmenté Adaptatif Incompressible (`"AAL"` / `"ADAPTATIVE_AUGMENTED_LAGRANGIAN"`)
Détermine cellule par cellule le paramètre optimal $r(t, M)$ (Vincent et al., *Computers & Fluids* 33, 2004, Éq. 14) en pondérant les sept échelles physiques locales par `solver.k_aal` ($K_{\text{aal}} \in [10, 1000]$, défaut `100.0`) :

$$r_{i,j} = K_{\text{aal}} \max\left( \alpha_{\text{time}}\rho \frac{L_0^2}{t_0}, \; \rho u_0 L_0, \; \rho \frac{L_0^2}{u_0}\|\mathbf{g}\|, \; \frac{p_0 L_0}{u_0}, \; \mu, \; \mu \frac{\Delta x \Delta y}{\mathcal{K}}, \; \frac{\sigma}{u_0} \right)$$

Les échelles de référence (`l0_ref`, `t0_ref`, `u0_ref`, `p0_ref`) sont calculées automatiquement si laissées à `0.0`.

### 3. Lagrangien Augmenté Compressible Bas-Mach (`"C-AAL"` / `"CAAL"`)
Active automatiquement `compressible = true` et `adaptative_al = true`. Régularise la pénalité dans les zones compressibles ($\chi_T > 0$) et actualise la pression selon la contrainte barotrope :

$$r_{\text{eff}, i,j} = \frac{r_{i,j}}{1 + r_{i,j} \frac{\chi_{T,i,j}}{\Delta t}}, \qquad p_{i,j}^{k+1} = p_{i,j}^k - r_{\text{eff}, i,j} \left[ (\nabla_h \cdot \mathbf{u}^{k+1})_{i,j} + \frac{\chi_{T,i,j}}{\Delta t}(p_{i,j}^k - p_{i,j}^n) \right]$$

### 4. Résolution Monolithique Tout-Couplé (`"MONOLITHIC"` / `"FULLY_COUPLED"` / `"FC"`)
Assemble simultanément la vitesse et la pression dans le système point-selle sans pénalité ($r = 0$) :

$$\begin{pmatrix} \mathbf{A}_u & \mathbf{A}_{uv} & \mathbf{G}_x \\ \mathbf{A}_{vu} & \mathbf{A}_v & \mathbf{G}_y \\ \mathbf{D}_x & \mathbf{D}_y & -\mathbf{C}_p \end{pmatrix} \begin{pmatrix} \mathbf{u}^{n+1} \\ \mathbf{v}^{n+1} \\ \mathbf{p}^{n+1} \end{pmatrix} = \begin{pmatrix} \mathbf{f}_u \\ \mathbf{f}_v \\ \mathbf{g}_p \end{pmatrix}$$

où $\mathbf{C}_p = \text{diag}(\chi_T \Delta x \Delta y / \Delta t)$ en compressible ($\mathbf{0}$ en incompressible). La ligne de continuité n'étant pas divisée par $\rho$, ce système garantit une divergence discrète à la tolérance linéaire `solver.tol` (typiquement $10^{-10}$) quel que soit le saut $\rho_1 / \rho_2$.

### 5. Méthode de Projection Standard (`"STANDARD_PROJECTION"` / `"PROJECTION"`)
Résout successivement la prédiction visqueuse $\mathbf{u}^*$, l'équation de Poisson en pression $\nabla \cdot (\frac{1}{\rho}\nabla \phi) = \frac{1}{\Delta t}\nabla \cdot \mathbf{u}^*$, puis corrige $\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho}\nabla \phi$.

---

## Solveurs Linéaires de Krylov et Préconditionneurs

COVOF 2D implémente une architecture **Matrix-Free** (l'opérateur complet $A \mathbf{x}$ est appliqué par parcours des pochoirs MAC avec synchronisation MPI des halos via `apply_system_matrix`) complétée par des blocs locaux pour le préconditionnement (`Fuu`, `Fvv`, `SchurLaplace` en CSR côté C++, stencils 5-points côté Fortran).

### Solveurs Linéaires (`solver.linear_solver`)

| Mot-Clé | Algorithme | Domaine d'Emploi Recommandé |
| :--- | :--- | :--- |
| `"BICGSTAB2"` (défaut) | BiCGSTAB(2) (Sleijpen & van der Vorst 1993) | Opérateurs asymétriques convectifs, AAL, C-AAL et Monolithique |
| `"BICGSTAB"` | BiCGSTAB standard (van der Vorst 1992) | Systèmes asymétriques modérément raides |
| `"CG"` / `"CONJUGATE_GRADIENT"` | Gradient Conjugué Préconditionné | Opérateurs symétriques définis positifs (Stokes pur, Poisson, diffusion) |
| `"GMRES"` | GMRES(30) redémarré (Arnoldi + Givens) | Systèmes fortement non normaux |
| `"MINRES"` | Minimal Residual (Paige & Saunders 1975) | Systèmes point-selle symétriques indéfinis |
| `"DIRECT"` / `"MUMPS"` | Factorisation LU LAPACK (`dgetrf`/`dgetrs`) | Petits maillages de référence ou validation exacte (avec mise en cache si $\Delta t$ constant) |

*(Remarque : C++ dispose également des passerelles `PetscSolverBridge` (`"PETSC"`) et `HypreSolverBridge` (`"HYPRE"`) qui basculent automatiquement sur `BiCGSTAB2Solver` si PETSc/HYPRE ne sont pas liés).*

### Préconditionneurs (`solver.precond_type` / `solver.preconditioner`)

| Mot-Clé (et Alias) | Formulation | Recommandation |
| :--- | :--- | :--- |
| `"JACOBI"` / `"DIAGONAL"` (défaut) | Inversion de la diagonale locale enrichie (`diag_u`, `diag_v`, `diag_p`) | Optimal et 100 % parallèle pour `AAL` et `C-AAL` |
| `"BLOCK_TRIANGULAR"` / `"BT_SBP"` / `"SCHUR"` | Factorisation par blocs triangulaires $P = \begin{pmatrix} D_u & G \\ 0 & S \end{pmatrix}$ avec 8 balayages SSOR sur le complément de Schur $S$ et 4 balayages sur $F_u, F_v$ | Indispensable pour le mode `"MONOLITHIC"` |
| `"SSOR"` / `"GAUSS_SEIDEL"` | Lissage avant/arrière symétrique ($\omega = 1.0$) | Écoulements modérément visqueux |
| `"ILU0"` / `"ILU"` | Factorisation incomplète ILU(0) locale par sous-domaine | Systèmes raides sur sous-domaines compacts |
| `"CUSTOM"` | Hook utilisateur (`CustomPreconditioner::set_custom_operator` / `custom_precond_hook`) | Prototypage de nouveaux préconditionneurs |
| `"NONE"` | Identité ($P = I$) | Diagnostic de conditionnement |
# Géométrie, Maillages et Conditions aux Limites

## Génération et Topologie du Maillage Cartésien MAC

Le solveur **COVOF 2D** opère sur des domaines rectangulaires $[0, L_x] \times [0, L_y]$ discrétisés par un maillage cartésien structuré uniforme selon la topologie décalée MAC (*Marker-and-Cell*).

### Définition Géométrique du Domaine (`mesh`)

Le domaine de calcul est entièrement spécifié dans le bloc `mesh` du fichier JSON par quatre paramètres exacts :
- `mesh.nx` (`integer`, défaut `64`) : nombre global de cellules de calcul selon l'axe horizontal $X$.
- `mesh.ny` (`integer`, défaut `64`) : nombre global de cellules de calcul selon l'axe vertical $Y$.
- `mesh.lx` (`double`, défaut `1.0` m) : longueur physique totale du domaine selon $X$.
- `mesh.ly` (`double`, défaut `1.0` m) : hauteur physique totale du domaine selon $Y$.

Les pas d'espace constants sont calculés à l'initialisation selon :

$$\Delta x = \frac{L_x}{N_x}, \qquad \Delta y = \frac{L_y}{N_y}$$

### Coordonnées des Centres, Faces et Sommets

Pour une maille repérée localement par $(i, j)$ avec décalage global `(i_offset, j_offset)` en C++ (base 0) ou `(i_start_global - 1, j_start_global - 1)` en Fortran (base 1), les fonctions d'accès physiques retournent exactement :

- **Centres de cellules (Pression $p$, Fraction VOF $C$, $C_s$, $\widetilde{C}$, $\rho$, $\mu$, $T$, $\kappa$, $\chi_T$) :**
  $$x_c(i) = \left(i_{\text{offset}} + i - \frac{1}{2}\right)\Delta x, \qquad y_c(j) = \left(j_{\text{offset}} + j - \frac{1}{2}\right)\Delta y$$
- **Faces verticales (Vitesse $u$, $\rho_u$, $F_{\sigma, x}$) :**
  $$x_u(i) = \left(i_{\text{offset}} + i - 1\right)\Delta x, \qquad y_c(j) = \left(j_{\text{offset}} + j - \frac{1}{2}\right)\Delta y$$
  La face `u(1, j)` coïncide avec le bord gauche $x = 0$ du sous-domaine (ou du domaine global), et la face `u(nx_local + 1, j)` avec le bord droit.
- **Faces horizontales (Vitesse $v$, $\rho_v$, $F_{\sigma, y}$) :**
  $$x_c(i) = \left(i_{\text{offset}} + i - \frac{1}{2}\right)\Delta x, \qquad y_v(j) = \left(j_{\text{offset}} + j - 1\right)\Delta y$$
  La face `v(i, 1)` coïncide avec le bord inférieur $y = 0$, et `v(i, ny_local + 1)` avec le bord supérieur.
- **Sommets / Nœuds (Viscosité harmonique nodale $\mu_{\text{node}}$) :**
  $$\big(x_u(i), \, y_v(j)\big)$$

---

## Décomposition de Domaine et Parallélisme Cartésien 2D (MPI pur)

Pour l'exécution distribuée, COVOF 2D partitionne la grille globale $(N_x, N_y)$ en $P = P_x \times P_y$ sous-domaines rectangulaires gérés par des processus MPI indépendants (**MPI pur**, sans OpenMP).

```text
       +-----------------------+-----------------------+
       |  Rang MPI (0, 1)      |  Rang MPI (1, 1)      |
       |  is_boundary_left     |  is_boundary_right    |
       |  is_boundary_top      |  is_boundary_top      |
Ny / 2 +-----------------------+-----------------------+
       |  Rang MPI (0, 0)      |  Rang MPI (1, 0)      |
       |  is_boundary_left     |  is_boundary_right    |
       |  is_boundary_bottom   |  is_boundary_bottom   |
       +-----------------------+-----------------------+
       0                    Nx / 2                     Nx
```

### Topologie Virtuelle Cartésienne (`MPI_Cart_create`)

1. **Factorisation automatique et gestion des restes :**
   `MPI_Dims_create(num_procs, 2, dims)` détermine la grille $(P_x, P_y)$ la plus équilibrée. Si $N_x$ n'est pas divisible par $P_x$, les `nx_global % px` premiers rangs reçoivent une maille supplémentaire (`Mesh2D::init` / `init_mesh`).
2. **Périodicité et voisins cardinaux :**
   Les drapeaux `periodic_x` et `periodic_y` sont activés automatiquement si `left`/`right` ou `bottom`/`top` ont pour type `"PERIODIC"`. `MPI_Cart_shift` identifie les 4 voisins directs (`rank_left`, `rank_right`, `rank_bottom`, `rank_top` en C++ ; `neighbor_west`, `neighbor_east`, `neighbor_south`, `neighbor_north` en Fortran), valant `MPI_PROC_NULL` sur les parois physiques non périodiques.

### Couche Unique de Cellules Fantômes (*Halos*)

Chaque sous-domaine alloue exactement **1 couche de cellules fantômes** de chaque côté :
- Dimensions allouées pour les scalaires centrés : `nx_ghost = nx_local + 2`, `ny_ghost = ny_local + 2`, soit des indices locaux valides $i \in [0, N_{x,\text{loc}}+1]$ et $j \in [0, N_{y,\text{loc}}+1]$.
- Zone intérieure utile : $i \in [1, N_{x,\text{loc}}]$, $j \in [1, N_{y,\text{loc}}]$.
- L'échange non bloquant (`MpiComm::sync_cell_scalar`, `sync_u`, `sync_v` en C++ ; `exchange_halos_scalar`, `exchange_halos_velocity` en Fortran) procède en deux passes séquentielles (`MPI_Irecv`/`MPI_Isend` + `MPI_Waitall`) : d'abord selon $X$, puis selon $Y$ sur toute la largeur `[0 .. nx_local + 1]` afin de **propager exactement les 4 coins diagonaux** sans nécessiter de communications diagonales.

---

## Catalogue Exhaustif des Conditions aux Limites (CL)

Les quatre frontières du domaine (`left`, `right`, `bottom`, `top` dans `boundary_conditions`) sont configurées indépendamment. L'analyseur JSON reconnaît **quatre types hydrodynamiques fondamentaux (`BCType`)**, complétés par le mot-clé `profile` pour les entrées spécialisées et `T_type` pour la thermique.

:::note
Contrairement à certaines conventions génériques, les mots-clés `NO_SLIP`, `PRESSURE`, `PRESSURE_OUTLET`, `OUTFLOW`, `SYMMETRY` ou `PISTON` **n'existent pas** dans les analyseurs `Config.cpp` et `module_config.f90` (toute chaîne inconnue rebascule silencieusement sur `DIRICHLET`). Les comportements correspondants s'expriment exactement à l'aide des quatre types ci-dessous et de `profile`.
:::

### 1. Condition de Dirichlet (`"DIRICHLET"`, défaut)

Impose le vecteur vitesse $\mathbf{u} = (u_{\text{bc}}, v_{\text{bc}})^T$ sur la frontière (`u`, `v`).

- **Paroi solide immobile (*No-Slip*) :** `"type": "DIRICHLET", "u": 0.0, "v": 0.0`.
- **Paroi mobile (ex. cavité entraînée) :** `"type": "DIRICHLET", "u": 1.0, "v": 0.0` sur `"top"`.
- **Piston mobile compressible :** `"type": "DIRICHLET", "u": 0.0, "v": v_piston` sur `"bottom"`.

#### Implémentation Discrète MAC (`MpiComm.cpp` / `module_mpi_comm.f90`)
- **Composante normale sur face frontière :** imposée directement sur la face et sa fantôme externe (ex. sur `left` : $u_{1,j} = u_{0,j} = u_{\text{inlet}}$ ; sur `top` : $v_{i, N_y+1} = v_{i, N_y+2} = v_{\text{top}}$).
- **Composante tangentielle :**
  - Sur les parois horizontales (`bottom` et `top`), **extrapolation quadratique d'ordre 2 à 3 points** (exacte pour les profils paraboliques de Poiseuille lorsque $N_{y,\text{loc}} \ge 2$) :
    $$u_{i, 0} = \frac{8 u_{\text{bot}} - 6 u_{i, 1} + u_{i, 2}}{3}, \qquad u_{i, N_y+1} = \frac{8 u_{\text{top}} - 6 u_{i, N_y} + u_{i, N_y-1}}{3}$$
    *(avec repli linéaire $2 u_{\text{wall}} - u_{\text{int}}$ si $N_{y,\text{loc}} < 2$)*.
  - Sur les parois verticales (`left` et `right`), réflexion linéaire d'ordre 2 :
    $$v_{0, j} = 2 v_{\text{left}} - v_{1, j}, \qquad v_{N_x+1, j} = 2 v_{\text{right}} - v_{N_x, j}$$
- **Pression et VOF :** Neumann homogène sur les parois ($p_{\text{ghost}} = p_{\text{int}}$, $C_{\text{ghost}} = C_{\text{int}}$).

#### Profils Spécialisés (`boundary_conditions.<face>.profile`)

| Valeur de `profile` | Face Supportée | Paramètres Associés | Formulation Mathématique Exacte |
| :--- | :--- | :--- | :--- |
| `"UNIFORM"` (défaut) | Toutes | `u`, `v` | $u = u_{\text{bc}}, \; v = v_{\text{bc}}$ uniforme sur toute la face |
| `"PARABOLIC"` | `left` | `u_mean` (C++), `u` (Fortran) | $u(y) = 6 \, u_{\text{mean}} \frac{y (L_y - y)}{L_y^2}$ (profil de Poiseuille établi, $\bar{u} = u_{\text{mean}}$) |
| `"POISEUILLE_STEP"` | `left` | `u_mean`, `y_step` (C++) / `obs_y1` (Fortran) | $u(y) = 0$ si $y < y_{\text{step}}$ ; $u(y) = 6 \, u_{\text{mean}} \frac{(y - y_{\text{step}})(L_y - y)}{(L_y - y_{\text{step}})^2}$ si $y \ge y_{\text{step}}$ |
| `"JET_INLET"` | `top` | `v`, `jet_center`, `jet_width` | Injecteur centré en $x_c$ sur $|x - x_c| \le W/2$ avec $v_{\text{in,eff}} = v_{\text{top}} \frac{W}{D_{\text{disc}}}$, $C_{\text{ghost}} = 0$ (Fluide 2), $u=0$ ; et **débit compensatoire d'évacuation du gaz** hors buse $v_{\text{out}} = -v_{\text{top}} \frac{W}{L_x - D_{\text{disc}}}$ assurant $\int_{\partial\Omega} \mathbf{u}\cdot\mathbf{n}\,ds = 0$ |

### 2. Condition de Neumann (`"NEUMANN"`)

Prescrit une dérivée normale nulle pour la vitesse, la pression et la fraction volumique :

$$\frac{\partial u}{\partial n} = 0, \qquad \frac{\partial v}{\partial n} = 0, \qquad \frac{\partial p}{\partial n} = 0, \qquad \frac{\partial C}{\partial n} = 0$$

- **Implémentation :** Copie directe de l'intérieur vers la couche fantôme (ex. sur `right` : $u_{N_x+2, j} = u_{N_x+1, j}$, $v_{N_x+1, j} = v_{N_x, j}$, $p_{N_x+1, j} = p_{N_x, j}$).
- **Usage :** Sortie libre (*outflow*) en aval de canal, marche descendante ou sillage de cylindre.

### 3. Glissement Parfait / Symétrie (`"FREE_SLIP"`)

Impose l'imperméabilité normale et l'absence de cisaillement tangentiel :

$$\mathbf{u} \cdot \mathbf{n} = 0, \qquad \frac{\partial (\mathbf{u} \cdot \boldsymbol{\tau})}{\partial n} = 0$$

- **Implémentation :** Sur `bottom` par exemple, $v_{i, 1} = v_{i, 0} = 0$ (imperméabilité stricte) et $u_{i, 0} = u_{i, 1}$ (cisaillement nul).
- **Usage :** Plans de symétrie axiale, parois sans couche limite pour la rupture de barrage ou l'ascension de bulle, surface supérieure libre idéalisée. Notez la graphie exacte `"FREE_SLIP"` (avec tiret bas).

### 4. Périodicité Spatiale (`"PERIODIC"`)

Relie la frontière Ouest à la frontière Est (`left` et `right` en `"PERIODIC"`) ou Sud à Nord (`bottom` et `top` en `"PERIODIC"`). Les communications MPI connectent directement les sous-domaines extrêmes via `periods` dans `MPI_Cart_create` (ou copient directement sur 1 rang), et la migration des marqueurs CVI applique le repliement géométrique $\pm L_x$ ou $\pm L_y$.
- **Usage :** Nappe liquide cisaillée de Kelvin-Helmholtz (`LIQUID_SHEET`), instabilités de Rayleigh-Taylor périodiques latéralement, canaux infinis.

---

## Conditions aux Limites Thermiques et Comportement des Marqueurs

### Conditions Thermiques (`T_type`, `T`)

Indépendamment de la condition hydrodynamique, chaque bord supporte `T_type` :

| `T_type` | Équation Continue | Implémentation Fantôme (`apply_thermal_boundary_conditions`) |
| :--- | :--- | :--- |
| `"DIRICHLET"` | $T = T_{\text{bc}}$ sur $\partial\Omega$ | $T_{\text{ghost}} = 2 \, T_{\text{bc}} - T_{\text{interior}}$ (précision ordre 2) |
| `"NEUMANN"` (défaut) | $\frac{\partial T}{\partial n} = 0$ (adiabatique) | $T_{\text{ghost}} = T_{\text{interior}}$ |
| `"PERIODIC"` (C++) | $T(0,y) = T(L_x, y)$ | Synchronisation directe MPI ou copie si `dims == 1` |

### Traitement des Frontières pour les Marqueurs CVI (`InterfaceMarkers`)

Lors de l'advection RK2-Heun (`advect_rk2`), les positions prédites $\mathbf{x}_p^*$ sont traitées de manière cohérente avec les frontières physiques :
- Sur bord `"PERIODIC"` : repliement modulo $L_x$ ou $L_y$ dans $[0, L_x[ \times [0, L_y[$.
- Sur bord physique non périodique : maintien à $\epsilon_b \Delta x$ à l'intérieur des parois lors de la sous-étape d'évaluation de vitesse (réflexion sur parois `FREE_SLIP`, suppression ou sortie en zone pure hors domaine, et ré-ensemencement continu à la buse pour `JET_INLET`).

---

## Modélisation Géométrique des Obstacles Immergés (`obstacles`)

COVOF 2D représente les corps solides ou milieux poreux fixes par le tableau JSON `obstacles` (C++ supporte un nombre quelconque d'obstacles dans `std::vector<Obstacle>` ; Fortran lit le premier élément du tableau `"obstacles": [...]` via `json_get_obstacle`). Deux formes géométriques exactes sont implémentées :

### 1. Obstacle Rectangulaire (`"BOX"`)

Défini par les bornes $[x_0, x_1] \times [y_0, y_1]$ (`x0`, `y0`, `x1`, `y1`) et sa perméabilité `permeability` :

$$\Omega_s = \left\{ (x, y) \in \Omega \;\middle|\; x_0 - 10^{-9} \le x \le x_1 + 10^{-9} \;\text{ et }\; y_0 - 10^{-9} \le y \le y_1 + 10^{-9} \right\}$$

Usage typique : marche descendante (*backward-facing step*), chicane, plaque mince verticale, bloc poreux rectangulaire.

### 2. Obstacle Cylindrique Circulaire (`"CYLINDER"`)

Défini par son centre $(x_c, y_c)$ (`xc`, `yc`), son rayon $R$ (`radius`) et sa perméabilité `permeability` :

$$\Omega_s = \left\{ (x, y) \in \Omega \;\middle|\; (x - x_c)^2 + (y - y_c)^2 \le R^2 \right\}$$

Usage typique : écoulement autour d'un cylindre rigide (allée de von Kármán à $Re = 100$) ou matrice cylindrique poreuse.

Dans les deux cas, l'appartenance est testée indépendamment au centre de chaque face $u$ (`is_u_face_in_obstacle`), au centre de chaque face $v$ (`is_v_face_in_obstacle`) et au centre de cellule pour l'échelle AAL (`is_cell_in_obstacle`). Si `permeability <= 1.0e-12`, la face est traitée comme un **solide rigide adhérent** ($u = 0$ ou $v = 0$) ; si `permeability > 1.0e-12`, le terme de traînée de Darcy-Brinkman $\frac{\mu}{\mathcal{K}}\Delta x \Delta y$ est ajouté à la matrice et à la diagonale du préconditionneur.
# Dictionnaire Encyclopédique des Paramètres d'Entrée

## Structure Générale du Fichier de Configuration JSON

COVOF 2D utilise un fichier d'entrée textuel au format **JSON** analysé nativement sans dépendance externe par `Config::load_from_json` (`COVOF_2D-C++/src/Config.cpp`) et `load_config_from_json` (`COVOF_2D-FORTRAN/src/module_config.f90`). Les commentaires sur une ligne débutant par `//` sont automatiquement supprimés avant l'extraction.

Le dictionnaire complet s'organise autour des **sections et blocs JSON réels** suivants :

1. **Clés globales (racine JSON)** et bloc optionnel `laplacian`
2. **`mesh`** : Dimensions et résolution du maillage cartésien MAC
3. **`physics`** : Régimes, schémas convectifs, lissage SVOF/diffusif, thermique, EOS et sous-blocs `fluid1` / `fluid2`
4. **`time`** : Discrétisation temporelle, pas de temps et critère CFL adaptatif
5. **`boundary_conditions`** : Conditions hydrodynamiques, profils d'entrée et conditions thermiques sur `left`, `right`, `bottom`, `top`
6. **`initial_condition`** (alias `initial_conditions`) : Géométrie initiale de la fraction volumique $C$ et champs initiaux
7. **`solver`** : Couplage vitesse-pression (AL, AAL, C-AAL, Monolithique, Projection), solveur de Krylov et préconditionneur
8. **`vof`** et **`markers`** : Advection VOF, reconstruction géométrique, limiteurs TVD et marqueurs lagrangiens CVI
9. **`obstacles`** : Tableau d'obstacles immergés `BOX` / `CYLINDER` (pénalisation de Darcy-Brinkman)
10. **`output`** : Répertoire, préfixe de cas, intervalles VTK/logs et diagnostics CSV

:::note
Les clés sont sensibles à la casse. Lorsqu'une valeur par défaut diffère légèrement entre `Config.hpp` (C++) et `module_config.f90` (Fortran), les deux valeurs exactes sont explicitement indiquées (`C++ / Fortran`) : il est donc recommandé de renseigner explicitement les paramètres critiques dans le fichier JSON pour garantir une parité stricte.
:::

---

## Dictionnaire Exhaustif des Paramètres par Section

### 1. Clés Globales à la Racine (et Bloc `laplacian`)

| Clé JSON | Type | Défaut (C++ / Fortran) | Options / Alias / Plage | Description Physique et Numérique |
| :--- | :--- | :--- | :--- | :--- |
| `mode` | `string` | `"NAVIER_STOKES"` | `"NAVIER_STOKES"`, `"LAPLACIAN"`, `"VORTEX_ADVECTION"` | Sélectionne le pilote d'exécution dans `Simulation::run()`. `"LAPLACIAN"` résout $-\Delta\psi = f$ ; `"VORTEX_ADVECTION"` exécute le test cinématique de Rider-Kothe ; toute autre valeur lance la boucle Navier-Stokes. |
| `test_case` | `string` | `"CAS_3_1"` | Texte libre (ex. `"CAS_3_2_NEUMANN"`, `"CAS_3_7"`, `"CAS_3_8"`, `"POISEUILLE"`) | Identifiant de cas utilisé pour activer les exports analytiques dédiés (`CAS_3_2_NEUMANN`, `CAS_3_7` Poiseuille stationnaire, `CAS_3_8` Poiseuille instationnaire). |
| `case_name` | `string` | — / `"covof_case"` | Chaîne valide | Nom de cas racine lu en repli par Fortran (priorité à `output.case_name`). |
| `time_scheme` | `string` | `"EULER"` | `"EULER"`, `"GEAR"`, `"BDF2"` | Peut être spécifié à la racine ou dans `time.time_scheme` (priorité au bloc `time`). |
| `convective_linearization` | `string` | `"PICARD"` | `"PICARD"`, `"ADAMS_BASHFORTH"`, `"AB2"`, `"NEWTON"` | Peut être spécifié à la racine ou dans `physics.convective_linearization` (priorité à `physics`). |
| `enable_svof_smoothing` | `boolean` | `false` | `true`, `false` | Peut être spécifié à la racine ou dans `physics.enable_svof_smoothing`. |
| `svof_filter_width` | `double` | `1.0` | $> 0.0$ | Peut être spécifié à la racine ou dans `physics.svof_filter_width`. |
| `steady_state_tol` | `double` | `0.0` | $\ge 0.0$ | Critère d'arrêt sur $\|\mathbf{u}^{n+1} - \mathbf{u}^n\|_\infty < \text{steady\_state\_tol}$ (après `step > 10`). Désactivé si `0.0` (ignoré en diphasique Fortran). Également lu dans `solver.steady_state_tol`. |
| `source_term` (alias `f`) | `double` | `1.0` | Réel quelconque | Terme source scalaire du mode `"LAPLACIAN"`. Également lu dans `laplacian.source_term` / `laplacian.f` et `physics.source_term` / `physics.f`. |

*(Remarque : en C++, les clés des marqueurs `use_markers`, `markers_per_cell`, `marker_advection_scheme`, `marker_curvature`, `marker_normals_plic`, `marker_max_per_cell`, `write_markers_vtk` sont aussi acceptées à la racine avant d'être surchargées par `vof.*` puis `markers.*`.)*

### 2. Section `mesh` (Maillage Cartésien MAC)

| Clé JSON | Type | Défaut | Plage Valide | Description Physique et Numérique |
| :--- | :--- | :--- | :--- | :--- |
| `mesh.nx` | `integer` | `64` | $\ge 4$ | Nombre global de cellules de pression selon l'axe $X$. |
| `mesh.ny` | `integer` | `64` | $\ge 4$ | Nombre global de cellules de pression selon l'axe $Y$. |
| `mesh.lx` | `double` | `1.0` | $> 0.0\,\text{m}$ | Longueur physique totale du domaine $\Omega$ selon $X$. |
| `mesh.ly` | `double` | `1.0` | $> 0.0\,\text{m}$ | Hauteur physique totale du domaine $\Omega$ selon $Y$. |

*(Note : Le domaine physique commence toujours à l'origine $(0, 0)$ : le maillage couvre exactement $[0, L_x] \times [0, L_y]$ avec $\Delta x = L_x / N_x$ et $\Delta y = L_y / N_y$.)*

### 3. Section `physics` (Modèles Physiques, Thermique, EOS et Fluides)

#### 3.1. Paramètres Généraux du Bloc `physics`

| Clé JSON | Type | Défaut (C++ / Fortran) | Options / Plage | Description Physique et Numérique |
| :--- | :--- | :--- | :--- | :--- |
| `physics.enable_two_phase` | `boolean` | `false` | `true`, `false` | Active l'advection VOF de $C$, la mise à jour des propriétés de mélange et éventuellement les marqueurs CVI. |
| `physics.enable_surface_tension` | `boolean` | `false` | `true`, `false` | Active la force capillaire CSF équilibrée $\mathbf{f}_\sigma = \sigma \kappa \nabla C$. |
| `physics.sigma` | `double` | `0.0` | $\ge 0.0\,\text{N}\cdot\text{m}^{-1}$ | Coefficient de tension superficielle entre Fluide 1 et Fluide 2. |
| `physics.gx` | `double` | `0.0` | Réel $[\text{m}\cdot\text{s}^{-2}]$ | Composante horizontale de la gravité. |
| `physics.gy` | `double` | `0.0` | Réel $[\text{m}\cdot\text{s}^{-2}]$ | Composante verticale de la gravité (ex. `-9.81`). |
| `physics.enable_convection` | `boolean` | `true` | `true`, `false` | Active le terme d'inertie $(\mathbf{u}\cdot\nabla)\mathbf{u}$ (`false` pour résoudre les équations de Stokes pures). |
| `physics.convective_scheme` | `string` | `"UPWIND"` | `"UPWIND"`, `"CENTERED"`, `"TVD_SUPERBEE"` | Schéma spatial des flux convectifs de quantité de mouvement. |
| `physics.convective_linearization` | `string` | `"PICARD"` | `"PICARD"`, `"ADAMS_BASHFORTH"`, `"AB2"`, `"NEWTON"` | Méthode de linéarisation temporelle de $(\mathbf{u}\cdot\nabla)\mathbf{u}$ dans l'opérateur implicite. |
| `physics.enable_svof_smoothing` | `boolean` | `false` | `true`, `false` | Active le filtre de Helmholtz $(I - D\Delta)C_s = C$ pour $\rho(C_s), \mu(C_s), \mu_{\text{node}}(C_s)$ tout en gardant $C$ raide. |
| `physics.svof_filter_width` | `double` | `1.0` | $> 0.0$ | Facteur $k$ de la diffusivité du filtre Helmholtz $D = k \, \Delta x \Delta y$. |
| `physics.enable_diffusion_smoothing` | `boolean` | `false` | `true`, `false` | Active le lissage par diffusion parabolique explicite $\widetilde{C}$ pour le calcul de la courbure CSF classique. |
| `physics.diffusion_smoothing_steps` | `integer` | `3` | $\ge 1$ | Nombre de pas de pseudo-temps $\Delta\tau = 0.2 \min(\Delta x^2, \Delta y^2)$ pour calculer $\widetilde{C}$. |
| `physics.compressible` | `boolean` | `false` | `true`, `false` | Active la formulation compressible barotrope / *All-Mach* ($\chi_T > 0$, équations d'état EOS). |
| `physics.flow_regime` | `string` | `"INCOMPRESSIBLE"` | `"INCOMPRESSIBLE"`, `"COMPRESSIBLE"`, `"ALL_MACH"` | Si `"COMPRESSIBLE"` ou `"ALL_MACH"`, force `compressible = true`. |
| `physics.enable_ssp_rk3_momentum` | `boolean` | `compressible` (`false` si incomp., `true` si comp.) | `true`, `false` | Active le transport inertiel explicite conservatif SSP-RK3 de $(\rho u, \rho v)$ avant le couplage vitesse-pression. |
| `physics.enable_energy` | `boolean` | `false` | `true`, `false` | Active l'équation thermique (Boussinesq si incompressible, énergie totale $\rho e$ si compressible). *(Fortran accepte aussi l'alias `thermal.enable_energy`.)* |
| `physics.Ra` | `double` | `1.0e5` / `1.0e4` | $> 0.0$ | Nombre de Rayleigh pour la convection naturelle thermique. |
| `physics.Pr` | `double` | `0.71` | $> 0.0$ | Nombre de Prandtl $\nu / \alpha$. |
| `physics.T_ref` | `double` | `0.0` | Réel $[\text{K}]$ ou adim. | Température de référence pour la flottabilité de Boussinesq (et initialisation si compressible). |
| `physics.beta` | `double` | `3.33e-3` / `3.41e-3` | $\ge 0.0\,\text{K}^{-1}$ | Coefficient de dilatation thermique volumique. |
| `physics.viscosity_model` | `string` | `"CONSTANT"` | `"CONSTANT"`, `"EXPONENTIAL"` (et `"ASYMMETRIC_EXPONENTIAL"` en C++) | Modèle de viscosité spatiale en monophasique ($\mu(y) = \mu_0 \exp(a((y-e)/e)^2)$ si `"EXPONENTIAL"`). |
| `physics.visc_mu0` | `double` | `1.85e-5` / `1.0e-3` | $> 0.0\,\text{Pa}\cdot\text{s}$ | Viscosité de référence $\mu_0$ du modèle exponentiel. |
| `physics.visc_a` | `double` | `0.02` / `2.0` | Réel | Coefficient d'amplitude $a$ du modèle exponentiel. |
| `physics.visc_e` | `double` | `0.01` | $> 0.0\,\text{m}$ | Échelle caractéristique $e$ du modèle exponentiel. |
| `physics.source_term` (alias `f`) | `double` | `1.0` | Réel quelconque | Alias du terme source scalaire en mode `"LAPLACIAN"`. |

#### 3.2. Sous-Blocs `physics.fluid1` et `physics.fluid2`

Chaque phase (`fluid1` correspondant à $C = 1$, et `fluid2` correspondant à $C = 0$) possède les paramètres suivants :

| Clé JSON (`physics.fluid{1,2}.*`) | Type | Défaut (C++ `f1`/`f2` \| Fortran `f1`/`f2`) | Options / Plage | Description Physique et Numérique |
| :--- | :--- | :--- | :--- | :--- |
| `rho` | `double` | `1.0` / `1.0` \| `1.0` / `rho1` | $> 0.0\,\text{kg}\cdot\text{m}^{-3}$ | Masse volumique de référence de la phase. |
| `mu` | `double` | `1.0` / `1.0` \| `1.0e-3` / `mu1` | $> 0.0\,\text{Pa}\cdot\text{s}$ | Viscosité dynamique de la phase. |
| `Cp` | `double` | `1006.0` \| `1.0` (`fluid1` uniquement en Fortran) | $> 0.0\,\text{J}\cdot\text{kg}^{-1}\cdot\text{K}^{-1}$ | Capacité thermique massique à pression constante. |
| `k_thermal` | `double` | `0.026` \| `0.014` (`fluid1` uniquement en Fortran) | $> 0.0\,\text{W}\cdot\text{m}^{-1}\cdot\text{K}^{-1}$ | Conductivité thermique du fluide. |
| `beta` | `double` | `3.33e-3` (C++) | $\ge 0.0\,\text{K}^{-1}$ | Dilatation thermique propre à la phase (C++ ; Fortran utilise `physics.beta`). |
| `eos` | `string` | `"INCOMPRESSIBLE"` | `"INCOMPRESSIBLE"`, `"LINEAR"` (`"ACOUSTIC"`), `"IDEAL_GAS"` (`"PERFECT_GAS"`), `"NASG"` (`"STIFFENED_GAS"`) | Loi d'état thermodynamique utilisée lorsque `physics.compressible = true`. |
| `chi_T` | `double` | `0.0` | $\ge 0.0\,\text{Pa}^{-1}$ | Compressibilité isotherme $\chi_T = \frac{1}{\rho c_s^2}$ pour l'EOS `"LINEAR"`. |
| `gamma` | `double` | `1.4` | $> 1.0$ | Rapport des chaleurs spécifiques $\gamma = C_p / C_v$ (`IDEAL_GAS`, `NASG`). |
| `r_gas` | `double` | `287.05` | $> 0.0\,\text{J}\cdot\text{kg}^{-1}\cdot\text{K}^{-1}$ | Constante massique des gaz parfaits (`IDEAL_GAS`). |
| `cv` | `double` | `717.5` | $> 0.0\,\text{J}\cdot\text{kg}^{-1}\cdot\text{K}^{-1}$ | Capacité calorifique à volume constant (`NASG`, énergie totale). |
| `p_inf` | `double` | `0.0` | $\ge 0.0\,\text{Pa}$ | Pression de cohésion $p_\infty$ de l'EOS `NASG`. |
| `b_nasg` | `double` | `0.0` | $\ge 0.0\,\text{m}^3\cdot\text{kg}^{-1}$ | Covolume $b$ de l'EOS `NASG`. |
| `q_nasg` | `double` | `0.0` | Réel $[\text{J}\cdot\text{kg}^{-1}]$ | Énergie de référence $q$ de l'EOS `NASG`. |

### 4. Section `time` (Intégration Temporelle)

| Clé JSON | Type | Défaut (C++ / Fortran) | Options / Plage | Description Physique et Numérique |
| :--- | :--- | :--- | :--- | :--- |
| `time.time_scheme` | `string` | `"EULER"` | `"EULER"`, `"GEAR"`, `"BDF2"` | Schéma temporel implicite pour Navier-Stokes (`"GEAR"` et `"BDF2"` activent BDF2 d'ordre 2). |
| `time.dt` | `double` | `0.001` / `0.01` | $> 0.0\,\text{s}$ | Pas de temps fixe (si `adaptive_dt = false`) ou pas de temps plafond initial (si `adaptive_dt = true`). |
| `time.t_final` | `double` | `1.0` | $> 0.0\,\text{s}$ | Temps physique final d'arrêt de la simulation. |
| `time.cfl_max` | `double` | `0.4` / `0.5` | $]0.0, 1.0]$ | Nombre de Courant maximal autorisé lors du calcul automatique de $\Delta t_{\text{CFL}}$. |
| `time.adaptive_dt` | `boolean` | `false` | `true`, `false` | Active la limitation dynamique $\Delta t = \min(\Delta t_{\text{init}}, \Delta t_{\text{CFL}}, \Delta t_\sigma)$. |

### 5. Section `boundary_conditions` (`left`, `right`, `bottom`, `top`)

Chaque sous-objet `boundary_conditions.<face>` (où `<face>` $\in \{\text{left}, \text{right}, \text{bottom}, \text{top}\}$) accepte les paramètres suivants :

| Clé JSON | Type | Défaut | Options / Alias / Plage | Description Physique et Numérique |
| :--- | :--- | :--- | :--- | :--- |
| `type` | `string` | `"DIRICHLET"` | `"DIRICHLET"`, `"NEUMANN"`, `"FREE_SLIP"`, `"PERIODIC"` | Condition aux limites hydrodynamique sur la vitesse et la pression. |
| `u` | `double` | `0.0` | Réel $[\text{m}\cdot\text{s}^{-1}]$ (Fortran accepte aussi `u_mean`, `u_val`) | Valeur de vitesse horizontale imposée si `type = "DIRICHLET"`. |
| `v` | `double` | `0.0` | Réel $[\text{m}\cdot\text{s}^{-1}]$ (Fortran accepte aussi `v_val`) | Valeur de vitesse verticale imposée si `type = "DIRICHLET"`. |
| `profile` | `string` | `"UNIFORM"` | `"UNIFORM"`, `"PARABOLIC"`, `"POISEUILLE_STEP"` (C++), `"JET_INLET"` | Profil spatial de vitesse à l'entrée (`left` pour `PARABOLIC`/`POISEUILLE_STEP`, `top` pour `JET_INLET`). |
| `u_mean` | `double` | `0.0` | Réel $[\text{m}\cdot\text{s}^{-1}]$ | Vitesse débitante moyenne pour les profils `"PARABOLIC"` et `"POISEUILLE_STEP"` en C++ (Fortran lit `u`). |
| `y_step` | `double` | `0.0` | $[0, L_y[\,\text{m}$ | Hauteur de la marche sur le bord gauche pour `"POISEUILLE_STEP"` en C++ (Fortran utilise `obs_y1` si `OBS_BOX`). |
| `jet_width` | `double` | `0.08` | $> 0.0\,\text{m}$ | Largeur de la buse d'injection pour `profile = "JET_INLET"` sur la face `top`. |
| `jet_center` | `double` | `0.5` | $]0, L_x[\,\text{m}$ | Abscisse centrale $x_c$ de l'injecteur pour `profile = "JET_INLET"`. |
| `T_type` | `string` | `"NEUMANN"` | `"DIRICHLET"`, `"NEUMANN"` (et `"PERIODIC"` en C++) | Type de condition aux limites thermique sur $T$. |
| `T` | `double` | `0.0` | Réel $[\text{K}]$ (Fortran accepte aussi `T_val`) | Température imposée si `T_type = "DIRICHLET"`. |

### 6. Section `initial_condition` (Alias `initial_conditions`)

| Clé JSON | Type | Défaut (C++ / Fortran) | Options / Plage | Description Physique et Numérique |
| :--- | :--- | :--- | :--- | :--- |
| `type` | `string` | `"SINGLE_PHASE"` / `"UNIFORM"` | `"SINGLE_PHASE"`, `"DROPLET"`, `"BUBBLE"`, `"CYLINDER"`, `"DAM_BREAK"`, `"LIQUID_SHEET"`, `"VISCOUS_JET"` (`"JET_BUCKLING"`), `"RAYLEIGH_TAYLOR"`, `"VORTEX_SHEAR"`, `"SHEAR_LAYER"` (C++), `"UNIFORM"` / `"POISEUILLE"` / `"LAPLACIAN_TEST"` (Fortran) | Forme géométrique initiale de la fraction volumique $C$ (et du champ de vitesse initial associé pour `LIQUID_SHEET` et `VISCOUS_JET`). |
| `u`, `v` | `double` | `0.0` | Réel $[\text{m}\cdot\text{s}^{-1}]$ | Champ de vitesse initial uniforme dans tout le domaine. |
| `center_x`, `center_y` | `double` | `0.5` | Coordonnées $[\text{m}]$ | Centre de la goutte/bulle/cylindre/vortex, ou hauteur médiane pour `LIQUID_SHEET` / `RAYLEIGH_TAYLOR` / `SHEAR_LAYER`, ou centre X de l'amorce pour `VISCOUS_JET`. |
| `radius` | `double` | `0.25` | $> 0.0\,\text{m}$ | Rayon initial pour `DROPLET`, `BUBBLE`, `CYLINDER`, `VORTEX_SHEAR`. |
| `axis_a`, `axis_b` | `double` | `0.0` | $\ge 0.0\,\text{m}$ | Demi-axes elliptiques $a$ (X) et $b$ (Y) pour bulle/goutte oscillante de Lamb. Si `axis_a <= 0` et `pert_amp > 0`, auto-calculés par $a = R(1+\varepsilon)$ et $b = R/(1+\varepsilon)$. |
| `pert_amp` | `double` | `0.0` (`0.01`/`0.05` si RT) | $\ge 0.0$ | Amplitude de perturbation géométrique (sinusoïde Rayleigh-Taylor/Shear-Layer ou ellipticité de Lamb). |
| `pert_k` | `double` | `1.0` | $> 0.0$ | Nombre de modes / nombre d'onde de la perturbation initiale. |
| `box_x0`, `box_y0`, `box_x1`, `box_y1` | `double` | `0,0,1,1` / `0,0,0.5,0.5` | Coordonnées $[\text{m}]$ | Emprise rectangulaire $[x_0, x_1] \times [y_0, y_1]$ de Fluide 1 pour `DAM_BREAK`. |
| `thickness` | `double` | `0.0003` | $> 0.0\,\text{m}$ | Épaisseur $\delta$ de la nappe liquide horizontale ($C=0$) pour `LIQUID_SHEET`. |
| `jet_width` | `double` | `0.08` | $> 0.0\,\text{m}$ | Largeur de l'amorce liquide initiale sous la paroi supérieure pour `VISCOUS_JET`. |
| `p_init` (alias `p0`) | `double` | `0.0` | Réel $[\text{Pa}]$ | Pression thermodynamique initiale uniforme (ex. `1.01325e5` en régime `compressible = true`). |

### 7. Section `solver` (Couplage Navier-Stokes, Krylov et Préconditionneurs)

| Clé JSON | Type | Défaut (C++ / Fortran) | Options / Alias / Plage | Description Physique et Numérique |
| :--- | :--- | :--- | :--- | :--- |
| `solver.type` | `string` | `"AUGMENTED_LAGRANGIAN"` | `"AUGMENTED_LAGRANGIAN"`, `"AAL"` (`"ADAPTATIVE_AUGMENTED_LAGRANGIAN"`, `"ADAPTIVE_AUGMENTED_LAGRANGIAN"`), `"C-AAL"` (`"CAAL"`, `"COMPRESSIBLE_AAL"`), `"MONOLITHIC"` (`"FULLY_COUPLED"`, `"FC"`), `"STANDARD_PROJECTION"` (`"PROJECTION"`), `"DIRECT"` (`"MUMPS"`), `"CONJUGATE_GRADIENT"` (Fortran) | Formulation de couplage vitesse-pression. `"AAL"` et `"C-AAL"` activent automatiquement `adaptative_al = true` (`"C-AAL"` active aussi `compressible = true`). |
| `solver.linear_solver` | `string` | `"BICGSTAB2"` | `"BICGSTAB2"`, `"BICGSTAB"`, `"CG"` (`"CONJUGATE_GRADIENT"`), `"GMRES"`, `"MINRES"`, `"DIRECT"` (`"MUMPS"`), et `"PETSC"`/`"HYPRE"` en C++ | Algorithme de résolution du système linéaire (par défaut `"DIRECT"` si `solver.type = "DIRECT"`). |
| `solver.precond_type` (alias `preconditioner`) | `string` | `"JACOBI"` | `"JACOBI"` (`"DIAGONAL"`), `"BLOCK_TRIANGULAR"` (`"SCHUR"`, `"BT_SBP"`, `"BT-SBP"`, `"BLOCK"`, `"BLOCK_PRECONDITIONER"`), `"SSOR"` (`"GAUSS_SEIDEL"`), `"ILU0"` (`"ILU"`), `"CUSTOM"`, `"NONE"` | Préconditionneur appliqué à gauche ($P^{-1} A \mathbf{x} = P^{-1} \mathbf{b}$). |
| `solver.tol` | `double` | `1.0e-7` / `1.0e-5` | $]0.0, 1.0[$ | Tolérance relative d'arrêt du solveur de Krylov $\|\mathbf{r}_k\|_2 / \|\mathbf{r}_0\|_2$. |
| `solver.max_iter` | `integer` | `1000` / `300` | $\ge 1$ | Nombre maximal d'itérations de Krylov par résolution linéaire. |
| `solver.r_al` | `double` | `500.0` / `100.0` | $> 0.0\,\text{Pa}\cdot\text{s}$ | Paramètre scalaire de pénalité AL (ou valeur de repli si AAL produit $0$). |
| `solver.adaptative_al` (alias `adaptative_r`, `adaptative_r_al`) | `boolean` | `false` | `true`, `false` | Active le calcul adaptatif spatio-temporel $r(t, M)$ (Vincent et al. 2004, Éq. 14). |
| `solver.k_aal` (alias `r_factor`) | `double` | `100.0` | $> 0.0$ (typ. $[10, 1000]$) | Facteur d'amplification adimensionnel $K_{\text{aal}}$ de l'Équation 14. |
| `solver.uzawa_iter` | `integer` | `1` | $\ge 1$ | Nombre maximal d'itérations externes d'Uzawa par pas de temps (en mode AL/AAL/C-AAL). |
| `solver.uzawa_tol` | `double` | `1.0e-12` | $\ge 0.0\,\text{s}^{-1}$ | Critère d'arrêt anticipé d'Uzawa sur $\|\nabla\cdot\mathbf{u} + \frac{\chi_T}{\Delta t}(p - p^n)\|_\infty \le \text{uzawa\_tol}$. |
| `solver.l0_ref`, `t0_ref`, `u0_ref`, `p0_ref` | `double` | `0.0` | $\ge 0.0$ (`0.0` = automatique) | Échelles physiques manuelles ($L_0, t_0, u_0, p_0$) pour l'adimensionnement AAL. |
| `solver.use_precond` | `boolean` | `true` | `true`, `false` | Active (`true`) ou court-circuite (`false`) le préconditionneur. |
| `solver.steady_state_tol` | `double` | `0.0` | $\ge 0.0$ | Alias du seuil d'arrêt stationnaire dans `solver`. |

### 8. Sections `vof` et `markers` (Suivi d'Interface et Marqueurs CVI)

:::note
Les clés de marqueurs peuvent être placées soit directement dans le bloc **`vof`** (préfixe `marker_*`), soit dans le bloc dédié **`markers`** (avec les alias indiqués ci-dessous), le bloc `markers` ayant la priorité finale s'il est présent.
:::

| Clé dans `vof` | Alias dans `markers` | Type | Défaut (C++ / Fortran) | Options / Plage | Description Physique et Numérique |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `vof.advection_method` | — | `string` | `"WEYMOUTH_YUE"` | `"WEYMOUTH_YUE"`, `"CIAM"`, `"VOF_TVD"` (`"TVD"`) | Algorithme d'advection de la fraction volumique $C$. |
| `vof.reconstruction_method` | — | `string` | `"PLIC"` | `"PLIC"`, `"ELVIRA"`, `"LVIRA"`, `"SLIC"`, `"NONE"` | Méthode de reconstruction affine (forcée à `"NONE"` si `VOF_TVD`, et à `"PLIC"` si géométrique). |
| `vof.tvd_limiter` | — | `string` | `"NONE"` (`"SUPERBEE"` si `VOF_TVD`) | `"SUPERBEE"`, `"VAN_LEER"`, `"MUSCL"`, `"HRIC"`, `"NONE"` | Limiteur algébrique (forcé à `"NONE"` pour `WEYMOUTH_YUE` et `CIAM`). |
| `vof.cfl_vof` | — | `double` | `0.5` | $]0.0, 0.5]$ | Nombre de CFL VOF de référence. |
| `vof.use_markers` | `markers.enabled`, `markers.use_markers` | `boolean` | `false` | `true`, `false` | Active le nuage de marqueurs lagrangiens non connectés (méthode FC-WY-CVI-Marker). |
| `vof.markers_per_cell` | `markers.markers_per_cell` | `integer` | `4` | $\ge 2$ | Nombre de marqueurs semés initialement (et lors du ré-ensemencement) par segment PLIC coupé. |
| `vof.marker_advection_scheme` | `markers.advection_scheme`, `markers.marker_advection_scheme` | `string` | `"CVI_RK2"` | `"CVI_RK2"`, `"BILINEAR_RK2"` (`"BILINEAR"`) | Interpolation conservative CVI (Variante B exacte) ou bilinéaire MAC pour l'advection RK2-Heun. |
| `vof.marker_curvature` | `markers.use_marker_curvature`, `markers.marker_curvature` | `boolean` | `true` / `false` | `true`, `false` | Active la courbure par régression parabolique pondérée dans le repère $(\tau, \eta)$ pour le modèle CSF. |
| `vof.marker_normals_plic` | `markers.use_marker_normals_plic`, `markers.marker_normals_plic` | `boolean` | `true` / `false` | `true`, `false` | Injecte les normales issues des marqueurs dans les segments PLIC lors de l'advection VOF. |
| `vof.marker_max_per_cell` | `markers.max_per_cell`, `markers.marker_max_per_cell` | `integer` | `8` | $\ge 4$ | Seuil d'élagage par sous-octants dans une cellule afin d'éviter la surpopulation locale. |
| `vof.write_markers_vtk` | `markers.write_vtk`, `markers.write_markers_vtk` *(et `output.write_markers_vtk` en C++)* | `boolean` | `false` | `true`, `false` | Exporte le nuage de points au format VTK PolyData (`<case_name>_markers_<step:06d>.vtk`). |

### 9. Section `obstacles` (Tableau d'Obstacles Immergés Brinkman)

Déclarée sous forme de tableau JSON `"obstacles": [ { ... } ]`. C++ supporte plusieurs objets dans le tableau ; Fortran lit le premier objet défini.

| Clé JSON | Type | Défaut (C++ / Fortran) | Options / Plage | Description Physique et Numérique |
| :--- | :--- | :--- | :--- | :--- |
| `type` | `string` | `"BOX"` | `"BOX"`, `"CYLINDER"` | Forme géométrique de l'obstacle immergé. |
| `x0`, `y0`, `x1`, `y1` | `double` | `0.0` | Coordonnées $[\text{m}]$ | Bornes du rectangle $[x_0, x_1] \times [y_0, y_1]$ (si `type = "BOX"`). |
| `xc`, `yc`, `radius` | `double` | `0.0` | Coordonnées & rayon $[\text{m}]$ | Centre $(x_c, y_c)$ et rayon $R$ du disque (si `type = "CYLINDER"`). |
| `permeability` | `double` | `1.0e-14` / `1.0e-40` | $> 0.0\,\text{m}^2$ | Perméabilité de Darcy $\mathcal{K}$ : $\mathcal{K} \le 10^{-12}$ impose un solide rigide adhérent ; $\mathcal{K} > 10^{-12}$ applique la traînée poreuse $-\frac{\mu}{\mathcal{K}}\mathbf{u}$. |

### 10. Section `output` (Sorties VTK, Journaux Console et CSV)

| Clé JSON | Type | Défaut (C++ / Fortran) | Options / Plage | Description Physique et Numérique |
| :--- | :--- | :--- | :--- | :--- |
| `output.output_dir` (alias Fortran `directory`) | `string` | `"output"` | Chemin de dossier | Répertoire de sortie (créé automatiquement sur le rang 0). |
| `output.case_name` | `string` | `"covof"` / `"covof2d_case"` | Préfixe de fichier | Préfixe commun des fichiers `.vtr`, `.pvtr`, `_markers_*.vtk` et `_history.csv`. |
| `output.vtk_interval` | `integer` | `10` / `50` | $\ge 1$ | Période d'écriture des fichiers `.vtr`/`.pvtr` (et `_markers_*.vtk` si activé), en nombre de pas de temps (le pas `0` et le pas final sont toujours écrits). |
| `output.log_interval` | `integer` | `1` / `10` | $\ge 1$ | Période d'affichage du tableau de suivi sur la sortie standard du rang 0. |
| `output.write_csv_history` | `boolean` | `true` | `true`, `false` | Active l'écriture à chaque pas de temps du fichier `<output_dir>/<case_name>_history.csv`. |
| `output.write_markers_vtk` | `boolean` | `false` | `true`, `false` | Alias C++ pour activer l'export VTK PolyData des marqueurs. |

---

## Fichier de Configuration Maître Intégral et Validé

Le fichier ci-dessous n'utilise **exclusivement que des clés reconnues à 100 %** par `Config.cpp` et `module_config.f90` (exemple : goutte capillaire ou bulle ascendante avec méthode hybride **FC-WY-CVI-Marker** et lissage SVOF) :

```json
{
  // =========================================================================
  // COVOF 2D — FICHIER DE CONFIGURATION MAÎTRE (C++17 & FORTRAN 2018)
  // =========================================================================
  "mode": "NAVIER_STOKES",
  "test_case": "MAITRE_FC_WY_CVI_MARKER",
  "steady_state_tol": 0.0,

  "mesh": {
    "nx": 64,
    "ny": 128,
    "lx": 1.0,
    "ly": 2.0
  },

  "physics": {
    "enable_two_phase": true,
    "enable_surface_tension": true,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_svof_smoothing": true,
    "svof_filter_width": 0.5,
    "enable_diffusion_smoothing": false,
    "diffusion_smoothing_steps": 3,
    "enable_energy": false,
    "compressible": false,
    "flow_regime": "INCOMPRESSIBLE",
    "enable_ssp_rk3_momentum": false,
    "sigma": 24.5,
    "gx": 0.0,
    "gy": -0.98,
    "Ra": 1.0e5,
    "Pr": 0.71,
    "T_ref": 300.0,
    "beta": 3.33e-3,
    "viscosity_model": "CONSTANT",
    "fluid1": {
      "rho": 1000.0,
      "mu": 10.0,
      "Cp": 4182.0,
      "k_thermal": 0.6,
      "beta": 2.1e-4,
      "eos": "INCOMPRESSIBLE",
      "chi_T": 0.0,
      "gamma": 1.4,
      "r_gas": 287.05,
      "cv": 717.5,
      "p_inf": 0.0,
      "b_nasg": 0.0,
      "q_nasg": 0.0
    },
    "fluid2": {
      "rho": 100.0,
      "mu": 1.0,
      "Cp": 1006.0,
      "k_thermal": 0.026,
      "beta": 3.33e-3,
      "eos": "INCOMPRESSIBLE",
      "chi_T": 0.0,
      "gamma": 1.4,
      "r_gas": 287.05,
      "cv": 717.5,
      "p_inf": 0.0,
      "b_nasg": 0.0,
      "q_nasg": 0.0
    }
  },

  "time": {
    "time_scheme": "GEAR",
    "dt": 1.0e-3,
    "t_final": 3.0,
    "cfl_max": 0.35,
    "adaptive_dt": true
  },

  "initial_condition": {
    "type": "BUBBLE",
    "u": 0.0,
    "v": 0.0,
    "center_x": 0.5,
    "center_y": 0.5,
    "radius": 0.25,
    "axis_a": 0.0,
    "axis_b": 0.0,
    "pert_amp": 0.0,
    "pert_k": 1.0,
    "p_init": 0.0
  },

  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0, "profile": "UNIFORM", "T_type": "NEUMANN", "T": 300.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0, "profile": "UNIFORM", "T_type": "NEUMANN", "T": 300.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0, "profile": "UNIFORM", "T_type": "NEUMANN", "T": 300.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0, "profile": "UNIFORM", "T_type": "NEUMANN", "T": 300.0 }
  },

  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "ELVIRA",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.35,
    "use_markers": true,
    "markers_per_cell": 4,
    "marker_advection_scheme": "CVI_RK2",
    "marker_curvature": true,
    "marker_normals_plic": true,
    "marker_max_per_cell": 8,
    "write_markers_vtk": true
  },

  "solver": {
    "type": "MONOLITHIC",
    "linear_solver": "BICGSTAB2",
    "precond_type": "BLOCK_TRIANGULAR",
    "tol": 1.0e-10,
    "max_iter": 1000,
    "r_al": 500.0,
    "adaptative_al": false,
    "k_aal": 100.0,
    "uzawa_tol": 1.0e-10,
    "uzawa_iter": 1,
    "l0_ref": 0.0,
    "t0_ref": 0.0,
    "u0_ref": 0.0,
    "p0_ref": 0.0,
    "use_precond": true
  },

  "obstacles": [],

  "output": {
    "output_dir": "output/bulle_cvi",
    "case_name": "rising_bubble",
    "vtk_interval": 50,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```
# Exécution, Post-Traitement et Dépannage

## Lancement d'une Simulation en Ligne de Commande

Les exécutables **C++17** (`covof2d_cpp` et son lien `covof2d`) et **Fortran 2018** (`covof2d_fortran`) acceptent exactement la même interface en ligne de commande et fonctionnent en **MPI pur** (aucun thread OpenMP).

### Syntaxe d'Exécution

```bash
# Exécution séquentielle (1 processus MPI) :
./COVOF_2D-C++/bin/covof2d_cpp chemin/vers/mon_cas.json
./COVOF_2D-FORTRAN/bin/covof2d_fortran chemin/vers/mon_cas.json

# Exécution parallèle distribuée MPI (exemple sur 4 rangs MPI) :
mpirun -np 4 ./COVOF_2D-C++/bin/covof2d_cpp chemin/vers/mon_cas.json
mpirun -np 4 ./COVOF_2D-FORTRAN/bin/covof2d_fortran chemin/vers/mon_cas.json

# Affichage de l'aide intégrée :
./COVOF_2D-C++/bin/covof2d_cpp --help
./COVOF_2D-FORTRAN/bin/covof2d_fortran --help
```

Si aucun argument n'est fourni sur la ligne de commande, les deux programmes tentent de charger le fichier par défaut `cases/cavity.json`. En cas d'échec d'ouverture du fichier JSON, le rang 0 affiche un message d'erreur et déclenche `MPI_Abort(MPI_COMM_WORLD, 1)`.

---

## Interprétation des Journaux Console (`main.cpp` / `Config::print()` / `Simulation::run()`)

Dès le lancement, seul le processus maître (`rank == 0`) affiche l'en-tête d'accueil, le résumé complet de la configuration (`Config::print()` en C++), puis l'en-tête de la boucle temporelle.

### Exemple Réel de Sortie Console (C++17)

```text
================================================================================
   COVOF 2D — Solveur CFD Navier-Stokes Incompressible Parallèle (C++17/MPI)   
================================================================================
  Nombre de processus MPI  : 4
  Fichier de configuration : rising_bubble.json
================================================================================
========================================================
       COVOF 2D PARALLELE - CONFIGURATION DU CAS         
========================================================
 Mode calcul : NAVIER_STOKES
 Maillage : 64 x 128 cellules (Lx = 1 m, Ly = 2 m)
 dx = 0.015625 m, dy = 0.015625 m
 Physique : Diphasique | Convection : Oui | Tension de surface : Oui
 Fluide 1 : rho=1000 kg/m3, mu=10 Pa.s
 Fluide 2 : rho=100 kg/m3, mu=1 Pa.s
 Tension superficielle sigma = 24.5 N/m
 VOF : Advection=WEYMOUTH_YUE | Reconstruction=ELVIRA | Marqueurs CVI=CVI_RK2 (N_init=4, Courbure=Oui, NormalesPLIC=Oui)
 Gravite : g = (0, -0.98) m/s^2
   Schema temporel : GEAR | Linearisation convective : ADAMS_BASHFORTH | Lissage SVOF : Oui (k=0.500000)
 Temps : dt = 0.001 s, t_final = 3 s
 Solveur : MONOLITHIC | BiCGSTAB(2), tol=1e-10, max_iter=1000, r_al=500, uzawa_iter=1
 Sorties : Dossier 'output/bulle_cvi', intervalle VTK=50
========================================================

>>> DEMARRAGE DU CALCUL COVOF 2D PARALLELE <<<
Pas     Temps [s]     dt [s]      Iter Krylov Vitesse Max     Div Max         
-------------------------------------------------------------------------------------
10      0.01000       1.000e-03   14          1.2410e-02      3.1842e-11      
20      0.02000       1.000e-03   15          2.4635e-02      2.9410e-11      
...
------------------------------------------------------------------------------
>>> SIMULATION COVOF 2D TERMINEE AVEC SUCCES A t = 3.00000 s <<<
Fichiers de visualisation VTK et CSV disponibles dans : output/bulle_cvi
```

*(Note : en Fortran, la bannière de démarrage de boucle affiche `>>> DEMARRAGE DU CALCUL COVOF 2D FORTRAN 2018 PARALLELE <<<` et inclut systématiquement la colonne `Nu Moyen`, qui vaut `0.0000` lorsque `enable_energy = false` ; en C++, la colonne `Nu Moyen` n'est affichée que si `enable_energy = true`.)*

### Signification des Colonnes de la Boucle Temporelle

| Colonne | Signification Physique et Numérique |
| :--- | :--- |
| `Pas` | Indice de l'itération temporelle courante (`step`). |
| `Temps [s]` | Temps physique cumulé $t = \sum \Delta t$ $[\text{s}]$. |
| `dt [s]` | Pas de temps effectif $\Delta t$ utilisé à l'itération courante $[\text{s}]$. |
| `Iter Krylov` | Nombre cumulé d'itérations du solveur linéaire (`total_iters`) sur toutes les itérations d'Uzawa du pas. |
| `Vitesse Max` | Norme maximale de vitesse observée sur la grille globale (`global_max_u`, $[\text{m}\cdot\text{s}^{-1}]$). |
| `Div Max` | Divergence discrète maximale $\max_{i,j} \left|\frac{u_{i+1,j}-u_{i,j}}{\Delta x} + \frac{v_{i,j+1}-v_{i,j}}{\Delta y}\right|$ $[\text{s}^{-1}]$. |
| `Nu Moyen` | Nombre de Nusselt moyen à la paroi chaude $x=0$ (si `enable_energy = true`, ou toujours affiché en Fortran). |

---

## Fichiers de Sortie, Visualisation ParaView et Post-Traitement

Tous les fichiers générés sont placés dans le répertoire `output.output_dir` en utilisant le préfixe `output.case_name` (noté `<case>` ci-dessous).

### 1. Historique Temporel CSV : `<case>_history.csv`

Écrit à chaque pas de temps lorsque `output.write_csv_history = true` (`VTKWriter.cpp` / `module_vtk_writer.f90`).

- **Format C++ (9 colonnes) :**
  `step,time,dt,center_of_mass_y,rise_velocity,bubble_volume,kinetic_energy,max_u,max_div`
- **Format Fortran (11 colonnes) :**
  `step,time,dt,center_of_mass_y,rise_velocity,bubble_volume,kinetic_energy,max_u,max_div,avg_rho_gas,avg_p_gas`

Définition exacte des métriques calculées sur l'indicateur de la phase 2 (`1.0 - C`) :
- `center_of_mass_y` : ordonnée du barycentre de la phase 2, $y_c = \frac{\int y (1-C)\,d\Omega}{\int (1-C)\,d\Omega}$.
- `rise_velocity` : vitesse verticale moyenne de la phase 2, $v_c = \frac{\int v (1-C)\,d\Omega}{\int (1-C)\,d\Omega}$.
- `bubble_volume` : aire totale (volume 2D) de la phase 2, $\int (1-C)\,d\Omega$.
- `kinetic_energy` : énergie cinétique discrète aux nœuds de pression, $\frac{1}{2}\sum_{i,j}\rho_{i,j}\left(\frac{u_{i,j}^2 + u_{i+1,j}^2}{2} + \frac{v_{i,j}^2 + v_{i,j+1}^2}{2}\right)\Delta x \Delta y$.
- `max_u`, `max_div` : vitesse maximale au centre des mailles et divergence discrète maximale.
- `avg_rho_gas`, `avg_p_gas` (Fortran) : masse volumique gaz et pression moyenne dans la phase $1-C$ (pour suivi compressible).

### 2. Champs Volumiques Parallèles XML VTK : `.vtr` et `.pvtr`

Écrits au pas `0`, tous les `output.vtk_interval` pas, et au pas final :

1. **Pièces locales par rang MPI :** `<case>_p<rank>_<step:06d>.vtr` (ex. `rising_bubble_p0_000050.vtr`). Chaque processus écrit son propre sous-domaine `RectilinearGrid` avec ses coordonnées `XCoordinates`, `YCoordinates`, `ZCoordinates`.
2. **Fichier maître parallèle :** `<case>_<step:06d>.pvtr` (ex. `rising_bubble_000050.pvtr`), écrit par le rang 0 et référençant l'ensemble des `<Piece Extent="..." Source="<case>_p<rank>_<step:06d>.vtr"/>`.

Les tableaux `<CellData>` exportés dans chaque `.vtr` et `.pvtr` sont exactement :
- `PhaseIndicator` : fraction volumique brute $C_{i,j}$
- `PhaseIndicator_Smoothed` : fraction volumique lissée SVOF $C_{s,i,j}$ (présente si `C_s` est alloué/non vide)
- `Pressure` : pression $p_{i,j}$
- `Density` : masse volumique $\rho_{i,j}$
- `Viscosity` : viscosité dynamique $\mu_{i,j}$
- `Divergence` : divergence discrète $(\nabla_h \cdot \mathbf{u})_{i,j}$
- `Temperature` : température $T_{i,j}$
- `Curvature` : courbure locale $\kappa_{i,j}$
- `Velocity` : vecteur vitesse $(u_c, v_c, 0.0)$ à 3 composantes interpolé aux centres de mailles

### 3. Nuage de Marqueurs Lagrangiens CVI : `<case>_markers_<step:06d>.vtk`

Lorsque `vof.use_markers = true` et `vof.write_markers_vtk = true`, les rangs MPI rassemblent par `MPI_Gatherv` vers le rang 0 l'ensemble des marqueurs actifs et écrivent un fichier VTK Legacy ASCII `DATASET POLYDATA` (`InterfaceMarkers::write_vtk` / `write_markers_vtk`) :
- `POINTS <N> double` : coordonnées $(x_p, y_p, 0.0)$
- `VERTICES <N> <2N>` : sommets ponctuels pour affichage direct dans ParaView
- `POINT_DATA <N>` contenant :
  - `VECTORS MarkerNormal double` : normale unitaire $(n_{x,p}, n_{y,p}, 0.0)$ orientée vers $C=1$
  - `SCALARS MarkerCurvature double 1` : courbure locale $\kappa_p$ portée par le marqueur

### 4. Sorties Analytiques Spécifiques aux Modes de Vérification

Selon la configuration (`mode` ou `test_case`), les routines dédiées de `Simulation.cpp` et `module_simulation.f90` exportent automatiquement des fichiers CSV de métrologie exacte dans `output_dir` :

| Cas / Mode | Fichiers Générés | Colonnes Principales |
| :--- | :--- | :--- |
| `mode = "LAPLACIAN"` (`CAS_3_1`, `CAS_3_2_NEUMANN`) | `laplacian_error.csv`, `solver_history.csv`, `solution_field_2d.csv`, `profile_x_mid.csv`, `profile_y_mid.csv`, `solution_field.vtk` | `nx,ny,dx,dy,iters,err_l1,err_l2,err_linf` et `iter,residual` |
| `test_case = "CAS_3_7"` / `"POISEUILLE"` | `poiseuille_error.csv`, `profile_u_y.csv`, `profile_p_x.csv`, `solution_field_2d.csv`, `solution_field.vtk` | `nx,ny,dx,dy,iters,err_l1,err_l2,err_linf,err_grad_p,cpu_time` |
| `test_case = "CAS_3_8"` / `"CAS_3_8_POISEUILLE_UNSTEADY"` | `poiseuille_unsteady_error.csv`, `profile_u_y.csv`, `solution_field_2d.csv`, `solution_field.vtk` | `nx,ny,dx,dy,total_step,final_time,err_l1,err_l2,err_linf,cpu_time` |
| `mode = "VORTEX_ADVECTION"` | `vortex_error.csv` (+ `.vtr`/`.pvtr` et `_markers_*.vtk`) | `nx,ny,dt,T,advection,reconstruction,limiter,err_l1,err_l2,err_linf,mass_init,mass_final,mass_conservation` |

### Exploitation sous ParaView

1. **Champs volumiques :** Dans ParaView, ouvrir le groupe temporel `<case>_*.pvtr` et cliquer sur **Apply**.
2. **Iso-contour d'interface VOF :** Appliquer le filtre **Cell Data to Point Data**, puis **Contour** sur `PhaseIndicator` avec la valeur `0.5`.
3. **Visualisation des marqueurs CVI :** Ouvrir simultanément le groupe `<case>_markers_*.vtk`, appliquer une taille de point visible (*Point Size = 4*) ou le filtre **Glyph** orienté selon `MarkerNormal` et coloré par `MarkerCurvature`.
4. **Vecteurs et lignes de courant :** Utiliser **Glyph** ou **Stream Tracer** sur le champ vectoriel `Velocity`.

### Scripts Python de Post-Traitement (`pyvista`, `numpy`, `matplotlib`)

#### Script 1 : Extraction de Profils et Superposition des Marqueurs depuis `.pvtr` et `_markers_*.vtk`

```python
#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

def plot_fields_and_markers(pvtr_file, markers_vtk=None):
    data = pv.get_reader(pvtr_file).read()
    bounds = data.bounds  # [xmin, xmax, ymin, ymax, zmin, zmax]
    x_mid = 0.5 * (bounds[0] + bounds[1])

    # 1. Extraction d'un profil vertical au centre X = Lx / 2
    line = data.sample_over_line([x_mid, bounds[2], 0.0],
                                 [x_mid, bounds[3], 0.0], resolution=400)
    y = line.points[:, 1]
    u = line["Velocity"][:, 0]
    p = line["Pressure"]

    fig, axes = plt.subplots(1, 2, figsize=(11, 5), dpi=200)
    axes[0].plot(u, y, "b-", lw=2, label="u(x_mid, y)")
    axes[0].set_xlabel("Vitesse horizontale u [m/s]")
    axes[0].set_ylabel("Ordonnée y [m]")
    axes[0].grid(True, ls="--", alpha=0.6)
    axes[0].legend()

    axes[1].plot(p, y, "r-", lw=2, label="p(x_mid, y)")
    axes[1].set_xlabel("Pression p [Pa]")
    axes[1].set_ylabel("Ordonnée y [m]")
    axes[1].grid(True, ls="--", alpha=0.6)
    axes[1].legend()

    if markers_vtk is not None:
        mk = pv.read(markers_vtk)
        print(f"Nombre de marqueurs CVI actifs : {mk.n_points}")
        print(f"Courbure moyenne |kappa|       : {np.mean(np.abs(mk['MarkerCurvature'])):.4f} m^-1")

    plt.tight_layout()
    plt.savefig("covof_postproc_profiles.png")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 postproc.py <case_000050.pvtr> [case_markers_000050.vtk]")
        sys.exit(1)
    mk_arg = sys.argv[2] if len(sys.argv) > 2 else None
    plot_fields_and_markers(sys.argv[1], mk_arg)
```

#### Script 2 : Suivi de l'Historique CSV (`<case>_history.csv`)

```python
#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

csv_path = sys.argv[1]
data = np.genfromtxt(csv_path, delimiter=",", names=True)

fig, ax = plt.subplots(1, 2, figsize=(10, 4), dpi=200)
ax[0].plot(data["time"], data["rise_velocity"], "k-", lw=1.8)
ax[0].set_xlabel("Temps [s]")
ax[0].set_ylabel("Vitesse d'ascension [m/s]")
ax[0].grid(True, ls="--", alpha=0.6)

ax[1].semilogy(data["time"], np.maximum(data["max_div"], 1e-16), "m-", lw=1.8)
ax[1].set_xlabel("Temps [s]")
ax[1].set_ylabel("Divergence max [s^-1]")
ax[1].grid(True, ls="--", alpha=0.6)

plt.tight_layout()
plt.savefig("history_diagnostics.png")
```

---

## Guide Complet de Dépannage (*Troubleshooting*)

| Symptôme Observé | Cause Physique / Numérique Probable | Action Corrective Recommandée |
| :--- | :--- | :--- |
| **Divergence de Krylov (`Iter Krylov` = `max_iter`, `NaN` dans les logs)** | 1. Pas de temps $\Delta t$ violant la CFL convective ou capillaire ($\Delta t_\sigma$).<br>2. Saut $\rho_1/\rho_2$ ou $\mu_1/\mu_2$ brutal sans lissage SVOF.<br>3. Préconditionneur inadapté au mode `MONOLITHIC`. | 1. Activer `time.adaptive_dt = true` avec `time.cfl_max = 0.35`.<br>2. Activer `physics.enable_svof_smoothing = true` (`svof_filter_width = 0.5`).<br>3. Avec `solver.type = "MONOLITHIC"`, utiliser impérativement `solver.precond_type = "BLOCK_TRIANGULAR"`. |
| **Divergence maximale `Div Max` élevée ($> 10^{-4}$)** | Trop peu d'itérations d'Uzawa en mode AL/AAL ou tolérance linéaire trop relâchée. | Passer `solver.type = "MONOLITHIC"` (ou `"AAL"` avec `uzawa_iter = 5`, `uzawa_tol = 1.0e-8`) et abaisser `solver.tol` à `1.0e-9` ou `1.0e-10`. |
| **Courants parasites intenses autour d'une goutte statique** | Courbure CSF bruitée par différentiation finie de $C$. | Activer la méthode hybride CVI-Marqueurs (`vof.use_markers = true`, `vof.marker_curvature = true`, `vof.marker_advection_scheme = "CVI_RK2"`) et la reconstruction `"ELVIRA"`. |
| **Fragmentation d'un filament mince ou d'une jupe de bulle** | Sous-résolution sous-maille ($\delta < \Delta x$) ou perte de cohérence de normale. | Activer `vof.use_markers = true`, `vof.marker_normals_plic = true` et `vof.advection_method = "WEYMOUTH_YUE"` avec `cfl_vof <= 0.35`. |
| **Oscillations spatiales en damier dans les zones à fort gradient** | Schéma convectif `CENTERED` utilisé alors que $Re_\Delta > 2$. | Basculer `physics.convective_scheme = "UPWIND"` ou `"TVD_SUPERBEE"`, ou raffiner `(mesh.nx, mesh.ny)`. |
| **Dépassement ou blocage MPI lors du lancement** | Nombre de processus MPI incompatible avec `nx`, `ny`, ou rangs bloqués sur barrière suite à une erreur JSON. | Vérifier que $N_x / P_x \ge 2$ et $N_y / P_y \ge 2$, et contrôler le JSON (clés exactes, pas de virgule traînante en fin de bloc). |
| **Clé JSON ignorée silencieusement** | Utilisation d'un ancien nom de clé (ex. `physics.rho1`, `time.t_end`, `thermal.*`, `convection.*`, `NO_SLIP`). | Consulter le Chapitre 6 : utiliser `physics.fluid1.rho`, `time.t_final`, `physics.enable_energy`, `physics.convective_*`, `DIRICHLET`, `FREE_SLIP`, `BOX`. |
