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
| **Solveurs & Préconditionneurs** | Solveurs linéaires | `BICGSTAB2`, `BICGSTAB`, `CG`, `GMRES`, `MINRES`, `DIRECT` | Implémentations Krylov *matrix-free* MPI (`MINRES` pour systèmes symétriques, `BICGSTAB2`/`GMRES` pour systèmes convectifs ou bloc-triangulaires) + LU directe LAPACK `dgetrf`/`dgetrs` (passerelles `MUMPS`/`PETSC`/`HYPRE` en C++) |
| | Préconditionneurs | `JACOBI`, `SSOR`, `ILU0`, `BLOCK_TRIANGULAR` (`BT_SBP`), `CUSTOM`, `NONE` | Préconditionnement diagonal enrichi, SSOR, ILU(0) et bloc triangulaire à complément de Schur CSR 100 % parallèle MPI |

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
  - `Run_Valide-CASES/` : 12 cas de validation rapides pilotés par `run_all_validations.sh` (`validation_cas_1_laplacian` à `validation_cas_7_vortex_shear`, puis `validation_cas_8_cavite_entrainee` à `validation_cas_12_convection_poreux_darcy_brinkman`).
  - `Run_CASES/` : 9 campagnes avancées ordonnées par complexité croissante (`validation_cas_1_obstacle_fin` à `validation_cas_9_compression_bulle`), incluant la chute libre du cylindre dense et les cas diphasiques/compressibles lourds, ainsi que le sous-dossier `RUN_FC-WY-CVI-Marker/` dédié aux 8 campagnes hybrides CVI-Marqueurs et au script `verify_cvi_markers.py`.
  - `Theoritical_Description/` : Référentiel théorique consolidé (`Validation_Description.pdf`) et spécifications des cas de `Run_Valide-CASES/`, `Run_CASES/` et `RUN_FC-WY-CVI-Marker/`.
- `Guides_Tutos/` : Suite documentaire complète (Volumes I, II et III, miroir `Run_FC-WY-CVI-Marker/` et `verify_cvi_markers.py`, synthèses théoriques `SYNTHESE_METHODE_FC_WY_CVI_MARQUEURS.pdf` et `FEUILLE_DE_ROUTE_IMPLEMENTATION_COVOF_CVI_MARQUEURS.pdf`).
- `output/` : Répertoire de sortie par défaut généré lors des simulations si `output.output_dir` vaut `"output"`.
