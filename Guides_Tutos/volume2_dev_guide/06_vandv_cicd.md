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
