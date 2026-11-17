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
