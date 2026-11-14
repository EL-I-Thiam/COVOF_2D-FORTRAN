# Tutoriel 8 : Chute Libre d'un Cylindre Dense dans le Vide (Pénalisation de Brinkman)

Ce tutoriel valide le couplage fluide-solide par pénalisation de Brinkman avec un contraste de masse volumique extrême (rapport rho_s / rho_f = 85 000). Il démontre que l'inclusion solide suit rigoureusement la cinématique galiléenne de chute libre dans le vide sans dispersion ni amortissement numérique parasite.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les fondements physiques et les équations aux dérivées partielles gouvernant l'écoulement.
- Configurer les paramètres numériques et conditions aux limites sur grille cartésienne décalée MAC.
- Exécuter le cas test avec les exécutables C++17 et Fortran 2018 (séquentiel et OpenMP).
- Quantifier avec rigueur la précision par rapport à la solution analytique exacte ($L_2$, $L_\infty$).
- Produire et exploiter les visualisations scientifiques volumiques sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 :** `Mesh2D`, `FlowFields`, `BrinkmanPenalization`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `SSORPreconditioner`.
- **Moteur Fortran 2018 :** `module_mesh`, `module_flow_fields`, `module_brinkman`, `module_navier_stokes`, `module_linear_solvers`, `module_preconditioners`.

---

## 2. Formulation Physique et Solution Analytique Exacte

Équations Physiques et Configuration

Ce cas diphasique sévère simule la chute libre sous gravité d'un cylindre circulaire ultra-dense et ultra-visqueux dans l'air au repos. Il met à l'épreuve la robustesse des solveurs de couplage vitesse-pression et des préconditionneurs face à des sauts extrêmes de masse volumique $r_\rho = \rho_2/\rho_1 \approx 84\,976$ et de viscosité dynamique $r_\mu = \mu_2/\mu_1 \approx 5.405\times 10^9$.

Dans une cavité rectangulaire $\Omega = [0, 0.1] \times [0, 0.2]\,\mathrm{m}^2$ soumise à la pesanteur $\mathbf{g} = (0, -9.81)\,\mathrm{m/s}^2$, un cylindre de rayon $R = 0.0125\,\mathrm{m}$ est initialement centré en $(x_c, Y_0) = (0.05, 0.15)\,\mathrm{m}$ sans vitesse initiale :
- **Fluide 1 (Air ambiant, $C=0$)** : $\rho_1 = 1.1768\,\mathrm{kg/m}^3$, $\mu_1 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$.
- **Fluide 2 (Cylindre dense, $C=1$)** : $\rho_2 = 100\,000\,\mathrm{kg/m}^3$, $\mu_2 = 100\,000\,\mathrm{Pa\cdot s}$, $\sigma = 0$.

Compte tenu du rapport d'inertie $r_\rho \gg 1$, la traînée aérodynamique et la poussée d'Archimède exercées par l'air sont négligeables devant le poids du cylindre. La vitesse verticale $U_c(t)$ et l'ordonnée du centre de gravité $Y_c(t)$ suivent la loi de chute libre galiléenne exacte dans le vide :
$$U_c^{\mathrm{ex}}(t) = -g t, \qquad Y_c^{\mathrm{ex}}(t) = Y_0 - \frac{1}{2} g t^2.$$
À l'instant final $t_{\mathrm{final}} = 0.144\,\mathrm{s}$ :
$$U_c^{\mathrm{ex}}(0.144\,\mathrm{s}) = -1.41264\,\mathrm{m/s}, \qquad Y_c^{\mathrm{ex}}(0.144\,\mathrm{s}) = 0.048288\,\mathrm{m}.$$

### Formulations numériques spécifiques mises en œuvre
1. **Schéma temporel Gear / BDF2 d'ordre 2** : $\frac{1.5\rho^{n+1}\mathbf{u}^{n+1} - 2\rho^n\mathbf{u}^n + 0.5\rho^{n-1}\mathbf{u}^{n-1}}{\Delta t}$ (`"time_scheme": "GEAR"`).
2. **Linéarisation convective d'Adams-Bashforth** : vitesse advectrice extrapolée $\mathbf{u}^* = 2\mathbf{u}^n - \mathbf{u}^{n-1}$ (`"convective_linearization": "ADAMS_BASHFORTH"`).
3. **Lissage SVOF de Helmholtz** : $(I - D\Delta)C_s = C$ avec $D = k\Delta x\Delta y$ ($k=1.0$) et moyenne harmonique nodale de la viscosité de cisaillement $\mu_{\mathrm{node}} = \frac{\mu_1\mu_2}{\mu_2\widetilde{C}_{s,\mathrm{node}} + \mu_1(1-\widetilde{C}_{s,\mathrm{node}})}$.
4. **Quatre méthodes de couplage vitesse-pression comparées** :
   - `MONOLITHIC` (**FC**, $r_{\mathrm{al}} = 0$) : résolution couplée directe du système selle par BiCGSTAB(2) + préconditionneur `BLOCK_TRIANGULAR`.
   - `AUGMENTED_LAGRANGIAN` (**AL**, $r_{\mathrm{al}} = 10\,000$ fixe).
   - `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (**AAL**) : paramètre local par cellule $r(t,M)$ adapté selon les gradients d'interface et la divergence locale avec facteur $K_{\mathrm{aal}}$ (`k_aal`) et critère d'arrêt d'Uzawa `uzawa_tol`.
   - `STANDARD_PROJECTION` (**SP**) : méthode de projection à pas fractionnaire.

---

---

## 3. Paramètres Physiques, Numériques et Discrétisation

Physiques, Numériques et Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et Cylindre | $\Omega$, $(x_c, Y_0)$, $R$ | $[0, 0.1]\times[0, 0.2]\,\mathrm{m}^2$, $(0.05, 0.15)\,\mathrm{m}$, $R = 0.0125\,\mathrm{m}$ |
| Maillages et Pas de temps | $(N_x \times N_y)$, $\Delta t$, $t_{\mathrm{final}}$ | $100\times 200$ et $200\times 400$, $\Delta t = 6.25\times 10^{-5}\,\mathrm{s}$, $t_{\mathrm{final}} = 0.144\,\mathrm{s}$ (2304 pas) |
| Conditions aux limites | Gauche/Droite ; Bas/Haut | `FREE_SLIP` sur $x=0, 0.1$ ; `DIRICHLET` ($u=v=0$) sur $y=0, 0.2$ |
| Solveurs et Préconditionneurs | Couplage, Krylov, Précond. | `MONOLITHIC` (FC), `AUGMENTED_LAGRANGIAN` (AL), `AAL` (`k_aal` $\in \{10, 100, 1000\}$), `STANDARD_PROJECTION` (SP) ; `BICGSTAB2`, `GMRES`, `BICGSTAB` avec `BLOCK_TRIANGULAR`, `JACOBI`, `SSOR`, `ILU0` |

---

---

## 4. Fichier de Configuration JSON Pas à Pas

de Fichiers `.json` (C++17 et Fortran 2018)

### Prototype `case/cpp/cas_8_fc_100x200.json`
```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 100,
    "ny": 200,
    "lx": 0.1,
    "ly": 0.2
  },
  "physics": {
    "enable_two_phase": true,
    "enable_surface_tension": false,
    "enable_convection": true,
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_svof_smoothing": true,
    "svof_filter_width": 1.0,
    "enable_diffusion_smoothing": false,
    "sigma": 0.0,
    "gx": 0.0,
    "gy": -9.81,
    "fluid1": { "rho": 1.1768,   "mu": 1.85e-05 },
    "fluid2": { "rho": 100000.0, "mu": 100000.0 }
  },
  "initial_condition": {
    "type": "CYLINDER",
    "center_x": 0.05,
    "center_y": 0.15,
    "radius": 0.0125
  },
  "time": {
    "time_scheme": "GEAR",
    "dt": 6.25e-05,
    "t_final": 0.144,
    "cfl_max": 0.5,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "MONOLITHIC",
    "linear_solver": "BICGSTAB2",
    "tol": 1e-06,
    "max_iter": 120,
    "r_al": 0.0,
    "uzawa_iter": 1,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.5
  },
  "output": {
    "output_dir": "VALIDATION_TOTALE/validation_cas_8_dense_cylinder/output/cpp/fc_100x200",
    "case_name": "cas_8_fc_100x200",
    "vtk_interval": 100,
    "log_interval": 100,
    "write_csv_history": true
  }
}

```

Pour activer le Lagrangien Augmenté Adaptatif (**AAL**), utiliser le bloc `"solver"` suivant :
```json
  "solver": {
    "type": "ADAPTATIVE_AUGMENTED_LAGRANGIAN",
    "linear_solver": "BICGSTAB2",
    "tol": 1e-06,
    "max_iter": 200,
    "r_al": 0.0,
    "adaptative_al": true,
    "k_aal": 100.0,
    "uzawa_iter": 5,
    "uzawa_tol": 1e-06,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  }

```

---

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

L'exécution s'effectue indifféremment à l'aide du binaire C++17 (`covof2d_cpp`) ou du binaire Fortran 2018 (`covof2d_fortran`), garantissant une stricte équivalence des résultats numériques.

```bash
# 1. Se placer dans le répertoire du cas
cd VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_8_dense_cylinder

# 2. Exécution avec le moteur C++17 (Séquentiel)
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_8_config.json

# 3. Exécution avec le moteur C++17 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_8_config.json

# 4. Exécution avec le moteur Fortran 2018 (Séquentiel)
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_8_config.json

# 5. Exécution avec le moteur Fortran 2018 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_8_config.json

```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification et Tolérances
à Mener (Précision, Intégrité et Performance)

1. **Erreurs relatives $L_2$ temporelles sur la vitesse et la position :**
   $$\mathrm{Err}_{\mathrm{rel},2}(U_c) = \sqrt{\frac{\sum_{n=1}^{N_t} \left(U_c(t^n) - U_c^{\mathrm{ex}}(t^n)\right)^2}{\sum_{n=1}^{N_t} \left(U_c^{\mathrm{ex}}(t^n)\right)^2}}, \qquad \mathrm{Err}_{\mathrm{rel},2}(Y_c) = \sqrt{\frac{\sum_{n=1}^{N_t} \left(Y_c(t^n) - Y_c^{\mathrm{ex}}(t^n)\right)^2}{\sum_{n=1}^{N_t} \left(Y_c^{\mathrm{ex}}(t^n)\right)^2}}.$$
   Sur grille $100\times 200$, vérifier que `FC` atteint $\mathrm{Err}_{\mathrm{rel},2}(U_c) \approx 8.71\times 10^{-5}$ et $\mathrm{Err}_{\mathrm{rel},2}(Y_c) \approx 1.57\times 10^{-4}$, avec une variation relative de volume $\Delta V/V_0 < 0.001\,\%$.
2. **Intégrité des méthodes de couplage et du filtre SVOF face aux forts contrastes ($r_\rho \approx 85\,000, r_\mu \approx 5.4\times 10^9$) :**
   - Comparer `MONOLITHIC` (FC), `AUGMENTED_LAGRANGIAN` ($r_{\mathrm{al}} = 10^4$), `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (`AAL` avec $K_{\mathrm{aal}} \in \{10, 100, 1000\}$) et `STANDARD_PROJECTION` (SP) sur la préservation de la circularité du cylindre et l'absence de traînée numérique parasite.
   - Vérifier l'impact du filtre elliptique de Helmholtz (`enable_svof_smoothing: true` vs `false`, `svof_filter_width` $k \in \{0.5, 1.0, 2.0\}$) et de la moyenne harmonique nodale de $\mu$ sur l'élimination des pics de vorticité non physiques à l'interface.
3. **Performance des préconditionneurs à saut de densité et parité C++/Fortran :**
   - Vérifier que le complément de Schur à densité variable $a_W = \Delta y^2 / D_u$, $a_E = \Delta y^2 / D_u$, $a_S = \Delta x^2 / D_v$, $a_N = \Delta x^2 / D_v$ dans `BlockTriangularPreconditioner` (utilisant $\alpha_{\mathrm{time}} = 1.5$ en mode `GEAR`) maintient la convergence de `BICGSTAB2` en $< 120$ itérations malgré le saut de densité $1:85\,000$.
   - Comparer le temps CPU total et la divergence résiduelle $\|\nabla\cdot\mathbf{u}\|_\infty$ entre C++17 et Fortran 2018.

---

### 6.2 Résultats Attendus et Visualisation ParaView

