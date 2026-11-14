# Tutoriel 7 : Advection d'une Tâche Circulaire dans un Tourbillon Cisaillé

Ce tutoriel isole le transport cinématique pur de l'interface fluide VOF sans résolution hydrodynamique. Une tâche circulaire est soumise à un champ de vitesse cisaillé instationnaire réversible qui l'étire en une fine nappe spirale avant d'inverser le flux à t = T/2 pour reformer le disque initial à t = T.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les fondements physiques et les équations aux dérivées partielles gouvernant l'écoulement.
- Configurer les paramètres numériques et conditions aux limites sur grille cartésienne décalée MAC.
- Exécuter le cas test avec les exécutables C++17 et Fortran 2018 (séquentiel et OpenMP).
- Quantifier avec rigueur la précision par rapport à la solution analytique exacte ($L_2$, $L_\infty$).
- Produire et exploiter les visualisations scientifiques volumiques sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 :** `Mesh2D`, `FlowFields`, `WeymouthYueVOF`, `InterfaceReconstruction`, `VTKWriter`.
- **Moteur Fortran 2018 :** `module_mesh`, `module_flow_fields`, `module_weymouth_yue`, `module_plic_reconstruction`, `module_vtk_writer`.

---

## 2. Formulation Physique et Solution Analytique Exacte

Équations Physiques et Configuration

Ce cas cinématique valide les schémas de transport et de reconstruction d'interface Volume-of-Fluid (VOF) sous étirement tourbillonnaire sévère et réversibilité temporelle exacte.

Dans la cavité carrée unitaire $\Omega = [0, 1] \times [0, 1]\,\mathrm{m}^2$, un disque fluide initial ($C = 1$) de rayon $R_0 = 0.125\,\mathrm{m}$ est centré en $(x_c, y_c) = (0.5, 0.75)\,\mathrm{m}$ à $t = 0$. L'indicatrice de phase $C(\mathbf{x}, t) \in [0, 1]$ obéit à l'équation d'advection pure sans terme source :
$$\frac{\partial C}{\partial t} + \nabla \cdot (C \mathbf{u}) = 0.$$
Le champ de vitesse solénoïdal ($\nabla \cdot \mathbf{u} \equiv 0$) dérive de la fonction de courant non stationnaire modulée par un cosinus temporel de période $T = 6.0\,\mathrm{s}$ :
$$\Psi(x, y, t) = \frac{1}{\pi} \sin^2(\pi x) \sin^2(\pi y) \cos\left(\frac{\pi t}{T}\right),$$
soit aux faces des mailles décalées MAC :
$$u(x, y, t) = -\frac{\partial \Psi}{\partial y} = -\sin^2(\pi x) \sin(2\pi y) \cos\left(\frac{\pi t}{T}\right),$$
$$v(x, y, t) = +\frac{\partial \Psi}{\partial x} = +\sin(2\pi x) \sin^2(\pi y) \cos\left(\frac{\pi t}{T}\right).$$
De $t = 0$ à $t = T/2 = 3.0\,\mathrm{s}$, le disque est étiré en une fine spirale tourbillonnaire ; à $t = 3.0\,\mathrm{s}$, le champ de vitesse change de signe et réenroule le filament jusqu'à $t = T = 6.0\,\mathrm{s}$ où la solution exacte coïncide identiquement avec l'état initial :
$$C_{\mathrm{ex}}(\mathbf{x}, T) \equiv C(\mathbf{x}, 0).$$

---

---

## 3. Paramètres Physiques, Numériques et Discrétisation

Physiques, Numériques et Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et Disque initial | $\Omega$, $(x_c, y_c)$, $R_0$ | $[0,1]\times[0,1]\,\mathrm{m}^2$, $(x_c, y_c) = (0.5, 0.75)\,\mathrm{m}$, $R_0 = 0.125\,\mathrm{m}$ |
| Période cinématique | $T = t_{\mathrm{final}}$ | $6.0\,\mathrm{s}$ (inversion à $T/2 = 3.0\,\mathrm{s}$) |
| Séquence de grilles et $\Delta t$ | $(N\times N, \Delta t)$ ($\mathrm{CFL} \le 0.32$) | $(64^2, 0.005\,\mathrm{s})$, $(128^2, 0.0025\,\mathrm{s})$, $(256^2, 0.00125\,\mathrm{s})$, $(512^2, 0.0005\,\mathrm{s})$ |
| Méthodes d'advection VOF | `vof.advection_method` | 1) `"WEYMOUTH_YUE"`, 2) `"CIAM"`, 3) `"VOF_TVD"` |
| Reconstructions et limiteurs | `reconstruction_method` / `tvd_limiter` | `"PLIC"`, `"ELVIRA"`, `"LVIRA"`, `"SLIC"` ; `"SUPERBEE"`, `"VAN_LEER"`, `"MUSCL"`, `"HRIC"` |
| Couplage Solveurs (Phase 2) | `solver.type` & `precond_type` | `AL+DIRECT`, `AL+BICGSTAB2+BT`, `AAL+BICGSTAB2+BT`, `MONOLITHIC+BICGSTAB2+BT` |

---

---

## 4. Fichier de Configuration JSON Pas à Pas

de Fichiers `.json` (C++17 et Fortran 2018)

### Prototype `case/cpp/cas_7_weymouth_yue_plic_n128.json`
```json
{
  "mode": "VORTEX_ADVECTION",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 1.0,
    "ly": 1.0
  },
  "physics": {
    "enable_two_phase": true,
    "enable_surface_tension": false,
    "enable_convection": true,
    "sigma": 0.0,
    "fluid1": { "rho": 1000.0, "mu": 0.001 },
    "fluid2": { "rho": 1.0,    "mu": 1.0e-05 }
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.5
  },
  "initial_condition": {
    "type": "VORTEX_SHEAR",
    "center_x": 0.5,
    "center_y": 0.75,
    "radius": 0.125
  },
  "time": {
    "dt": 0.0025,
    "t_final": 6.0,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "right":  { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "linear_solver": "DIRECT",
    "tol": 1e-12,
    "max_iter": 500,
    "r_al": 100.0,
    "uzawa_iter": 1,
    "use_precond": false,
    "precond_type": "NONE"
  },
  "output": {
    "output_dir": "VALIDATION_TOTALE/validation_cas_7_vortex_shear/output/cpp/weymouth_yue_plic/n128",
    "case_name": "cas_7_weymouth_yue_plic_n128",
    "vtk_interval": 240,
    "log_interval": 200,
    "write_csv_history": true
  }
}

```

---

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

L'exécution s'effectue indifféremment à l'aide du binaire C++17 (`covof2d_cpp`) ou du binaire Fortran 2018 (`covof2d_fortran`), garantissant une stricte équivalence des résultats numériques.

```bash
# 1. Se placer dans le répertoire du cas
cd VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_7_vortex_shear

# 2. Exécution avec le moteur C++17 (Séquentiel)
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_7_config.json

# 3. Exécution avec le moteur C++17 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_7_config.json

# 4. Exécution avec le moteur Fortran 2018 (Séquentiel)
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_7_config.json

# 5. Exécution avec le moteur Fortran 2018 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_7_config.json

```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification et Tolérances
à Mener (Précision, Intégrité et Performance)

1. **Erreur géométrique de retour $L_1$ à $t = T = 6.0\,\mathrm{s}$ et ordre de convergence :**
   $$E_{L_1} = \frac{1}{N^2}\sum_{i=1}^{N}\sum_{j=1}^{N} |C_{i,j}(T) - C_{i,j}(0)|.$$
   Vérifier que `WEYMOUTH_YUE + PLIC` décroît de $7.35\times 10^{-3}$ ($N=64$) à $2.32\times 10^{-3}$ ($N=128$, $p=1.66$), $9.42\times 10^{-4}$ ($N=256$, $p=1.30$) et $4.29\times 10^{-4}$ ($N=512$, $p=1.14$), sans aucune saturation sur grille fine.
2. **Conservation de la masse, bornitude et exploration des options VOF du code :**
   - Vérifier que $\frac{\sum C_{i,j}(t)}{\sum C_{i,j}(0)} = 100.0000\,\%$ pour `WEYMOUTH_YUE` et `VOF_TVD` (et $> 98.9\,\%$ sur grille fine pour `CIAM`), et que $0 \le C_{i,j}(t) \le 1$ en tout point et à tout instant.
   - Vérifier la règle d'exclusion mutuelle implémentée dans `Config.cpp` / `module_config.f90` : lorsque `advection_method = "VOF_TVD"`, `reconstruction_method` est automatiquement forcé à `"NONE"` et le limiteur par défaut à `"SUPERBEE"` ; lorsque `advection_method` est géométrique (`"WEYMOUTH_YUE"` ou `"CIAM"`), `tvd_limiter` est forcé à `"NONE"`.
   - Comparer la netteté de l'interface et l'erreur $E_{L_1}$ entre les limiteurs algébriques `SUPERBEE`, `VAN_LEER`, `MUSCL` et `HRIC`, ainsi qu'entre les reconstructions `PLIC`, `ELVIRA`, `LVIRA` et `SLIC`.
3. **Invariance vis-à-vis du solveur de couplage et parité C++17 / Fortran 2018 :**
   - Vérifier en Phase 2 que les erreurs $E_{L_1}, E_{L_2}, E_{L_\infty}$ sont rigoureusement identiques entre `AL+DIRECT`, `AL+BICGSTAB2+BLOCK_TRIANGULAR`, `AAL+BICGSTAB2+BLOCK_TRIANGULAR` et `MONOLITHIC+BICGSTAB2+BLOCK_TRIANGULAR`.
   - Vérifier la parité machine C++17 vs Fortran 2018 : $|E_{L_1}^{\mathrm{C++}} - E_{L_1}^{\mathrm{Fortran}}| < 10^{-14}$.

---

### 6.2 Résultats Attendus et Visualisation ParaView

