# Tutoriel 3 : Écoulement de Poiseuille Plan Laminaire Stationnaire

Ce tutoriel traite de l'écoulement laminaire cisaillé 1D stationnaire en canal plan entre deux parois immobiles, entraîné par un gradient de pression longitudinal imposé. Il valide la discrétisation de Navier-Stokes incompressible et démontre l'intérêt capital de l'extrapolation quadratique à 3 points dans les mailles fantômes pour atteindre une précision machine sur la vitesse tangentielle.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les fondements physiques et les équations aux dérivées partielles gouvernant l'écoulement.
- Configurer les paramètres numériques et conditions aux limites sur grille cartésienne décalée MAC.
- Exécuter le cas test avec les exécutables C++17 et Fortran 2018 (séquentiel et OpenMP).
- Quantifier avec rigueur la précision par rapport à la solution analytique exacte ($L_2$, $L_\infty$).
- Produire et exploiter les visualisations scientifiques volumiques sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 :** `Mesh2D`, `FlowFields`, `Config`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `SSORPreconditioner`, `VTKWriter`.
- **Moteur Fortran 2018 :** `module_mesh`, `module_flow_fields`, `module_config`, `module_navier_stokes`, `module_linear_solvers`, `module_preconditioners`, `module_vtk_writer`.

---

## 2. Formulation Physique et Solution Analytique Exacte

Équations Physiques et Configuration

Ce cas valide la résolution stationnaire des équations de Navier-Stokes incompressibles en canal plan laminaire et vérifie la propriété d'exactitude polynomiale d'ordre 2 du stencil visqueux couplé à l'extrapolation quadratique à 3 points dans les mailles fantômes de paroi.

Dans un canal plan $\Omega = [0, L] \times [0, H]$ de longueur $L = 1.0\,\mathrm{m}$ et de hauteur $H = 2e = 0.1\,\mathrm{m}$ ($e = 0.05\,\mathrm{m}$), un fluide newtonien de masse volumique $\rho = 1.0\,\mathrm{kg/m}^3$ et de viscosité $\mu = 1.0\times 10^{-3}\,\mathrm{Pa\cdot s}$ s'écoule avec une vitesse maximale d'entrée $u_{\max} = 1.5\times 10^{-3}\,\mathrm{m/s}$ ($v = 0$, $\mathrm{Re}_H = \rho U_{\mathrm{moy}} H / \mu = 0.1$).

### Solution analytique exacte
Le profil de vitesse longitudinale est purement parabolique en $y$ et le champ de pression décroît linéairement en $x$ :
$$u_{\mathrm{ex}}(y) = \frac{4 u_{\max}}{H^2} y (H - y) = 0.6\,y\,(0.1 - y), \qquad v_{\mathrm{ex}}(x,y) = 0,$$
$$-\frac{\mathrm{d}p_{\mathrm{ex}}}{\mathrm{d}x} = \frac{8 \mu u_{\max}}{H^2} = \frac{2 \mu u_{\max}}{e^2} = 1.2\times 10^{-3}\,\mathrm{Pa/m}, \qquad \Delta P_{\mathrm{ex}} = 1.2\times 10^{-3}\,\mathrm{Pa}.$$

### Extrapolation quadratique à 3 points dans les mailles fantômes
Sur le maillage décalé MAC, les vitesses $u_{i,j}$ sont décalées de $\Delta y/2$ des parois $y=0$ et $y=H$. COVOF 2D implémente l'extrapolation quadratique de Lagrange aux abscisses $(- \Delta y/2, 0, +\Delta y/2, +3\Delta y/2)$ :
$$u_{i,0} = \frac{8}{3} u_{\mathrm{wall}} - 2 u_{i,1} + \frac{1}{3} u_{i,2}, \qquad u_{i,N_y+1} = \frac{8}{3} u_{\mathrm{wall}} - 2 u_{i,N_y} + \frac{1}{3} u_{i,N_y-1}.$$
En injectant cette relation dans le stencil visqueux centré au premier nœud intérieur $j=1$, on obtient :
$$\frac{u_{i,0} - 2 u_{i,1} + u_{i,2}}{\Delta y^2} = \frac{\frac{8}{3}u_{\mathrm{wall}} - 4 u_{i,1} + \frac{4}{3} u_{i,2}}{\Delta y^2},$$
qui est **exact** pour tout polynôme du second degré $u(y) \in \mathbb{P}_2$, garantissant une erreur spatiale nulle à la précision machine sur tout maillage.

---

---

## 3. Paramètres Physiques, Numériques et Discrétisation

Physiques, Numériques et Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine | $[0, L_x] \times [0, L_y]$ | $L_x = 1.0\,\mathrm{m}$, $L_y = H = 0.1\,\mathrm{m}$ ($e = 0.05\,\mathrm{m}$) |
| Propriétés du fluide | $\rho, \mu, u_{\max}$ | $\rho = 1.0\,\mathrm{kg/m}^3$, $\mu = 10^{-3}\,\mathrm{Pa\cdot s}$, $u_{\max} = 1.5\times 10^{-3}\,\mathrm{m/s}$ |
| Séquence de grilles | $(N_x \times N_y)$ | $16\times 16$, $32\times 32$, $64\times 64$, $128\times 128$ |
| Conditions aux limites | Gauche, Droite, Bas/Haut | `INLET` (`"profile": "PARABOLIC"`), `OUTLET` ($p=0$), `DIRICHLET` ($u=v=0$) |
| Solveurs de couplage | `solver.type` | `DIRECT`, `AUGMENTED_LAGRANGIAN` ($r_{\mathrm{al}} \in \{0, 1, 100\}$), `ADAPTATIVE_AUGMENTED_LAGRANGIAN` ($K_{\mathrm{aal}}=100$), `MONOLITHIC`, `STANDARD_PROJECTION` |
| Solveurs linéaires & Précond. | `linear_solver` / `precond_type` | `DIRECT`, `BICGSTAB2`, `BICGSTAB`, `GMRES`, `MINRES` avec `NONE`, `JACOBI`, `SSOR`, `ILU0`, `BLOCK_TRIANGULAR` |

---

---

## 4. Fichier de Configuration JSON Pas à Pas

de Fichiers `.json` (C++17 et Fortran 2018)

### Prototype `case/cpp/cas_3_poiseuille_grid64_al.json`
```json
{
  "case_name": "cas_3_poiseuille_grid64_al",
  "mode": "NAVIER_STOKES",
  "test_case": "POISEUILLE",
  "u_max_poiseuille": 0.0015,
  "steady_state_tol": 1e-13,
  "mesh": {
    "nx": 64,
    "ny": 64,
    "lx": 1.0,
    "ly": 0.1
  },
  "physics": {
    "enable_two_phase": false,
    "enable_surface_tension": false,
    "enable_convection": true,
    "convective_scheme": "UPWIND",
    "gx": 0.0,
    "gy": 0.0,
    "fluid1": {
      "rho": 1.0,
      "mu": 0.001
    }
  },
  "time": {
    "dt": 0.1,
    "t_final": 0.5,
    "cfl_max": 0.5,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "INLET",     "u": 0.0015, "v": 0.0, "profile": "PARABOLIC" },
    "right":  { "type": "OUTLET",    "u": 0.0,    "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0,    "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0,    "v": 0.0 }
  },
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "linear_solver": "BICGSTAB2",
    "precond_type": "BLOCK_TRIANGULAR",
    "use_precond": true,
    "tol": 1e-12,
    "max_iter": 1000,
    "r_al": 1.0,
    "uzawa_iter": 1
  },
  "output": {
    "output_dir": "VALIDATION_TOTALE/validation_cas_3_poiseuille/output/cpp/cas_3_poiseuille_grid64_al",
    "case_name": "cas_3_poiseuille_grid64_al",
    "vtk_interval": 5,
    "log_interval": 1,
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
cd VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_3_poiseuille

# 2. Exécution avec le moteur C++17 (Séquentiel)
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_3_config.json

# 3. Exécution avec le moteur C++17 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_3_config.json

# 4. Exécution avec le moteur Fortran 2018 (Séquentiel)
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_3_config.json

# 5. Exécution avec le moteur Fortran 2018 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_3_config.json

```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification et Tolérances
à Mener (Précision, Intégrité et Performance)

1. **Exactitude machine du profil de vitesse et du gradient de pression :**
   - Vérifier que sur toutes les grilles ($16\times 16$ à $128\times 128$), l'erreur sur $u$ satisfait $E_{L_1}, E_{L_2}, E_{L_\infty} < 10^{-12}\,\mathrm{m/s}$ (et $\sim 10^{-18}\,\mathrm{m/s}$ avec `DIRECT`), et que la vitesse transversale vérifie $\|v\|_\infty < 10^{-15}\,\mathrm{m/s}$.
   - Vérifier que le gradient axial de pression restitué par régression linéaire sur $p(x, H/2)$ vaut $-\mathrm{d}p/\mathrm{d}x = 1.200000\times 10^{-3}\,\mathrm{Pa/m}$ avec une erreur relative $< 10^{-11}$.
2. **Intégrité des schémas convectifs, des couplages et des préconditionneurs :**
   - Vérifier l'invariance vis-à-vis du schéma convectif (`UPWIND`, `CENTERED`, `TVD_SUPERBEE`) et de la linéarisation (`PICARD`, `ADAMS_BASHFORTH`, `NEWTON`) puisque $\partial u_{\mathrm{ex}}/\partial x = 0$ et $v_{\mathrm{ex}} = 0$.
   - Comparer la vitesse de convergence entre `DIRECT`, `MONOLITHIC`, `AUGMENTED_LAGRANGIAN` ($r_{\mathrm{al}} \in \{0, 1, 100\}$), `ADAPTATIVE_AUGMENTED_LAGRANGIAN` et `STANDARD_PROJECTION`.
   - Vérifier que le préconditionneur `BLOCK_TRIANGULAR` surpasse `JACOBI`, `SSOR` et `ILU0` sur le système point-selle complet $(u,v,p)$.
3. **Parité C++17 / Fortran 2018 :** Vérifier que $\|u^{\mathrm{C++}} - u^{\mathrm{Fortran}}\|_\infty < 10^{-14}\,\mathrm{m/s}$ et $\|p^{\mathrm{C++}} - p^{\mathrm{Fortran}}\|_\infty < 10^{-14}\,\mathrm{Pa}$.

---

### 6.2 Résultats Attendus et Visualisation ParaView

