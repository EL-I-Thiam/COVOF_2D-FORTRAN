# Tutoriel 11 : Compression Isotrope d'une Bulle d'Air en Cavité Carrée

Ce tutoriel évalue la capacité de l'algorithme C-AAL à traiter la compression axisymétrique ou uniforme d'une bulle de gaz immergée. L'évolution temporelle du rayon de la bulle et de la densité de la phase compressible est confrontée à la solution exacte déduite du bilan de masse global.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les fondements physiques et les équations aux dérivées partielles gouvernant l'écoulement.
- Configurer les paramètres numériques et conditions aux limites sur grille cartésienne décalée MAC.
- Exécuter le cas test avec les exécutables C++17 et Fortran 2018 (séquentiel et OpenMP).
- Quantifier avec rigueur la précision par rapport à la solution analytique exacte ($L_2$, $L_\infty$).
- Produire et exploiter les visualisations scientifiques volumiques sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 :** `Mesh2D`, `FlowFields`, `CAALSolver`, `VOFCompressible`, `BiCGSTAB2Solver`.
- **Moteur Fortran 2018 :** `module_mesh`, `module_flow_fields`, `module_caal`, `module_vof_compressible`, `module_linear_solvers`.

---

## 2. Formulation Physique et Solution Analytique Exacte

Équations Physiques et Configuration

Ce treizième cas valide la compression multidimensionnelle centripète d'une bulle d'air circulaire immergée dans de l'eau au moyen du modèle 1-fluide faiblement compressible (`compressible: true`).

Dans une cavité carrée $\Omega = [0, 0.1] \times [0, 0.1]\,\mathrm{m}^2$ sans gravité ($\mathbf{g}=\mathbf{0}$), une bulle d'air de rayon initial $R_0 = 0.03\,\mathrm{m}$ ($C=0$, $\rho_{2,0} = 1.176829268\,\mathrm{kg/m}^3$, $\mu_2 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$, $\chi_{T,2} = 9.8692327\times 10^{-6}\,\mathrm{Pa}^{-1}$) est centrée en $(0.05, 0.05)\,\mathrm{m}$ dans de l'eau ($C=1$, $\rho_{1,0} = 1000\,\mathrm{kg/m}^3$, $\mu_1 = 10^{-3}\,\mathrm{Pa\cdot s}$, $\chi_{T,1} = 0.444\times 10^{-9}\,\mathrm{Pa}^{-1}$). Sur chacune des 4 parois, de l'eau est injectée normalement vers l'intérieur avec une vitesse $V_0 = 0.0025\,\mathrm{m/s}$, soit un débit entrant total $q = 4 L V_0 = 1.0\times 10^{-3}\,\mathrm{m}^2/\mathrm{s}$.

#### Loi analytique d'évolution du rayon, de la masse volumique et de la pression
L'aire de la bulle décroît linéairement selon $\pi R^2(t) = \pi R_0^2 - q\,t$ (temps d'effondrement complet $t_{\mathrm{collapse}} = \pi R_0^2 / q = 2.827433388\,\mathrm{s}$), d'où :
$$R(t) = R_0\sqrt{1 - \frac{q\,t}{\pi R_0^2}}, \qquad \rho_{\mathrm{air}}(t) = \frac{\rho_{\mathrm{air},0}}{1 - \dfrac{q\,t}{\pi R_0^2}}, \qquad p(t) = \frac{p_0}{1 - \dfrac{q\,t}{\pi R_0^2}}.$$
À $t = 1.5\,\mathrm{s}$ ($\Delta t = 0.00125\,\mathrm{s}$), on obtient $\rho_{\mathrm{th}}(1.5\,\mathrm{s}) = 2.50664658\,\mathrm{kg/m}^3$, $R_{\mathrm{th}}(1.5\,\mathrm{s}) = 0.0205543\,\mathrm{m}$ et $p_{\mathrm{th}}(1.5\,\mathrm{s}) = 215\,821.8\,\mathrm{Pa}$.

---

---

## 3. Paramètres Physiques, Numériques et Discrétisation

Physiques, Numériques et Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et bulle initiale | $[0, L]^2$, $(x_c, y_c)$, $R_0$ | $L = 0.1\,\mathrm{m}$, centre $(0.05, 0.05)\,\mathrm{m}$, $R_0 = 0.03\,\mathrm{m}$, $p_0 = 101\,325\,\mathrm{Pa}$ |
| Fluides 1 (Eau) et 2 (Air) | $\rho_{i,0}, \mu_i, \chi_{T,i}$ | Eau : $\rho_{1,0}=1000$, $\chi_{T,1}=4.44\times 10^{-10}$ ; Air : $\rho_{2,0}=1.176829268$, $\chi_{T,2}=9.8692327\times 10^{-6}\,\mathrm{Pa}^{-1}$ |
| Conditions d'injection (4 bords) | Gauche, Droite, Bas, Haut | `DIRICHLET` entrant : $u_{\mathrm{left}}=+V_0$, $u_{\mathrm{right}}=-V_0$, $v_{\mathrm{bottom}}=+V_0$, $v_{\mathrm{top}}=-V_0$ ($V_0=0.0025\,\mathrm{m/s})$ |
| Séquence de grilles et $\Delta t$ | $N \in \{32, 64, 128, 256, 512\}$ | $\Delta t = 0.00125\,\mathrm{s}$, $t_{\mathrm{final}} = 1.5\,\mathrm{s}$ ($1200$ pas) et test prolongé jusqu'à $t = 2.5\,\mathrm{s}$ |
| Solveurs et schémas VOF | `solver` \& `vof` | `MONOLITHIC` avec `BICGSTAB2` + `BLOCK_TRIANGULAR` / `DIRECT` ; `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC`, `VOF_TVD` |

---

---

## 4. Fichier de Configuration JSON Pas à Pas

de Fichiers `.json` (C++17 et Fortran 2018)

### Prototype `case/cpp/cas_11_compression_bulle_n128.json`
*(Pour le moteur Fortran, remplacer `output/cpp/` par `output/fortran/` dans `output_dir`)*
```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 0.1,
    "ly": 0.1
  },
  "physics": {
    "enable_two_phase": true,
    "compressible": true,
    "enable_convection": true,
    "enable_surface_tension": false,
    "gx": 0.0,
    "gy": 0.0,
    "fluid1": { "rho": 1000.0, "mu": 1.0e-3, "chi_t": 4.44e-10 },
    "fluid2": { "rho": 1.176829268, "mu": 1.85e-5, "chi_t": 9.8692327e-6 }
  },
  "initial_condition": {
    "type": "BUBBLE",
    "center_x": 0.05,
    "center_y": 0.05,
    "radius": 0.03,
    "p_init": 101325.0
  },
  "time": {
    "dt": 0.00125,
    "t_final": 1.5,
    "time_scheme": "EULER",
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u":  0.0025, "v":  0.0 },
    "right":  { "type": "DIRICHLET", "u": -0.0025, "v":  0.0 },
    "bottom": { "type": "DIRICHLET", "u":  0.0,    "v":  0.0025 },
    "top":    { "type": "DIRICHLET", "u":  0.0,    "v": -0.0025 }
  },
  "solver": {
    "type": "MONOLITHIC",
    "linear_solver": "BICGSTAB2",
    "tol": 1e-10,
    "max_iter": 1000,
    "r_al": 0.0,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "SUPERBEE",
    "cfl_vof": 0.5
  },
  "output": {
    "output_dir": "VALIDATION_TOTALE/validation_cas_11_compression_bulle/output/cpp/n128",
    "case_name": "cas_11_compression_bulle_n128",
    "vtk_interval": 100,
    "log_interval": 20,
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
cd VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_11_compression_bulle

# 2. Exécution avec le moteur C++17 (Séquentiel)
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_11_config.json

# 3. Exécution avec le moteur C++17 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_11_config.json

# 4. Exécution avec le moteur Fortran 2018 (Séquentiel)
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_11_config.json

# 5. Exécution avec le moteur Fortran 2018 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_11_config.json

```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification et Tolérances
à Mener (Précision, Intégrité et Performance)

1. **Convergence vers $\rho_{\mathrm{th}}(1.5\,\mathrm{s}) = 2.50664658\,\mathrm{kg/m}^3$ :** Évaluer $\overline{\rho}_{\mathrm{air}}(1.5\,\mathrm{s})$ sur $32^2, 64^2, 128^2, 256^2, 512^2$. Vérifier que l'erreur relative décroît de $1.96\,\%$ ($32^2$) et $0.45\,\%$ ($64^2$) à $0.015\,\%$ ($512^2$), avec un ordre global $p \approx 1.76 - 2.00$.
2. **Circularité de la bulle comprimée et bilan acoustique global :** Vérifier que $\int_\Omega \chi_T(C)\frac{p^{n+1}-p^n}{\Delta t}\,\mathrm{d}\Omega = q$ absorbe exactement le débit injecté sur les 4 parois, et comparer la circularité du contour $C=0.5$ entre `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC` et `VOF_TVD`.
3. **Performance des préconditionneurs et parité dual-code :** Comparer `BLOCK_TRIANGULAR`, `ILU0` et `SSOR` avec `BICGSTAB2` et `GMRES`, et vérifier la parité C++17 / Fortran 2018.

---

### 6.2 Résultats Attendus et Visualisation ParaView

