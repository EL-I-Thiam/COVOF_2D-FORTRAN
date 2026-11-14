# Tutoriel 9 : Sédimentation Confinée 2D d'une Inclusion Cylindrique et Loi de Faxén

Ce tutoriel modélise la sédimentation d'un cylindre solide en régime rampant de Stokes (Re << 1) entre deux parois planes verticales. Il confronte la vitesse limite de chute numérique stationnaire avec la loi analytique de Faxén intégrant la correction de confinement pariétal.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les fondements physiques et les équations aux dérivées partielles gouvernant l'écoulement.
- Configurer les paramètres numériques et conditions aux limites sur grille cartésienne décalée MAC.
- Exécuter le cas test avec les exécutables C++17 et Fortran 2018 (séquentiel et OpenMP).
- Quantifier avec rigueur la précision par rapport à la solution analytique exacte ($L_2$, $L_\infty$).
- Produire et exploiter les visualisations scientifiques volumiques sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 :** `Mesh2D`, `FlowFields`, `BrinkmanPenalization`, `StokesSolver`, `BiCGSTAB2Solver`, `SSORPreconditioner`.
- **Moteur Fortran 2018 :** `module_mesh`, `module_flow_fields`, `module_brinkman`, `module_stokes`, `module_linear_solvers`, `module_preconditioners`.

---

## 2. Formulation Physique et Solution Analytique Exacte

Équations Physiques et Configuration

Ce vingt-deuxième cas valide la méthode de pénalisation de viscosité combinée au Lagrangien Augmenté Adaptatif (`ADAPTATIVE_AUGMENTED_LAGRANGIAN`) pour la sédimentation lente d'un cylindre rigide confiné entre deux parois verticales, comparée à la loi analytique de Faxén ($Re_p \ll 1, Ga < 1$).

Dans une cavité carrée fermée $\Omega = [0, 0.1] \times [0, 0.1]\,\mathrm{m}^2$ ($2l = 0.1\,\mathrm{m}$, $l = 0.05\,\mathrm{m}$) à parois adhérentes ($\mathbf{u} = \mathbf{0}$), un cylindre de rayon $a = 0.00625\,\mathrm{m}$ (rapport de confinement $k_x = a/l = 1/8 = 0.125$) est centré en $(0.05, 0.05)\,\mathrm{m}$. Le fluide porteur ($C=0$) a $\rho_f = 1000\,\mathrm{kg/m}^3, \mu_f = 0.153125\,\mathrm{Pa\cdot s}$, et la particule ($C=1$) a $\rho_p = 1500\,\mathrm{kg/m}^3, \mu_p = 1000\,\mathrm{Pa\cdot s}$ sous $g_y = -9.81\,\mathrm{m/s}^2$.

#### Loi analytique de vitesse terminale de Faxén
En régime de Stokes 2D confiné ($k_x = a/l = 0.125$), la vitesse terminale analytique vaut :
$$U_{c,\mathrm{th}} = \frac{(\rho_f - \rho_p)\,g\,a^2}{4\,\mu_f}\,\left(-0.9157 - \ln(k_x) + 1.7244\,k_x^2 - 1.7302\,k_x^4\right) = -0.05060\,\mathrm{m/s},$$
$$E_{c,\mathrm{th}} \approx 2.7554\times 10^{-3}\,\mathrm{J}.$$

---

---

## 3. Paramètres Physiques, Numériques et Discrétisation

Physiques, Numériques et Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et confinement $k_x$ | $[0, 2l]^2$, $a$, $k_x = a/l$ | $2l = 0.1\,\mathrm{m}$ ($l=0.05\,\mathrm{m}$), $a = 0.00625\,\mathrm{m}$, $k_x = 1/8 = 0.125$, centre $(0.05, 0.05)\,\mathrm{m}$ |
| Particule ($C=1$) et Fluide ($C=0$) | $\rho_p, \mu_p$ ; $\rho_f, \mu_f$ | $\rho_1 = \rho_p = 1500\,\mathrm{kg/m}^3, \mu_1 = 1000\,\mathrm{Pa\cdot s}$ ; $\rho_2 = \rho_f = 1000\,\mathrm{kg/m}^3, \mu_2 = 0.153125\,\mathrm{Pa\cdot s}$ |
| Séquence de grilles et $\Delta t$ | $64\times 128$, $128\times 128$, $256\times 256$ | $\Delta t = 5\times 10^{-4}\,\mathrm{s}$, $t_{\mathrm{final}} = 0.15\,\mathrm{s}$ (atteinte du plateau quasi-stationnaire dès $t \approx 0.04\,\mathrm{s}$) |
| Solveurs et filtrage SVOF | `solver` \& `svof` | `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (`k_aal`=100), `AUGMENTED_LAGRANGIAN`, `MONOLITHIC` ; `enable_svof_smoothing: true` |

---

---

## 4. Fichier de Configuration JSON Pas à Pas

de Fichiers `.json` (C++17 et Fortran 2018)

### Prototype `case/cpp/cas_9_sedimentation_faxen_n128.json`
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
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "PICARD",
    "enable_surface_tension": false,
    "enable_svof_smoothing": true,
    "svof_filter_width": 0.5,
    "gx": 0.0,
    "gy": -9.81,
    "fluid1": { "rho": 1500.0, "mu": 1000.0 },
    "fluid2": { "rho": 1000.0, "mu": 0.153125 }
  },
  "initial_condition": {
    "type": "DROPLET",
    "center_x": 0.05,
    "center_y": 0.05,
    "radius": 0.00625
  },
  "time": {
    "dt": 0.0005,
    "t_final": 0.15,
    "time_scheme": "GEAR",
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "right":  { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "ADAPTATIVE_AUGMENTED_LAGRANGIAN",
    "adaptative_al": true,
    "k_aal": 100.0,
    "uzawa_tol": 1e-7,
    "uzawa_iter": 6,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-08,
    "max_iter": 800,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.3
  },
  "output": {
    "output_dir": "VALIDATION_TOTALE/validation_cas_9_sedimentation_faxen/output/cpp/n128",
    "case_name": "cas_9_sedimentation_faxen_n128",
    "vtk_interval": 30,
    "log_interval": 10,
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
cd VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_9_sedimentation_faxen

# 2. Exécution avec le moteur C++17 (Séquentiel)
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_9_config.json

# 3. Exécution avec le moteur C++17 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_9_config.json

# 4. Exécution avec le moteur Fortran 2018 (Séquentiel)
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_9_config.json

# 5. Exécution avec le moteur Fortran 2018 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_9_config.json

```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification et Tolérances
à Mener (Précision, Intégrité et Performance)

1. **Précision de la vitesse terminale vis-à-vis de la loi analytique de Faxén ($U_{c,\mathrm{th}} = -0.05060\,\mathrm{m/s}$) :** Extraire $V_{\mathrm{drop}}(t)$ et $E_k(t)$ sur $64\times 128, 128\times 128, 256\times 256$ et vérifier que le plateau de vitesse concorde à moins de $3.5\,\%$ avec $-0.05060\,\mathrm{m/s}$ ($E_{k,\max} \approx 2.755\times 10^{-3}\,\mathrm{J}$).
2. **Indéformabilité du cylindre ($\mu_p/\mu_f \approx 6530$) et rôle de `AAL` :** Vérifier que la pénalité adaptative $r_{i,j}(t)$ verrouille l'incompressibilité et la rigidité sur l'interface du cylindre tout en accélérant `BICGSTAB2`.
3. **Parité C++17 / Fortran 2018 :** Vérifier l'identité des courbes $V_{\mathrm{drop}}(t)$ et $E_k(t)$ entre C++ et Fortran.

---

### 6.2 Résultats Attendus et Visualisation ParaView

