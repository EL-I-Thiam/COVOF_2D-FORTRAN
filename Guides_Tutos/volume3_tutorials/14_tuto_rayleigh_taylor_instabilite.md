# Tutoriel 14 : Instabilité Interfaciale Diphasique de Rayleigh-Taylor (Dispersion Linéaire)

Ce dernier tutoriel simule l'instabilité hydrodynamique se développant à l'interface entre un fluide lourd superposé à un fluide léger soumis à la gravité. Le taux de croissance exponentielle en régime linéaire est validé face à la relation de dispersion analytique exacte hydro-capillaire.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les fondements physiques et les équations aux dérivées partielles gouvernant l'écoulement.
- Configurer les paramètres numériques et conditions aux limites sur grille cartésienne décalée MAC.
- Exécuter le cas test avec les exécutables C++17 et Fortran 2018 (séquentiel et OpenMP).
- Quantifier avec rigueur la précision par rapport à la solution analytique exacte ($L_2$, $L_\infty$).
- Produire et exploiter les visualisations scientifiques volumiques sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 :** `Mesh2D`, `FlowFields`, `WeymouthYueVOF`, `CSFModel`, `NavierStokesSolver`, `UzawaSolver`.
- **Moteur Fortran 2018 :** `module_mesh`, `module_flow_fields`, `module_vof`, `module_csf`, `module_navier_stokes`, `module_uzawa`.

---

## 2. Formulation Physique et Solution Analytique Exacte

Équations Physiques et Configuration

Ce vingt-troisième cas valide la compétition entre l'instabilité gravitationnelle de Rayleigh-Taylor (fluide lourd au-dessus d'un fluide léger, $At = (\rho_1-\rho_2)/(\rho_1+\rho_2) = 0.5$) et la stabilisation capillaire CSF, du régime linéaire exponentiel jusqu'à l'enroulement non-linéaire en champignon.

Dans un domaine $\Omega = [0, 1.0] \times [0, 4.0]\,\mathrm{m}^2$ avec parois latérales `FREE_SLIP` (ou `PERIODIC`) et parois haute/basse `DIRICHLET` adhérentes, le fluide lourd ($\rho_1 = 3.0\,\mathrm{kg/m}^3, \mu_1 = 10^{-3}\,\mathrm{Pa\cdot s}$) surmonte initialement le fluide léger ($\rho_2 = 1.0\,\mathrm{kg/m}^3, \mu_2 = 10^{-3}\,\mathrm{Pa\cdot s}$) sous $g_y = -9.81\,\mathrm{m/s}^2$. L'interface initiale est perturbée selon (`type: "RAYLEIGH_TAYLOR"`, `center_y: 2.0`, `pert_amp: 0.05`, `pert_k: 1.0`) :
$$y_{\mathrm{int}}(x, 0) = 2.0 + A_0 \cos(k x), \qquad k = 2\pi\,\mathrm{m}^{-1}, \quad A_0 = 0.05\,\mathrm{m}.$$

#### Relation de dispersion hydro-capillaire linéaire
Le taux de croissance théorique $\omega_{\mathrm{th}}(k)$ et le seuil de coupure capillaire $\sigma_c$ valent :
$$\omega_{\mathrm{th}}^2(k) = \frac{(\rho_1 - \rho_2)\,g\,k - \sigma\,k^3}{\rho_1 + \rho_2}, \qquad \sigma_c = \frac{(\rho_1 - \rho_2)\,g}{k^2} = 0.4970\,\mathrm{N/m}.$$
Pour $\sigma = 0$, $\omega_{\mathrm{th}} = 5.5515\,\mathrm{s}^{-1}$ ; pour $\sigma = 0.20\,\mathrm{N/m}$, $\omega_{\mathrm{th}} = 4.2915\,\mathrm{s}^{-1}$ ; et pour $\sigma = 0.60\,\mathrm{N/m} > \sigma_c$, l'interface est stabilisée et oscille sans croître.

---

---

## 3. Paramètres Physiques, Numériques et Discrétisation

Physiques, Numériques et Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et perturbation initiale | $[0, 1]\times[0, 4]\,\mathrm{m}^2$, $y_c, A_0, k$ | $y_c = 2.0\,\mathrm{m}$, $A_0 = 0.05\,\mathrm{m}$, $m=1$ ($k = 2\pi\,\mathrm{m}^{-1}$, `RAYLEIGH_TAYLOR`) |
| Fluides ($At = 0.5$) et gravité | $\rho_1, \rho_2, \mu_1, \mu_2, g_y$ | $\rho_1 = 3.0\,\mathrm{kg/m}^3$, $\rho_2 = 1.0\,\mathrm{kg/m}^3$, $\mu_1 = \mu_2 = 10^{-3}\,\mathrm{Pa\cdot s}$, $g_y = -9.81\,\mathrm{m/s}^2$ |
| Trois régimes capillaires $\sigma$ | $\sigma \in \{0.0, 0.20, 0.60\}\,\mathrm{N/m}$ | Instable pur ($\omega_{\mathrm{th}}=5.55\,\mathrm{s}^{-1}$), ralenti ($\omega_{\mathrm{th}}=4.29\,\mathrm{s}^{-1}$), stabilisé ($\sigma > \sigma_c=0.497\,\mathrm{N/m}$) |
| Séquence de grilles et $\Delta t$ | $64\times 256$, $128\times 512$, $256\times 1024$ | $\Delta t = 10^{-3}\,\mathrm{s}$, $t_{\mathrm{final}} = 1.0\,\mathrm{s}$ |
| Schémas VOF et solveurs | `vof` \& `solver` | `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC`, `VOF_TVD` ; `ADAPTATIVE_AUGMENTED_LAGRANGIAN`, `AUGMENTED_LAGRANGIAN`, `MONOLITHIC` |

---

---

## 4. Fichier de Configuration JSON Pas à Pas

de Fichiers `.json` (C++17 et Fortran 2018)

### Prototype `case/cpp/cas_14_rayleigh_taylor_sigma0_128x512.json`
*(Pour le moteur Fortran, remplacer `output/cpp/` par `output/fortran/` dans `output_dir`)*
```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 512,
    "lx": 1.0,
    "ly": 4.0
  },
  "physics": {
    "enable_two_phase": true,
    "enable_convection": true,
    "convective_scheme": "TVD_SUPERBEE",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_surface_tension": false,
    "sigma": 0.0,
    "gx": 0.0,
    "gy": -9.81,
    "fluid1": { "rho": 3.0, "mu": 1.0e-3 },
    "fluid2": { "rho": 1.0, "mu": 1.0e-3 }
  },
  "initial_condition": {
    "type": "RAYLEIGH_TAYLOR",
    "center_y": 2.0,
    "pert_amp": 0.05,
    "pert_k": 1.0
  },
  "time": {
    "dt": 0.001,
    "t_final": 1.0,
    "time_scheme": "GEAR",
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "ADAPTATIVE_AUGMENTED_LAGRANGIAN",
    "adaptative_al": true,
    "k_aal": 100.0,
    "uzawa_tol": 1e-6,
    "uzawa_iter": 4,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-08,
    "max_iter": 800,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "SUPERBEE",
    "cfl_vof": 0.3
  },
  "output": {
    "output_dir": "VALIDATION_TOTALE/validation_cas_14_instabilite_interface_gravitaire/output/cpp/sigma0_128x512",
    "case_name": "cas_14_rayleigh_taylor_sigma0_128x512",
    "vtk_interval": 50,
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
cd VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_14_instabilite_interface_gravitaire

# 2. Exécution avec le moteur C++17 (Séquentiel)
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_14_config.json

# 3. Exécution avec le moteur C++17 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_14_config.json

# 4. Exécution avec le moteur Fortran 2018 (Séquentiel)
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_14_config.json

# 5. Exécution avec le moteur Fortran 2018 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_14_config.json

```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification et Tolérances
à Mener (Précision, Intégrité et Performance)

1. **Taux de croissance linéaire $\omega_{\mathrm{num}}$ et seuil de stabilisation capillaire $\sigma_c$ :** Ajuster $A(t) = A_0 \cosh(\omega t)$ sur $t \in [0.05, 0.30]\,\mathrm{s}$ pour $\sigma \in \{0.0, 0.20, 0.60\}\,\mathrm{N/m}$ et vérifier l'accord avec $\omega_{\mathrm{th}}(0) = 5.5515\,\mathrm{s}^{-1}$, $\omega_{\mathrm{th}}(0.2) = 4.2915\,\mathrm{s}^{-1}$ et la stabilisation capillaire pour $\sigma = 0.60 > \sigma_c$.
2. **Enroulement non-linéaire en champignon ($t \in [0.5, 1.0]\,\mathrm{s}$), symétrie et conservation de masse :** Vérifier la symétrie $C(L_x-x,y,t) = C(x,y,t)$, $\Delta V / V_0 < 10^{-5}$, et comparer `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC` et `VOF_TVD`.
3. **Performance des solveurs (`AAL`, `AL`, `MONOLITHIC`) et parité C++17 / Fortran 2018 :** Comparer les temps CPU et vérifier la superposition des interfaces C++ et Fortran.

---

### 6.2 Résultats Attendus et Visualisation ParaView

