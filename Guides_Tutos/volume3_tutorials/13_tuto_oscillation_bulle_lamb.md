# Tutoriel 13 : Oscillations Capillaires Amorties d'une Bulle de Gaz (Théorie de Lamb)

Ce tutoriel met en œuvre une bulle initialement déformée en ellipse selon le mode n=2 sous l'action combinée de la tension superficielle et de la dissipation visqueuse. La fréquence propre d'oscillation et le taux d'amortissement sont confrontés à la théorie linéaire exacte de Lamb.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les fondements physiques et les équations aux dérivées partielles gouvernant l'écoulement.
- Configurer les paramètres numériques et conditions aux limites sur grille cartésienne décalée MAC.
- Exécuter le cas test avec les exécutables C++17 et Fortran 2018 (séquentiel et OpenMP).
- Quantifier avec rigueur la précision par rapport à la solution analytique exacte ($L_2$, $L_\infty$).
- Produire et exploiter les visualisations scientifiques volumiques sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 :** `Mesh2D`, `FlowFields`, `WeymouthYueVOF`, `CSFModel`, `NavierStokesSolver`.
- **Moteur Fortran 2018 :** `module_mesh`, `module_flow_fields`, `module_vof`, `module_csf`, `module_navier_stokes`.

---

## 2. Formulation Physique et Solution Analytique Exacte

Équations Physiques et Configuration

Ce douzième cas valide la dynamique instationnaire du couplage entre l'inertie diphasique à fort rapport de densité ($\rho_1/\rho_2 \approx 5950$), la dissipation visqueuse et la force de tension superficielle CSF lors des oscillations libres d'une bulle de gaz initialement elliptique.

Dans une cavité carrée $\Omega = [0, 0.1] \times [0, 0.1]\,\mathrm{m}^2$ sans gravité ($\mathbf{g} = \mathbf{0}$) avec parois à glissement libre (`FREE_SLIP`), une bulle d'air ($\rho_2 = 1.17683\,\mathrm{kg/m}^3$, $\mu_2 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$) est centrée en $(0.05, 0.05)\,\mathrm{m}$ dans un métal liquide dense ($\rho_1 = 7000\,\mathrm{kg/m}^3$, $\mu_1 = 4.0\times 10^{-3}\,\mathrm{Pa\cdot s}$) avec $\sigma = 1.8\,\mathrm{N/m}$. À $t=0$, la bulle de rayon équivalent $R_0 = 0.02\,\mathrm{m}$ est déformée sur le mode quadrupolaire $n=2$ ($a_0 = R_0(1+\varepsilon)$, $b_0 = R_0/(1+\varepsilon)$, avec $\varepsilon = 0.20$).

#### Relation de dispersion capillaire linéaire
La pulsation propre $\omega_n$ et la fréquence $f_{\mathrm{th}}$ du mode $n=2$ en milieu infini valent :
$$\omega_n^2 = \frac{(n^3 - n)\,\sigma}{(\rho_1 + \rho_2)\,R_0^3} \implies \omega_2 = \sqrt{\frac{6\,\sigma}{(\rho_1 + \rho_2)\,R_0^3}}, \qquad f_{\mathrm{th}} = \frac{\omega_2}{2\pi} = 5.70632\,\mathrm{s}^{-1} \quad (T_{\mathrm{th}} = 0.17524\,\mathrm{s}).$$
Sur les grilles $32^2 \dots 512^2$ en cavité confinée, la fréquence numérique croît de $4.368\,\mathrm{Hz}$ ($32^2$), $4.792\,\mathrm{Hz}$ ($64^2$), $4.952\,\mathrm{Hz}$ ($128^2$), $4.990\,\mathrm{Hz}$ ($256^2$) à $5.038\,\mathrm{Hz}$ ($512^2$).

---

---

## 3. Paramètres Physiques, Numériques et Discrétisation

Physiques, Numériques et Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et bulle initiale | $[0, L]^2$, $(x_c, y_c)$, $R_0$, $\varepsilon$ | $L = 0.1\,\mathrm{m}$, centre $(0.05, 0.05)\,\mathrm{m}$, $R_0 = 0.02\,\mathrm{m}$, $\varepsilon \in \{0.05, 0.20\}$ |
| Fluide 1 (Métal liquide) | $\rho_1, \mu_1, \chi_{T,1}$ | $\rho_1 = 7000.0\,\mathrm{kg/m}^3$, $\mu_1 = 4.0\times 10^{-3}\,\mathrm{Pa\cdot s}$, $\chi_{T,1} = 10^{-9}\,\mathrm{Pa}^{-1}$ |
| Fluide 2 (Gaz / Air) | $\rho_2, \mu_2, \chi_{T,2}$ | $\rho_2 = 1.17683\,\mathrm{kg/m}^3$, $\mu_2 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$, $\chi_{T,2} = 9.8692327\times 10^{-6}\,\mathrm{Pa}^{-1}$ |
| Tension superficielle et lissage | $\sigma$, `enable_diffusion_smoothing` | $\sigma = 1.8\,\mathrm{N/m}$, `enable_diffusion_smoothing: true`, `diffusion_smoothing_steps` $\in \{2, 3\}$ |
| Séquence de grilles et $\Delta t$ | $N \in \{32, 64, 128, 256, 512\}$ | $\Delta t = 2.5\times 10^{-4}\,\mathrm{s}$, $t_{\mathrm{final}} = 0.6\,\mathrm{s}$ (plus de 3 périodes complètes) |
| Solveurs et schémas VOF | `solver` \& `vof` | `AUGMENTED_LAGRANGIAN`, `ADAPTATIVE_AUGMENTED_LAGRANGIAN`, `MONOLITHIC` ; `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC` |

---

---

## 4. Fichier de Configuration JSON Pas à Pas

de Fichiers `.json` (C++17 et Fortran 2018)

### Prototype `case/cpp/cas_13_oscillation_bulle_n128.json`
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
    "compressible": false,
    "enable_convection": true,
    "enable_surface_tension": true,
    "enable_diffusion_smoothing": true,
    "diffusion_smoothing_steps": 3,
    "sigma": 1.8,
    "gx": 0.0,
    "gy": 0.0,
    "fluid1": { "rho": 7000.0, "mu": 4.0e-3 },
    "fluid2": { "rho": 1.17683, "mu": 1.85e-5 }
  },
  "initial_condition": {
    "type": "BUBBLE",
    "center_x": 0.05,
    "center_y": 0.05,
    "radius": 0.02,
    "axis_a": 0.024,
    "axis_b": 0.0166666667
  },
  "time": {
    "dt": 0.00025,
    "t_final": 0.6,
    "time_scheme": "GEAR",
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "top":    { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "ADAPTATIVE_AUGMENTED_LAGRANGIAN",
    "adaptative_al": true,
    "k_aal": 100.0,
    "uzawa_tol": 1e-6,
    "uzawa_iter": 6,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-08,
    "max_iter": 600,
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
    "output_dir": "VALIDATION_TOTALE/validation_cas_13_oscillation_bulle_capillaire/output/cpp/n128",
    "case_name": "cas_13_oscillation_bulle_n128",
    "vtk_interval": 20,
    "log_interval": 5,
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
cd VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_13_oscillation_bulle_capillaire

# 2. Exécution avec le moteur C++17 (Séquentiel)
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_13_config.json

# 3. Exécution avec le moteur C++17 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_13_config.json

# 4. Exécution avec le moteur Fortran 2018 (Séquentiel)
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_13_config.json

# 5. Exécution avec le moteur Fortran 2018 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_13_config.json

```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification et Tolérances
à Mener (Précision, Intégrité et Performance)

1. **Fréquence d'oscillation quadrupolaire $f_{\mathrm{num}}$ et convergence en maillage :** Extraire les demi-axes horizontal $a(t)$ et vertical $b(t)$ du contour $C=0.5$ ainsi que l'énergie cinétique $E_k(t)$ dans `*_history.csv`. Vérifier la convergence de $f_{\mathrm{num}}$ sur $32^2 \dots 512^2$ ($4.368 \to 4.792 \to 4.952 \to 4.990 \to 5.038\,\mathrm{Hz}$ pour $\varepsilon=0.20$).
2. **Robustesse du couplage vitesse-pression à haut rapport de densité ($\rho_1/\rho_2 \approx 5950$) :** Comparer `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (`k_aal=100`), `AUGMENTED_LAGRANGIAN` ($r_{\mathrm{al}}=1000$) et `MONOLITHIC` avec `GEAR` vs `EULER`. Vérifier le gain d'itérations procuré par `AAL` et `BLOCK_TRIANGULAR`.
3. **Conservation du volume de la bulle et parité C++17 / Fortran 2018 :** Vérifier que l'erreur relative sur l'aire de la bulle reste $< 10^{-4}$ avec `WEYMOUTH_YUE+PLIC` et que les courbes $a(t), b(t), E_k(t)$ C++ et Fortran se superposent.

---

### 6.2 Résultats Attendus et Visualisation ParaView

