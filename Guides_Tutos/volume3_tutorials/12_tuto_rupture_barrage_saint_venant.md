# Tutoriel 12 : Rupture de Barrage 2D sur Fond Sec (Solution de Saint-Venant)

Ce tutoriel aborde l'effondrement gravitationnel instantané d'une colonne d'eau sur fond sec. La cinématique d'avancée du front sec et le profil parabolique de la surface libre sont confrontés à la solution exacte issue de la théorie des caractéristiques de Saint-Venant (solution de Ritter).

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les fondements physiques et les équations aux dérivées partielles gouvernant l'écoulement.
- Configurer les paramètres numériques et conditions aux limites sur grille cartésienne décalée MAC.
- Exécuter le cas test avec les exécutables C++17 et Fortran 2018 (séquentiel et OpenMP).
- Quantifier avec rigueur la précision par rapport à la solution analytique exacte ($L_2$, $L_\infty$).
- Produire et exploiter les visualisations scientifiques volumiques sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 :** `Mesh2D`, `FlowFields`, `WeymouthYueVOF`, `NavierStokesSolver`, `UzawaSolver`.
- **Moteur Fortran 2018 :** `module_mesh`, `module_flow_fields`, `module_vof`, `module_navier_stokes`, `module_uzawa`.

---

## 2. Formulation Physique et Solution Analytique Exacte

Équations Physiques et Configuration

Ce vingtième cas valide le couplage Navier-Stokes diphasique à surface libre ($\rho_1/\rho_2 \approx 850$) sous gravité sur le problème de la rupture de barrage 2D sur fond sec (*Dam Break*).

Dans un canal $\Omega = [0, 1.2] \times [0, 0.2]\,\mathrm{m}^2$ à parois glissantes (`FREE_SLIP`), une colonne d'eau occupe initialement $[0, 0.6] \times [0, 0.1]\,\mathrm{m}^2$ ($x_0 = 0.6\,\mathrm{m}, h_g = 0.1\,\mathrm{m}$) avec $\rho_1 = 1000\,\mathrm{kg/m}^3, \mu_1 = 10^{-3}\,\mathrm{Pa\cdot s}$, $\rho_2 = 1.1768\,\mathrm{kg/m}^3, \mu_2 = 10^{-5}\,\mathrm{Pa\cdot s}$, $\sigma = 0.075\,\mathrm{N/m}$ et $g_y = -9.81\,\mathrm{m/s}^2$.

#### Solution analytique des caractéristiques de Saint-Venant
Avec la célérité d'onde $c_0 = \sqrt{g h_g} = 0.990454\,\mathrm{m/s}$, la hauteur d'eau analytique en eau peu profonde vaut :
$$h_{\mathrm{SV}}(x,t) = \begin{cases} h_g & \text{si } x \le x_0 - c_0 t, \\ \dfrac{1}{9g}\left(2c_0 - \dfrac{x - x_0}{t}\right)^2 & \text{si } x_0 - c_0 t < x < x_0 + 2c_0 t, \\ 0 & \text{si } x \ge x_0 + 2c_0 t. \end{cases}$$
Les profils de surface libre $h_{\mathrm{num}}(x,t) = \int_0^{L_y} C(x,y,t)\,\mathrm{d}y$ sont comparés à $h_{\mathrm{SV}}(x,t)$ aux instants $t = 0.24\,\mathrm{s}$ et $t = 0.40\,\mathrm{s}$.

---

---

## 3. Paramètres Physiques, Numériques et Discrétisation

Physiques, Numériques et Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et retenue initiale | $[0, 1.2]\times[0, 0.2]\,\mathrm{m}^2$, $x_0, h_g$ | $x_0 = 0.6\,\mathrm{m}$, $h_g = 0.1\,\mathrm{m}$ (`DAM_BREAK`), $g_y = -9.81\,\mathrm{m/s}^2$ |
| Fluides 1 (Eau) et 2 (Air) | $\rho_1, \mu_1, \rho_2, \mu_2, \sigma$ | $\rho_1 = 1000\,\mathrm{kg/m}^3$, $\mu_1 = 10^{-3}\,\mathrm{Pa\cdot s}$ ; $\rho_2 = 1.1768\,\mathrm{kg/m}^3$, $\mu_2 = 10^{-5}\,\mathrm{Pa\cdot s}$ ; $\sigma = 0.075\,\mathrm{N/m}$ |
| Séquence de grilles et $\Delta t$ | $120\times 20$, $240\times 40$, $480\times 80$ | $\Delta t \in \{10^{-3}, 5\times 10^{-4}, 2.5\times 10^{-4}\}\,\mathrm{s}$, $t_{\mathrm{final}} = 0.50\,\mathrm{s}$ (instants clés $0.24\,\mathrm{s}$ et $0.40\,\mathrm{s}$) |
| Schémas VOF et solveurs | `vof` \& `solver` | `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC`, `VOF_TVD` (`SUPERBEE`) ; `ADAPTATIVE_AUGMENTED_LAGRANGIAN`, `AUGMENTED_LAGRANGIAN`, `MONOLITHIC` |

---

---

## 4. Fichier de Configuration JSON Pas à Pas

de Fichiers `.json` (C++17 et Fortran 2018)

### Prototype `case/cpp/cas_12_rupture_barrage_240x40.json`
*(Pour le moteur Fortran, remplacer `output/cpp/` par `output/fortran/` dans `output_dir`)*
```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 240,
    "ny": 40,
    "lx": 1.2,
    "ly": 0.2
  },
  "physics": {
    "enable_two_phase": true,
    "enable_convection": true,
    "convective_scheme": "TVD_SUPERBEE",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_surface_tension": true,
    "enable_diffusion_smoothing": true,
    "diffusion_smoothing_steps": 2,
    "sigma": 0.075,
    "gx": 0.0,
    "gy": -9.81,
    "fluid1": { "rho": 1000.0, "mu": 1.0e-3 },
    "fluid2": { "rho": 1.1768, "mu": 1.0e-5 }
  },
  "initial_condition": {
    "type": "DAM_BREAK",
    "box_x0": 0.0,
    "box_y0": 0.0,
    "box_x1": 0.6,
    "box_y1": 0.1
  },
  "time": {
    "dt": 0.0005,
    "t_final": 0.50,
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
    "uzawa_iter": 5,
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
    "output_dir": "VALIDATION_TOTALE/validation_cas_12_rupture_barrage/output/cpp/grid_240x40",
    "case_name": "cas_12_rupture_barrage_240x40",
    "vtk_interval": 40,
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
cd VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_12_rupture_barrage

# 2. Exécution avec le moteur C++17 (Séquentiel)
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_12_config.json

# 3. Exécution avec le moteur C++17 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_12_config.json

# 4. Exécution avec le moteur Fortran 2018 (Séquentiel)
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_12_config.json

# 5. Exécution avec le moteur Fortran 2018 (Parallèle OpenMP)
export OMP_NUM_THREADS=4
../../../COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_12_config.json

```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification et Tolérances
à Mener (Précision, Intégrité et Performance)

1. **Profils de surface libre $h(x,t)$ à $t=0.24\,\mathrm{s}$ et $t=0.40\,\mathrm{s}$ et front d'onde $x_{\mathrm{front}}(t)$ :** Extraire $h(x,t) = \int_0^{L_y} C(x,y,t)\,\mathrm{d}y$ sur $120\times 20, 240\times 40, 480\times 80$ et comparer à la solution analytique de Saint-Venant $h_{\mathrm{SV}}(x,t)$.
2. **Conservation du volume d'eau et comparaison des 3 schémas VOF :** Vérifier que $\Delta V / V_0 < 10^{-5}$ et comparer `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC` et `VOF_TVD`.
3. **Performance de `AAL` vs `AL` et parité C++17 / Fortran 2018 :** Comparer les performances des solveurs et vérifier la parité C++/Fortran.

---

### 6.2 Résultats Attendus et Visualisation ParaView

