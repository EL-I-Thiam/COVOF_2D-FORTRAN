# Tutoriel 8 : Chute Libre d'un Cylindre Dense dans le Vide (Forts Contrastes)

Ce tutoriel valide le couplage fluide-solide par approche 1-fluide à fort contraste de masse volumique et de viscosité ($\rho_2 / \rho_1 \approx 85\,000$, $\mu_2 / \mu_1 \approx 5.4\times 10^9$). Il démontre que l'inclusion cylindrique dense suit rigoureusement la cinématique galiléenne de chute libre dans le vide sans dispersion ni amortissement numérique parasite.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Maîtriser le comportement du couplage 1-fluide Navier-Stokes/VOF face à des sauts extrêmes de densité et de viscosité.
- Configurer le schéma temporel `GEAR` (BDF2), la linéarisation `ADAMS_BASHFORTH` et le lissage SVOF de Helmholtz (`enable_svof_smoothing`).
- Comparer les formulations de couplage `MONOLITHIC` (FC), `AUGMENTED_LAGRANGIAN` (AL), `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (AAL) et `STANDARD_PROJECTION` (SP).
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Évaluer dynamiquement les erreurs $L_2$ temporelles sur $U_c(t)$ et $Y_c(t)$ depuis `*_history.csv`.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker`, `InterfaceReconstruction`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_interface_reconstruction`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas diphasique sévère simule la chute libre sous gravité d'un cylindre circulaire ultra-dense et ultra-visqueux dans l'air au repos. Il met à l'épreuve les solveurs de couplage vitesse-pression et le préconditionneur `BLOCK_TRIANGULAR` à densité variable face à des sauts $r_\rho = \rho_2/\rho_1 \approx 84\,976$ et $r_\mu = \mu_2/\mu_1 \approx 5.405\times 10^9$.

Dans une cavité rectangulaire $\Omega = [0, 0.1] \times [0, 0.2]\,\mathrm{m}^2$ soumise à la pesanteur $\mathbf{g} = (0, -9.81)\,\mathrm{m/s}^2$, un cylindre de rayon $R = 0.0125\,\mathrm{m}$ est initialement centré en $(x_c, Y_0) = (0.05, 0.15)\,\mathrm{m}$ sans vitesse initiale (`"type": "CYLINDER"`, qui affecte $C = 0$ à l'intérieur du cylindre et $C = 1$ à l'extérieur) :
- **Fluide 1 (Air ambiant, $C=1$)** : $\rho_1 = 1.1768\,\mathrm{kg/m}^3$, $\mu_1 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$.
- **Fluide 2 (Cylindre dense, $C=0$)** : $\rho_2 = 100\,000\,\mathrm{kg/m}^3$, $\mu_2 = 100\,000\,\mathrm{Pa\cdot s}$, $\sigma = 0$.

### 2.2 Solution Analytique Galiléenne et Formulations Numériques
Compte tenu du rapport d'inertie $r_\rho \gg 1$, la traînée aérodynamique et la poussée d'Archimède exercées par l'air sont négligeables devant le poids du cylindre. La vitesse verticale $U_c(t)$ et l'ordonnée du centre de masse $Y_c(t)$ obéissent à la chute libre galiléenne :
$$U_c^{\mathrm{ex}}(t) = g_y t = -9.81\,t, \qquad Y_c^{\mathrm{ex}}(t) = Y_0 + \frac{1}{2} g_y t^2 = 0.15 - 4.905\,t^2.$$
À l'instant final $t_{\mathrm{final}} = 0.144\,\mathrm{s}$, la cinématique théorique vaut :
$$U_c^{\mathrm{ex}}(0.144\,\mathrm{s}) = -1.41264\,\mathrm{m/s}, \qquad Y_c^{\mathrm{ex}}(0.144\,\mathrm{s}) = 0.048288\,\mathrm{m}.$$

Les briques numériques clés activées pour ce cas sont :
1. **Schéma temporel Gear / BDF2** (`"time_scheme": "GEAR"`) et linéarisation convective d'Adams-Bashforth (`"convective_linearization": "ADAMS_BASHFORTH"`).
2. **Lissage SVOF de Helmholtz** (`(I - k\Delta x\Delta y \Delta)C_s = C`) avec moyenne harmonique nodale de la viscosité de cisaillement.
3. **Couplages comparés :** `MONOLITHIC` (FC), `AUGMENTED_LAGRANGIAN` (AL), `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (AAL avec facteur `k_aal` et critère `uzawa_tol`) et `STANDARD_PROJECTION` (SP).

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et Cylindre | $\Omega$, $(x_c, Y_0)$, $R$ | $[0, 0.1]\times[0, 0.2]\,\mathrm{m}^2$, $(0.05, 0.15)\,\mathrm{m}$, $R = 0.0125\,\mathrm{m}$ (`axis_a = axis_b = 0.0125`, `pert_amp = 0.0`) |
| Maillages et Pas de temps | $(N_x \times N_y)$, $\Delta t$, $t_{\mathrm{final}}$ | $100\times 200$ et $200\times 400$, $\Delta t = 6.25\times 10^{-5}\,\mathrm{s}$, $t_{\mathrm{final}} = 0.144\,\mathrm{s}$ (2304 pas) |
| Conditions aux limites | Gauche/Droite ; Bas/Haut | `FREE_SLIP` sur $x=0, 0.1$ ; `DIRICHLET` ($u=v=0$) sur $y=0, 0.2$ |
| Solveurs et Préconditionneurs | Couplage, Krylov, Précond. | `MONOLITHIC` (FC), `AUGMENTED_LAGRANGIAN` (AL), `AAL` (`k_aal` $\in \{10, 100, 1000\}$), `STANDARD_PROJECTION` ; `BICGSTAB2`, `GMRES`, `BICGSTAB` ; `BLOCK_TRIANGULAR`, `JACOBI`, `SSOR`, `ILU0` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_8_fc_100x200.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/fc_100x200"`.)*

```json
{
  "case_name": "cas_8_fc_100x200",
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
    "convective_scheme": "UPWIND",
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
    "radius": 0.0125,
    "axis_a": 0.0125,
    "axis_b": 0.0125,
    "pert_amp": 0.0
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
    "output_dir": "output/cpp/fc_100x200",
    "case_name": "cas_8_fc_100x200",
    "vtk_interval": 100,
    "log_interval": 100,
    "write_csv_history": true
  }
}
```

Pour activer le Lagrangien Augmenté Adaptatif (`AAL`), remplacer le bloc `"solver"` par :

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

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

L'exécution utilise **MPI pur** (`mpirun -np 4` pour `MONOLITHIC`, `AL` et `AAL` avec solveurs Krylov ; `mpirun -np 1` si `STANDARD_PROJECTION` est sélectionné).

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/fc_100x200 output/fortran/fc_100x200

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_8_fc_100x200.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_8_fc_100x200.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Erreurs relatives $L_2$ temporelles sur $U_c(t)$ et $Y_c(t)$ :**
   À partir des colonnes `time`, `center_of_mass_y` ($Y_c$) et `rise_velocity` ($U_c$) de `<case_name>_history.csv` (car pour `IC_CYLINDER`, $1 - C$ représente le cylindre dense), calculer dynamiquement :
   $$\mathrm{Err}_{\mathrm{rel},2}(U_c) = \sqrt{\frac{\sum_{n=1}^{N_t} \big(U_c(t^n) - U_c^{\mathrm{ex}}(t^n)\big)^2}{\sum_{n=1}^{N_t} \big(U_c^{\mathrm{ex}}(t^n)\big)^2}}, \qquad \mathrm{Err}_{\mathrm{rel},2}(Y_c) = \sqrt{\frac{\sum_{n=1}^{N_t} \big(Y_c(t^n) - Y_c^{\mathrm{ex}}(t^n)\big)^2}{\sum_{n=1}^{N_t} \big(Y_c^{\mathrm{ex}}(t^n)\big)^2}}.$$
   Critères cibles : $\mathrm{Err}_{\mathrm{rel},2}(U_c) < 10^{-2}$, $\mathrm{Err}_{\mathrm{rel},2}(Y_c) < 10^{-2}$, et variation relative de volume $|V(t_{\mathrm{final}}) - V(0)| / V(0) < 10^{-4}$ mesurée depuis `PhaseIndicator` dans les `.vtr`.
2. **Comparaison des couplages (`FC`, `AL`, `AAL`, `SP`) et du filtre SVOF :**
   Comparer la préservation de la circularité du cylindre, l'absence de traînée parasite et le nombre moyen d'itérations Krylov relevé dans les journaux d'exécution.
3. **Parité C++17 / Fortran 2018 :**
   Vérifier la superposition des trajectoires $Y_c(t)$ et des vitesses $U_c(t)$ entre `covof2d_cpp` et `covof2d_fortran` ainsi que le maintien de $\|\nabla\cdot\mathbf{u}\|_\infty$ (`max_div`).

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` : `step,time,dt,center_of_mass_y,rise_velocity,bubble_volume,kinetic_energy,max_u,max_div` (`9` colonnes C++ / `11` colonnes Fortran avec `avg_rho_gas,avg_p_gas`).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : `PhaseIndicator`, `PhaseIndicator_Smoothed`, `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Velocity`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 8 (Chute libre cylindre dense) depuis *_history.csv."""
import csv
import math
from pathlib import Path

def verify_free_fall(csv_path: str, y0: float = 0.15, gy: float = -9.81) -> None:
    times, yc, uc = [], [], []
    with Path(csv_path).open(newline="") as f:
        for r in csv.DictReader(f):
            t = float(r["time"])
            if t > 0.0:
                times.append(t)
                yc.append(float(r["center_of_mass_y"]))
                uc.append(float(r["rise_velocity"]))

    u_ex = [gy * t for t in times]
    y_ex = [y0 + 0.5 * gy * t * t for t in times]
    err_u = math.sqrt(sum((a - b) ** 2 for a, b in zip(uc, u_ex)) / sum(b ** 2 for b in u_ex))
    err_y = math.sqrt(sum((a - b) ** 2 for a, b in zip(yc, y_ex)) / sum(b ** 2 for b in y_ex))
    print(f"[{csv_path}] t_end={times[-1]:.4f} s | Uc={uc[-1]:.6f} (ex={u_ex[-1]:.6f}) "
          f"| Yc={yc[-1]:.6f} (ex={y_ex[-1]:.6f}) | Err_L2(U)={err_u:.3e} | Err_L2(Y)={err_y:.3e}")

if __name__ == "__main__":
    verify_free_fall("output/cpp/fc_100x200/cas_8_fc_100x200_history.csv")
    verify_free_fall("output/fortran/fc_100x200/cas_8_fc_100x200_history.csv")
```

### 6.4 Visualisation sous ParaView
1. Ouvrir la série temporelle `<case_name>_..pvtr` dans ParaView (**File -> Open**).
2. Afficher le champ `Density` en échelle logarithmique pour distinguer le cylindre dense ($\rho_2 = 10^5\,\mathrm{kg/m}^3$) de l'air ambiant ($\rho_1 = 1.1768\,\mathrm{kg/m}^3$), et tracer l'isocontour `PhaseIndicator = 0.5` afin de vérifier l'absence de déformation du solide rigide durant sa chute.
3. Utiliser le filtre **Glyph** sur `Velocity` pour visualiser l'écoulement d'air déplacé autour du cylindre.
