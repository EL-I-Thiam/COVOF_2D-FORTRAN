# Tutoriel 14 : Campagne Dual-Code Newton, AAL et SVOF sur Cylindre Dense ($50 \times 100$)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_CASES/validation_cas_2_cylindre_dense_newton_aal/` *(campagne avancée `Run_CASES`, cas 2)*.

Ce tutoriel constitue une campagne paramétrique rapide ($50\times 100$, $\Delta x = \Delta y = 0.002\,\mathrm{m}$) dédiée à l'analyse croisée de la linéarisation convective de Newton (`NEWTON` vs `ADAMS_BASHFORTH` vs `PICARD`), du Lagrangien Augmenté Adaptatif (`AAL` vs `MONOLITHIC`), du schéma temporel (`GEAR` vs `EULER`) et du lissage SVOF ($\sigma = 0$) entre `covof2d_cpp` et `covof2d_fortran`.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comparer les trois linéarisations convectives (`NEWTON`, `ADAMS_BASHFORTH`, `PICARD`) et les trois schémas spatiaux (`CENTERED`, `UPWIND`, `TVD_SUPERBEE`) sur un saut $r_\rho \approx 85\,000$, $r_\mu \approx 5.4\times 10^9$.
- Quantifier le comportement de `AAL` ($K_{\mathrm{aal}} \in \{10, 100\}$, `uzawa_iter` $\in \{1, 3, 5\}$) et `MONOLITHIC` avec le préconditionneur `BLOCK_TRIANGULAR` 100 % parallèle MPI.
- Mesurer l'impact rhéologique pur du lissage SVOF de Helmholtz (`enable_svof_smoothing`) lorsque $\sigma = 0$ (`enable_diffusion_smoothing: false`).
- Vérifier la parité stricte C++17 $\leftrightarrow$ Fortran 2018 sur la trajectoire galiléenne $U_c(t) = -gt$, $Y_c(t) = Y_0 - \frac{1}{2}gt^2$.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker`, `InterfaceReconstruction`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_interface_reconstruction`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Matrice des Paramètres

Sur $\Omega = [0, 0.1] \times [0, 0.2]\,\mathrm{m}^2$ discrétisé en $50\times 100$ mailles ($\Delta t = 6.25\times 10^{-5}\,\mathrm{s}$, $t_{\mathrm{final}} = 0.144\,\mathrm{s}$), le cylindre de rayon $R = 0.0125\,\mathrm{m}$ centré en $(0.05, 0.15)\,\mathrm{m}$ (`axis_a = axis_b = 0.0125`, `pert_amp = 0.0`) chute librement dans l'air ($\rho_1 = 1.1768$, $\mu_1 = 1.85\times 10^{-5}$, $\rho_2 = 10^5$, $\mu_2 = 10^5$, $\sigma = 0$). La solution analytique à $t_{\mathrm{final}} = 0.144\,\mathrm{s}$ vaut $U_c^{\mathrm{ex}} = -1.41264\,\mathrm{m/s}$ et $Y_c^{\mathrm{ex}} = 0.048288\,\mathrm{m}$.

| Catégorie | Clé JSON | Plage / Options étudiées |
| :--- | :--- | :--- |
| Convection | `physics.convective_scheme` / `convective_linearization` | `CENTERED`, `UPWIND`, `TVD_SUPERBEE` ; `NEWTON`, `ADAMS_BASHFORTH`, `PICARD` |
| Temps & SVOF | `time.time_scheme` ; `physics.enable_svof_smoothing` | `GEAR` vs `EULER` ; `true` ($k = 1.0$) vs `false` |
| Solveurs & Précond. | `solver` | `AAL` ($K_{\mathrm{aal}} \in \{10, 100\}$, `uzawa_iter` $\in \{1, 3, 5\}$), `MONOLITHIC` ; `BLOCK_TRIANGULAR`, `JACOBI`, `SSOR`, `ILU0` ($\mathrm{tol}=10^{-6}$, `max_iter=500`) |

---

## 3. Fichier de Configuration JSON Pas à Pas (`case/cpp/cas_2_aal_gear_newton_svof.json`)
*(Pour Fortran 2018, utiliser `"output_dir": "output/fortran/aal_gear_newton_svof"`.)*

```json
{
  "mode": "NAVIER_STOKES",
  "mesh": { "nx": 50, "ny": 100, "lx": 0.1, "ly": 0.2 },
  "physics": {
    "enable_two_phase": true,
    "enable_surface_tension": false,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "NEWTON",
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
  "time": { "time_scheme": "GEAR", "dt": 6.25e-05, "t_final": 0.144, "cfl_max": 0.5, "adaptive_dt": false },
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
    "uzawa_tol": 1e-06,
    "uzawa_iter": 3,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-06,
    "max_iter": 500,
    "r_al": 0.0,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "vof": { "advection_method": "WEYMOUTH_YUE", "reconstruction_method": "PLIC", "tvd_limiter": "NONE", "cfl_vof": 0.5, "use_markers": false },
  "output": { "output_dir": "output/cpp/aal_gear_newton_svof", "case_name": "cas_2_aal_gear_newton_svof", "vtk_interval": 100, "log_interval": 100, "write_csv_history": true }
}
```

---

## 4. Procédure d'Exécution et Vérification Dynamique Dual-Code

```bash
mkdir -p output/cpp/aal_gear_newton_svof output/fortran/aal_gear_newton_svof
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_2_aal_gear_newton_aal.json
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_2_aal_gear_newton_aal.json
```

```python
#!/usr/bin/env python3
"""Vérification de parité C++/Fortran pour le Tutoriel 14 (Cylindre dense 50x100)."""
import csv, math
from pathlib import Path

def eval_csv(path: str, y0: float = 0.15, gy: float = -9.81):
    with Path(path).open(newline="") as f:
        rows = [r for r in csv.DictReader(f) if float(r["time"]) > 0.0]
    t = [float(r["time"]) for r in rows]
    uc = [float(r["rise_velocity"]) for r in rows]
    yc = [float(r["center_of_mass_y"]) for r in rows]
    err_u = math.sqrt(sum((a - gy*ti)**2 for a, ti in zip(uc, t)) / sum((gy*ti)**2 for ti in t))
    err_y = math.sqrt(sum((a - (y0 + 0.5*gy*ti*ti))**2 for a, ti in zip(yc, t)) / sum((y0 + 0.5*gy*ti*ti)**2 for ti in t))
    print(f"[{path}] Err_L2(U)={err_u:.3e} | Err_L2(Y)={err_y:.3e} | Uc_end={uc[-1]:.6f}")

if __name__ == "__main__":
    eval_csv("output/cpp/aal_gear_newton_svof/cas_2_aal_gear_newton_svof_history.csv")
    eval_csv("output/fortran/aal_gear_newton_svof/cas_2_aal_gear_newton_svof_history.csv")
```
