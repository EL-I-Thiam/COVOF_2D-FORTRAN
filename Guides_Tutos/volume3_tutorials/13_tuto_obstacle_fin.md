# Tutoriel 13 : Écoulement Laminaire autour d'un Obstacle Rectangulaire Fin en Canal ($Re = 82.5$)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_CASES/validation_cas_1_obstacle_fin/` *(campagne avancée `Run_CASES`, cas 1)*.

Ce tutoriel met en œuvre la pénalisation rigide de Darcy-Brinkman sur un obstacle mince vertical fixe en canal (`BOX`, $e = 2\,\mathrm{mm}$, $S/H = 0.514$) soumis à une entrée parabolique de Poiseuille (`profile: "PARABOLIC"`).

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Configurer un obstacle mince rigide de Brinkman (`obstacles` avec `permeability: 1e-14`) et analyser la sensibilité à $\mathcal{K} \in \{10^{-10}, 10^{-12}, 10^{-14}, 10^{-16}\}\,\mathrm{m}^2$.
- Paramétrer le profil d'entrée parabolique plein canal (`"profile": "PARABOLIC"`, `u_mean`, `u`) couplé à une sortie `NEUMANN`.
- Extraire dynamiquement les longueurs de recirculation aval $x_1/S$ (paroi inférieure) et $(x_2/S, x_3/S)$ (paroi supérieure), la perte de charge $\Delta P$ et les 9 profils axiaux $u(x^*, y/H)$.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Grandeurs de Référence

Dans le canal $\Omega = [0, 0.50] \times [0, 0.025]\,\mathrm{m}^2$ ($H = 0.025\,\mathrm{m}$), un obstacle mince de hauteur $S = 0.01285\,\mathrm{m}$ et d'épaisseur $e = 0.002\,\mathrm{m}$ est posé en $x_{\mathrm{obs}} = 0.15\,\mathrm{m}$ (`[0.149, 0.151] x [0.0, 0.01285]`). Le fluide ($\rho = 1.0\,\mathrm{kg/m}^3$, $\mu = 1.557576\times 10^{-4}\,\mathrm{Pa\cdot s}$) entre avec $u_{\mathrm{mean}} = 1.0\,\mathrm{m/s}$ ($u_{\max} = 1.5\,\mathrm{m/s}$), donnant $Re = \rho \bar{U} S / \mu = 82.5$. La littérature de référence rapporte un tourbillon primaire aval en paroi inférieure ($x_1/S \approx 4.99$) et un décollement supérieur entre $x_2/S \approx 4.01$ et $x_3/S \approx 9.19$.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et maillages | `mesh.{nx,ny,lx,ly}` | $L_x = 0.50\,\mathrm{m}$, $L_y = 0.025\,\mathrm{m}$ ; $250\times 50$, $500\times 100$, $1000\times 200$ |
| Obstacle mince | `obstacles[0]` | `{"type": "BOX", "x0": 0.149, "y0": 0.0, "x1": 0.151, "y1": 0.01285, "permeability": 1e-14}` |
| Fluide et convection | `physics` | `enable_two_phase: false`, $\rho=1.0$, $\mu = 1.557576\times 10^{-4}$, `CENTERED` vs `TVD_SUPERBEE` vs `UPWIND`, `NEWTON` vs `AB2` |
| Temps | `time.{time_scheme,dt,t_final}` | `"GEAR"`, $\Delta t = 5\times 10^{-4}\,\mathrm{s}$, $t_{\mathrm{final}} = 2.0\,\mathrm{s}$ |
| Conditions aux limites | `boundary_conditions` | `left`: `DIRICHLET` (`PARABOLIC`, `u_mean: 1.0`, `u: 1.5`) ; `right`: `NEUMANN` ; `bottom`/`top`: `DIRICHLET (0,0)` |
| Solveurs et précond. | `solver` | `AL` ($r_{\mathrm{al}}=10$), `AAL`, `MONOLITHIC` ; `BICGSTAB2` ; `BLOCK_TRIANGULAR` ($\mathrm{tol}=10^{-8}$, `max_iter=800`) |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_1_obstacle_fin_500x100.json`
*(Pour Fortran 2018, utiliser `"output_dir": "output/fortran/grid_500x100"`.)*

```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 500,
    "ny": 100,
    "lx": 0.50,
    "ly": 0.025
  },
  "physics": {
    "enable_two_phase": false,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "NEWTON",
    "fluid1": { "rho": 1.0, "mu": 0.0001557576 },
    "gx": 0.0,
    "gy": 0.0
  },
  "obstacles": [
    {
      "type": "BOX",
      "x0": 0.149,
      "y0": 0.0,
      "x1": 0.151,
      "y1": 0.01285,
      "permeability": 1.0e-14
    }
  ],
  "time": {
    "time_scheme": "GEAR",
    "dt": 0.0005,
    "t_final": 2.0,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "profile": "PARABOLIC", "u_mean": 1.0, "u": 1.5, "v": 0.0 },
    "right":  { "type": "NEUMANN",   "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "r_al": 10.0,
    "uzawa_iter": 3,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-08,
    "max_iter": 800,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "output": {
    "output_dir": "output/cpp/grid_500x100",
    "case_name": "cas_1_obstacle_fin_500x100",
    "vtk_interval": 200,
    "log_interval": 20,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

```bash
mkdir -p output/cpp/grid_500x100 output/fortran/grid_500x100
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_1_obstacle_fin_500x100.json
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_1_obstacle_fin_500x100.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

1. **Longueurs $x_1/S, x_2/S, x_3/S$ :** Détecter dynamiquement l'inversion de $u(x, \Delta y/2)$ et $u(x, H - \Delta y/2)$ pour $x > 0.151\,\mathrm{m}$ à partir du dernier `.pvtr`.
2. **Rigidité Brinkman et perte de charge :** Vérifier $\max |u|$ dans l'obstacle et $\Delta P = \bar{p}(0.05) - \bar{p}(0.45)$.

```python
#!/usr/bin/env python3
"""Extraction dynamique de x1/S sur le Tutoriel 13 (Obstacle fin en canal)."""
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

def check_thin_obstacle(pvtr_path: str, lx: float = 0.50, x_obs: float = 0.151, S: float = 0.01285):
    p = Path(pvtr_path)
    root = ET.parse(p).getroot().find("PRectilinearGrid")
    ext = [int(x) for x in root.attrib["WholeExtent"].split()]
    nx, ny = ext[1] - ext[0], ext[3] - ext[2]
    u = np.zeros((ny, nx))
    for piece in root.findall("Piece"):
        pe = [int(x) for x in piece.attrib["Extent"].split()]
        vtr = ET.parse(p.parent / piece.attrib["Source"]).getroot().find("RectilinearGrid").find("Piece")
        for da in vtr.find("CellData").findall("DataArray"):
            if da.attrib.get("Name") == "Velocity":
                u[pe[2]:pe[3], pe[0]:pe[1]] = np.fromstring(da.text, sep=" ").reshape((pe[3]-pe[2], pe[1]-pe[0], 3))[:, :, 0]
    xc = (np.arange(nx) + 0.5) * (lx / nx)
    m = xc > x_obs + 0.002
    sub_x, sub_u = xc[m], u[0, m]
    idx = np.where((sub_u[:-1] <= 0.0) & (sub_u[1:] > 0.0))[0]
    x1_s = float("nan") if len(idx) == 0 else (sub_x[idx[0]] - sub_u[idx[0]]*(sub_x[idx[0]+1]-sub_x[idx[0]])/(sub_u[idx[0]+1]-sub_u[idx[0]]) - 0.15) / S
    print(f"[{pvtr_path}] x1/S = {x1_s:.4f}")

if __name__ == "__main__":
    check_thin_obstacle("output/cpp/grid_500x100/cas_1_obstacle_fin_500x100_004000.pvtr")
```
