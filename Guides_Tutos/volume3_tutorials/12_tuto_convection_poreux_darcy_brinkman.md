# Tutoriel 12 : Convection Naturelle en Cavité Poreuse de Darcy-Brinkman ($Ra_f = 10^7, Da = 10^{-6}, Pr = 1$)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_12_convection_poreux_darcy_brinkman/` *(validation rapide `Run_Valide-CASES`, cas 12)*.

Ce tutoriel valide le couplage entre la convection naturelle de Boussinesq et la perméabilité volumique finie de Darcy-Brinkman (`obstacles` avec $\mathcal{K} > 10^{-12}\,\mathrm{m}^2$ sur l'ensemble du domaine), permettant d'explorer la transition du régime de Darcy pur ($Da \ll 1$) vers le régime de Navier-Stokes libre.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Configurer un milieu poreux isotrope saturé par une boîte poreuse recouvrant $\Omega$ (`"permeability": 1.0e-8`, $Da = \mathcal{K}/H^2 = 10^{-6}$).
- Vérifier l'activation simultanée de la résistance de Darcy $-(\mu/\mathcal{K})\mathbf{u}$, de la diffusion visqueuse de Brinkman $\mu \Delta \mathbf{u}$ et de l'équation d'énergie thermique.
- Mesurer l'ordre de convergence spatiale 2 sur $\overline{Nu}(N)$ ainsi que l'effet du nombre de Darcy $Da \in \{10^{-6}, 10^{-4}, 10^{-2}\}$.
- Valider la parité C++17 / Fortran 2018 en parallèle MPI pur (`covof2d_cpp` et `covof2d_fortran`).

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `NavierStokesSolver` (`get_*_obstacle_perm`, `solve_energy`), `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Grandeurs de Référence

Dans la cavité carrée $\Omega = [0, 0.1] \times [0, 0.1]\,\mathrm{m}^2$ ($H = 0.1\,\mathrm{m}$), la matrice poreuse uniforme possède une perméabilité $\mathcal{K} = 10^{-8}\,\mathrm{m}^2$. Avec $\rho = 1.0\,\mathrm{kg/m}^3$, $\mu = 10^{-5}\,\mathrm{Pa\cdot s}$, $C_p = 1000\,\mathrm{J/(kg\cdot K)}$, $k_{\mathrm{th}} = 0.01\,\mathrm{W/(m\cdot K)}$, $\beta = 10^{-3}\,\mathrm{K}^{-1}$, $g_y = -10.0\,\mathrm{m/s}^2$, $T_G = +0.5$ et $T_D = -0.5$ ($\Delta T = 1.0\,\mathrm{K}$, $T_{\mathrm{ref}} = 0$) :
$$Pr = \frac{\mu C_p}{k_{\mathrm{th}}} = 1.0, \qquad Da = \frac{\mathcal{K}}{H^2} = 10^{-6}, \qquad Ra_f = \frac{\rho^2 g \beta \Delta T H^3 C_p}{\mu k_{\mathrm{th}}} = 10^7, \qquad Ra^* = Ra_f \, Da = 10.$$

Référence de la littérature (Lauriat & Prasad 1989 ; Nithiarasu et al. 1997) : en régime de Darcy pur ($Da \to 0$) à $Ra^* = 10$, $\overline{Nu}_{\mathrm{th}} \approx 1.079$, valeur légèrement réduite par l'effet visqueux de paroi de Brinkman à $Da = 10^{-6}$.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et maillages | `mesh.{nx,ny,lx,ly}` | $H = 0.1\,\mathrm{m}$, grilles $32^2, 64^2, 128^2, 256^2, 512^2$ |
| Matrice poreuse | `obstacles[0]` | `{"type": "BOX", "x0": 0.0, "y0": 0.0, "x1": 0.1, "y1": 0.1, "permeability": 1e-8}` ($Da \in \{10^{-6}, 10^{-4}, 10^{-2}\}$) |
| Thermophysique | `physics` | `enable_energy: true`, `T_ref: 0.0`, $\rho=1.0, \mu=10^{-5}, C_p=1000, k_{\mathrm{th}}=0.01, \beta=10^{-3}, g_y=-10.0$ |
| Temps | `time.{time_scheme,dt,t_final}` | `"GEAR"`, $\Delta t = 0.05\,\mathrm{s}$, $t_{\mathrm{final}} = 10.0\,\mathrm{s}$ |
| Conditions aux limites | `boundary_conditions` | `u=v=0` (`DIRICHLET`) ; `left`: $T=+0.5$ ; `right`: $T=-0.5$ ; `bottom`/`top`: `T_type: "NEUMANN"` |
| Solveurs et précond. | `solver` | `AL` ($r_{\mathrm{al}}=100$), `MONOLITHIC`, `DIRECT` ; `BICGSTAB2` ; `BLOCK_TRIANGULAR` ($\mathrm{tol}=10^{-10}$, `max_iter=600`) |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_12_poreux_n128.json`
*(Pour Fortran 2018, utiliser une copie dans `case/fortran/` avec `"output_dir": "output/fortran/n128"`.)*

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
    "enable_two_phase": false,
    "enable_energy": true,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "PICARD",
    "T_ref": 0.0,
    "gx": 0.0,
    "gy": -10.0,
    "fluid1": {
      "rho": 1.0,
      "mu": 1.0e-05,
      "Cp": 1000.0,
      "k_thermal": 0.01,
      "beta": 0.001
    }
  },
  "obstacles": [
    {
      "type": "BOX",
      "x0": 0.0,
      "y0": 0.0,
      "x1": 0.1,
      "y1": 0.1,
      "permeability": 1.0e-08
    }
  ],
  "initial_condition": {
    "t_init": 0.0,
    "u0": 0.0,
    "v0": 0.0
  },
  "time": {
    "time_scheme": "GEAR",
    "dt": 0.05,
    "t_final": 10.0,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u": 0.0, "v": 0.0, "T_type": "DIRICHLET", "T":  0.5 },
    "right":  { "type": "DIRICHLET", "u": 0.0, "v": 0.0, "T_type": "DIRICHLET", "T": -0.5 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0, "T_type": "NEUMANN",   "T":  0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0, "T_type": "NEUMANN",   "T":  0.0 }
  },
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "r_al": 100.0,
    "uzawa_iter": 3,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-10,
    "max_iter": 600,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "output": {
    "output_dir": "output/cpp/n128",
    "case_name": "cas_12_poreux_n128",
    "vtk_interval": 50,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/n128 output/fortran/n128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_12_poreux_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_12_poreux_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification Dynamique
1. **Convergence spatiale d'ordre 2 sur $\overline{Nu}$ :** Extraire $\overline{Nu}(N)$ de manière identique pour C++ et Fortran depuis `Temperature` dans `.vtr`/`.pvtr`, et calculer les ordres de Richardson $p$ sur $32^2 \to 512^2$.
2. **Transition Darcy $\to$ Brinkman :** Varier `permeability` ($\mathcal{K} \in \{10^{-8}, 10^{-6}, 10^{-4}\}\,\mathrm{m}^2$, soit $Da \in \{10^{-6}, 10^{-4}, 10^{-2}\}$) et observer l'épaississement de la couche limite de Brinkman et l'augmentation de $\overline{Nu}$.
3. **Parité C++17 / Fortran 2018 :** Vérifier $|\overline{Nu}_{\mathrm{C++}} - \overline{Nu}_{\mathrm{F90}}| \ll 10^{-6}$.

### 6.2 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique de Nu moyen pour le Tutoriel 12 (Cavité poreuse Darcy-Brinkman)."""
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

def porous_nusselt(pvtr_path: str, H: float = 0.1) -> float:
    p = Path(pvtr_path)
    root = ET.parse(p).getroot().find("PRectilinearGrid")
    ext = [int(x) for x in root.attrib["WholeExtent"].split()]
    nx, ny = ext[1] - ext[0], ext[3] - ext[2]
    T = np.zeros((ny, nx))
    for piece in root.findall("Piece"):
        pe = [int(x) for x in piece.attrib["Extent"].split()]
        vtr = ET.parse(p.parent / piece.attrib["Source"]).getroot().find("RectilinearGrid").find("Piece")
        for da in vtr.find("CellData").findall("DataArray"):
            if da.attrib.get("Name") == "Temperature":
                T[pe[2]:pe[3], pe[0]:pe[1]] = np.fromstring(da.text, sep=" ").reshape((pe[3]-pe[2], pe[1]-pe[0]))
    dx = H / nx
    nu_hot = float(np.mean(H * (8.0 * 0.5 - 9.0 * T[:, 0] + T[:, 1]) / (3.0 * dx)))
    print(f"[{pvtr_path}] Nu_hot = {nu_hot:.6f}")
    return nu_hot

if __name__ == "__main__":
    porous_nusselt("output/cpp/n128/cas_12_poreux_n128_000200.pvtr")
```

### 6.3 Visualisation sous ParaView
1. Ouvrir `output/cpp/n128/cas_12_poreux_n128_000200.pvtr` dans ParaView.
2. Tracer `Temperature` (isothermes quasiment conductives à $Ra^*=10$ puis fortement convectives lorsque $Da$ augmente) et `Velocity` le long des parois verticales.
