# Tutoriel 9 : Écoulement Laminaire derrière une Marche Descendante ($Re = 100$ à $1000$)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_9_marche_descendante/` *(validation rapide `Run_Valide-CASES`, cas 9)*.

Ce tutoriel valide la pénalisation volumique de Darcy-Brinkman pour obstacle solide rectangulaire (`type: "BOX"`), la condition d'entrée parabolique décalée `POISEUILLE_STEP`, la sortie ouverte `NEUMANN` et l'évolution des zones de recirculation en fonction du nombre de Reynolds.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Modéliser une géométrie à marche descendante (*Backward-Facing Step*) sur maillage cartésien par pénalisation de Brinkman (`obstacles` avec $K = 10^{-14}\,\mathrm{m}^2$).
- Configurer le profil d'entrée dédié `"profile": "POISEUILLE_STEP"` (`y_step`, `u_mean`).
- Mesurer dynamiquement la longueur de recollement primaire $X_r/H$ sur la paroi inférieure et l'apparition du décollement secondaire supérieur $(X_s/H, X_{rs}/H)$ pour $Re \ge 400$.
- Comparer les schémas convectifs `CENTERED`, `TVD_SUPERBEE` et `UPWIND` en parallèle MPI pur.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `NavierStokesSolver` (`is_*_in_obstacle`), `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Grandeurs de Référence

Dans le canal $\Omega = [0, 40.0] \times [0, 2.0612244898]\,\mathrm{m}^2$, la marche solide occupe $\Omega_{\mathrm{step}} = [0, 4.0] \times [0, 1.0]\,\mathrm{m}^2$ (hauteur de marche $H = 1.0\,\mathrm{m}$, rapport d'expansion $L_y / (L_y - H) \approx 1.942$). En entrée ($x=0$), le fluide ($\rho = 1.0\,\mathrm{kg/m}^3$, $\bar{U} = 1.0\,\mathrm{m/s}$) entre avec le profil établi sur $y \in [H, L_y]$ :
$$u(0, y) = 6\,\bar{U}\,\frac{(y - H)(L_y - y)}{(L_y - H)^2} \quad (y \in [H, L_y]), \qquad u(0, y) = 0 \quad (y \in [0, H]).$$

Le nombre de Reynolds est défini par $Re = \frac{2 \rho \bar{U} H}{\mu} \in [100, 1000]$, soit $\mu = \frac{2}{Re}\,\mathrm{Pa\cdot s}$. Les grandeurs extraites par changement de signe du cisaillement pariétal $\tau_w(x)$ sont comparées aux résultats de référence d'Armaly et al. (1983), Gartling (1990) et Erturk (2008).

---

## 3. Paramètres Physiques, Numériques et Discrétisation

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et maillages | `mesh.{nx,ny,lx,ly}` | $L_x = 40.0\,\mathrm{m}$, $L_y = 2.0612244898\,\mathrm{m}$ ; $128\times 32$, $256\times 64$, $512\times 128$ |
| Marche immergée | `obstacles[0]` | `{"type": "BOX", "x0": 0.0, "y0": 0.0, "x1": 4.0, "y1": 1.0, "permeability": 1e-14}` |
| Fluide et convection | `physics` | `enable_two_phase: false`, $\rho=1.0$, $\mu = 2/Re$, `CENTERED`/`TVD_SUPERBEE`/`UPWIND`, `NEWTON`/`AB2`/`PICARD` |
| Temps | `time.{time_scheme,dt,t_final}` | `"GEAR"`, $\Delta t = 0.05\,\mathrm{s}$, $t_{\mathrm{final}} = 100.0\,\mathrm{s}$ |
| Conditions aux limites | `boundary_conditions` | `left`: `DIRICHLET` (`POISEUILLE_STEP`, `y_step: 1.0`, `u_mean: 1.0`) ; `right`: `NEUMANN` ; `bottom`/`top`: `DIRICHLET (0,0)` |
| Solveurs et précond. | `solver` | `AL` ($r_{\mathrm{al}}=10$), `AAL`, `MONOLITHIC`, `DIRECT` ; `BICGSTAB2` / `GMRES` ; `BLOCK_TRIANGULAR` ($\mathrm{tol}=10^{-8}$, `max_iter=800`) |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_9_marche_re400_256x64.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/re400_256x64"`.)*

```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 256,
    "ny": 64,
    "lx": 40.0,
    "ly": 2.0612244898
  },
  "physics": {
    "enable_two_phase": false,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "NEWTON",
    "fluid1": { "rho": 1.0, "mu": 0.005 },
    "gx": 0.0,
    "gy": 0.0
  },
  "obstacles": [
    {
      "type": "BOX",
      "x0": 0.0,
      "y0": 0.0,
      "x1": 4.0,
      "y1": 1.0,
      "permeability": 1.0e-14
    }
  ],
  "time": {
    "time_scheme": "GEAR",
    "dt": 0.05,
    "t_final": 100.0,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "profile": "POISEUILLE_STEP", "y_step": 1.0, "u_mean": 1.0, "u": 1.0, "v": 0.0 },
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
    "output_dir": "output/cpp/re400_256x64",
    "case_name": "cas_9_marche_re400_256x64",
    "vtk_interval": 200,
    "log_interval": 20,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/re400_256x64 output/fortran/re400_256x64

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_9_marche_re400_256x64.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_9_marche_re400_256x64.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification Dynamique
1. **Longueurs de recollement :** Déterminer par interpolation linéaire sur la première couche de mailles près de $y=0$ (pour $x > 4.0$) le zéro de $u(x, \Delta y/2)$ donnant $X_r/H = (x_{\mathrm{recol}} - 4.0)/H$, ainsi que les points de séparation/recollement sur $y = L_y$ pour $Re \ge 400$.
2. **Conservation globale du débit et rigidité de l'obstacle :** Contrôler $\int u(0, y)\,\mathrm{d}y \approx \int u(L_x, y)\,\mathrm{d}y$ et vérifier que la vitesse résiduelle dans $\Omega_{\mathrm{step}}$ est $\mathcal{O}(K) \ll 10^{-8}\,\mathrm{m/s}$.
3. **Parité C++17 / Fortran 2018 et solveurs :** Comparer $X_r/H$, les itérations Krylov et la divergence maximale entre les deux moteurs.

### 6.2 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Extraction dynamique de Xr/H pour le Tutoriel 9 (Marche descendante)."""
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

def extract_reattachment(pvtr_path: str, lx: float = 40.0, x_step: float = 4.0) -> float:
    p = Path(pvtr_path)
    root = ET.parse(p).getroot().find("PRectilinearGrid")
    ext = [int(x) for x in root.attrib["WholeExtent"].split()]
    nx = ext[1] - ext[0]
    u_bot = np.zeros(nx)
    for piece in root.findall("Piece"):
        pe = [int(x) for x in piece.attrib["Extent"].split()]
        if pe[2] == 0:
            vtr = ET.parse(p.parent / piece.attrib["Source"]).getroot().find("RectilinearGrid").find("Piece")
            for da in vtr.find("CellData").findall("DataArray"):
                if da.attrib.get("Name") == "Velocity":
                    arr = np.fromstring(da.text, sep=" ").reshape((pe[3]-pe[2], pe[1]-pe[0], 3))
                    u_bot[pe[0]:pe[1]] = arr[0, :, 0]
    dx = lx / nx
    xc = (np.arange(nx) + 0.5) * dx
    mask = xc > x_step + 0.25
    sub_x, sub_u = xc[mask], u_bot[mask]
    idx = np.where((sub_u[:-1] <= 0.0) & (sub_u[1:] > 0.0))[0]
    if len(idx) == 0:
        return float("nan")
    i = idx[0]
    xr = sub_x[i] - sub_u[i] * (sub_x[i+1] - sub_x[i]) / (sub_u[i+1] - sub_u[i])
    return (xr - x_step)

if __name__ == "__main__":
    xr = extract_reattachment("output/cpp/re400_256x64/cas_9_marche_re400_256x64_002000.pvtr")
    print(f"Xr/H = {xr:.4f}")
```

### 6.3 Visualisation sous ParaView
1. Ouvrir `output/cpp/re400_256x64/cas_9_marche_re400_256x64_002000.pvtr` dans ParaView.
2. Filtrer la marche solide par **Threshold** ou **Clip** et tracer les lignes de courant (**Stream Tracer**) derrière l'arête $(4.0, 1.0)$ pour visualiser la bulle primaire inférieure et la zone de séparation supérieure.
