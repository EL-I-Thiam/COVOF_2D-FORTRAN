# Tutoriel 10 : Sillage Laminaire Stationnaire autour d'un Cylindre Immergé ($Re = 20$)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_10_sillage_cylindre/` *(validation rapide `Run_Valide-CASES`, cas 10)*.

Ce tutoriel valide la représentation d'obstacles courbes immergés sur grille cartésienne par pénalisation volumique de Darcy-Brinkman (`"type": "CYLINDER"` dans `"obstacles"`) sur l'écoulement laminaire stationnaire autour d'un cylindre circulaire à $Re = 20$.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Distinguer un obstacle rigide fixe de Darcy-Brinkman (`obstacles` avec `type: "CYLINDER"`) d'une phase fluide mobile diphasique (`initial_condition.type: "CYLINDER"`).
- Calculer les efforts hydrodynamiques par intégration volumique de la force de Brinkman : $F_x = \int_{\Omega_{\mathrm{cyl}}} (\mu / \mathcal{K}) u\,\mathrm{d}\Omega$.
- Mesurer dynamiquement la longueur de recirculation symétrique $L_r/D$, l'angle de décollement $\theta_s$ et le coefficient de traînée $C_d$.
- Valider la symétrie transverse et la convergence spatiale entre les deux moteurs MPI `covof2d_cpp` et `covof2d_fortran`.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Grandeurs de Référence

Dans le domaine $\Omega = [0, 30] \times [0, 15]\,\mathrm{m}^2$, un cylindre circulaire fixe de diamètre $D = 1.0\,\mathrm{m}$ ($R = 0.5\,\mathrm{m}$) est placé en $(x_c, y_c) = (7.5, 7.5)\,\mathrm{m}$ (`permeability: 1e-14`). Une vitesse uniforme $U_\infty = 1.0\,\mathrm{m/s}$ est prescrite à l'entrée $x=0$, avec sortie `NEUMANN` à $x=30$ et glissement libre `FREE_SLIP` sur $y=0$ et $y=15$. Avec $\rho = 1.0\,\mathrm{kg/m}^3$ et $\mu = 0.05\,\mathrm{Pa\cdot s}$ :
$$Re = \frac{\rho U_\infty D}{\mu} = 20, \qquad C_d = \frac{2 F_x}{\rho U_\infty^2 D}, \qquad C_l = \frac{2 F_y}{\rho U_\infty^2 D}.$$

Références externes de la littérature (Dennis & Chang 1970 ; Fornberg 1980 ; Linnick & Fasel 2005) : $L_r/D \approx 0.90\text{--}0.93$, $C_d \approx 2.04\text{--}2.09$ et $C_l = 0$ par symétrie géométrique.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et maillages | `mesh.{nx,ny,lx,ly}` | $L_x = 30.0\,\mathrm{m}$, $L_y = 15.0\,\mathrm{m}$ ; $200\times 100$, $400\times 200$, $800\times 400$ |
| Cylindre immergé | `obstacles[0]` | `{"type": "CYLINDER", "xc": 7.5, "yc": 7.5, "radius": 0.5, "permeability": 1e-14}` |
| Fluide et convection | `physics` | `enable_two_phase: false`, $\rho=1.0$, $\mu=0.05$, `CENTERED` vs `TVD_SUPERBEE`, `NEWTON` vs `AB2` |
| Temps | `time.{time_scheme,dt,t_final}` | `"GEAR"`, $\Delta t = 0.05\,\mathrm{s}$, $t_{\mathrm{final}} = 50.0\,\mathrm{s}$ |
| Conditions aux limites | `boundary_conditions` | `left`: `DIRICHLET (1,0)` ; `right`: `NEUMANN` ; `bottom`/`top`: `FREE_SLIP` |
| Solveurs et précond. | `solver` | `AL` ($r_{\mathrm{al}}=10$), `AAL`, `MONOLITHIC`, `SP` ; `BICGSTAB2` ; `BLOCK_TRIANGULAR` ($\mathrm{tol}=10^{-8}$, `max_iter=800`) |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_10_sillage_re20_400x200.json`
*(Pour Fortran 2018, utiliser une copie dans `case/fortran/` avec `"output_dir": "output/fortran/grid_400x200"`.)*

```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 400,
    "ny": 200,
    "lx": 30.0,
    "ly": 15.0
  },
  "physics": {
    "enable_two_phase": false,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "NEWTON",
    "fluid1": { "rho": 1.0, "mu": 0.05 },
    "gx": 0.0,
    "gy": 0.0
  },
  "obstacles": [
    {
      "type": "CYLINDER",
      "xc": 7.5,
      "yc": 7.5,
      "radius": 0.5,
      "permeability": 1.0e-14
    }
  ],
  "time": {
    "time_scheme": "GEAR",
    "dt": 0.05,
    "t_final": 50.0,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u": 1.0, "v": 0.0 },
    "right":  { "type": "NEUMANN",   "u": 0.0, "v": 0.0 },
    "bottom": { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "top":    { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 }
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
    "output_dir": "output/cpp/grid_400x200",
    "case_name": "cas_10_sillage_re20_400x200",
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
mkdir -p output/cpp/grid_400x200 output/fortran/grid_400x200

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_10_sillage_re20_400x200.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_10_sillage_re20_400x200.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification Dynamique
1. **Longueur du sillage $L_r/D$ :** Identifier sur la ligne médiane $y = y_c = 7.5\,\mathrm{m}$ en aval du bord arrière ($x > 8.0\,\mathrm{m}$) l'abscisse $x_{\mathrm{stag}}$ où $u(x, y_c)$ s'annule et calculer $L_r/D = (x_{\mathrm{stag}} - 8.0) / 1.0$.
2. **Traînée de Brinkman $C_d$ et symétrie $|C_l|$ :** Intégrer $(\mu/\mathcal{K}) u\,\Delta x\Delta y$ et $(\mu/\mathcal{K}) v\,\Delta x\Delta y$ sur les cellules contenues dans $(x-x_c)^2+(y-y_c)^2 \le R^2$ et contrôler la symétrie transverse $|u(x, y_c+\eta) - u(x, y_c-\eta)|$.
3. **Convergence et parité C++17 / Fortran 2018 :** Étudier la convergence sur $200\times 100 \to 400\times 200 \to 800\times 400$ et vérifier la cohérence inter-codes et inter-solveurs.

### 6.2 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Extraction dynamique de Lr/D et Cd de Brinkman pour le Tutoriel 10 (Sillage cylindre Re=20)."""
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

def analyze_wake(pvtr_path: str, lx: float = 30.0, ly: float = 15.0,
                 mu: float = 0.05, perm: float = 1e-14):
    p = Path(pvtr_path)
    root = ET.parse(p).getroot().find("PRectilinearGrid")
    ext = [int(x) for x in root.attrib["WholeExtent"].split()]
    nx, ny = ext[1] - ext[0], ext[3] - ext[2]
    u, v = np.zeros((ny, nx)), np.zeros((ny, nx))
    for piece in root.findall("Piece"):
        pe = [int(x) for x in piece.attrib["Extent"].split()]
        vtr = ET.parse(p.parent / piece.attrib["Source"]).getroot().find("RectilinearGrid").find("Piece")
        for da in vtr.find("CellData").findall("DataArray"):
            if da.attrib.get("Name") == "Velocity":
                arr = np.fromstring(da.text, sep=" ").reshape((pe[3]-pe[2], pe[1]-pe[0], 3))
                u[pe[2]:pe[3], pe[0]:pe[1]] = arr[:, :, 0]
                v[pe[2]:pe[3], pe[0]:pe[1]] = arr[:, :, 1]
    dx, dy = lx / nx, ly / ny
    xc = (np.arange(nx) + 0.5) * dx
    yc = (np.arange(ny) + 0.5) * dy
    X, Y = np.meshgrid(xc, yc)
    inside = (X - 7.5)**2 + (Y - 7.5)**2 <= 0.5**2
    fx = float(np.sum((mu / perm) * u[inside]) * dx * dy)
    cd = 2.0 * fx / (1.0 * 1.0 * 1.0)
    u_mid = 0.5 * (u[ny//2 - 1, :] + u[ny//2, :])
    mask = xc > 8.0
    sub_x, sub_u = xc[mask], u_mid[mask]
    idx = np.where((sub_u[:-1] <= 0.0) & (sub_u[1:] > 0.0))[0]
    lr = float("nan") if len(idx) == 0 else (sub_x[idx[0]] - sub_u[idx[0]]*(sub_x[idx[0]+1]-sub_x[idx[0]])/(sub_u[idx[0]+1]-sub_u[idx[0]]) - 8.0)
    print(f"[{pvtr_path}] Cd = {cd:.4f} | Lr/D = {lr:.4f}")

if __name__ == "__main__":
    analyze_wake("output/cpp/grid_400x200/cas_10_sillage_re20_400x200_001000.pvtr")
```

### 6.3 Visualisation sous ParaView
1. Ouvrir `output/cpp/grid_400x200/cas_10_sillage_re20_400x200_001000.pvtr` dans ParaView.
2. Superposer le champ `Pressure` et les lignes de courant (**Stream Tracer**) autour de $(7.5, 7.5)$ pour observer les deux lobes tourbillonnaires symétriques attachés en aval du cylindre.
