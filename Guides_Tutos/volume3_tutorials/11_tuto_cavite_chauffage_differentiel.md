# Tutoriel 11 : Convection Naturelle Laminaire en Cavité Différentiellement Chauffée ($Pr = 0.71, Ra = 10^3$ à $10^6$)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_11_cavite_chauffage_differentiel/` *(validation rapide `Run_Valide-CASES`, cas 11)*.

Ce tutoriel valide le solveur thermique implicite (`physics.enable_energy = true`), les conditions aux limites thermiques mixtes (Dirichlet isothermes et Neumann adiabatiques) et le terme de flottabilité de Boussinesq sur le benchmark de De Vahl Davis (1983) / Le Quéré (1991).

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Renseigner avec précision les propriétés thermiques requises (`Cp`, `k_thermal`, `beta` dans `physics.fluid1`, et `T_ref` avec `T` majuscule dans `physics`).
- Paramétrer indépendamment les conditions aux limites hydrodynamiques (`DIRICHLET` adhérent) et thermiques (`T_type: "DIRICHLET"` ou `"NEUMANN"`).
- Extraire dynamiquement le nombre de Nusselt moyen $\overline{Nu}(x=0)$ et $\overline{Nu}(x=H)$ depuis `Temperature` dans `.vtr`/`.pvtr`.
- Vérifier l'ordre de convergence spatiale 2 et la centro-symétrie exacte du champ thermo-convectif pour $Ra \in \{10^3, 10^4, 10^5, 10^6\}$.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `NavierStokesSolver` (`solve_energy`, flottabilité), `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm` (`apply_thermal_bc`), `module_navier_stokes` (`ns_solve_energy`), `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Références Bibliographiques

Sur la cavité carrée $\Omega = [0, H] \times [0, H]$ ($H = 1.0\,\mathrm{m}$) à parois rigides adhérentes ($\mathbf{u} = \mathbf{0}$), le système couple Navier-Stokes incompressible avec poussée de Boussinesq $-\rho_0 \beta (T - T_{\mathrm{ref}})\mathbf{g}$ et l'équation de l'énergie :
$$\rho_0 C_p \left(\frac{\partial T}{\partial t} + \mathbf{u} \cdot \nabla T\right) = k_{\mathrm{th}} \Delta T.$$

Avec $\rho_0 = 1.0$, $\mu = 1.0$, $C_p = 0.71$, $k_{\mathrm{th}} = 1.0$, $\beta = 1.0$, $T_{\mathrm{ref}} = 0$, $T_G(0,y) = +0.5$, $T_D(H,y) = -0.5$ ($\Delta T = 1.0$) et $\partial_y T = 0$ en $y=0, H$, le nombre de Rayleigh est fixé par $g_y = -Ra / Pr$ ($Pr = \mu C_p / k_{\mathrm{th}} = 0.71$). Les valeurs de référence (De Vahl Davis 1983 ; Le Quéré 1991) sont :
$$\overline{Nu}_{\mathrm{r\acute{e}f}}(10^3) = 1.11780, \qquad \overline{Nu}_{\mathrm{r\acute{e}f}}(10^4) = 2.24475, \qquad \overline{Nu}_{\mathrm{r\acute{e}f}}(10^5) = 4.52164, \qquad \overline{Nu}_{\mathrm{r\acute{e}f}}(10^6) = 8.82513.$$

---

## 3. Paramètres Physiques, Numériques et Discrétisation

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et maillages | `mesh.{nx,ny,lx,ly}` | $H = 1.0\,\mathrm{m}$, $N \in \{16, 32, 64, 128, 256, 512\}$ |
| Thermophysique | `physics` | `enable_energy: true`, `T_ref: 0.0`, $\rho_0=1.0, \mu=1.0, C_p=0.71, k_{\mathrm{th}}=1.0, \beta=1.0$, $g_y = -Ra/0.71$ |
| Convection et temps | `physics` & `time` | `CENTERED` vs `UPWIND`, `PICARD`/`AB2`, `"GEAR"`, $\Delta t = 0.002\,\mathrm{s}$ (pour $Ra=10^4$), $t_{\mathrm{final}} = 2.0\,\mathrm{s}$ |
| Conditions aux limites | `boundary_conditions` | `u=v=0` partout ; `left`: `T_type: "DIRICHLET", T: 0.5` ; `right`: `T_type: "DIRICHLET", T: -0.5` ; `bottom`/`top`: `T_type: "NEUMANN"` |
| Solveurs et précond. | `solver` | `AL` ($r_{\mathrm{al}}=100$), `AAL`, `MONOLITHIC` ; `BICGSTAB2` / `GMRES` ; `BLOCK_TRIANGULAR` ($\mathrm{tol}=10^{-9}$, `max_iter=800`) |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_11_thermique_ra1e4_n128.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/ra1e4_n128"`.)*

```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 1.0,
    "ly": 1.0
  },
  "physics": {
    "enable_two_phase": false,
    "enable_energy": true,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "PICARD",
    "T_ref": 0.0,
    "gx": 0.0,
    "gy": -14084.507042253521,
    "fluid1": {
      "rho": 1.0,
      "mu": 1.0,
      "Cp": 0.71,
      "k_thermal": 1.0,
      "beta": 1.0
    }
  },
  "initial_condition": {
    "t_init": 0.0,
    "u0": 0.0,
    "v0": 0.0
  },
  "time": {
    "time_scheme": "GEAR",
    "dt": 0.002,
    "t_final": 2.0,
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
    "tol": 1e-09,
    "max_iter": 800,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "output": {
    "output_dir": "output/cpp/ra1e4_n128",
    "case_name": "cas_11_thermique_ra1e4_n128",
    "vtk_interval": 100,
    "log_interval": 20,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/ra1e4_n128 output/fortran/ra1e4_n128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_11_thermique_ra1e4_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_11_thermique_ra1e4_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification Dynamique
1. **Nombre de Nusselt moyen :** Calculer $\overline{Nu}(x=0) = \frac{1}{H}\int_0^H -\partial_x T(0,y)\,\mathrm{d}y$ à l'ordre 2 sur la paroi chaude (`T_0 = 0.5`) et comparer avec $\overline{Nu}(x=H)$ et la référence de De Vahl Davis / Le Quéré.
2. **Centro-symétrie et bilan énergétique :** Vérifier $\max |T(H-x, H-y) + T(x, y)| \ll 10^{-4}$ et $|\overline{Nu}_{x=0} - \overline{Nu}_{x=H}| \ll 10^{-3}$.
3. **Convergence de Richardson et parité C++17 / Fortran 2018 :** Vérifier l'ordre 2 sur $\overline{Nu}(N)$ pour $N \in \{16, 32, 64, 128, 256, 512\}$.

### 6.2 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Calcul dynamique du Nusselt moyen et de la centro-symétrie (Tutoriel 11)."""
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

def check_thermal_cavity(pvtr_path: str):
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
    dx = 1.0 / nx
    # Dérivée quadratique à la paroi gauche T_wall = 0.5 (centres en dx/2 et 3dx/2)
    nu_hot = float(np.mean((8.0 * 0.5 - 9.0 * T[:, 0] + T[:, 1]) / (3.0 * dx)))
    sym_err = float(np.max(np.abs(T + T[::-1, ::-1])))
    print(f"[{pvtr_path}] Nu_hot = {nu_hot:.6f} | sym_err(T) = {sym_err:.3e}")

if __name__ == "__main__":
    check_thermal_cavity("output/cpp/ra1e4_n128/cas_11_thermique_ra1e4_n128_001000.pvtr")
```

### 6.3 Visualisation sous ParaView
1. Ouvrir `output/cpp/ra1e4_n128/cas_11_thermique_ra1e4_n128_001000.pvtr` dans ParaView.
2. Appliquer le filtre **Contour** sur `Temperature` (10 isothermes entre `-0.5` et `+0.5`) et **Stream Tracer** sur `Velocity` pour observer la stratification thermique centrale et les couches limites thermiques pariétales.
