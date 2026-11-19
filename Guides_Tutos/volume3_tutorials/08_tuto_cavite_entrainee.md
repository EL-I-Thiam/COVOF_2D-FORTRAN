# Tutoriel 8 : Écoulement de la Cavité Carrée Entraînée 2D à $Re = 1000$

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_8_cavite_entrainee/` *(validation rapide `Run_Valide-CASES`, cas 8)*.

Ce tutoriel valide la précision spatiale d'ordre 2 des opérateurs convectifs et visqueux incompressibles de COVOF 2D sur le benchmark classique de la cavité carrée entraînée (*Lid-Driven Cavity*) à $Re = 1000$. Il permet de comparer les schémas spatiaux (`CENTERED`, `TVD_SUPERBEE`, `UPWIND`), les linéarisations (`NEWTON`, `ADAMS_BASHFORTH`, `PICARD`) et les couplages vitesse-pression (`AL`, `AAL`, `MONOLITHIC`).

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Configurer un écoulement tourbillonnaire confiné avec cisaillement pariétal entraîné sur la paroi supérieure.
- Comparer quantitativement la dissipation numérique des schémas convectifs `CENTERED`, `TVD_SUPERBEE` et `UPWIND` câblés dans `apply_system_matrix`.
- Évaluer dynamiquement l'ordre de convergence spatiale de Richardson sur la fonction de courant $\psi(x, y)$, la vorticité $\omega_z$ et les profils médians $u(0.5, y)$, $v(x, 0.5)$.
- Exécuter le cas en parallèle MPI pur (`covof2d_cpp` et `covof2d_fortran`) et exploiter les sorties `.pvtr`/`.vtr`.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Références Bibliographiques

Sur le domaine carré $\Omega = [0, 1] \times [0, 1]\,\mathrm{m}^2$, un fluide newtonien incompressible ($\rho = 1.0\,\mathrm{kg/m}^3$, $\mu = 10^{-3}\,\mathrm{Pa\cdot s}$) est entraîné par la paroi supérieure mobile à vitesse tangentielle $u_{\mathrm{top}} = -1.0\,\mathrm{m/s}$, les parois gauche, droite et inférieure étant fixes et adhérentes ($\mathbf{u} = \mathbf{0}$) :
$$Re = \frac{\rho |u_{\mathrm{top}}| L}{\mu} = 1000.$$

La littérature de référence (Botella & Peyret 1998 ; Ghia et al. 1982, adaptée au sens $u_{\mathrm{top}} = -1$) fournit les repères suivants :
- **Tourbillon primaire central :** $\psi_{\max} \approx 0.1189366$, $\omega_z \approx 2.067753$ en $(0.4692, 0.5652)$.
- **Tourbillon secondaire inférieur gauche (aval) :** $\psi_{\min,1} \approx -1.729717\times 10^{-3}$, $\omega_z \approx -1.109789$ en $(0.1360, 0.1118)$.
- **Tourbillon secondaire inférieur droit (amont) :** $\psi_{\min,2} \approx -2.334528\times 10^{-4}$, $\omega_z \approx -0.35214$ en $(0.9167, 0.0781)$.
- **Extrema des profils médians :** $\max_y u(0.5, y) \approx 0.3885698$ vers $y \approx 0.1717$ ; sur $v(x, 0.5)$, $v_{\max} \approx 0.5270771$ vers $x \approx 0.0909$ et $v_{\min} \approx -0.3769447$ vers $x \approx 0.8422$.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et maillages | `mesh.{nx,ny,lx,ly}` | $L_x = L_y = 1.0\,\mathrm{m}$, grilles $32^2, 64^2, 128^2, 256^2, 512^2$ |
| Fluide et convection | `physics` | `enable_two_phase: false`, $\rho=1.0$, $\mu=10^{-3}$, `CENTERED` / `TVD_SUPERBEE` / `UPWIND`, `NEWTON` / `AB2` / `PICARD` |
| Temps | `time.{time_scheme,dt,t_final}` | `"GEAR"`, $\Delta t = 0.05\,\mathrm{s}$ ($0.02\,\mathrm{s}$ sur $256^2, 512^2$), $t_{\mathrm{final}} = 50.0\,\mathrm{s}$ |
| Conditions aux limites | `boundary_conditions` | `top`: `DIRICHLET` ($u = -1.0, v = 0$) ; `left`, `right`, `bottom`: `DIRICHLET` ($u=v=0$) |
| Solveurs et précond. | `solver` | `AL` ($r_{\mathrm{al}}=100$), `AAL`, `MONOLITHIC` ; `BICGSTAB2` / `GMRES` ; `BLOCK_TRIANGULAR` ($\mathrm{tol}=10^{-9}$, `max_iter=800`) |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_8_cavite_n128.json`
*(Pour Fortran 2018, utiliser une copie dans `case/fortran/` avec `"output_dir": "output/fortran/n128"`.)*

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
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "NEWTON",
    "fluid1": { "rho": 1.0, "mu": 0.001 },
    "gx": 0.0,
    "gy": 0.0
  },
  "time": {
    "time_scheme": "GEAR",
    "dt": 0.05,
    "t_final": 50.0,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u":  0.0, "v": 0.0 },
    "right":  { "type": "DIRICHLET", "u":  0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u":  0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": -1.0, "v": 0.0 }
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
    "output_dir": "output/cpp/n128",
    "case_name": "cas_8_cavite_n128",
    "vtk_interval": 200,
    "log_interval": 20,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Exécuter impérativement avec `cwd` placé dans le répertoire du cas de validation et répertoires de sortie relatifs courts pré-créés :

```bash
# 1. Créer les répertoires de sortie
mkdir -p output/cpp/n128 output/fortran/n128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_8_cavite_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_8_cavite_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification Dynamique
1. **Fonction de courant et vorticité :** Intégrer $\partial_y \psi = u$ (avec $\psi = 0$ sur la paroi inférieure) et calculer $\omega_z = \partial_x v - \partial_y u$ à partir de `Velocity` lu dans le dernier `.pvtr`/`.vtr`.
2. **Convergence spatiale d'ordre 2 :** Calculer les ordres de Richardson $p$ entre $32^2 \to 64^2 \to 128^2 \to 256^2 \to 512^2$ sur les extrema $\psi_{\max}$, $\psi_{\min,1}$ et les profils médians.
3. **Invariance et comparaison des schémas :** Quantifier l'écart entre `CENTERED`, `TVD_SUPERBEE` et `UPWIND` (notamment sur la capture des tourbillons secondaires de coin) et la parité C++17 / Fortran 2018.

### 6.2 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 8 (Cavité entraînée Re=1000) depuis .pvtr/.vtr."""
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

def read_pvtr_velocity(pvtr_path: str):
    p = Path(pvtr_path)
    root = ET.parse(p).getroot().find("PRectilinearGrid")
    ext = [int(x) for x in root.attrib["WholeExtent"].split()]
    nx, ny = ext[1] - ext[0], ext[3] - ext[2]
    u = np.zeros((ny, nx))
    v = np.zeros((ny, nx))
    for piece in root.findall("Piece"):
        pe = [int(x) for x in piece.attrib["Extent"].split()]
        vtr = ET.parse(p.parent / piece.attrib["Source"]).getroot().find("RectilinearGrid").find("Piece")
        for da in vtr.find("CellData").findall("DataArray"):
            if da.attrib.get("Name") == "Velocity":
                arr = np.fromstring(da.text, sep=" ").reshape((pe[3]-pe[2], pe[1]-pe[0], 3))
                u[pe[2]:pe[3], pe[0]:pe[1]] = arr[:, :, 0]
                v[pe[2]:pe[3], pe[0]:pe[1]] = arr[:, :, 1]
    return u, v

if __name__ == "__main__":
    u, v = read_pvtr_velocity("output/cpp/n128/cas_8_cavite_n128_001000.pvtr")
    ny, nx = u.shape
    dy = 1.0 / ny
    psi = np.cumsum(u, axis=0) * dy
    print(f"psi_max = {psi.max():.6f} | psi_min = {psi.min():.6e} | max|u(0.5,y)| = {np.max(np.abs(u[:, nx//2])):.6f}")
```

### 6.3 Visualisation sous ParaView
1. Ouvrir `output/cpp/n128/cas_8_cavite_n128_001000.pvtr` avec **File -> Open**.
2. Tracer les lignes de courant avec le filtre **Stream Tracer** (semences diagonales) et calculer la vorticité via **Compute Derivatives** pour observer le tourbillon central et les deux tourbillons secondaires de coin.
