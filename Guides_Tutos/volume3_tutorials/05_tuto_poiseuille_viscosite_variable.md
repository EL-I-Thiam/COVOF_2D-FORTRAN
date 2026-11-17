# Tutoriel 5 : Écoulement de Poiseuille Plan à Viscosité Variable Exponentielle

Ce tutoriel explore l'effet d'une hétérogénéité spatiale continue de la viscosité dynamique $\mu(y) = \mu_0 \exp\big(a (y - e)^2 / e^2\big)$ sur un écoulement de Poiseuille. Il valide l'évaluation conservative des contraintes visqueuses aux nœuds de vorticité sur grille MAC ainsi que le comportement des solveurs et préconditionneurs à coefficients variables.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre l'influence d'une viscosité variable spatialement sur le profil d'équilibre de quantité de mouvement.
- Configurer le modèle `"viscosity_model": "EXPONENTIAL"` (`visc_mu0`, `visc_a`, `visc_e`) dans COVOF 2D.
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Quantifier dynamiquement l'ordre de convergence spatiale $\mathcal{O}(\Delta y^2)$ et la linéarité de la contrainte de cisaillement $\tau_{xz}(z) = -G_x z$.
- Visualiser les champs de viscosité, vitesse et cisaillement sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker` (`update_physical_properties`), `NavierStokesSolver`, `BiCGSTAB2Solver`, `DirectSolver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas vérifie la discrétisation conservative du terme de diffusion visqueuse $\nabla \cdot \big[\mu(\mathbf{x})(\nabla \mathbf{u} + \nabla \mathbf{u}^T)\big]$ lorsque la viscosité dynamique varie fortement dans l'espace.

Dans un canal plan $\Omega = [0, L_x] \times [0, L_y]$ avec $L_x = 0.05\,\mathrm{m}$ et $L_y = 2e = 0.02\,\mathrm{m}$ (demi-largeur $e = 0.01\,\mathrm{m}$, coordonnée transversale centrée $z = y - e \in [-e, +e]$), la viscosité dynamique suit une loi exponentielle symétrique :
$$\mu(z) = \mu_0 \exp\left(a \frac{z^2}{e^2}\right), \qquad \text{avec } \mu_0 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}, \quad a = 2.0, \quad e = 0.01\,\mathrm{m}.$$
Le rapport de viscosité entre les parois ($|z|=e$) et l'axe central ($z=0$) vaut $\mu(\pm e)/\mu_0 = \mathrm{e}^a = \mathrm{e}^2 \approx 7.3891$. Sous un forçage axial constant $G_x = \rho g_x = 0.294207317\,\mathrm{Pa/m}$ ($\rho = 1.0\,\mathrm{kg/m}^3$), l'équation d'équilibre stationnaire s'écrit :
$$\frac{\mathrm{d}}{\mathrm{d}z}\left(\mu(z)\frac{\mathrm{d}u}{\mathrm{d}z}\right) = -G_x, \qquad u(-e) = u(+e) = 0.$$

### 2.2 Solution Analytique Exacte
Une première intégration avec la condition de symétrie $\frac{\mathrm{d}u}{\mathrm{d}z}(0) = 0$ fournit la contrainte de cisaillement linéaire :
$$\tau_{xz}(z) = \mu(z)\frac{\mathrm{d}u}{\mathrm{d}z} = -G_x z.$$
En intégrant une seconde fois avec $u(\pm e) = 0$, on obtient le profil analytique gaussien tronqué :
$$u_{\mathrm{ex}}(z) = \frac{G_x e^2}{2 a \mu_0}\left[\exp\left(-a\frac{z^2}{e^2}\right) - \exp(-a)\right],$$
de vitesse maximale au centre $u_{\max} = u_{\mathrm{ex}}(0) = \frac{G_x e^2}{2 a \mu_0}(1 - \mathrm{e}^{-a}) \approx 0.343771\,\mathrm{m/s}$ pour $a = 2.0$. Sur le maillage décalé MAC, les flux de cisaillement aux faces Nord et Sud des volumes de contrôle de $u_{i,j}$ utilisent la viscosité évaluée aux nœuds de vorticité (`mu_node`), assurant une précision formelle d'ordre 2 ($\mathcal{O}(\Delta y^2)$).

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------------- | :----------------------- | :------------------------------------------------------------------ |
| Domaine | $[0, L_x] \times [0, L_y]$ | $L_x = 0.05\,\mathrm{m}$, $L_y = 0.02\,\mathrm{m}$ ($e = 0.01\,\mathrm{m}$) |
| Loi de viscosité | `"viscosity_model": "EXPONENTIAL"` | $\mu_0 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$, $a = 2.0$, $e = 0.01\,\mathrm{m}$, $\rho = 1.0\,\mathrm{kg/m}^3$ |
| Forçage axial | $g_x, g_y$ | $g_x = 0.294207317\,\mathrm{m/s}^2$, $g_y = 0.0$ |
| Séquence de grilles | $N_x \times N_y$ | $N_x = 16$, $N_y \in \{16, 32, 64, 128, 256, 512\}$ |
| Intégration temporelle | $\Delta t$, $t_{\mathrm{final}}$, `steady_state_tol` | $\Delta t = 0.5\,\mathrm{s}$, $t_{\mathrm{final}} = 50.0\,\mathrm{s}$, $\mathrm{steady\_state\_tol} = 10^{-12}$ |
| Solveurs testés | Formulations & Krylov | `AUGMENTED_LAGRANGIAN`, `AAL`, `MONOLITHIC` avec `BICGSTAB2` (`BLOCK_TRIANGULAR`, `SSOR`, `ILU0`) ou `DIRECT` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_5_var_visc_al_ny128.json`
*(Pour Fortran 2018, utiliser une copie dans `case/fortran/` avec `"output_dir": "output/fortran/al_ny128"`.)*

```json
{
  "case_name": "cas_5_var_visc_al_ny128",
  "mode": "NAVIER_STOKES",
  "test_case": "CAS_5_VAR_MU",
  "steady_state_tol": 1e-12,
  "mesh": {
    "nx": 16,
    "ny": 128,
    "lx": 0.05,
    "ly": 0.02
  },
  "physics": {
    "enable_two_phase": false,
    "enable_surface_tension": false,
    "enable_convection": false,
    "viscosity_model": "EXPONENTIAL",
    "visc_a": 2.0,
    "visc_e": 0.01,
    "visc_mu0": 1.85e-05,
    "gx": 0.294207317,
    "gy": 0.0,
    "fluid1": {
      "rho": 1.0,
      "mu": 1.85e-05
    }
  },
  "time": {
    "time_scheme": "EULER",
    "dt": 0.5,
    "t_final": 50.0,
    "cfl_max": 0.5,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "PERIODIC",  "u": 0.0, "v": 0.0 },
    "right":  { "type": "PERIODIC",  "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "linear_solver": "BICGSTAB2",
    "tol": 1e-12,
    "max_iter": 5000,
    "r_al": 100.0,
    "uzawa_iter": 1,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR",
    "steady_state_tol": 1e-12
  },
  "output": {
    "output_dir": "output/cpp/al_ny128",
    "case_name": "cas_5_var_visc_al_ny128",
    "vtk_interval": 200,
    "log_interval": 50,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs utilisent **MPI pur** (`mpirun -np 4`) pour les solveurs itératifs de Krylov, et `mpirun -np 1` lorsque `"linear_solver": "DIRECT"` est sélectionné.

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/al_ny128 output/fortran/al_ny128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_5_var_visc_al_ny128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_5_var_visc_al_ny128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Convergence spatiale d'ordre 2 face à $u_{\mathrm{ex}}(z)$ :**
   Extraire le profil stationnaire $u(y_j)$ depuis le dernier fichier `.pvtr`/`.vtr` et calculer dynamiquement $E_{L_1}$, $E_{L_2}$, $E_{L_\infty}$ et l'ordre de Richardson sur la séquence $N_y \in \{16, 32, 64, 128, 256, 512\}$. Vérifier $p_{L_1}, p_{L_2}, p_{L_\infty} \in [1.95, 2.05]$.
2. **Linéarité de la contrainte de cisaillement et sensibilité à $a$ (`visc_a`) :**
   Vérifier que $\tau_{xz}(y_{j+1/2}) = \mu(y_{j+1/2})\frac{u_{j+1}-u_j}{\Delta y}$ coïncide avec $-G_x(y_{j+1/2} - e)$. Pour $a \to 0$, vérifier que la solution rejoint le profil parabolique classique à viscosité constante.
3. **Parité C++17 / Fortran 2018 :**
   Vérifier depuis `<case_name>_history.csv` et le dernier `.pvtr` que $|u_{\max}^{\mathrm{C++}} - u_{\max}^{\mathrm{Fortran}}| < 10^{-12}\,\mathrm{m/s}$ et $\|\mathbf{u}^{\mathrm{C++}} - \mathbf{u}^{\mathrm{Fortran}}\|_\infty < 10^{-12}\,\mathrm{m/s}$.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` : `step,time,dt,center_of_mass_y,rise_velocity,bubble_volume,kinetic_energy,max_u,max_div` (`9` colonnes C++ / `11` colonnes Fortran avec `avg_rho_gas,avg_p_gas`).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : contiennent `PhaseIndicator`, `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Velocity`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 5 (Viscosité variable exponentielle) depuis .vtr."""
import math
import xml.etree.ElementTree as ET
from pathlib import Path

def load_pvtr_velocity_x(pvtr_path: str):
    pvtr = Path(pvtr_path)
    root = ET.parse(pvtr).getroot()
    ys_all, ux_all = [], []
    for piece in root.findall(".//Piece"):
        vtr = ET.parse(pvtr.parent / piece.attrib["Source"]).getroot()
        ext = [int(v) for v in piece.attrib["Extent"].split()]
        nx_l, ny_l = ext[1] - ext[0], ext[3] - ext[2]
        y_coords = [float(v) for v in vtr.find(".//DataArray[@Name='YCoordinates']").text.split()]
        vel = [float(v) for v in vtr.find(".//DataArray[@Name='Velocity']").text.split()]
        for j in range(ny_l):
            yc = 0.5 * (y_coords[j] + y_coords[j + 1])
            u_row = [vel[3 * (j * nx_l + i)] for i in range(nx_l)]
            ys_all.append(yc)
            ux_all.append(sum(u_row) / nx_l)
    return ys_all, ux_all

def check_exact(pvtr_path: str, gx: float = 0.294207317, mu0: float = 1.85e-5, a: float = 2.0, e: float = 0.01):
    ys, u_num = load_pvtr_velocity_x(pvtr_path)
    errs = []
    for y, un in zip(ys, u_num):
        z = y - e
        u_ex = (gx * e * e / (2.0 * a * mu0)) * (math.exp(-a * (z / e) ** 2) - math.exp(-a))
        errs.append(abs(un - u_ex))
    l2 = math.sqrt(sum(v * v for v in errs) / len(errs))
    print(f"[{pvtr_path}] N_samples={len(errs)} | L1={sum(errs)/len(errs):.4e} | L2={l2:.4e} | Linf={max(errs):.4e}")

if __name__ == "__main__":
    p = sorted(Path("output/cpp/al_ny128").glob("*.pvtr"))[-1]
    check_exact(str(p))
```

### 6.4 Visualisation sous ParaView
1. Ouvrir le dernier fichier `<case_name>_<step:06d>.pvtr` dans ParaView (**File -> Open**).
2. Inspecter le champ `Viscosity` pour visualiser le gradient exponentiel symétrique entre le cœur du canal ($\mu_0$) et les parois ($\mu_0 \mathrm{e}^a$).
3. Appliquer le filtre **Plot Over Line** transversalement en $x = 0.025\,\mathrm{m}$ sur `Velocity` (composante $X$) pour observer l'aplatissement au centre et le resserrement du profil gaussien tronqué près des parois.
