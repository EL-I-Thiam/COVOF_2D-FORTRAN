# Tutoriel 6 : Goutte Statique en Équilibre Capillaire et Loi de Young-Laplace

> **Référence dépôt :** `VALIDATION_TOTALE/Run_Valide-CASES/validation_cas_6_laplace_droplet/` *(validation rapide `Run_Valide-CASES`, cas 6)*.

Ce tutoriel constitue le banc d'essai de référence de l'hydrostatique diphasique capillaire. Il valide le modèle de force superficielle continue CSF (*Continuum Surface Force*) à équilibre de forces (*Balanced-Force*), les lissages d'interface (diffusif et SVOF de Helmholtz), ainsi que l'option hybride par marqueurs lagrangiens CVI pour la maîtrise des courants parasites autour d'une goutte circulaire à l'équilibre.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre le modèle CSF et la loi hydrostatique de Young-Laplace sur maillage décalé MAC.
- Respecter la contrainte de pas de temps capillaire explicite de Brackbill lors du raffinement du maillage.
- Configurer un fichier `.json` conforme avec spécification explicite d'une goutte circulaire (`axis_a`, `axis_b`, `pert_amp`).
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Quantifier dynamiquement le saut de pression $\Delta P$, l'intensité des courants parasites $\|\mathbf{u}\|_\infty$ et la conservation du volume depuis les fichiers `.vtr`.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker`, `InterfaceReconstruction`, `InterfaceMarkers` (optionnel CVI), `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_interface_reconstruction`, `module_markers` (optionnel CVI), `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas valide le couplage hydro-capillaire entre les équations de Navier-Stokes diphasiques et la force de tension superficielle CSF discrétisée aux faces des vitesses exactement aux mêmes points que $-\nabla p$.

Dans un domaine carré fermé $\Omega = [0, 1] \times [0, 1]\,\mathrm{m}^2$ en apesanteur ($\mathbf{g} = \mathbf{0}$) avec conditions de glissement libre (`FREE_SLIP`) sur les 4 parois, une goutte circulaire de rayon $R_0 = 0.25\,\mathrm{m}$ ($C = 1$) est centrée en $(x_c, y_c) = (0.5, 0.5)\,\mathrm{m}$ au repos initial $\mathbf{u}(\mathbf{x}, 0) = \mathbf{0}$. Les deux fluides ont des propriétés identiques $\rho_1 = \rho_2 = 1.0\,\mathrm{kg/m}^3$, $\mu_1 = \mu_2 = 1.0\,\mathrm{Pa\cdot s}$, et un coefficient de tension superficielle $\sigma = 1.0\,\mathrm{N/m}$.

### 2.2 Loi Analytique de Young-Laplace et Contrainte de Stabilité
À l'équilibre exact, la vitesse est identiquement nulle et la courbure géométrique vaut $\kappa_{\mathrm{ex}} = 1/R_0 = 4.0000\,\mathrm{m}^{-1}$. Le saut de pression à travers l'interface vaut :
$$\mathbf{u}_{\mathrm{ex}}(\mathbf{x}, t) \equiv \mathbf{0}, \qquad \Delta P_{\mathrm{ex}} = P_{\mathrm{int}} - P_{\mathrm{ext}} = \frac{\sigma}{R_0} = \frac{1.0}{0.25} = 4.0000\,\mathrm{Pa}.$$

Deux mécanismes numériques gouvernent la précision et la stabilité :
1. **Contrainte de pas de temps capillaire (Brackbill) :** La force capillaire $\mathbf{F}_\sigma = \sigma \kappa \nabla C$ étant traitée explicitement dans le second membre, le pas de temps doit respecter :
   $$\Delta t < \Delta t_{\mathrm{cap}} = \sqrt{\frac{(\rho_1 + \rho_2)\,\Delta x^3}{4\pi \sigma}}.$$
2. **Régularisation et marqueurs :** L'activation de `enable_diffusion_smoothing` (ou de `enable_svof_smoothing`) régularise le calcul de la courbure $\kappa = -\nabla \cdot (\nabla \widetilde{C}/\|\nabla \widetilde{C}\|)$, tandis que l'option `use_markers: true` (`CVI_RK2`) calcule une courbure lagrangienne régulière sans diffusion de l'interface.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et Goutte | $\Omega$, $(x_c, y_c)$, $R_0$ | $[0,1]\times[0,1]\,\mathrm{m}^2$, $(0.5, 0.5)\,\mathrm{m}$, $R_0 = 0.25\,\mathrm{m}$ (`axis_a = axis_b = 0.25`, `pert_amp = 0.0`) |
| Propriétés physiques | $\rho_1=\rho_2, \mu_1=\mu_2, \sigma, \mathbf{g}$ | $\rho = 1.0\,\mathrm{kg/m}^3$, $\mu = 1.0\,\mathrm{Pa\cdot s}$, $\sigma = 1.0\,\mathrm{N/m}$, $\mathbf{g} = \mathbf{0}$ |
| Séquence de grilles et $\Delta t \le 0.5\Delta t_{\mathrm{cap}}$ | $N \in \{32, 64, 128, 256\}$, $t_{\mathrm{final}}=0.02\,\mathrm{s}$ | $\Delta t \in \{10^{-3}, 3\times 10^{-4}, 10^{-4}, 4\times 10^{-5}\}\,\mathrm{s}$ |
| Schémas VOF et régularisation | `vof` & `physics` | `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC`, `VOF_TVD` ; `enable_diffusion_smoothing`, `enable_svof_smoothing`, `use_markers` |
| Couplage & Solveurs | `solver.type`, `linear_solver`, `precond_type` | `AUGMENTED_LAGRANGIAN`, `AAL`, `MONOLITHIC`, `STANDARD_PROJECTION` ; `BICGSTAB2`, `GMRES`, `MINRES` ; `BLOCK_TRIANGULAR`, `JACOBI`, `SSOR`, `ILU0` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_6_wy_plic_n128.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/wy_plic_n128"`.)*

```json
{
  "case_name": "cas_6_wy_plic_n128",
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 1.0,
    "ly": 1.0
  },
  "physics": {
    "enable_two_phase": true,
    "enable_surface_tension": true,
    "enable_convection": false,
    "enable_svof_smoothing": false,
    "svof_filter_width": 0.5,
    "enable_diffusion_smoothing": true,
    "diffusion_smoothing_steps": 3,
    "sigma": 1.0,
    "gx": 0.0,
    "gy": 0.0,
    "fluid1": { "rho": 1.0, "mu": 1.0 },
    "fluid2": { "rho": 1.0, "mu": 1.0 }
  },
  "initial_condition": {
    "type": "DROPLET",
    "center_x": 0.5,
    "center_y": 0.5,
    "radius": 0.25,
    "axis_a": 0.25,
    "axis_b": 0.25,
    "pert_amp": 0.0
  },
  "time": {
    "time_scheme": "GEAR",
    "dt": 0.0001,
    "t_final": 0.02,
    "cfl_max": 0.3,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "top":    { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "linear_solver": "BICGSTAB2",
    "tol": 1e-08,
    "max_iter": 400,
    "r_al": 100.0,
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
    "output_dir": "output/cpp/wy_plic_n128",
    "case_name": "cas_6_wy_plic_n128",
    "vtk_interval": 20,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs utilisent **MPI pur** (`mpirun -np 4` avec les solveurs de Krylov, ou `mpirun -np 1` si `STANDARD_PROJECTION` ou `DIRECT` est choisi).

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/wy_plic_n128 output/fortran/wy_plic_n128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_6_wy_plic_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_6_wy_plic_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Saut de pression de Young-Laplace et convergence en maillage :**
   À $t = t_{\mathrm{final}}$, extraire dynamiquement depuis les fichiers `.vtr` la pression moyenne interne $P_{\mathrm{int}}$ sur $\{r < 0.15\,\mathrm{m},\, C > 0.99\}$ et externe $P_{\mathrm{ext}}$ sur $\{r > 0.38\,\mathrm{m},\, C < 0.01\}$, puis évaluer :
   $$\Delta P_{\mathrm{num}} = P_{\mathrm{int}} - P_{\mathrm{ext}}, \qquad E_{\mathrm{rel}}(\Delta P) = \frac{|\Delta P_{\mathrm{num}} - \sigma/R_0|}{\sigma/R_0}.$$
   Vérifier la décroissance de $E_{\mathrm{rel}}(\Delta P)$ avec le raffinement du maillage ($N = 32 \to 256$).
2. **Maîtrise des courants parasites et conservation du volume :**
   Calculer à $t_{\mathrm{final}}$ le nombre capillaire parasite $\mathrm{Ca}_{\mathrm{parasite}} = \mu \|\mathbf{u}\|_\infty / \sigma$ ainsi que la variation de volume $|\sum C_{i,j}(t_{\mathrm{final}}) - \sum C_{i,j}(0)| / \sum C_{i,j}(0)$ directement à partir de `PhaseIndicator` dans les `.vtr`. Comparer les variantes sans lissage, avec lissage diffusif, avec lissage SVOF et avec marqueurs CVI (`marker_curvature: true`).
3. **Invariance des solveurs, préconditionnement MPI et parité C++17 / Fortran 2018 :**
   Vérifier l'absence de dérive de pression moyenne (projection globale à moyenne nulle dans `BlockTriangularPreconditioner` en MPI pur) et comparer $\Delta P_{\mathrm{num}}$ ainsi que $\|\mathbf{u}\|_\infty$ entre `covof2d_cpp` et `covof2d_fortran` *(note : utiliser `BICGSTAB2` ou `GMRES` avec `BLOCK_TRIANGULAR`, et réserver `MINRES` aux formulations strictement symétriques)*.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` : `step,time,dt,center_of_mass_y,rise_velocity,bubble_volume,kinetic_energy,max_u,max_div` (`9` colonnes C++ / `11` colonnes Fortran avec `avg_rho_gas,avg_p_gas`).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : `PhaseIndicator`, `PhaseIndicator_Smoothed` (si SVOF activé), `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Curvature`, `Velocity`.
- `<case_name>_markers_<step:06d>.vtk` (si `use_markers: true` et `write_markers_vtk: true`).

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 6 (Young-Laplace) sur fichiers .pvtr/.vtr réels."""
import math
import xml.etree.ElementTree as ET
from pathlib import Path

def analyze_young_laplace(pvtr_path: str, sigma: float = 1.0, r0: float = 0.25) -> None:
    pvtr = Path(pvtr_path)
    root = ET.parse(pvtr).getroot()
    p_in, p_out, u_max, vol_drop = [], [], 0.0, 0.0

    for piece in root.findall(".//Piece"):
        vtr = ET.parse(pvtr.parent / piece.attrib["Source"]).getroot()
        ext = [int(v) for v in piece.attrib["Extent"].split()]
        nx_l, ny_l = ext[1] - ext[0], ext[3] - ext[2]
        xc = [float(v) for v in vtr.find(".//DataArray[@Name='XCoordinates']").text.split()]
        yc = [float(v) for v in vtr.find(".//DataArray[@Name='YCoordinates']").text.split()]
        C = [float(v) for v in vtr.find(".//DataArray[@Name='PhaseIndicator']").text.split()]
        P = [float(v) for v in vtr.find(".//DataArray[@Name='Pressure']").text.split()]
        V = [float(v) for v in vtr.find(".//DataArray[@Name='Velocity']").text.split()]
        for j in range(ny_l):
            y = 0.5 * (yc[j] + yc[j + 1])
            dy = yc[j + 1] - yc[j]
            for i in range(nx_l):
                x = 0.5 * (xc[i] + xc[i + 1])
                dx = xc[i + 1] - xc[i]
                k = j * nx_l + i
                vol_drop += C[k] * dx * dy
                r = math.hypot(x - 0.5, y - 0.5)
                if r < 0.15 and C[k] > 0.99:
                    p_in.append(P[k])
                elif r > 0.38 and C[k] < 0.01:
                    p_out.append(P[k])
                spd = math.hypot(V[3 * k], V[3 * k + 1])
                if spd > u_max:
                    u_max = spd

    dp_num = (sum(p_in) / len(p_in)) - (sum(p_out) / len(p_out))
    dp_ex = sigma / r0
    print(f"[{pvtr.name}] dP_num={dp_num:.6f} Pa (ex={dp_ex:.6f} Pa, err={abs(dp_num-dp_ex)/dp_ex:.3e}) "
          f"| max|u|={u_max:.4e} m/s | V_drop={vol_drop:.6f} m^2")

if __name__ == "__main__":
    last_pvtr = sorted(Path("output/cpp/wy_plic_n128").glob("*.pvtr"))[-1]
    analyze_young_laplace(str(last_pvtr))
```

### 6.4 Visualisation sous ParaView
1. Ouvrir le dernier fichier `<case_name>_<step:06d>.pvtr` dans ParaView (**File -> Open**).
2. Tracer l'isocontour `PhaseIndicator = 0.5` avec le filtre **Contour**, puis afficher `Pressure` en fond pour observer le saut de pression uniforme à l'intérieur de la goutte.
3. Appliquer le filtre **Plot Over Line** le long de $y = 0.5\,\mathrm{m}$ sur `Pressure` pour vérifier la raideur de la marche de Young-Laplace et la planéité des plateaux interne et externe.
4. Superposer des vecteurs **Glyph** sur `Velocity` pour localiser les éventuels courants parasites résiduels près de l'interface.
