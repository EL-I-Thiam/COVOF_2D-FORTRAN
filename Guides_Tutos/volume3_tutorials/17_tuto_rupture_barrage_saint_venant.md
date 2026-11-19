# Tutoriel 17 : Rupture de Barrage 2D sur Fond Sec (Solution de Saint-Venant)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_CASES/validation_cas_5_rupture_barrage/` *(campagne avancée `Run_CASES`)*.

Ce tutoriel aborde l'effondrement gravitationnel instantané d'une colonne d'eau sur fond sec (*Dam Break*). La cinématique d'avancée du front d'onde et le profil parabolique de la surface libre sont confrontés à la solution analytique des caractéristiques de Saint-Venant (solution de Ritter).

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre la physique de l'effondrement d'une colonne liquide à grand rapport de densité ($\rho_1/\rho_2 \approx 850$) et le modèle analytique en eau peu profonde de Saint-Venant.
- Configurer la condition initiale rectangulaire `"DAM_BREAK"` (`box_x0`, `box_y0`, `box_x1`, `box_y1`) et les schémas `GEAR` + `ADAMS_BASHFORTH` + `TVD_SUPERBEE`.
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Extraire dynamiquement la hauteur d'eau intégrée $h_{\mathrm{num}}(x,t)$ et l'abscisse du front $x_{\mathrm{front}}(t)$ depuis les fichiers `.vtr`.
- Visualiser la propagation de la vague et les tourbillons d'air sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker`, `InterfaceReconstruction`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_interface_reconstruction`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas valide le couplage Navier-Stokes diphasique à surface libre ($\rho_1/\rho_2 \approx 850$, $\mu_1/\mu_2 = 100$) sous gravité sur le problème classique de rupture de barrage 2D sur fond sec.

Dans un canal $\Omega = [0, 1.2] \times [0, 0.2]\,\mathrm{m}^2$ à parois glissantes (`FREE_SLIP`), une colonne d'eau occupe initialement $[0, 0.6] \times [0, 0.1]\,\mathrm{m}^2$ ($x_0 = 0.6\,\mathrm{m}, h_g = 0.1\,\mathrm{m}$) avec $\rho_1 = 1000\,\mathrm{kg/m}^3, \mu_1 = 10^{-3}\,\mathrm{Pa\cdot s}$, $\rho_2 = 1.1768\,\mathrm{kg/m}^3, \mu_2 = 10^{-5}\,\mathrm{Pa\cdot s}$, $\sigma = 0.075\,\mathrm{N/m}$ et $g_y = -9.81\,\mathrm{m/s}^2$.

### 2.2 Solution Analytique des Caractéristiques de Saint-Venant (Ritter)
Avec la célérité d'onde initiale $c_0 = \sqrt{g h_g} \approx 0.990454\,\mathrm{m/s}$ ($g = |g_y| = 9.81\,\mathrm{m/s}^2$), la hauteur d'eau analytique en théorie d'eau peu profonde (avant réflexion sur le mur droit $x = 1.2\,\mathrm{m}$ à $t_{\mathrm{refl}} = (L_x - x_0)/(2c_0) \approx 0.303\,\mathrm{s}$ pour le front et sans effet de bord amont) vaut :
$$h_{\mathrm{SV}}(x,t) = \begin{cases} h_g & \text{si } x \le x_0 - c_0 t, \\ \dfrac{1}{9g}\left(2c_0 - \dfrac{x - x_0}{t}\right)^2 & \text{si } x_0 - c_0 t < x < x_0 + 2c_0 t, \\ 0 & \text{si } x \ge x_0 + 2c_0 t. \end{cases}$$
L'abscisse théorique du front d'onde sec est $x_{\mathrm{front,SV}}(t) = x_0 + 2 c_0 t$. Numériquement, le profil de surface libre se calcule directement par intégration verticale de la fraction volumique d'eau :
$$h_{\mathrm{num}}(x_i, t) = \sum_{j=1}^{N_y} C_{i,j}(t)\,\Delta y.$$

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et retenue initiale | $[0, 1.2]\times[0, 0.2]\,\mathrm{m}^2$, $x_0, h_g$ | $x_0 = 0.6\,\mathrm{m}$, $h_g = 0.1\,\mathrm{m}$ (`DAM_BREAK`), $g_y = -9.81\,\mathrm{m/s}^2$ |
| Fluides 1 (Eau) et 2 (Air) | $\rho_1, \mu_1, \rho_2, \mu_2, \sigma$ | $\rho_1 = 1000\,\mathrm{kg/m}^3$, $\mu_1 = 10^{-3}\,\mathrm{Pa\cdot s}$ ; $\rho_2 = 1.1768\,\mathrm{kg/m}^3$, $\mu_2 = 10^{-5}\,\mathrm{Pa\cdot s}$ ; $\sigma = 0.075\,\mathrm{N/m}$ |
| Séquence de grilles et $\Delta t$ | $120\times 20$, $240\times 40$, $480\times 80$ | $\Delta t \in \{10^{-3}, 5\times 10^{-4}, 2.5\times 10^{-4}\}\,\mathrm{s}$, $t_{\mathrm{final}} = 0.50\,\mathrm{s}$ (contrôles à $0.24\,\mathrm{s}$ et $0.40\,\mathrm{s}$) |
| Schémas VOF et solveurs | `vof` & `solver` | `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC`, `VOF_TVD` (`SUPERBEE`) ; `ADAPTATIVE_AUGMENTED_LAGRANGIAN`, `AUGMENTED_LAGRANGIAN`, `MONOLITHIC` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_12_rupture_barrage_240x40.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/grid_240x40"`.)*

```json
{
  "case_name": "cas_12_rupture_barrage_240x40",
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 240,
    "ny": 40,
    "lx": 1.2,
    "ly": 0.2
  },
  "physics": {
    "enable_two_phase": true,
    "enable_convection": true,
    "convective_scheme": "TVD_SUPERBEE",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_surface_tension": true,
    "enable_diffusion_smoothing": true,
    "diffusion_smoothing_steps": 2,
    "sigma": 0.075,
    "gx": 0.0,
    "gy": -9.81,
    "fluid1": { "rho": 1000.0, "mu": 1.0e-3 },
    "fluid2": { "rho": 1.1768, "mu": 1.0e-5 }
  },
  "initial_condition": {
    "type": "DAM_BREAK",
    "box_x0": 0.0,
    "box_y0": 0.0,
    "box_x1": 0.6,
    "box_y1": 0.1
  },
  "time": {
    "dt": 0.0005,
    "t_final": 0.50,
    "time_scheme": "GEAR",
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "top":    { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "ADAPTATIVE_AUGMENTED_LAGRANGIAN",
    "adaptative_al": true,
    "k_aal": 100.0,
    "uzawa_tol": 1e-6,
    "uzawa_iter": 5,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-08,
    "max_iter": 800,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.3
  },
  "output": {
    "output_dir": "output/cpp/grid_240x40",
    "case_name": "cas_12_rupture_barrage_240x40",
    "vtk_interval": 40,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs utilisent **MPI pur** (`mpirun -np 4`).

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/grid_240x40 output/fortran/grid_240x40

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_12_rupture_barrage_240x40.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_12_rupture_barrage_240x40.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Profil de surface libre $h_{\mathrm{num}}(x,t)$ et front d'onde $x_{\mathrm{front}}(t)$ :**
   À $t = 0.24\,\mathrm{s}$ (avant réflexion sur la paroi droite), extraire $h_{\mathrm{num}}(x_i, t) = \sum_{j=1}^{N_y} C_{i,j}(t)\Delta y$ depuis les `.vtr` et comparer à $h_{\mathrm{SV}}(x_i, t)$ en norme $L_2$. À $t = 0.40\,\mathrm{s}$, observer la remontée de vague sur la paroi droite.
2. **Conservation du volume d'eau et comparaison des schémas VOF :**
   Vérifier par sommation de `PhaseIndicator * dx * dy` que la variation relative de volume d'eau $|\Delta V / V_0| < 10^{-4}$ avec `WEYMOUTH_YUE+PLIC` et `VOF_TVD`.
3. **Performance de `AAL` vs `AL` et parité C++17 / Fortran 2018 :**
   Comparer la convergence Krylov et vérifier l'accord de $h_{\mathrm{num}}(x,t)$ et $E_k(t)$ entre `covof2d_cpp` et `covof2d_fortran`.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` (`9` colonnes en C++ / `11` colonnes en Fortran).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : `PhaseIndicator`, `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Curvature`, `Velocity`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 17 (Rupture de barrage / Saint-Venant)."""
import math
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

def h_saint_venant(x: float, t: float, x0: float = 0.6, hg: float = 0.1, g: float = 9.81) -> float:
    c0 = math.sqrt(g * hg)
    if x <= x0 - c0 * t:
        return hg
    if x < x0 + 2.0 * c0 * t:
        return ((2.0 * c0 - (x - x0) / t) ** 2) / (9.0 * g)
    return 0.0

def extract_height_profile(pvtr_path: str):
    pvtr = Path(pvtr_path)
    root = ET.parse(pvtr).getroot()
    cols = defaultdict(float)
    vol_water = 0.0
    for piece in root.findall(".//Piece"):
        vtr = ET.parse(pvtr.parent / piece.attrib["Source"]).getroot()
        ext = [int(v) for v in piece.attrib["Extent"].split()]
        nx_l, ny_l = ext[1] - ext[0], ext[3] - ext[2]
        xc = [float(v) for v in vtr.find(".//DataArray[@Name='XCoordinates']").text.split()]
        yc = [float(v) for v in vtr.find(".//DataArray[@Name='YCoordinates']").text.split()]
        C = [float(v) for v in vtr.find(".//DataArray[@Name='PhaseIndicator']").text.split()]
        for j in range(ny_l):
            dy = yc[j + 1] - yc[j]
            for i in range(nx_l):
                dx = xc[i + 1] - xc[i]
                x = 0.5 * (xc[i] + xc[i + 1])
                cval = C[j * nx_l + i]
                cols[round(x, 8)] += cval * dy
                vol_water += cval * dx * dy
    return sorted(cols.items()), vol_water

if __name__ == "__main__":
    p = Path("output/cpp/grid_240x40/cas_12_rupture_barrage_240x40_000480.pvtr")
    if p.exists():
        prof, vol = extract_height_profile(str(p))
        l2 = math.sqrt(sum((h - h_saint_venant(x, 0.24)) ** 2 for x, h in prof) / len(prof))
        x_front = max((x for x, h in prof if h > 1e-3), default=0.6)
        print(f"t=0.24 s | L2(h)={l2:.4e} m | x_front={x_front:.4f} m (SV={0.6+2*math.sqrt(9.81*0.1)*0.24:.4f}) | V={vol:.6f}")
```

### 6.4 Visualisation sous ParaView
1. Ouvrir la série `<case_name>_..pvtr` dans ParaView (**File -> Open**).
2. Tracer l'isocontour `PhaseIndicator = 0.5` sur fond coloré par la norme de `Velocity` afin de suivre le collapse gravitationnel, l'avancée du front sec et la réflexion contre la paroi aval.
3. Utiliser **Glyph** sur `Velocity` pour observer le contre-écoulement de l'air au-dessus de la nappe liquide.
