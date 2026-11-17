# Tutoriel 7 : Advection d'une Tâche Circulaire dans un Tourbillon Cisaillé

Ce tutoriel isole le transport cinématique pur de l'interface fluide VOF sans résolution hydrodynamique (mode `"VORTEX_ADVECTION"`). Une tâche circulaire est soumise à un champ de vitesse cisaillé instationnaire réversible qui l'étire en une fine spirale avant d'inverser le flux à $t = T/2$ pour reformer le disque initial à $t = T$.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les schémas d'advection VOF géométriques (`WEYMOUTH_YUE`, `CIAM`) et algébriques (`VOF_TVD`), ainsi que la règle de cohérence avec les méthodes de reconstruction et les limiteurs TVD.
- Configurer le mode cinématique `"VORTEX_ADVECTION"` et l'option marqueurs lagrangiens CVI.
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Quantifier dynamiquement l'erreur géométrique de retour $E_{L_1}$, $E_{L_2}$, $E_{L_\infty}$ et la conservation de masse depuis `vortex_error.csv`.
- Visualiser l'étirement maximal à $t = T/2$ et le retour à $t = T$ sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker`, `InterfaceReconstruction`, `InterfaceMarkers` (optionnel CVI), `VTKWriter`, `Simulation` (`Simulation::run_vortex_advection`).
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_interface_reconstruction`, `module_markers` (optionnel CVI), `module_vtk_writer`, `module_simulation` (`sim_run_vortex_advection`).

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas cinématique (benchmark de Rider & Kothe) valide les schémas de transport et de reconstruction d'interface Volume-of-Fluid (VOF) sous étirement tourbillonnaire sévère et réversibilité temporelle exacte.

Dans la cavité carrée unitaire $\Omega = [0, 1] \times [0, 1]\,\mathrm{m}^2$, un disque fluide initial ($C = 1$) de rayon $R_0 = 0.125\,\mathrm{m}$ est centré en $(x_c, y_c) = (0.5, 0.75)\,\mathrm{m}$ à $t = 0$. L'indicatrice de phase $C(\mathbf{x}, t) \in [0, 1]$ obéit à l'équation d'advection pure sans terme source :
$$\frac{\partial C}{\partial t} + \nabla \cdot (C \mathbf{u}) = 0.$$
Le champ de vitesse solénoïdal ($\nabla \cdot \mathbf{u} \equiv 0$) dérive de la fonction de courant instationnaire modulée par un cosinus temporel de période $T = 6.0\,\mathrm{s}$ :
$$\Psi(x, y, t) = \frac{1}{\pi} \sin^2(\pi x) \sin^2(\pi y) \cos\left(\frac{\pi t}{T}\right),$$
évalué au demi-pas de temps $t + \Delta t/2$ aux faces des mailles décalées MAC pour assurer une précision temporelle d'ordre 2 :
$$u(x, y, t) = -\frac{\partial \Psi}{\partial y} = -\sin^2(\pi x) \sin(2\pi y) \cos\left(\frac{\pi t}{T}\right),$$
$$v(x, y, t) = +\frac{\partial \Psi}{\partial x} = +\sin(2\pi x) \sin^2(\pi y) \cos\left(\frac{\pi t}{T}\right).$$

### 2.2 Réversibilité Exacte et Règle de Cohérence VOF
De $t = 0$ à $t = T/2 = 3.0\,\mathrm{s}$, le disque est étiré en une fine spirale ; à $t = 3.0\,\mathrm{s}$, le champ de vitesse change de signe et réenroule le filament jusqu'à $t = T = 6.0\,\mathrm{s}$ où la solution exacte coïncide identiquement avec l'état initial :
$$C_{\mathrm{ex}}(\mathbf{x}, T) \equiv C(\mathbf{x}, 0), \qquad M(T) = \int_\Omega C(\mathbf{x}, T)\,\mathrm{d}\Omega = M(0).$$
Dans `Config.cpp` et `module_config.f90`, une règle stricte de cohérence s'applique :
- Si `advection_method` vaut `"VOF_TVD"`, `reconstruction_method` est forcé à `"NONE"` et un limiteur (`"SUPERBEE"`, `"VAN_LEER"`, `"MUSCL"` ou `"HRIC"`) est utilisé.
- Si `advection_method` vaut `"WEYMOUTH_YUE"` ou `"CIAM"`, `tvd_limiter` est forcé à `"NONE"` et la reconstruction géométrique (`"PLIC"`, `"ELVIRA"`, `"LVIRA"` ou `"SLIC"`) est activée.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et Disque initial | $\Omega$, $(x_c, y_c)$, $R_0$ | $[0,1]\times[0,1]\,\mathrm{m}^2$, $(0.5, 0.75)\,\mathrm{m}$, $R_0 = 0.125\,\mathrm{m}$ (`axis_a = axis_b = 0.125`, `pert_amp = 0.0`) |
| Période cinématique | $T = t_{\mathrm{final}}$ | $6.0\,\mathrm{s}$ (étirement maximal à $T/2 = 3.0\,\mathrm{s}$) |
| Séquence de grilles et $\Delta t$ | $(N\times N, \Delta t)$ | $(64^2, 0.005\,\mathrm{s})$, $(128^2, 0.0025\,\mathrm{s})$, $(256^2, 0.00125\,\mathrm{s})$, $(512^2, 0.0005\,\mathrm{s})$ |
| Méthodes d'advection VOF | `vof.advection_method` | `"WEYMOUTH_YUE"`, `"CIAM"`, `"VOF_TVD"` |
| Reconstructions et limiteurs | `reconstruction_method` / `tvd_limiter` | `"PLIC"`, `"ELVIRA"`, `"LVIRA"`, `"SLIC"` ; `"SUPERBEE"`, `"VAN_LEER"`, `"MUSCL"`, `"HRIC"` |
| Option marqueurs lagrangiens | `vof.use_markers` | `false` (VOF pur) ou `true` (`"CVI_RK2"`, `markers_per_cell: 4`, `write_markers_vtk: true`) |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_7_weymouth_yue_plic_n128.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/wy_plic_n128"`.)*

```json
{
  "case_name": "cas_7_weymouth_yue_plic_n128",
  "mode": "VORTEX_ADVECTION",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 1.0,
    "ly": 1.0
  },
  "physics": {
    "enable_two_phase": true,
    "enable_surface_tension": false,
    "enable_convection": true,
    "sigma": 0.0,
    "fluid1": { "rho": 1000.0, "mu": 0.001 },
    "fluid2": { "rho": 1.0,    "mu": 1.0e-05 }
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.5
  },
  "initial_condition": {
    "type": "VORTEX_SHEAR",
    "center_x": 0.5,
    "center_y": 0.75,
    "radius": 0.125,
    "axis_a": 0.125,
    "axis_b": 0.125,
    "pert_amp": 0.0
  },
  "time": {
    "dt": 0.0025,
    "t_final": 6.0,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "right":  { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "AUGMENTED_LAGRANGIAN",
    "linear_solver": "BICGSTAB2",
    "tol": 1e-12,
    "max_iter": 500,
    "r_al": 100.0,
    "uzawa_iter": 1,
    "use_precond": true,
    "precond_type": "JACOBI"
  },
  "output": {
    "output_dir": "output/cpp/wy_plic_n128",
    "case_name": "cas_7_weymouth_yue_plic_n128",
    "vtk_interval": 240,
    "log_interval": 200,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Le mode `"VORTEX_ADVECTION"` est entièrement parallélisé via **MPI pur** (échange des halos VOF et réductions globales des erreurs $L_1, L_2, L_\infty$ et des masses).

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/wy_plic_n128 output/fortran/wy_plic_n128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_7_weymouth_yue_plic_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_7_weymouth_yue_plic_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Erreur géométrique de retour à $t = T = 6.0\,\mathrm{s}$ :**
   Extraire directement depuis `vortex_error.csv` les normes calculées par les binaires :
   $$E_{L_1} = \frac{1}{N_x N_y}\sum_{i=1}^{N_x}\sum_{j=1}^{N_y} |C_{i,j}(T) - C_{i,j}(0)|, \qquad E_{L_2} = \sqrt{\frac{1}{N_x N_y}\sum_{i=1}^{N_x}\sum_{j=1}^{N_y} |C_{i,j}(T) - C_{i,j}(0)|^2}.$$
   Calculer l'ordre de convergence spatiale $p = \ln(E_{L_1}(2h)/E_{L_1}(h)) / \ln(2)$ entre grilles successives ($p \in [1.0, 2.0]$ attendu selon la résolution de la queue du filament).
2. **Conservation de la masse et bornitude :**
   Vérifier dans `vortex_error.csv` que `mass_conservation = mass_final / mass_init` reste à $1.00000 \pm 10^{-5}$ pour `WEYMOUTH_YUE` et `VOF_TVD`, et que $0 \le C_{i,j}(t) \le 1$ dans chaque fichier `.vtr`.
3. **Parité C++17 / Fortran 2018 :**
   Vérifier que $|E_{L_1}^{\mathrm{C++}} - E_{L_1}^{\mathrm{Fortran}}| < 10^{-12}$ et que les erreurs $E_{L_2}$ et `mass_conservation` concordent.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `vortex_error.csv` : `nx,ny,dt,T,advection,reconstruction,limiter,err_l1,err_l2,err_linf,mass_init,mass_final,mass_conservation`.
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` (exportés à $t=0$, à $t \approx T/2$, tous les `vtk_interval` pas, et à $t = T$).
- `<case_name>_markers_<step:06d>.vtk` (si `use_markers: true` et `write_markers_vtk: true`).

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 7 (Vortex cisaillé) sur vortex_error.csv."""
import csv
import math
from pathlib import Path

def read_vortex_error(out_dir: str) -> dict:
    with (Path(out_dir) / "vortex_error.csv").open(newline="") as f:
        return next(csv.DictReader(f))

if __name__ == "__main__":
    prev = None
    for n in [64, 128, 256, 512]:
        d = f"output/cpp/wy_plic_n{n}"
        if not (Path(d) / "vortex_error.csv").exists():
            continue
        r = read_vortex_error(d)
        l1 = float(r["err_l1"])
        mc = float(r["mass_conservation"]) * 100.0
        order = math.log(float(prev["err_l1"]) / l1) / math.log(n / int(prev["nx"])) if prev else float("nan")
        print(f"N={n:3d} | {r['advection']}+{r['reconstruction']} | L1={l1:.4e} "
              f"| L2={float(r['err_l2']):.4e} | Masse={mc:.5f}% | p={order:.2f}")
        prev = r
```

### 6.4 Visualisation sous ParaView
1. Ouvrir la série temporelle `<case_name>_..pvtr` dans ParaView (**File -> Open**).
2. Afficher `PhaseIndicator` et appliquer **Contour** sur la valeur $0.5$ aux instants clés $t = 0$, $t = T/2 = 3.0\,\mathrm{s}$ (étirement maximal en spirale) et $t = T = 6.0\,\mathrm{s}$.
3. Comparer visuellement la superposition de l'isocontour final $t = T$ avec l'isocontour initial $t = 0$ (et, si `use_markers: true`, charger les fichiers `*_markers_*.vtk` avec rendu `Point Gaussian`).
