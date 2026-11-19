# Tutoriel 18 : Oscillations Capillaires Amorties d'une Bulle de Gaz (Théorie de Lamb)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_CASES/validation_cas_6_oscillation_bulle_capillaire/` *(campagne avancée `Run_CASES`)*.

Ce tutoriel met en œuvre une bulle de gaz initialement déformée en ellipse selon le mode quadrupolaire $n = 2$ dans un métal liquide dense, sous l'action combinée de la tension superficielle et de la dissipation visqueuse. La fréquence propre d'oscillation et l'amortissement sont confrontés à la théorie linéaire de Rayleigh-Lamb.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre la dynamique inertio-capillaire d'une bulle oscillante à fort rapport de masse volumique ($\rho_1/\rho_2 \approx 5950$).
- Configurer une bulle elliptique explicite (`"axis_a"`, `"axis_b"`, `"pert_amp": 0.0`) et la régularisation CSF (`enable_diffusion_smoothing`).
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Extraire dynamiquement les demi-axes $a(t), b(t)$, la fréquence d'oscillation $f_{\mathrm{num}}$ et le bilan énergétique depuis les `.vtr` et `*_history.csv`.
- Visualiser les cycles d'oscillation et la pression capillaire sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker`, `InterfaceReconstruction`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_interface_reconstruction`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas valide la dynamique instationnaire du couplage entre l'inertie diphasique à fort rapport de densité ($\rho_1/\rho_2 \approx 5950$), la dissipation visqueuse et la force de tension superficielle CSF lors des oscillations libres d'une bulle initialement elliptique.

Dans une cavité carrée $\Omega = [0, 0.1] \times [0, 0.1]\,\mathrm{m}^2$ sans gravité ($\mathbf{g} = \mathbf{0}$) avec parois à glissement libre (`FREE_SLIP`), une bulle d'air ($\rho_2 = 1.17683\,\mathrm{kg/m}^3$, $\mu_2 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$) est centrée en $(0.05, 0.05)\,\mathrm{m}$ dans un métal liquide dense ($\rho_1 = 7000\,\mathrm{kg/m}^3$, $\mu_1 = 4.0\times 10^{-3}\,\mathrm{Pa\cdot s}$) avec $\sigma = 1.8\,\mathrm{N/m}$. À $t=0$, la bulle de rayon équivalent $R_0 = 0.02\,\mathrm{m}$ est déformée sur le mode $n=2$ avec une excentricité $\varepsilon = 0.20$ :
$$a_0 = R_0 (1 + \varepsilon) = 0.024\,\mathrm{m}, \qquad b_0 = \frac{R_0}{1 + \varepsilon} \approx 0.0166666667\,\mathrm{m} \qquad (a_0 b_0 = R_0^2).$$

### 2.2 Relation de Dispersion Capillaire Linéaire
En milieu infini non visqueux, la pulsation propre $\omega_n$ et la fréquence $f_{\mathrm{th}}$ du mode $n = 2$ de Rayleigh-Lamb cylindrique valent :
$$\omega_n^2 = \frac{(n^3 - n)\,\sigma}{(\rho_1 + \rho_2)\,R_0^3} \implies \omega_2 = \sqrt{\frac{6\,\sigma}{(\rho_1 + \rho_2)\,R_0^3}}, \qquad f_{\mathrm{th}} = \frac{\omega_2}{2\pi} \approx 5.70632\,\mathrm{Hz} \quad (T_{\mathrm{th}} \approx 0.17524\,\mathrm{s}).$$
En cavité carrée confinée ($L / R_0 = 5$) et pour une amplitude finie $\varepsilon = 0.20$, l'effet de confinement pariétal abaisse légèrement la fréquence d'oscillation effective par rapport à la valeur en domaine infini, tandis que l'énergie cinétique $E_k(t)$ oscille à la fréquence double $2 f_{\mathrm{num}}$.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et bulle initiale | $[0, L]^2$, $(x_c, y_c)$, $a_0, b_0$ | $L = 0.1\,\mathrm{m}$, centre $(0.05, 0.05)\,\mathrm{m}$, $R_0 = 0.02\,\mathrm{m}$, $a_0 = 0.024\,\mathrm{m}$, $b_0 = 0.0166666667\,\mathrm{m}$ (`pert_amp = 0.0`) |
| Fluide 1 (Métal liquide) | $\rho_1, \mu_1$ | $\rho_1 = 7000.0\,\mathrm{kg/m}^3$, $\mu_1 = 4.0\times 10^{-3}\,\mathrm{Pa\cdot s}$ |
| Fluide 2 (Air) | $\rho_2, \mu_2$ | $\rho_2 = 1.17683\,\mathrm{kg/m}^3$, $\mu_2 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$ |
| Tension superficielle et lissage | $\sigma$, `enable_diffusion_smoothing` | $\sigma = 1.8\,\mathrm{N/m}$, `enable_diffusion_smoothing: true`, `diffusion_smoothing_steps: 3` |
| Séquence de grilles et $\Delta t$ | $N \in \{32, 64, 128, 256, 512\}$ | $\Delta t = 2.5\times 10^{-4}\,\mathrm{s}$ (ou adapté selon $\Delta t_{\mathrm{cap}}$ sur grille fine), $t_{\mathrm{final}} = 0.6\,\mathrm{s}$ |
| Solveurs et schémas VOF | `solver` & `vof` | `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (`AAL`), `AUGMENTED_LAGRANGIAN`, `MONOLITHIC` ; `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_13_oscillation_bulle_n128.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/n128"`.)*

```json
{
  "case_name": "cas_13_oscillation_bulle_n128",
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 0.1,
    "ly": 0.1
  },
  "physics": {
    "enable_two_phase": true,
    "compressible": false,
    "enable_convection": true,
    "enable_surface_tension": true,
    "enable_diffusion_smoothing": true,
    "diffusion_smoothing_steps": 3,
    "sigma": 1.8,
    "gx": 0.0,
    "gy": 0.0,
    "fluid1": { "rho": 7000.0,  "mu": 4.0e-3 },
    "fluid2": { "rho": 1.17683, "mu": 1.85e-5 }
  },
  "initial_condition": {
    "type": "BUBBLE",
    "center_x": 0.05,
    "center_y": 0.05,
    "radius": 0.02,
    "axis_a": 0.024,
    "axis_b": 0.0166666667,
    "pert_amp": 0.0
  },
  "time": {
    "dt": 0.00025,
    "t_final": 0.6,
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
    "uzawa_iter": 6,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-08,
    "max_iter": 600,
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
    "output_dir": "output/cpp/n128",
    "case_name": "cas_13_oscillation_bulle_n128",
    "vtk_interval": 20,
    "log_interval": 5,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs utilisent **MPI pur** (`mpirun -np 4`).

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/n128 output/fortran/n128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_13_oscillation_bulle_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_13_oscillation_bulle_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Fréquence d'oscillation $f_{\mathrm{num}}$ et amortissement visqueux :**
   Extraire la série temporelle d'énergie cinétique $E_k(t)$ depuis `*_history.csv` et les demi-axes $a(t), b(t)$ depuis les fichiers `.vtr` (à partir des moments géométriques d'ordre 2 de $1 - C$). Déterminer la période d'oscillation $T_{\mathrm{num}}$ par détection des maxima de $a(t)$ (ou à partir de la demi-période des pics de $E_k(t)$) et vérifier la convergence de $f_{\mathrm{num}} = 1/T_{\mathrm{num}}$ avec le raffinement du maillage.
2. **Conservation du volume de la bulle et stabilité à haut rapport de densité ($\rho_1/\rho_2 \approx 5950$) :**
   Vérifier par sommation de $(1-C)\Delta x\Delta y$ dans les `.vtr` que $|V(t_{\mathrm{final}}) - V(0)| / V(0) < 10^{-4}$ avec `WEYMOUTH_YUE+PLIC`.
3. **Parité C++17 / Fortran 2018 :**
   Vérifier la superposition des courbes $E_k(t)$, $a(t)$ et $b(t)$ entre `covof2d_cpp` et `covof2d_fortran`.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` (`9` colonnes en C++ / `11` colonnes en Fortran).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : `PhaseIndicator`, `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Curvature`, `Velocity`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 18 (Oscillation de bulle / Lamb)."""
import csv
import math
import xml.etree.ElementTree as ET
from pathlib import Path

def estimate_freq_from_ke(csv_path: str) -> float:
    t_vals, ke_vals = [], []
    with Path(csv_path).open(newline="") as f:
        for r in csv.DictReader(f):
            t_vals.append(float(r["time"]))
            ke_vals.append(float(r["kinetic_energy"]))
    peaks = [t_vals[i] for i in range(1, len(ke_vals) - 1)
             if ke_vals[i] > ke_vals[i - 1] and ke_vals[i] >= ke_vals[i + 1] and ke_vals[i] > 1e-8]
    if len(peaks) < 3:
        return float("nan")
    dt_peaks = [(peaks[i + 1] - peaks[i]) for i in range(len(peaks) - 1)]
    t_half = sum(dt_peaks) / len(dt_peaks)
    return 1.0 / (2.0 * t_half)  # E_k oscille à 2*f

def bubble_axes_from_pvtr(pvtr_path: str):
    pvtr = Path(pvtr_path)
    root = ET.parse(pvtr).getroot()
    m0, mxx, myy = 0.0, 0.0, 0.0
    for piece in root.findall(".//Piece"):
        vtr = ET.parse(pvtr.parent / piece.attrib["Source"]).getroot()
        ext = [int(v) for v in piece.attrib["Extent"].split()]
        nx_l, ny_l = ext[1] - ext[0], ext[3] - ext[2]
        xc = [float(v) for v in vtr.find(".//DataArray[@Name='XCoordinates']").text.split()]
        yc = [float(v) for v in vtr.find(".//DataArray[@Name='YCoordinates']").text.split()]
        C = [float(v) for v in vtr.find(".//DataArray[@Name='PhaseIndicator']").text.split()]
        for j in range(ny_l):
            y = 0.5 * (yc[j] + yc[j + 1])
            dy = yc[j + 1] - yc[j]
            for i in range(nx_l):
                x = 0.5 * (xc[i] + xc[i + 1])
                dx = xc[i + 1] - xc[i]
                w = (1.0 - C[j * nx_l + i]) * dx * dy
                m0 += w
                mxx += w * (x - 0.05) ** 2
                myy += w * (y - 0.05) ** 2
    return 2.0 * math.sqrt(mxx / m0), 2.0 * math.sqrt(myy / m0), m0

if __name__ == "__main__":
    f_th = math.sqrt(6.0 * 1.8 / ((7000.0 + 1.17683) * 0.02 ** 3)) / (2.0 * math.pi)
    csv_f = Path("output/cpp/n128/cas_13_oscillation_bulle_n128_history.csv")
    if csv_f.exists():
        f_num = estimate_freq_from_ke(str(csv_f))
        print(f"Fréquence numérique f_num = {f_num:.4f} Hz (théorie milieu infini = {f_th:.4f} Hz)")
```

### 6.4 Visualisation sous ParaView
1. Ouvrir la série `<case_name>_..pvtr` dans ParaView (**File -> Open**).
2. Appliquer le filtre **Contour** sur `PhaseIndicator = 0.5` et colorer le fond par `Pressure` pour observer l'inversion alternée de la surpression aux pôles horizontaux puis verticaux.
3. Afficher les vecteurs **Glyph** de `Velocity` pour visualiser l'écoulement quadripolaire interne et externe.
