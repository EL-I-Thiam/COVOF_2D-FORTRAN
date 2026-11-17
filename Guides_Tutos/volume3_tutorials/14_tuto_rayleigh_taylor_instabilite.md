# Tutoriel 14 : Instabilité Interfaciale Diphasique de Rayleigh-Taylor (Dispersion Linéaire)

Ce dernier tutoriel simule l'instabilité hydrodynamique de Rayleigh-Taylor se développant à l'interface entre un fluide lourd superposé à un fluide léger soumis à la gravité. Le taux de croissance exponentiel en régime linéaire et le seuil de stabilisation capillaire sont confrontés à la relation de dispersion analytique hydro-capillaire.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre la compétition entre déstabilisation gravitationnelle (nombre d'Atwood $\mathrm{At}$) et stabilisation capillaire (seuil critique $\sigma_c$).
- Configurer la condition initiale `"RAYLEIGH_TAYLOR"` (`center_y`, `pert_amp`, `pert_k`) et les schémas d'ordre 2 (`GEAR`, `ADAMS_BASHFORTH`, `WEYMOUTH_YUE+PLIC`).
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Extraire dynamiquement l'amplitude de la perturbation $A(t)$ depuis les `.vtr` et estimer $\omega_{\mathrm{num}}$ par ajustement sur $A_0 \cosh(\omega t)$.
- Visualiser le développement du champignon non-linéaire sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker`, `InterfaceReconstruction`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_interface_reconstruction`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas valide la compétition entre l'instabilité gravitationnelle de Rayleigh-Taylor (nombre d'Atwood $\mathrm{At} = (\rho_1 - \rho_2)/(\rho_1 + \rho_2) = 0.5$) et la stabilisation capillaire CSF, du régime linéaire exponentiel jusqu'à l'enroulement non-linéaire en champignon.

Dans un domaine $\Omega = [0, 1.0] \times [0, 4.0]\,\mathrm{m}^2$ avec parois latérales `FREE_SLIP` (ou `PERIODIC`) et parois haute/basse `DIRICHLET` adhérentes ($\mathbf{u} = \mathbf{0}$), le fluide lourd ($\rho_1 = 3.0\,\mathrm{kg/m}^3, \mu_1 = 10^{-3}\,\mathrm{Pa\cdot s}$, $C=1$ pour $y \ge y_{\mathrm{int}}$) surmonte initialement le fluide léger ($\rho_2 = 1.0\,\mathrm{kg/m}^3, \mu_2 = 10^{-3}\,\mathrm{Pa\cdot s}$, $C=0$ pour $y < y_{\mathrm{int}}$) sous $g_y = -9.81\,\mathrm{m/s}^2$. L'interface initiale est perturbée selon (`"type": "RAYLEIGH_TAYLOR"`, `"center_y": 2.0`, `"pert_amp": 0.05`, `"pert_k": 1.0`) :
$$y_{\mathrm{int}}(x, 0) = 2.0 + A_0 \cos\left(\frac{2\pi m x}{L_x}\right), \qquad k = \frac{2\pi m}{L_x} = 2\pi\,\mathrm{m}^{-1}, \quad A_0 = 0.05\,\mathrm{m}.$$

### 2.2 Relation de Dispersion Hydro-Capillaire Linéaire
En régime linéaire non visqueux, le taux de croissance théorique $\omega_{\mathrm{th}}(k)$ et le seuil de coupure capillaire $\sigma_c$ ($g = |g_y| = 9.81\,\mathrm{m/s}^2$) valent :
$$\omega_{\mathrm{th}}^2(k) = \frac{(\rho_1 - \rho_2)\,g\,k - \sigma\,k^3}{\rho_1 + \rho_2}, \qquad \sigma_c = \frac{(\rho_1 - \rho_2)\,g}{k^2} \approx 0.4970\,\mathrm{N/m}.$$
- Pour $\sigma = 0\,\mathrm{N/m}$ : régime instable pur avec $\omega_{\mathrm{th}}(0) \approx 5.5515\,\mathrm{s}^{-1}$.
- Pour $\sigma = 0.20\,\mathrm{N/m}$ : croissance ralentie par la capillarité avec $\omega_{\mathrm{th}}(0.2) \approx 4.2915\,\mathrm{s}^{-1}$.
- Pour $\sigma = 0.60\,\mathrm{N/m} > \sigma_c$ : $\omega_{\mathrm{th}}^2 < 0$, l'interface est stabilisée et oscille de façon dispersive sans croître.

En régime linéaire ($k A(t) \ll 1$), l'amplitude de la perturbation depuis le repos évolue selon $A(t) = A_0 \cosh(\omega t)$.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et perturbation initiale | $[0, 1]\times[0, 4]\,\mathrm{m}^2$, $y_c, A_0, m$ | $y_c = 2.0\,\mathrm{m}$, $A_0 = 0.05\,\mathrm{m}$, $m = 1$ ($k = 2\pi\,\mathrm{m}^{-1}$, `RAYLEIGH_TAYLOR`) |
| Fluides ($\mathrm{At} = 0.5$) et gravité | $\rho_1, \rho_2, \mu_1, \mu_2, g_y$ | $\rho_1 = 3.0\,\mathrm{kg/m}^3$, $\rho_2 = 1.0\,\mathrm{kg/m}^3$, $\mu_1 = \mu_2 = 10^{-3}\,\mathrm{Pa\cdot s}$, $g_y = -9.81\,\mathrm{m/s}^2$ |
| Trois régimes capillaires $\sigma$ | $\sigma \in \{0.0, 0.20, 0.60\}\,\mathrm{N/m}$ | Instable pur ($\sigma=0$), ralenti ($\sigma=0.20$), stabilisé ($\sigma=0.60 > \sigma_c$) |
| Séquence de grilles et $\Delta t$ | $64\times 256$, $128\times 512$, $256\times 1024$ | $\Delta t = 10^{-3}\,\mathrm{s}$, $t_{\mathrm{final}} = 1.0\,\mathrm{s}$ |
| Schémas VOF et solveurs | `vof` & `solver` | `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC`, `VOF_TVD` ; `ADAPTATIVE_AUGMENTED_LAGRANGIAN`, `AUGMENTED_LAGRANGIAN`, `MONOLITHIC` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_14_rayleigh_taylor_sigma0_128x512.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/sigma0_128x512"`. Pour les cas $\sigma > 0$, passer `"enable_surface_tension": true` et `"enable_diffusion_smoothing": true`.)*

```json
{
  "case_name": "cas_14_rayleigh_taylor_sigma0_128x512",
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 512,
    "lx": 1.0,
    "ly": 4.0
  },
  "physics": {
    "enable_two_phase": true,
    "enable_convection": true,
    "convective_scheme": "TVD_SUPERBEE",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_surface_tension": false,
    "sigma": 0.0,
    "gx": 0.0,
    "gy": -9.81,
    "fluid1": { "rho": 3.0, "mu": 1.0e-3 },
    "fluid2": { "rho": 1.0, "mu": 1.0e-3 }
  },
  "initial_condition": {
    "type": "RAYLEIGH_TAYLOR",
    "center_y": 2.0,
    "pert_amp": 0.05,
    "pert_k": 1.0
  },
  "time": {
    "dt": 0.001,
    "t_final": 1.0,
    "time_scheme": "GEAR",
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "ADAPTATIVE_AUGMENTED_LAGRANGIAN",
    "adaptative_al": true,
    "k_aal": 100.0,
    "uzawa_tol": 1e-6,
    "uzawa_iter": 4,
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
    "output_dir": "output/cpp/sigma0_128x512",
    "case_name": "cas_14_rayleigh_taylor_sigma0_128x512",
    "vtk_interval": 50,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs utilisent **MPI pur** (`mpirun -np 4` ou plus sur maillage fin $256\times 1024$).

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/sigma0_128x512 output/fortran/sigma0_128x512

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_14_rayleigh_taylor_sigma0_128x512.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_14_rayleigh_taylor_sigma0_128x512.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Taux de croissance linéaire $\omega_{\mathrm{num}}$ et seuil capillaire $\sigma_c$ :**
   À partir des fichiers `.pvtr`/`.vtr` sur la fenêtre linéaire $t \in [0.05, 0.30]\,\mathrm{s}$, calculer la demi-amplitude crête-à-creux $A(t) = \frac{1}{2}\big(y_{\mathrm{int},\max}(t) - y_{\mathrm{int},\min}(t)\big)$ de l'isocontour $C = 0.5$ et estimer $\omega_{\mathrm{num}}$ par ajustement sur $A(t) = A_0 \cosh(\omega_{\mathrm{num}} t)$. Comparer avec $\omega_{\mathrm{th}}(\sigma)$ pour $\sigma \in \{0.0, 0.20\}\,\mathrm{N/m}$, et vérifier la stabilisation oscillatoire pour $\sigma = 0.60 > \sigma_c$.
2. **Symétrie bilatérale et conservation de masse :**
   Vérifier la symétrie axiale $C(L_x - x, y, t) = C(x, y, t)$ tout au long de l'enroulement non-linéaire et vérifier que la variation relative de volume $|\Delta V / V_0| < 10^{-4}$ depuis `PhaseIndicator` dans les `.vtr`.
3. **Parité C++17 / Fortran 2018 :**
   Comparer $A(t)$, `kinetic_energy` (`*_history.csv`) et la morphologie des panaches entre `covof2d_cpp` et `covof2d_fortran`.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` (`9` colonnes en C++ / `11` colonnes en Fortran).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : `PhaseIndicator`, `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Curvature`, `Velocity`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 14 (Rayleigh-Taylor) depuis les .pvtr."""
import math
import xml.etree.ElementTree as ET
from pathlib import Path

def interface_amplitude_and_volume(pvtr_path: str):
    pvtr = Path(pvtr_path)
    root = ET.parse(pvtr).getroot()
    y_mix, vol_heavy = [], 0.0
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
                dx = xc[i + 1] - xc[i]
                cval = C[j * nx_l + i]
                vol_heavy += cval * dx * dy
                if 0.05 < cval < 0.95:
                    y_mix.append(y)
    amp = 0.5 * (max(y_mix) - min(y_mix)) if y_mix else 0.0
    return amp, vol_heavy

if __name__ == "__main__":
    k = 2.0 * math.pi
    omega_th = math.sqrt((3.0 - 1.0) * 9.81 * k / (3.0 + 1.0))
    dt = 0.001
    for pvtr in sorted(Path("output/cpp/sigma0_128x512").glob("*.pvtr")):
        t = int(pvtr.stem.split("_")[-1]) * dt
        amp, vol = interface_amplitude_and_volume(str(pvtr))
        omega_est = math.acosh(max(1.0, amp / 0.05)) / t if 0.05 <= t <= 0.30 else float("nan")
        print(f"t={t:.3f} s | A={amp:.4f} m | omega_est={omega_est:.3f} s^-1 (th={omega_th:.3f}) | V1={vol:.6f}")
```

### 6.4 Visualisation sous ParaView
1. Ouvrir la série `<case_name>_..pvtr` dans ParaView (**File -> Open**).
2. Colorer par `PhaseIndicator` et tracer l'isocontour `PhaseIndicator = 0.5` pour observer le développement symétrique de la pointe descendante (*spike*) et de la bulle ascendante (*bubble*) avec enroulement de Kelvin-Helmholtz aux bords.
3. Utiliser le filtre **Compute Derivatives** (vorticité) ou **Glyph** sur `Velocity` pour visualiser les deux tourbillons contrarotatifs qui alimentent le champignon.
