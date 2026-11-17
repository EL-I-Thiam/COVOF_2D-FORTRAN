# Tutoriel 9 : Sédimentation Confinée 2D d'une Inclusion Cylindrique et Loi de Faxén

Ce tutoriel modélise la sédimentation gravitaire d'un cylindre solide en régime rampant de Stokes ($\mathrm{Re}_p \ll 1$) entre deux parois planes verticales. La vitesse terminale quasi-stationnaire obtenue numériquement est confrontée à la loi analytique de Faxén intégrant l'effet de confinement pariétal.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre la physique du régime de Stokes confiné et la formule analytique de correction pariétale de Faxén.
- Configurer la pénalisation de viscosité ($\mu_p \gg \mu_f$) avec le solveur `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (`AAL`) et le lissage SVOF.
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Extraire dynamiquement la vitesse verticale du centre de masse de la particule depuis les sorties `.vtr` et la comparer à la valeur analytique de Faxén.
- Visualiser le champ de vitesse autour du cylindre sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `InterfaceTracker`, `InterfaceReconstruction`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_interface_vof`, `module_interface_reconstruction`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas valide le couplage entre la pénalisation visqueuse d'une particule rigide ($\mu_p / \mu_f \approx 6530$) et la méthode du Lagrangien Augmenté Adaptatif (`ADAPTATIVE_AUGMENTED_LAGRANGIAN`) pour la sédimentation lente d'un cylindre confiné entre deux parois verticales ($\mathrm{Re}_p \ll 1$).

Dans une cavité carrée fermée $\Omega = [0, 0.1] \times [0, 0.1]\,\mathrm{m}^2$ ($2l = 0.1\,\mathrm{m}$, demi-largeur $l = 0.05\,\mathrm{m}$) à parois adhérentes (`DIRICHLET`, $\mathbf{u} = \mathbf{0}$), un cylindre de rayon $a = 0.00625\,\mathrm{m}$ (rapport de confinement $k_x = a/l = 1/8 = 0.125$) est centré en $(0.05, 0.05)\,\mathrm{m}$ avec `"type": "DROPLET"` ($C = 1$ dans la particule). Les propriétés des deux phases sont :
- **Fluide 1 (Particule solide, $C=1$)** : $\rho_p = 1500\,\mathrm{kg/m}^3$, $\mu_p = 1000\,\mathrm{Pa\cdot s}$.
- **Fluide 2 (Fluide porteur, $C=0$)** : $\rho_f = 1000\,\mathrm{kg/m}^3$, $\mu_f = 0.153125\,\mathrm{Pa\cdot s}$, sous gravité $g_y = -9.81\,\mathrm{m/s}^2$.

### 2.2 Loi Analytique de Vitesse Terminale de Faxén
En régime de Stokes 2D confiné ($k_x = a/l = 0.125$), le développement asymptotique de Faxén donne la vitesse terminale analytique du cylindre ($g = |g_y| = 9.81\,\mathrm{m/s}^2$) :
$$U_{c,\mathrm{th}} = \frac{(\rho_f - \rho_p)\,g\,a^2}{4\,\mu_f}\,\left(-0.9157 - \ln(k_x) + 1.7244\,k_x^2 - 1.7302\,k_x^4\right) \approx -0.05060\,\mathrm{m/s}.$$
*(Note technique sur les diagnostics : dans `VTKWriter::write_diagnostics`, les colonnes `center_of_mass_y` et `rise_velocity` de `*_history.csv` intègrent l'indicatrice $1 - C$, qui correspond au fluide porteur lorsque la particule est définie par `"type": "DROPLET"`. La vitesse verticale de la particule $V_p(t)$ s'extrait donc directement par intégration de $C \cdot v$ sur `PhaseIndicator` et `Velocity` dans les fichiers `.vtr`.)*

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et confinement $k_x$ | $[0, 2l]^2$, $a$, $k_x = a/l$ | $2l = 0.1\,\mathrm{m}$ ($l=0.05\,\mathrm{m}$), $a = 0.00625\,\mathrm{m}$, $k_x = 0.125$, centre $(0.05, 0.05)\,\mathrm{m}$ (`axis_a = axis_b = 0.00625`, `pert_amp = 0.0`) |
| Particule ($C=1$) et Fluide ($C=0$) | $\rho_1, \mu_1$ ; $\rho_2, \mu_2$ | $\rho_1 = \rho_p = 1500\,\mathrm{kg/m}^3, \mu_1 = 1000\,\mathrm{Pa\cdot s}$ ; $\rho_2 = \rho_f = 1000\,\mathrm{kg/m}^3, \mu_2 = 0.153125\,\mathrm{Pa\cdot s}$ |
| Séquence de grilles et $\Delta t$ | $64\times 128$, $128\times 128$, $256\times 256$ | $\Delta t = 5\times 10^{-4}\,\mathrm{s}$, $t_{\mathrm{final}} = 0.15\,\mathrm{s}$ |
| Solveurs et filtrage SVOF | `solver` & `physics` | `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (`k_aal=100`), `AUGMENTED_LAGRANGIAN`, `MONOLITHIC` ; `enable_svof_smoothing: true` (`svof_filter_width: 0.5`) |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_9_sedimentation_faxen_n128.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/n128"`.)*

```json
{
  "case_name": "cas_9_sedimentation_faxen_n128",
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 0.1,
    "ly": 0.1
  },
  "physics": {
    "enable_two_phase": true,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "PICARD",
    "enable_surface_tension": false,
    "enable_svof_smoothing": true,
    "svof_filter_width": 0.5,
    "gx": 0.0,
    "gy": -9.81,
    "fluid1": { "rho": 1500.0, "mu": 1000.0 },
    "fluid2": { "rho": 1000.0, "mu": 0.153125 }
  },
  "initial_condition": {
    "type": "DROPLET",
    "center_x": 0.05,
    "center_y": 0.05,
    "radius": 0.00625,
    "axis_a": 0.00625,
    "axis_b": 0.00625,
    "pert_amp": 0.0
  },
  "time": {
    "dt": 0.0005,
    "t_final": 0.15,
    "time_scheme": "GEAR",
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "right":  { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "DIRICHLET", "u": 0.0, "v": 0.0 }
  },
  "solver": {
    "type": "ADAPTATIVE_AUGMENTED_LAGRANGIAN",
    "adaptative_al": true,
    "k_aal": 100.0,
    "uzawa_tol": 1e-7,
    "uzawa_iter": 6,
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
    "output_dir": "output/cpp/n128",
    "case_name": "cas_9_sedimentation_faxen_n128",
    "vtk_interval": 30,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs fonctionnent en **MPI pur** (`mpirun -np 4`).

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/n128 output/fortran/n128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_9_sedimentation_faxen_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_9_sedimentation_faxen_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Précision sur la vitesse terminale de Faxén :**
   À partir des fichiers `.pvtr`/`.vtr` sur le plateau quasi-stationnaire ($t \in [0.05, 0.10]\,\mathrm{s}$, avant l'approche de la paroi inférieure), calculer la vitesse moyenne verticale du cylindre :
   $$V_p(t) = \frac{\int_\Omega C(\mathbf{x}, t)\,v(\mathbf{x}, t)\,\mathrm{d}\Omega}{\int_\Omega C(\mathbf{x}, t)\,\mathrm{d}\Omega}, \qquad \mathrm{Err}_{\mathrm{rel}}(V_p) = \frac{|V_p - U_{c,\mathrm{th}}|}{|U_{c,\mathrm{th}}|}.$$
   Vérifier la convergence de $V_p$ vers $U_{c,\mathrm{th}} \approx -0.05060\,\mathrm{m/s}$ à mesure que le maillage est raffiné ($64\times 128 \to 256\times 256$).
2. **Rigidité de la particule et énergie cinétique :**
   Vérifier la conservation de la forme circulaire et du volume $V_p(t) = \int_\Omega C\,\mathrm{d}\Omega$ ($|\Delta V_p / V_p(0)| < 10^{-4}$) ainsi que la stabilisation de `kinetic_energy` dans `*_history.csv`.
3. **Parité C++17 / Fortran 2018 :**
   Vérifier la concordance des courbes $V_p(t)$ et `kinetic_energy(t)` entre les deux exécutables.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` : `step,time,dt,center_of_mass_y,rise_velocity,bubble_volume,kinetic_energy,max_u,max_div` (`9` colonnes C++ / `11` colonnes Fortran avec `avg_rho_gas,avg_p_gas`).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : `PhaseIndicator`, `PhaseIndicator_Smoothed`, `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Velocity`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 9 (Sédimentation de Faxén) depuis les .pvtr."""
import math
import xml.etree.ElementTree as ET
from pathlib import Path

def particle_velocity_from_pvtr(pvtr_path: str):
    pvtr = Path(pvtr_path)
    root = ET.parse(pvtr).getroot()
    num, den = 0.0, 0.0
    for piece in root.findall(".//Piece"):
        vtr = ET.parse(pvtr.parent / piece.attrib["Source"]).getroot()
        ext = [int(v) for v in piece.attrib["Extent"].split()]
        nx_l, ny_l = ext[1] - ext[0], ext[3] - ext[2]
        xc = [float(v) for v in vtr.find(".//DataArray[@Name='XCoordinates']").text.split()]
        yc = [float(v) for v in vtr.find(".//DataArray[@Name='YCoordinates']").text.split()]
        C = [float(v) for v in vtr.find(".//DataArray[@Name='PhaseIndicator']").text.split()]
        V = [float(v) for v in vtr.find(".//DataArray[@Name='Velocity']").text.split()]
        for j in range(ny_l):
            dy = yc[j + 1] - yc[j]
            for i in range(nx_l):
                dx = xc[i + 1] - xc[i]
                k = j * nx_l + i
                w = C[k] * dx * dy
                num += w * V[3 * k + 1]
                den += w
    return num / den, den

if __name__ == "__main__":
    kx = 0.00625 / 0.05
    u_faxen = ((1000.0 - 1500.0) * 9.81 * (0.00625 ** 2) / (4.0 * 0.153125)) * (
        -0.9157 - math.log(kx) + 1.7244 * kx**2 - 1.7302 * kx**4
    )
    for p in sorted(Path("output/cpp/n128").glob("*.pvtr")):
        vp, area = particle_velocity_from_pvtr(str(p))
        print(f"{p.name} | Vp={vp:.6f} m/s (Faxen={u_faxen:.6f}, err={abs(vp-u_faxen)/abs(u_faxen):.2%}) | Aire={area:.6e}")
```

### 6.4 Visualisation sous ParaView
1. Ouvrir la série `<case_name>_..pvtr` dans ParaView (**File -> Open**).
2. Tracer l'isocontour `PhaseIndicator = 0.5` pour suivre le cylindre solide et appliquer le filtre **Stream Tracer** (ou **Glyph**) sur `Velocity` pour visualiser les deux tourbillons de recirculation générés entre la particule et les parois verticales.
