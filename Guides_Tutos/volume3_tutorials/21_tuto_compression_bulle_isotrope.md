# Tutoriel 21 : Compression Isotrope d'une Bulle d'Air en Cavité Carrée

> **Référence dépôt :** `VALIDATION_TOTALE/Run_CASES/validation_cas_9_compression_bulle/` *(campagne avancée `Run_CASES`)*.

Ce tutoriel évalue la capacité de l'algorithme bas-Mach compressible (`C-AAL` / `MONOLITHIC` compressible) à traiter la compression multidimensionnelle d'une bulle de gaz immergée sous injection symétrique sur les quatre parois. L'évolution temporelle du rayon équivalent de la bulle, de sa masse volumique et de la pression est confrontée à la solution analytique exacte déduite du bilan de masse global.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre la dynamique de compression d'une bulle gazeuse dans un liquide faiblement compressible à bas nombre de Mach.
- Configurer correctement les lois d'état (`"eos": "LINEAR"`, `"chi_T"`), la géométrie circulaire explicite (`axis_a`, `axis_b`, `pert_amp: 0.0`) et les injections pariétales.
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Quantifier dynamiquement $R(t)$, $\overline{\rho}_{\mathrm{air}}(t)$, $\overline{p}_{\mathrm{air}}(t)$ et la circularité de la bulle depuis les fichiers `.vtr`.
- Visualiser le rétrécissement isotrope de la bulle sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `EOS`, `InterfaceTracker`, `InterfaceReconstruction`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_eos`, `module_interface_vof`, `module_interface_reconstruction`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas valide la compression centripète d'une bulle d'air circulaire immergée dans de l'eau au moyen du modèle 1-fluide faiblement compressible (`"compressible": true`).

Dans une cavité carrée $\Omega = [0, 0.1] \times [0, 0.1]\,\mathrm{m}^2$ sans gravité ($\mathbf{g}=\mathbf{0}$) ni tension superficielle ($\sigma = 0$), une bulle d'air de rayon initial $R_0 = 0.03\,\mathrm{m}$ ($C=0$, $\rho_{2,0} = 1.176829268\,\mathrm{kg/m}^3$, $\mu_2 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$, $\chi_{T,2} = 9.8692327\times 10^{-6}\,\mathrm{Pa}^{-1}$) est centrée en $(0.05, 0.05)\,\mathrm{m}$ dans de l'eau ($C=1$, $\rho_{1,0} = 1000\,\mathrm{kg/m}^3$, $\mu_1 = 10^{-3}\,\mathrm{Pa\cdot s}$, $\chi_{T,1} = 4.44\times 10^{-10}\,\mathrm{Pa}^{-1}$), à la pression initiale $p_0 = 101\,325\,\mathrm{Pa}$. Sur chacune des 4 parois, de l'eau est injectée normalement vers l'intérieur à la vitesse $V_0 = 0.0025\,\mathrm{m/s}$, soit un débit surfacique entrant total :
$$q = 4 L V_0 = 4 \times 0.1 \times 0.0025 = 1.0\times 10^{-3}\,\mathrm{m}^2/\mathrm{s}.$$

### 2.2 Loi Analytique d'Évolution du Rayon, de la Densité et de la Pression
L'eau étant quasi-incompressible ($\chi_{T,1} \ll \chi_{T,2}$), le débit entrant $q$ comprime quasi intégralement le volume gazeux : $\pi R^2(t) = \pi R_0^2 - q\,t$, d'où un temps d'effondrement théorique $t_{\mathrm{collapse}} = \pi R_0^2 / q \approx 2.827433\,\mathrm{s}$ et les lois exactes pour $t < t_{\mathrm{collapse}}$ :
$$R_{\mathrm{th}}(t) = R_0\sqrt{1 - \frac{q\,t}{\pi R_0^2}}, \qquad \rho_{\mathrm{th}}(t) = \frac{\rho_{2,0}}{1 - \dfrac{q\,t}{\pi R_0^2}}, \qquad p_{\mathrm{th}}(t) = \frac{p_0}{1 - \dfrac{q\,t}{\pi R_0^2}}.$$
À $t = 1.5\,\mathrm{s}$, la solution analytique donne $R_{\mathrm{th}}(1.5\,\mathrm{s}) \approx 0.020554\,\mathrm{m}$, $\rho_{\mathrm{th}}(1.5\,\mathrm{s}) \approx 2.506647\,\mathrm{kg/m}^3$ et $p_{\mathrm{th}}(1.5\,\mathrm{s}) \approx 215\,821.8\,\mathrm{Pa}$.

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et bulle initiale | $[0, L]^2$, $(x_c, y_c)$, $R_0$ | $L = 0.1\,\mathrm{m}$, centre $(0.05, 0.05)\,\mathrm{m}$, $R_0 = 0.03\,\mathrm{m}$ (`axis_a = axis_b = 0.03`, `pert_amp = 0.0`), $p_0 = 101\,325\,\mathrm{Pa}$ |
| Fluides 1 (Eau) et 2 (Air) | $\rho_{i,0}, \mu_i, \chi_{T,i}$, `eos` | Eau : $\rho_{1,0}=1000$, $\chi_{T,1}=4.44\times 10^{-10}$ ; Air : $\rho_{2,0}=1.176829268$, $\chi_{T,2}=9.8692327\times 10^{-6}\,\mathrm{Pa}^{-1}$ (`"LINEAR"`) |
| Conditions d'injection (4 bords) | Gauche, Droite, Bas, Haut | `DIRICHLET` entrant : $u_{\mathrm{left}}=+V_0$, $u_{\mathrm{right}}=-V_0$, $v_{\mathrm{bottom}}=+V_0$, $v_{\mathrm{top}}=-V_0$ ($V_0=0.0025\,\mathrm{m/s}$) |
| Séquence de grilles et $\Delta t$ | $N \in \{32, 64, 128, 256, 512\}$ | $\Delta t = 0.00125\,\mathrm{s}$, $t_{\mathrm{final}} = 1.5\,\mathrm{s}$ (1200 pas ; extension possible jusqu'à $t = 2.5\,\mathrm{s}$) |
| Solveurs et schémas VOF | `solver` & `vof` | `MONOLITHIC` ($r_{\mathrm{al}}=0$) ou `C-AAL` avec `BICGSTAB2` / `GMRES` (`BLOCK_TRIANGULAR`, `ILU0`, `SSOR`) ; `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC`, `VOF_TVD` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_11_compression_bulle_n128.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/n128"`.)*

```json
{
  "case_name": "cas_11_compression_bulle_n128",
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 0.1,
    "ly": 0.1
  },
  "physics": {
    "enable_two_phase": true,
    "compressible": true,
    "flow_regime": "COMPRESSIBLE",
    "enable_convection": true,
    "enable_surface_tension": false,
    "gx": 0.0,
    "gy": 0.0,
    "fluid1": { "rho": 1000.0,      "mu": 1.0e-3, "eos": "LINEAR", "chi_T": 4.44e-10 },
    "fluid2": { "rho": 1.176829268, "mu": 1.85e-5, "eos": "LINEAR", "chi_T": 9.8692327e-6 }
  },
  "initial_condition": {
    "type": "BUBBLE",
    "center_x": 0.05,
    "center_y": 0.05,
    "radius": 0.03,
    "axis_a": 0.03,
    "axis_b": 0.03,
    "pert_amp": 0.0,
    "p_init": 101325.0
  },
  "time": {
    "dt": 0.00125,
    "t_final": 1.5,
    "time_scheme": "EULER",
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "DIRICHLET", "u":  0.0025, "v":  0.0 },
    "right":  { "type": "DIRICHLET", "u": -0.0025, "v":  0.0 },
    "bottom": { "type": "DIRICHLET", "u":  0.0,    "v":  0.0025 },
    "top":    { "type": "DIRICHLET", "u":  0.0,    "v": -0.0025 }
  },
  "solver": {
    "type": "MONOLITHIC",
    "linear_solver": "BICGSTAB2",
    "tol": 1e-10,
    "max_iter": 1000,
    "r_al": 0.0,
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
    "output_dir": "output/cpp/n128",
    "case_name": "cas_11_compression_bulle_n128",
    "vtk_interval": 100,
    "log_interval": 20,
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
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_11_compression_bulle_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_11_compression_bulle_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Convergence vers $\rho_{\mathrm{th}}(t)$, $p_{\mathrm{th}}(t)$ et $R_{\mathrm{th}}(t)$ :**
   À partir des fichiers `.pvtr`/`.vtr`, évaluer le rayon équivalent $R_{\mathrm{num}}(t) = \sqrt{\frac{1}{\pi}\int_\Omega (1-C)\,\mathrm{d}\Omega}$, la densité moyenne $\overline{\rho}_{\mathrm{air}}(t)$ et la pression moyenne $\overline{p}_{\mathrm{air}}(t)$ dans $\{C < 0.01\}$. Vérifier la décroissance des erreurs relatives lors du raffinement du maillage ($32^2 \to 512^2$).
2. **Conservation de la masse de gaz et isotropie :**
   Vérifier que la masse totale d'air $M_{\mathrm{air}}(t) = \int_\Omega (1-C)\rho\,\mathrm{d}\Omega$ reste conservée ($|M_{\mathrm{air}}(t_{\mathrm{final}}) - M_{\mathrm{air}}(0)|/M_{\mathrm{air}}(0) < 5\times 10^{-3}$) et que la bulle demeure symétrique par rapport aux axes $x = 0.05\,\mathrm{m}$ et $y = 0.05\,\mathrm{m}$.
3. **Parité C++17 / Fortran 2018 :**
   Vérifier la concordance de $R_{\mathrm{num}}(t)$, $\overline{\rho}_{\mathrm{air}}(t)$ et $\overline{p}_{\mathrm{air}}(t)$ entre les deux moteurs.

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` (`9` colonnes en C++ dont `bubble_volume` $= \pi R_{\mathrm{num}}^2(t)$ ; `11` colonnes en Fortran incluant `avg_rho_gas,avg_p_gas`).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : `PhaseIndicator`, `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Velocity`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 21 (Compression isotrope de bulle)."""
import math
import xml.etree.ElementTree as ET
from pathlib import Path

def eval_bubble(pvtr_path: str):
    pvtr = Path(pvtr_path)
    root = ET.parse(pvtr).getroot()
    vol_b, sum_rho, sum_p, sum_v_pure = 0.0, 0.0, 0.0, 0.0
    for piece in root.findall(".//Piece"):
        vtr = ET.parse(pvtr.parent / piece.attrib["Source"]).getroot()
        ext = [int(v) for v in piece.attrib["Extent"].split()]
        nx_l, ny_l = ext[1] - ext[0], ext[3] - ext[2]
        xc = [float(v) for v in vtr.find(".//DataArray[@Name='XCoordinates']").text.split()]
        yc = [float(v) for v in vtr.find(".//DataArray[@Name='YCoordinates']").text.split()]
        C = [float(v) for v in vtr.find(".//DataArray[@Name='PhaseIndicator']").text.split()]
        P = [float(v) for v in vtr.find(".//DataArray[@Name='Pressure']").text.split()]
        R = [float(v) for v in vtr.find(".//DataArray[@Name='Density']").text.split()]
        for j in range(ny_l):
            dy = yc[j + 1] - yc[j]
            for i in range(nx_l):
                dx = xc[i + 1] - xc[i]
                dvol = dx * dy
                k = j * nx_l + i
                vol_b += (1.0 - C[k]) * dvol
                if C[k] < 0.01:
                    sum_rho += R[k] * dvol
                    sum_p += P[k] * dvol
                    sum_v_pure += dvol
    r_num = math.sqrt(vol_b / math.pi)
    return r_num, sum_rho / sum_v_pure, sum_p / sum_v_pure

if __name__ == "__main__":
    r0, rho0, p0, q, dt = 0.03, 1.176829268, 101325.0, 1.0e-3, 0.00125
    for pvtr in sorted(Path("output/cpp/n128").glob("*.pvtr")):
        t = int(pvtr.stem.split("_")[-1]) * dt
        rn, rhon, pn = eval_bubble(str(pvtr))
        f = 1.0 - q * t / (math.pi * r0**2)
        print(f"t={t:.3f} s | R={rn:.6f} (th={r0*math.sqrt(f):.6f}) "
              f"| rho={rhon:.5f} (th={rho0/f:.5f}) | p={pn:.1f} (th={p0/f:.1f})")
```

### 6.4 Visualisation sous ParaView
1. Ouvrir la série `<case_name>_..pvtr` dans ParaView (**File -> Open**).
2. Appliquer le filtre **Contour** sur `PhaseIndicator = 0.5` et **Glyph** sur `Velocity` pour visualiser l'écoulement radial convergent qui comprime la bulle de manière isotrope.
3. Afficher `Pressure` et vérifier l'uniformité spatiale du champ de pression à chaque instant.
