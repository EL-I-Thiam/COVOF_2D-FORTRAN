# Tutoriel 20 : Compression d'Air par Injection d'Eau (Piston Compressible Bas Mach)

> **Référence dépôt :** `VALIDATION_TOTALE/Run_CASES/validation_cas_8_piston_compressible/` *(campagne avancée `Run_CASES`)*.

Ce tutoriel valide la formulation compressible à bas nombre de Mach (Lagrangien Augmenté Compressible `C-AAL` / `MONOLITHIC` compressible) sur un cas d'injection d'un liquide quasi-incompressible comprimant une poche de gaz dans une enceinte étanche. La montée en pression et en masse volumique est confrontée à la loi analytique exacte issue de la conservation de la masse.

---

## 1. Objectifs Pédagogiques et Modules Logiciels Activés

### 1.1 Compétences Acquises
- Comprendre les équations 1-fluide faiblement compressibles à bas nombre de Mach et l'équation de continuité avec compressibilité isotherme $\chi_T(C)$.
- Configurer correctement les paramètres EOS (`"eos"`, `"chi_T"` avec majuscule) et l'état initial compressible (`"p_init"`, `"u"`, `"v"`).
- Exécuter le cas test en parallèle MPI pur avec `covof2d_cpp` et `covof2d_fortran`.
- Quantifier dynamiquement $\overline{\rho}_{\mathrm{air}}(t)$ et $\overline{p}_{\mathrm{air}}(t)$ face à la loi analytique de compression.
- Visualiser le mouvement de l'interface plane et la montée en pression homogène sous ParaView.

### 1.2 Modules Logiciels Sollicités
- **Moteur C++17 (`COVOF_2D-C++`) :** `Config`, `Mesh2D`, `FlowFields`, `MpiComm`, `EOS`, `InterfaceTracker`, `InterfaceReconstruction`, `NavierStokesSolver`, `BiCGSTAB2Solver`, `BlockTriangularPreconditioner`, `VTKWriter`, `Simulation`.
- **Moteur Fortran 2018 (`COVOF_2D-FORTRAN`) :** `module_config`, `module_json`, `module_mesh`, `module_fields`, `module_mpi_comm`, `module_eos`, `module_interface_vof`, `module_interface_reconstruction`, `module_navier_stokes`, `module_linear_solvers`, `module_vtk_writer`, `module_simulation`.

---

## 2. Formulation Physique et Solution Analytique Exacte

### 2.1 Équations Physiques et Configuration
Ce cas valide l'extension du modèle diphasique 1-fluide aux écoulements faiblement compressibles à bas nombre de Mach (`"compressible": true`), dans lesquels un liquide quasi-incompressible comprime un gaz compressible au sein d'une cavité fermée.

Le système couple l'équation de continuité compressible isotherme $\chi_T(C) = C\,\chi_{T,1} + (1-C)\,\chi_{T,2}$, l'équation de quantité de mouvement, l'advection de la fonction couleur avec terme de dilatation $C\,\nabla\cdot\mathbf{u}$ et la loi d'état linéaire (`"eos": "LINEAR"`) :
$$\chi_T(C)\frac{\partial p}{\partial t} + \nabla \cdot \mathbf{u} = 0, \qquad \rho(C)\left(\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u}\right) = -\nabla p + \nabla\cdot\left[\mu(C)\left(\nabla\mathbf{u} + \nabla\mathbf{u}^T - \frac{2}{3}(\nabla\cdot\mathbf{u})\mathbf{I}\right)\right],$$
$$\frac{\partial C}{\partial t} + \nabla\cdot(\mathbf{u}C) = C\,\nabla\cdot\mathbf{u}, \qquad \rho_i(p) = \rho_{i,0}\left[1 + \chi_{T,i}\,(p - p_0)\right].$$

### 2.2 Configuration Géométrique et Loi Analytique Exacte
Dans une cavité carrée $\Omega = [0, 0.1] \times [0, 0.1]\,\mathrm{m}^2$ sans gravité ($\mathbf{g}=\mathbf{0}$) ni capillarité ($\sigma=0$), l'eau ($C=1$) occupe initialement $y \in [0, 0.025]\,\mathrm{m}$, surmontée d'une hauteur $H_0 = 0.075\,\mathrm{m}$ d'air ($C=0$) à $p_0 = 101\,325\,\mathrm{Pa}$. De l'eau est injectée par la face inférieure $y=0$ à la vitesse verticale constante $V_0 = 0.1\,\mathrm{m/s}$, les trois autres parois étant à glissement libre (`FREE_SLIP`). L'interface plane monte selon $y_{\mathrm{int}}(t) = 0.025 + V_0 t$, et la masse volumique ainsi que la pression de l'air vérifient la loi exacte :
$$\rho_{\mathrm{air}}(t) = \frac{\rho_{\mathrm{air},0}}{1 - \dfrac{V_0\,t}{H_0}}, \qquad p_{\mathrm{isotherme}}(t) = \frac{p_0}{1 - \dfrac{V_0\,t}{H_0}}.$$
*(Pour la loi d'état linéaire $\rho_2(p) = \rho_{2,0}[1 + \chi_{T,2}(p - p_0)]$ avec $\chi_{T,2} = 1/p_0 = 9.8692327\times 10^{-6}\,\mathrm{Pa}^{-1}$, on retrouve exactement $\rho_{\mathrm{air}}(t)$ et $p(t) = p_0 + \frac{1}{\chi_{T,2}}\frac{V_0 t / H_0}{1 - V_0 t / H_0} = \frac{p_0}{1 - V_0 t / H_0}$.)*

---

## 3. Paramètres Physiques, Numériques et Discrétisation

### 3.1 Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :----------------------- | :----------------------- | :------------------------------------------------------- |
| Domaine et interface initiale | $[0, L_x]\times[0, L_y]$, $h_0$, $H_0$ | $L_x = L_y = 0.1\,\mathrm{m}$, eau sur $y\in[0, 0.025]\,\mathrm{m}$, air sur $H_0 = 0.075\,\mathrm{m}$ |
| Fluide 1 (Eau liquide) | $\rho_{1,0}, \mu_1, \chi_{T,1}$, `eos` | $\rho_{1,0} = 1000.0\,\mathrm{kg/m}^3$, $\mu_1 = 10^{-3}\,\mathrm{Pa\cdot s}$, $\chi_{T,1} = 4.44\times 10^{-10}\,\mathrm{Pa}^{-1}$, `"LINEAR"` |
| Fluide 2 (Air compressible) | $\rho_{2,0}, \mu_2, \chi_{T,2}$, `eos` | $\rho_{2,0} = 1.176829268\,\mathrm{kg/m}^3$, $\mu_2 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$, $\chi_{T,2} = 9.8692327\times 10^{-6}\,\mathrm{Pa}^{-1}$, `"LINEAR"` |
| Conditions aux limites | Bas ($y=0$), Haut/Gauche/Droite | Injection `DIRICHLET` $(u=0, v=V_0=0.1\,\mathrm{m/s})$ en bas ; `FREE_SLIP` sur les 3 autres bords |
| Séquence de grilles et $\Delta t$ | $N \in \{32, 64, 128, 256, 512\}$ | $\Delta t = 10^{-3}\,\mathrm{s}$, $t_{\mathrm{final}} = 0.5\,\mathrm{s}$ (500 pas ; contrôles à $t \in \{0.1, 0.2, 0.3, 0.4, 0.5\}\,\mathrm{s}$) |
| Solveurs et schémas VOF | `solver` & `vof` | `MONOLITHIC` ($r_{\mathrm{al}}=0$) ou `C-AAL` avec `BICGSTAB2` / `GMRES` / `DIRECT` ; `WEYMOUTH_YUE+PLIC`, `CIAM+PLIC`, `VOF_TVD` |

---

## 4. Fichier de Configuration JSON Pas à Pas

### 4.1 Prototype `case/cpp/cas_10_piston_compressible_n128.json`
*(Pour Fortran 2018, placer une copie dans `case/fortran/` avec `"output_dir": "output/fortran/n128"`. Respecter strictement la casse `"chi_T"`.)*

```json
{
  "case_name": "cas_10_piston_compressible_n128",
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
    "type": "DAM_BREAK",
    "box_x0": 0.0,
    "box_y0": 0.0,
    "box_x1": 0.1,
    "box_y1": 0.025,
    "p_init": 101325.0,
    "u": 0.0,
    "v": 0.0
  },
  "time": {
    "dt": 0.001,
    "t_final": 0.5,
    "time_scheme": "EULER",
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.1 },
    "top":    { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 }
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
    "case_name": "cas_10_piston_compressible_n128",
    "vtk_interval": 50,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 5. Procédure d'Exécution (Moteurs C++17 et Fortran 2018)

Les deux moteurs utilisent **MPI pur** (`mpirun -np 4` avec solveurs itératifs Krylov ; `mpirun -np 1` si `DIRECT` est sélectionné).

```bash
# 1. Créer les répertoires de sortie relatifs courts
mkdir -p output/cpp/n128 output/fortran/n128

# 2. Exécution parallèle MPI avec le moteur C++17
mpirun -np 4 /chemin/vers/COVOF_2D-C++/bin/covof2d_cpp case/cpp/cas_10_piston_compressible_n128.json

# 3. Exécution parallèle MPI avec le moteur Fortran 2018
mpirun -np 4 /chemin/vers/COVOF_2D-FORTRAN/bin/covof2d_fortran case/fortran/cas_10_piston_compressible_n128.json
```

---

## 6. Analyse Quantitative, Script Python et Visualisation ParaView

### 6.1 Critères de Vérification, Formules Exactes et Tolérances
1. **Précision sur $\overline{\rho}_{\mathrm{air}}(t)$ et $\overline{p}_{\mathrm{air}}(t)$ :**
   À partir des fichiers `.pvtr`/`.vtr` aux instants $t \in \{0.1, 0.2, 0.3, 0.4, 0.5\}\,\mathrm{s}$, extraire la masse volumique et la pression moyennes dans la zone gazeuse pure $\{C < 0.01\}$ et calculer l'écart relatif par rapport à $\rho_{\mathrm{th}}(t)$ et $p_{\mathrm{th}}(t)$ :
   $$E_{\mathrm{rel}}(\rho, t) = \frac{|\overline{\rho}_{\mathrm{air}}(t) - \rho_{\mathrm{th}}(t)|}{\rho_{\mathrm{th}}(t)}, \qquad E_{\mathrm{rel}}(p, t) = \frac{|\overline{p}_{\mathrm{air}}(t) - p_{\mathrm{th}}(t)|}{p_{\mathrm{th}}(t)}.$$
   Vérifier la convergence de ces erreurs vers zéro lors du raffinement du maillage ($32^2 \to 512^2$).
2. **Conservation de la masse totale de gaz et planéité de l'interface :**
   Vérifier que la masse totale d'air $M_{\mathrm{air}}(t) = \int_\Omega (1 - C)\,\rho\,\mathrm{d}\Omega$ reste constante ($|M_{\mathrm{air}}(t) - M_{\mathrm{air}}(0)| / M_{\mathrm{air}}(0) < 10^{-3}$) et que le front liquide-gaz reste horizontal à $y_{\mathrm{int}}(t) = 0.025 + V_0 t$.
3. **Parité C++17 / Fortran 2018 :**
   Vérifier la concordance des courbes $\overline{\rho}_{\mathrm{air}}(t)$ et $\overline{p}_{\mathrm{air}}(t)$ entre `covof2d_cpp` et `covof2d_fortran` ($|\overline{p}^{\mathrm{C++}} - \overline{p}^{\mathrm{Fortran}}| / p_{\mathrm{th}} < 10^{-6}$).

### 6.2 Fichiers Réellement Générés dans `output_dir`
- `<case_name>_history.csv` (`9` colonnes en C++ ; `11` colonnes en Fortran incluant `avg_rho_gas,avg_p_gas`).
- `<case_name>_<step:06d>.pvtr` et `<case_name>_p<rank>_<step:06d>.vtr` : `PhaseIndicator`, `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Velocity`.

### 6.3 Script Python de Vérification Dynamique

```python
#!/usr/bin/env python3
"""Vérification dynamique du Tutoriel 20 (Piston compressible) sur fichiers .pvtr."""
import xml.etree.ElementTree as ET
from pathlib import Path

def eval_gas_state(pvtr_path: str):
    pvtr = Path(pvtr_path)
    root = ET.parse(pvtr).getroot()
    sum_rho, sum_p, sum_vol, sum_mass_air = 0.0, 0.0, 0.0, 0.0
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
                sum_mass_air += (1.0 - C[k]) * R[k] * dvol
                if C[k] < 0.01:
                    sum_rho += R[k] * dvol
                    sum_p += P[k] * dvol
                    sum_vol += dvol
    return sum_rho / sum_vol, sum_p / sum_vol, sum_mass_air

if __name__ == "__main__":
    rho0, p0, v0, h0 = 1.176829268, 101325.0, 0.1, 0.075
    dt = 0.001
    for pvtr in sorted(Path("output/cpp/n128").glob("*.pvtr")):
        step = int(pvtr.stem.split("_")[-1])
        t = step * dt
        rho_m, p_m, m_air = eval_gas_state(str(pvtr))
        factor = 1.0 / (1.0 - v0 * t / h0)
        rho_th, p_th = rho0 * factor, p0 * factor
        print(f"t={t:.3f} s | rho={rho_m:.6f} (th={rho_th:.6f}, err={abs(rho_m-rho_th)/rho_th:.3e}) "
              f"| p={p_m:.1f} (th={p_th:.1f}) | M_air={m_air:.6e}")
```

### 6.4 Visualisation sous ParaView
1. Ouvrir la série `<case_name>_..pvtr` dans ParaView (**File -> Open**).
2. Afficher `PhaseIndicator` avec un contour à $0.5$ pour vérifier la planéité parfaite de la surface libre ascendante.
3. Tracer un profil vertical via **Plot Over Line** ($x = 0.05\,\mathrm{m}$) sur `Pressure` et `Density` pour observer la pression uniforme dans la poche d'air et le saut de masse volumique à l'interface.
