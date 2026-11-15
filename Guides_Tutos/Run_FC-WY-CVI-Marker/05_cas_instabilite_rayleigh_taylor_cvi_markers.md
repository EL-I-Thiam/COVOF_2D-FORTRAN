# INSTRUCTIONS D'EXÉCUTION POUR AGENT / MODÈLE IA — CAS 5 (FC-WY-CVI-MARKER) : INSTABILITÉ GRAVITAIRE DE RAYLEIGH-TAYLOR DIPHASIQUE ET VOLUTES SECONDAIRES (`05_cas_instabilite_rayleigh_taylor_cvi_markers`)

## 1. Mission et Rôle de l'Agent

Vous êtes un agent autonome chargé d'exécuter, de vérifier et de post-traiter de bout en bout le **Cas de Validation 5 de la méthode hybride FC-WY-CVI-Marker (`05_cas_instabilite_rayleigh_taylor_cvi_markers`)** sur les moteurs **C++17** (`covof2d_cpp`) et **Fortran 2018** (`covof2d_fortran`) de la plateforme **COVOF 2D**.

Ce cas évalue l'apport du couplage **Solveur Monolithique (FC) + VOF-PLIC (Weymouth & Yue) + Interpolation CVI + Marqueurs Non Connectés** sur :
1. Le taux de croissance linéaire initial $\gamma$ de l'instabilité de Rayleigh-Taylor (avec et sans tension superficielle $\sigma$) ;
2. La résolution géométrique fine des **volutes secondaires de Kelvin-Helmholtz** qui s'enroulent aux extrémités du champignon lors de la phase non linéaire tardive ($t \in [0.5, 1.0]\,\mathrm{s}$), là où les schémas VOF purement eulériens brisent les spires minces.

---

## 2. Théorie : Relation de Dispersion de Rayleigh-Taylor et Régime Non Linéaire

Dans une cavité rectangulaire $\Omega = [0, L_x] \times [0, L_y]$ ($L_x = 0.5\,\mathrm{m}, L_y = 1.5\,\mathrm{m}$), un fluide lourd ($\rho_1 = 1.225\,\mathrm{kg/m}^3$, en haut) surmonte un fluide léger ($\rho_2 = 0.1694\,\mathrm{kg/m}^3$, en bas) sous la gravité $g_y = -9.81\,\mathrm{m/s}^2$, soit un nombre d'Atwood :
$$At = \frac{\rho_1 - \rho_2}{\rho_1 + \rho_2}.$$
L'interface initiale perturbée par un mode cosinusoïdal de nombre d'onde $k = 2\pi / L_x$ et de faible amplitude $\eta_0$ croît exponentiellement en phase linéaire ($\eta(t) = \eta_0 \cosh(\gamma t)$) avec le taux de croissance théorique :
$$\gamma_{\mathrm{th}} = \sqrt{At \, g \, k - \frac{\sigma \, k^3}{\rho_1 + \rho_2}}.$$

En phase non linéaire avancée, le panache descendant forme un champignon dont les bords latéraux s'enroulent en spirales très minces ($\delta \sim \Delta x$). Les marqueurs lagrangiens advectés par **CVI** (avec filtrage directionnel $\mathbf{n}_p \cdot \mathbf{n}_{i,j}^{\mathrm{ref}} \ge 0.5$ pour ne pas mélanger les faces opposées du bras spiralé) maintiennent la définition géométrique des volutes et stabilisent la courbure capillaire $\kappa^{\mathrm{marker}}$ lorsque $\sigma > 0$.

---

## 3. Matrice des Jeux de Paramètres et Prototype `.json`

| Axe d'étude | Grilles $N_x \times N_y$ | Tension superficielle $\sigma$ | Comparaison |
| :--- | :---: | :---: | :--- |
| **Axe A : Croissance linéaire et volutes non linéaires ($\sigma = 0$)** | $64\times 192$, $128\times 384$ | $\sigma = 0.0\,\mathrm{N/m}$ | `WY + PLIC` seul vs `WY + PLIC + CVI_RK2 (marker_normals_plic: true)` |
| **Axe B : Stabilisation capillaire ($\sigma > 0$)** | $64\times 192$, $128\times 384$ | $\sigma \in \{0.005, 0.02\}\,\mathrm{N/m}$ | Courbure CSF standard vs Courbure parabolique sur marqueurs CVI (`marker_curvature: true`) |
| **Axe C : Étude du Lissage SVOF de Helmholtz ($C_s$) vs Lissage Diffusif ($\widetilde{C}$)** | $64\times 192$, $128\times 384$ | $\sigma \in \{0.0, 0.005\}\,\mathrm{N/m}$ | Comparaison entre : (1) Sans lissage (`enable_svof_smoothing: false`), (2) Lissage diffusif (`enable_diffusion_smoothing: true, diffusion_smoothing_steps: 3`), (3) **Lissage SVOF de Helmholtz seul (`enable_svof_smoothing: true`, `svof_filter_width` $k \in \{0.5, 1.0\}$, `use_markers: false`)**, et (4) **Synergie SVOF de Helmholtz ($k = 0.5$) + Marqueurs CVI (`use_markers: true, marker_curvature: true`)** sur la finesse des spires de Kelvin-Helmholtz |

### Prototype `case/cpp/cas_05_rt_fc_wy_cvi_markers_128x384.json`
```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 384,
    "lx": 0.5,
    "ly": 1.5
  },
  "physics": {
    "enable_two_phase": true,
    "enable_surface_tension": true,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_svof_smoothing": true,
    "svof_filter_width": 0.5,
    "enable_diffusion_smoothing": false,
    "sigma": 0.005,
    "gx": 0.0,
    "gy": -9.81,
    "fluid1": { "rho": 1.225,  "mu": 3.13e-3 },
    "fluid2": { "rho": 0.1694, "mu": 3.13e-3 }
  },
  "initial_condition": {
    "type": "RAYLEIGH_TAYLOR",
    "center_y": 0.75,
    "pert_amp": 0.01,
    "pert_k": 1.0
  },
  "time": {
    "dt": 0.001,
    "t_final": 0.9,
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
    "type": "MONOLITHIC",
    "linear_solver": "BICGSTAB2",
    "tol": 1e-10,
    "max_iter": 1000,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.3,
    "use_markers": true,
    "markers_per_cell": 4,
    "marker_advection_scheme": "CVI_RK2",
    "marker_curvature": true,
    "marker_normals_plic": true,
    "marker_max_per_cell": 8,
    "write_markers_vtk": true
  },
  "output": {
    "output_dir": "Validation_Totale/Run_CASES/Run_FC-WY-CVI-Marker/05_cas_instabilite_rayleigh_taylor_cvi_markers/output/cpp/cvi_128x384",
    "case_name": "cas_05_rt_cvi_128x384",
    "vtk_interval": 25,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 4. Livrables et Post-Traitement Attendus

1. **Ajustement du taux de croissance linéaire $\gamma_{\mathrm{num}}$** comparé à $\gamma_{\mathrm{th}}$ (`figures/croissance_lineaire_rt.pdf`).
2. **Comparaison visuelle des volutes de Kelvin-Helmholtz à $t \in \{0.6, 0.8, 0.9\}\,\mathrm{s}$** (`figures/volutes_champignon_vof_vs_cvi_markers.png`) superposant le champ $C$ et les marqueurs `*_markers_*.vtk`.
3. **Rapport PDF complet (`rapport/rapport_validation_05_instabilite_rayleigh_taylor_cvi_markers.pdf`)**.
