# INSTRUCTIONS D'EXÉCUTION POUR AGENT / MODÈLE IA — CAS 8 (FC-WY-CVI-MARKER) : FLAMBEMENT D'UN JET LIQUIDE VISQUEUX 2D (CAS 3.3 *FLUIDS* 2021) AVEC MARQUEURS CVI (`08_cas_flambement_jet_visqueux_cvi_markers`)

## 1. Mission et Rôle de l'Agent

Vous êtes un agent autonome chargé d'exécuter, de vérifier et de post-traiter de bout en bout le **Cas de Validation 8 de la méthode hybride FC-WY-CVI-Marker (`08_cas_flambement_jet_visqueux_cvi_markers`)** sur les moteurs **C++17** (`covof2d_cpp`) et **Fortran 2018** (`covof2d_fortran`) de la plateforme **COVOF 2D**.

Ce cas correspond au **Cas 3.3 de [El Ouafa, Vincent & Le Chenadec (*Fluids* 2021)](https://doi.org/10.3390/fluids6010023)** couplé à l'architecture **Solveur Monolithique (FC) + VOF-PLIC (Weymouth & Yue) + Interpolation CVI + Marqueurs Non Connectés**. Il constitue le test ultime de robustesse rhéologique et géométrique du code ($\rho_2/\rho_1 = 1529.6$, $\mu_2/\mu_1 = 2.70\times 10^5$) :
1. Validation du diagramme de stabilité du flambement visqueux ($Re = 0.1728 \le 0.52$ instable vs $Re = 0.60$ et $Re = 3.60$ symétriques stables) ;
2. Suivi lagrangien sans perte de divergence (`CVI_RK2` sur le champ monolithique `FC`) des multiples replis en accordéon du cordon visqueux et des fines lames d'air piégées lors de l'empilement sur la paroi inférieure ($t \in [1.0, 2.0]\,\mathrm{s}$).

---

## 2. Théorie : Instabilité de Flambement Visqueux et Suivi des Replis par Marqueurs CVI

Dans une cavité carrée $\Omega = [0, 1] \times [0, 1]\,\mathrm{m}^2$ initialement remplie d'air ($\rho_1 = 1.1768\,\mathrm{kg/m}^3$, $\mu_1 = 1.85\times 10^{-3}\,\mathrm{Pa\cdot s}$, $C = 1$), un jet plan très visqueux ($\rho_2 = 1800.0\,\mathrm{kg/m}^3$, $\mu_2 = 500.0\,\mathrm{Pa\cdot s}$, $C = 0$) de largeur $D = 0.08\,\mathrm{m}$ est injecté verticalement depuis le centre de la paroi supérieure ($x_c = 0.5\,\mathrm{m}, y = 1.0\,\mathrm{m}$) à la vitesse $v_{\mathrm{in}} = -V_0 = -0.60\,\mathrm{m/s}$ sous la gravité $g_y = -9.81\,\mathrm{m/s}^2$.
Le flambement périodique se déclenche lorsque :
$$Re = \frac{\rho_2 D V_0}{\mu_2} \le Re_c \approx 0.52 \qquad \text{et} \qquad \frac{H}{D} > 3\pi \approx 9.4248.$$

### 2.1 Rôle du couplage FC + WY + CVI + Marqueurs lors de l'empilement des plis
1. **Solveur Monolithique FC + SVOF Helmholtz** : À $\mu_2/\mu_1 = 2.7\times 10^5$, la méthode de projection standard diverge vers $t \approx 1.4\,\mathrm{s}$ lors de l'impact des plis, tandis que le solveur monolithique **FC** maintient $\|B_u \mathbf{u}^{n+1}\|_\infty \le 10^{-10}$ jusqu'à $t = 2.0\,\mathrm{s}$.
2. **Ré-ensemencement automatique à l'injecteur `JET_INLET` et filtrage angulaire dans les plis** : À chaque pas de temps, les mailles d'interface au bord de l'injecteur supérieur ($y = 1.0\,\mathrm{m}$) sont automatiquement ré-ensemencées en marqueurs, et le filtrage directionnel $\mathbf{n}_p \cdot \mathbf{n}_{i,j}^{\mathrm{ref}} \ge 0.5$ préserve la séparation des lames d'air entre deux boucles repliées du jet.

---

## 3. Matrice des Trois Sous-Cas et Prototype `.json`

| Sous-cas | $D$ ($\mathrm{m}$) | $\mu_2$ ($\mathrm{Pa\cdot s}$) | $V_0$ ($\mathrm{m/s}$) | $Re$ | Régime physique & Étude du Lissage SVOF |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Sous-cas 1** | $0.08$ | $500.0$ | $0.60$ | $0.1728$ | **Flambement périodique instable** ($t_{\mathrm{final}} = 2.0\,\mathrm{s}$) — Tester : (1) Sans lissage (`enable_svof_smoothing: false`), (2) Lissage diffusif (`enable_diffusion_smoothing: true`), et (3) **Lissage elliptique SVOF de Helmholtz (`enable_svof_smoothing: true`, `svof_filter_width` $k \in \{0.5, 1.0, 1.5\}$)** couplé à `CVI_RK2` |
| **Sous-cas 2** | $0.10$ | $300.0$ | $1.00$ | $0.6000$ | **Étalement symétrique stable** sans flambement ($t_{\mathrm{final}} = 2.0\,\mathrm{s}$) avec `enable_svof_smoothing: true` ($k = 1.0$) + `CVI_RK2` |
| **Sous-cas 3** | $0.10$ | $300.0$ | $6.00$ | $3.6000$ | **Impact inertiel symétrique** et remontée pariétale ($t_{\mathrm{final}} = 0.9\,\mathrm{s}$) avec `enable_svof_smoothing: true` ($k = 1.0$) + `CVI_RK2` |

### Prototype `case/cpp/cas_08_re01728_fc_wy_cvi_markers_160x160.json`
```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 160,
    "ny": 160,
    "lx": 1.0,
    "ly": 1.0
  },
  "physics": {
    "enable_two_phase": true,
    "fluid1": { "rho": 1.1768, "mu": 0.00185 },
    "fluid2": { "rho": 1800.0, "mu": 500.0 },
    "gx": 0.0,
    "gy": -9.81,
    "sigma": 0.0,
    "enable_surface_tension": false,
    "convective_scheme": "CENTERED",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_svof_smoothing": true,
    "svof_filter_width": 1.0
  },
  "time": {
    "dt": 0.0005,
    "t_final": 2.0,
    "time_scheme": "GEAR",
    "adaptive_dt": false,
    "cfl_max": 0.3
  },
  "initial_condition": {
    "type": "VISCOUS_JET",
    "center_x": 0.5,
    "center_y": 1.0,
    "jet_width": 0.08
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "top": {
      "type": "DIRICHLET",
      "profile": "JET_INLET",
      "jet_center": 0.5,
      "jet_width": 0.08,
      "u": 0.0,
      "v": -0.6
    }
  },
  "solver": {
    "type": "MONOLITHIC",
    "r_al": 0.0,
    "linear_solver": "BICGSTAB2",
    "tol": 1e-08,
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
    "marker_curvature": false,
    "marker_normals_plic": true,
    "marker_max_per_cell": 8,
    "write_markers_vtk": true
  },
  "output": {
    "output_dir": "Validation_Totale/Run_CASES/Run_FC-WY-CVI-Marker/08_cas_flambement_jet_visqueux_cvi_markers/output/cpp/re01728_fc_cvi_160x160",
    "case_name": "cas_08_re01728_fc_cvi_160x160",
    "vtk_interval": 100,
    "log_interval": 20,
    "write_csv_history": true
  }
}
```

---

## 4. Livrables et Post-Traitement Attendus

1. **Séquence temporelle du flambement ($t \in \{0.25, 0.75, 1.25, 2.0\}\,\mathrm{s}$)** (`figures/sequence_flambement_re01728_cvi_markers.png`) superposant le champ $C$, les lignes de courant et les marqueurs `*_markers_*.vtk`.
2. **Vérification de la loi d'injection linéaire du volume $V_{\mathrm{jet}}(t) = V_{\mathrm{jet}}(0) + D V_0 t$ et de la divergence $\|\nabla\cdot\mathbf{u}\|_\infty(t)$** (`figures/volume_injecte_et_divergence_fc_cvi.pdf`).
3. **Rapport PDF complet (`rapport/rapport_validation_08_flambement_jet_visqueux_cvi_markers.pdf`)**.
