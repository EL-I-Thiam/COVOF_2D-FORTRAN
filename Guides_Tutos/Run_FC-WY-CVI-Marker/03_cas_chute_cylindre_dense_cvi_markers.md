# INSTRUCTIONS D'EXÉCUTION POUR AGENT / MODÈLE IA — CAS 3 (FC-WY-CVI-MARKER) : CHUTE LIBRE ET ADVECTION D'UN CYLINDRE DENSE (CAS 3.1 *FLUIDS* 2021) AVEC MARQUEURS CVI (`03_cas_chute_cylindre_dense_cvi_markers`)

## 1. Mission et Rôle de l'Agent

Vous êtes un agent autonome chargé d'exécuter, de vérifier et de post-traiter de bout en bout le **Cas de Validation 3 de la méthode hybride FC-WY-CVI-Marker (`03_cas_chute_cylindre_dense_cvi_markers`)** sur les moteurs **C++17** (`covof2d_cpp`) et **Fortran 2018** (`covof2d_fortran`) de la plateforme **COVOF 2D**.

Ce cas correspond au **Cas 3.1 de [El Ouafa, Vincent & Le Chenadec (*Fluids* 2021)](https://doi.org/10.3390/fluids6010023)** enrichi par le couplage **Solveur Monolithique (FC) + VOF-PLIC (Weymouth & Yue) + Interpolation CVI + Marqueurs Non Connectés**. Il valide simultanément :
1. La préservation de la circularité géométrique d'une inclusion fluide très dense ($\rho_2/\rho_1 \in \{10, 1000, 10^6\}$) en chute libre gravitaire ou en translation uniforme sans déformation facettée ;
2. La synergie entre le résidu de divergence nul du solveur **Monolithique FC** ($\|B_u \mathbf{u}^{n+1}\|_\infty \le 10^{-10}$), le lissage rhéologique **SVOF Helmholtz** ($k = 0.5$) et l'advection lagrangienne **CVI** ($\nabla \cdot \mathbf{u}_{\mathrm{CVI}} \equiv 0$).

---

## 2. Théorie : Cinématique de Chute Libre et Circularité Sous-Maille

Dans un canal vertical $\Omega = [0, 0.05] \times [0, 0.15]\,\mathrm{m}^2$, un cylindre de fluide 2 ($C = 0$ dans la convention `CYLINDER` de `COVOF_2D`, entouré du fluide ambiant $C = 1$) de rayon $R_0 = 0.01\,\mathrm{m}$ est lâché sans vitesse initiale sous la gravité $\mathbf{g} = (0, -9.81)\,\mathrm{m/s}^2$.
Lorsque le contraste de densité et de viscosité est très grand ($\rho_{\mathrm{cyl}} / \rho_{\mathrm{amb}} \ge 1000$, $\mu_{\mathrm{cyl}} / \mu_{\mathrm{amb}} \ge 1000$) et en régime inertiel à court temps :
$$v_c(t) \approx - \left(1 - \frac{\rho_{\mathrm{amb}}}{\rho_{\mathrm{cyl}}}\right) g \, t, \qquad y_c(t) \approx y_0 - \frac{1}{2}\left(1 - \frac{\rho_{\mathrm{amb}}}{\rho_{\mathrm{cyl}}}\right) g \, t^2.$$

### 2.1 Pourquoi la méthode FC + WY + CVI + Marqueurs améliore le Cas 3.1
Sur des grilles modérées ($32 \times 96$ ou $64 \times 192$), le transport VOF-PLIC standard d'un cylindre soumis aux contraintes visqueuses d'un fluide extérieur provoque un aplatissement frontal et un facettage angulaire liés aux erreurs d'estimation de la normale $\nabla C$ près des sauts de viscosité.
Avec le couplage **FC + SVOF + WY + CVI + Marqueurs** :
- Le lissage elliptique SVOF ($C_s - D\nabla^2 C_s = C$) et la viscosité harmonique aux nœuds (`mu_node`) stabilisent le solveur monolithique FC sans diffuser $C$ ;
- Les marqueurs advectés par **CVI** maintiennent l'enveloppe circulaire du cylindre et fournissent des normales sous-maille précises (`marker_normals_plic: true`) tout au long de la chute.

---

## 3. Matrice des Paramètres et Prototype `.json`

| Paramètre | Spécification |
| :--- | :--- |
| Géométrie et Maillages | $\Omega = [0, 0.05] \times [0, 0.15]\,\mathrm{m}^2$, $R_0 = 0.01\,\mathrm{m}$, $(x_0, y_0) = (0.025, 0.12)\,\mathrm{m}$ ; Grilles $32\times 96$, $64\times 192$, $128\times 384$ |
| Rapports rhéologiques | $r_\rho = \rho_2/\rho_1 \in \{10, 1000\}$, $r_\mu = \mu_2/\mu_1 \in \{10, 1000\}$ |
| Solveurs comparés | `MONOLITHIC` (FC) vs `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (AAL) vs `STANDARD_PROJECTION` |
| Options Marqueurs | `use_markers: true`, `marker_advection_scheme: "CVI_RK2"`, `marker_normals_plic: true` vs `use_markers: false` |

### Prototype `case/cpp/cas_03_fc_wy_cvi_markers_64x192.json`
```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 64,
    "ny": 192,
    "lx": 0.05,
    "ly": 0.15
  },
  "physics": {
    "enable_two_phase": true,
    "enable_surface_tension": false,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_svof_smoothing": true,
    "svof_filter_width": 0.5,
    "sigma": 0.0,
    "gx": 0.0,
    "gy": -9.81,
    "fluid1": { "rho": 1.0,    "mu": 1.0e-3 },
    "fluid2": { "rho": 1000.0, "mu": 1.0 }
  },
  "initial_condition": {
    "type": "CYLINDER",
    "center_x": 0.025,
    "center_y": 0.12,
    "radius": 0.01
  },
  "time": {
    "dt": 0.0002,
    "t_final": 0.10,
    "time_scheme": "GEAR",
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "DIRICHLET", "u": 0.0, "v": 0.0 },
    "top":    { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 }
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
    "marker_curvature": false,
    "marker_normals_plic": true,
    "marker_max_per_cell": 8,
    "write_markers_vtk": true
  },
  "output": {
    "output_dir": "Validation_Totale/Run_CASES/Run_FC-WY-CVI-Marker/03_cas_chute_cylindre_dense_cvi_markers/output/cpp/fc_cvi_64x192",
    "case_name": "cas_03_cylinder_cvi_64x192",
    "vtk_interval": 25,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 4. Livrables et Diagnostics Attendus

1. **Suivi cinématique ($y_c(t), v_c(t)$) et défaut de circularité $\Delta R_{\max}(t) = \max_p \left|\|\mathbf{x}_p(t) - \mathbf{x}_c(t)\| - R_0\right| / R_0$** (`figures/cinematique_et_circularite_cylindre.pdf`).
2. **Bilan de conservation du volume** $|\Delta V(t)| / V_0$ confirmant une erreur $< 10^{-12}$ avec le solveur `MONOLITHIC`.
3. **Rapport PDF complet (`rapport/rapport_validation_03_chute_cylindre_dense_cvi_markers.pdf`)**.
