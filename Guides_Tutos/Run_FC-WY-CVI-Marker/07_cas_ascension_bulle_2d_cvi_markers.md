# INSTRUCTIONS D'EXÉCUTION POUR AGENT / MODÈLE IA — CAS 7 (FC-WY-CVI-MARKER) : ASCENSION D'UNE BULLE 2D DÉFORMABLE ($Bo = 10$ ET $Bo = 125$ À JUPES FINES) AVEC MARQUEURS CVI (`07_cas_ascension_bulle_2d_cvi_markers`)

## 1. Mission et Rôle de l'Agent

Vous êtes un agent autonome chargé d'exécuter, de vérifier et de post-traiter de bout en bout le **Cas de Validation 7 de la méthode hybride FC-WY-CVI-Marker (`07_cas_ascension_bulle_2d_cvi_markers`)** sur les moteurs **C++17** (`covof2d_cpp`) et **Fortran 2018** (`covof2d_fortran`) de la plateforme **COVOF 2D**.

Ce cas applique le couplage **Solveur Monolithique (FC) + VOF-PLIC (Weymouth & Yue) + Interpolation CVI + Marqueurs Non Connectés** au benchmark de référence de Hysing et al. (2009) pour l'ascension d'une bulle 2D sous capillarité :
- **Régime (i) — Bulle ellipsoïdale ($Re = 35, Bo = 10, \rho_1/\rho_2 = 10, \mu_1/\mu_2 = 10$)** ;
- **Régime (ii) — Bulle en calotte à jupes fines ($Re = 35, Bo = 125, \rho_1/\rho_2 = 1000, \mu_1/\mu_2 = 100$)**.

---

## 2. Théorie : Dynamique Capillaire et Préservation des Jupes Gazeuses ($Bo = 125$)

Dans la colonne fermée $\Omega = [0, 1.0] \times [0, 2.0]\,\mathrm{m}^2$ ($g_y = -0.98\,\mathrm{m/s}^2$, $t_{\mathrm{final}} = 3.0\,\mathrm{s}$), une bulle initialement circulaire de rayon $R_0 = 0.25\,\mathrm{m}$ est centrée en $(0.5, 0.5)\,\mathrm{m}$.
Dans le Régime (ii) ($Bo = 125, \sigma = 1.96\,\mathrm{N/m}$), le sillage tourbillonnaire étire l'arrière de la bulle en deux **jupes gazeuses latérales extrêmement minces** à partir de $t = 1.8\,\mathrm{s}$.
- Avec un schéma VOF eulérien standard, dès que l'épaisseur des jupes devient inférieure à $1.5\,\Delta x$, la courbure eulérienne $\kappa = -\nabla \cdot (\nabla \widetilde{C}/\|\nabla \widetilde{C}\|)$ dégénère et fragmente prématurément les jupes en gouttelettes satellites.
- Avec le couplage **FC + WY + CVI + Marqueurs** (`marker_curvature: true`, `marker_normals_plic: true`), l'ajustement parabolique local dans le repère tourné $(\tau, \eta)$ avec filtrage angulaire $\mathbf{n}_p \cdot \mathbf{n}_{i,j}^{\mathrm{ref}} \ge 0.5$ fournit une courbure régulière aux pointes des jupes et maintient leur continuité lagrangienne.

---

## 3. Matrice des Paramètres et Prototype `.json`

| Catégorie | Régime (i) ($Bo = 10$) | Régime (ii) ($Bo = 125$) |
| :--- | :---: | :---: |
| Propriétés physiques | $\rho_1=1000, \rho_2=100\,\mathrm{kg/m}^3$, $\mu_1=10, \mu_2=1\,\mathrm{Pa\cdot s}$, $\sigma=24.5\,\mathrm{N/m}$ | $\rho_1=1000, \rho_2=1\,\mathrm{kg/m}^3$, $\mu_1=10, \mu_2=0.1\,\mathrm{Pa\cdot s}$, $\sigma=1.96\,\mathrm{N/m}$ |
| Séquence de grilles | $40\times 80$, $80\times 160$, $160\times 320$ | $80\times 160$, $160\times 320$ |
| Configurations comparées | `VOF + CSF` standard vs `FC + WY + CVI + Marqueurs` | `VOF + CSF` standard vs `FC + WY + CVI + Marqueurs` |

### Prototype `case/cpp/cas_07_bulle_bo125_fc_cvi_markers_160x320.json`
```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 160,
    "ny": 320,
    "lx": 1.0,
    "ly": 2.0
  },
  "physics": {
    "enable_two_phase": true,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_surface_tension": true,
    "enable_svof_smoothing": true,
    "svof_filter_width": 0.5,
    "enable_diffusion_smoothing": false,
    "diffusion_smoothing_steps": 0,
    "sigma": 1.96,
    "gx": 0.0,
    "gy": -0.98,
    "fluid1": { "rho": 1000.0, "mu": 10.0 },
    "fluid2": { "rho": 1.0,    "mu": 0.1 }
  },
  "initial_condition": {
    "type": "BUBBLE",
    "center_x": 0.5,
    "center_y": 0.5,
    "radius": 0.25
  },
  "time": {
    "dt": 0.00125,
    "t_final": 3.0,
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
    "output_dir": "Validation_Totale/Run_CASES/Run_FC-WY-CVI-Marker/07_cas_ascension_bulle_2d_cvi_markers/output/cpp/bo125_cvi_160x320",
    "case_name": "cas_07_bulle_bo125_cvi_160x320",
    "vtk_interval": 80,
    "log_interval": 20,
    "write_csv_history": true
  }
}
```

---

## 4. Livrables et Post-Traitement Attendus

1. **Évolution temporelle de la vitesse d'ascension $V_c(t)$ et du centroïde $y_c(t)$** (`figures/cinematique_bulle_bo10_bo125_cvi.pdf`) comparée aux valeurs de référence de Hysing et al. (2009).
2. **Comparaison de la géométrie des jupes fines à $t = 3.0\,\mathrm{s}$ ($Bo = 125$)** (`figures/jupes_bo125_vof_vs_cvi_markers.png`) superposant le contour $C=0.5$ et le nuage de marqueurs `*_markers_*.vtk`.
3. **Rapport PDF complet (`rapport/rapport_validation_07_ascension_bulle_2d_cvi_markers.pdf`)**.
