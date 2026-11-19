# INSTRUCTIONS D'EXÉCUTION POUR AGENT / MODÈLE IA — CAS 7 (FC-WY-CVI-MARKER) : ASCENSION D'UNE BULLE 2D DÉFORMABLE ($Bo = 10$ ET $Bo = 125$) (`07_cas_ascension_bulle_2d_cvi_markers`)

## 1. Mission, Contraintes Critiques et Interdictions Absolues

Vous êtes un agent autonome chargé de préparer, d'exécuter, de vérifier et de post-traiter de bout en bout le **Cas de Validation 7 (`07_cas_ascension_bulle_2d_cvi_markers`)** sur les moteurs précompilés **C++17** (`COVOF_2D-C++/bin/covof2d_cpp`) et **Fortran 2018** (`COVOF_2D-FORTRAN/bin/covof2d_fortran`).

### Interdictions absolues et règles de fidélité
1. **INTERDICTION DE MODIFIER LE CODE SOURCE OU LES BINAIRES** : Ne modifiez aucun fichier dans `COVOF_2D-C++/` ni `COVOF_2D-FORTRAN/` et ne recompilez pas les binaires.
2. **CONFORMITÉ DES CHEMINS (`output_dir <= 512` car.)** : Exécutez avec `cwd = VALIDATION_TOTALE/Run_CASES/RUN_FC-WY-CVI-Marker/07_cas_ascension_bulle_2d_cvi_markers`, pré-créez `output/{cpp,fortran}/<cas>`, et utilisez des chemins courts relatifs `"output/cpp/<cas>"` et `"output/fortran/<cas>"`.
3. **SPÉCIFICATION EXPLICITE DE LA GÉOMÉTRIE (`axis_a`, `axis_b`, `pert_amp`)** : Spécifiez explicitement **`"axis_a": 0.25, "axis_b": 0.25, "pert_amp": 0.0`** dans `"initial_condition"` de tous les JSON afin de garantir une bulle strictement circulaire ($R_0 = 0.25\,\mathrm{m}$) avec sous-échantillonnage net `16x16` dans toutes les configurations.
4. **STABILITÉ CAPILLAIRE ET BUDGET CPU** : Sur $Bo = 125$ ($\rho_1=1000, \rho_2=1.0, \sigma=1.96$), vérifiez si $\Delta t$ respecte le critère capillaire côté gaz ($\Delta t_{\sigma,\mathrm{gaz}} = \sqrt{2\rho_2 h^3/(4\pi\sigma)} \approx 1.26\times 10^{-3}\,\mathrm{s}$ sur $80\times 160$ et $4.46\times 10^{-4}\,\mathrm{s}$ sur $160\times 320$). Sur `160x320`, privilégiez `np=4` en MPI (pour `MONOLITHIC` / `AAL`, jamais `STANDARD_PROJECTION`) si le temps séquentiel menace de dépasser le timeout, et rapportez tout timeout ou instabilité sans rien masquer.
5. **ZÉRO RÉSULTAT EN DUR, ZÉRO FALLBACK** : Comparez toujours les méthodes au même instant exact ; relevez les itérations Krylov depuis `.log` et la conservation volumique depuis `.vtr`.

---

## 2. Arborescence et Théorie (Benchmark de Hysing et al., 2009)

```text
VALIDATION_TOTALE/Run_CASES/RUN_FC-WY-CVI-Marker/07_cas_ascension_bulle_2d_cvi_markers/
├── case/{cpp,fortran}/
├── output/{cpp,fortran}/
├── figures/
├── scripts/{create_cases.py,run_all_simulations.py,postprocess_analysis.py}
└── rapport/
```

Dans la colonne $\Omega = [0, 1.0] \times [0, 2.0]\,\mathrm{m}^2$ ($g_y = -0.98\,\mathrm{m/s}^2$, $t_{\mathrm{final}} = 3.0\,\mathrm{s}$, parois latérales `FREE_SLIP`, parois inférieure/supérieure `DIRICHLET` $u=v=0$), une bulle de gaz circulaire (`C = 0`, fluide 2) de rayon $R_0 = 0.25\,\mathrm{m}$ est centrée en $(0.5, 0.5)\,\mathrm{m}$ :
- **Régime (i) ($Bo = 10, Re = 35$)** : $\rho_1=1000, \rho_2=100\,\mathrm{kg/m}^3$, $\mu_1=10, \mu_2=1\,\mathrm{Pa\cdot s}$, $\sigma=24.5\,\mathrm{N/m}$. Références Hysing (groupe 1) : $V_{\max} \approx 0.2417\,\mathrm{m/s}$ à $t \approx 0.921\,\mathrm{s}$, $y_c(3.0\,\mathrm{s}) \approx 1.0813\,\mathrm{m}$.
- **Régime (ii) ($Bo = 125, Re = 35$)** : $\rho_1=1000, \rho_2=1\,\mathrm{kg/m}^3$, $\mu_1=10, \mu_2=0.1\,\mathrm{Pa\cdot s}$, $\sigma=1.96\,\mathrm{N/m}$ (formation de jupes latérales minces après $t \approx 1.8\,\mathrm{s}$).
*(Note sur `enable_diffusion_smoothing: true` : `C_tilde` est lissé dès $t=0$ avec conditions de Neumann aux bords puis à chaque pas $n \ge 1$.)*

---

## 3. Matrice Complète des Paramètres Physiques et Numériques

| Régime | Grilles $N_x \times N_y$ | Pas de temps $\Delta t$ [s] (`t_final=3.0`) | Configurations comparées (`MONOLITHIC` + `WY` + `PLIC`) |
| :--- | :---: | :---: | :--- |
| **Régime (i) : $Bo = 10$** | $40\times 80$, $80\times 160$, $160\times 320$ | $0.005,\, 0.0025,\, 0.00125$ | (1) `CSF Diffusif` (`diffusion_smoothing_steps: 2`, `use_markers: false`), (2) `CSF SVOF` (`k=0.5`, `use_markers: false`), (3) `Marqueurs CVI` sans SVOF, (4) `Marqueurs CVI + SVOF (k=0.5)` (`marker_curvature: true, marker_normals_plic: true`) |
| **Régime (ii) : $Bo = 125$** | $80\times 160$, $160\times 320$ | $0.00125,\, 0.0005$ ($\le \Delta t_{\sigma,\mathrm{gaz}}$) | Même comparaison à 4 méthodes + étude `CVI_RK2` vs `BILINEAR_RK2` (`80x160`, même `marker_normals_plic: true`) + parité C++/Fortran (`80x160`) |

### Prototype `case/cpp/cas_07_bo125_cvi_svof_80x160.json`
```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 80,
    "ny": 160,
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
    "radius": 0.25,
    "axis_a": 0.25,
    "axis_b": 0.25
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
    "output_dir": "output/cpp/bo125_cvi_svof_80x160",
    "case_name": "cas_07_bo125_cvi_80",
    "vtk_interval": 160,
    "log_interval": 40,
    "write_csv_history": true
  }
}
```

---

## 4. Livrables, Post-Traitement et Rédaction du Rapport

1. **Grandeurs extraites depuis `*_history.csv`, `.vtr`, `*_markers_*.vtk` et `.log`** :
   - $y_c(t)$, $V_c(t)$, $V_{\max}$ et $t(V_{\max})$ (comparés à Hysing 2009 pour $Bo=10$) ; circularité $\psi = P_{\mathrm{cercle}} / P_{\mathrm{bulle}}$ sur le contour $C=0.5$ du `.vtr` ; conservation $|V(t)-V(0)|/V(0)$ sur `.vtr` ; `n_markers(t)` ; itérations Krylov et `max_div`.
2. **Livrables** : `figures/cinematique_bulle_bo10_bo125.{png,pdf}`, `figures/jupes_bo125_contour_et_marqueurs.{png,pdf}`, `figures/dual_code_parity_cas07.{png,pdf}`, `rapport/rapport_validation_07_ascension_bulle_2d_cvi_markers.{md,pdf}` avec table complète des paramètres JSON et analyse objective (incluant tout éventuel pincement ou différence entre méthodes/codes).
