# INSTRUCTIONS D'EXÉCUTION POUR AGENT / MODÈLE IA — CAS 6 (FC-WY-CVI-MARKER) : CISAILLEMENT PÉRIODIQUE DE KELVIN-HELMHOLTZ D'UNE NAPPE LIQUIDE (*FLUIDS* 2021, CAS 3.2) (`06_cas_nappe_liquide_kelvin_helmholtz_cvi_markers`)

## 1. Mission, Contraintes Critiques et Interdictions Absolues

Vous êtes un agent autonome chargé de préparer, d'exécuter, de vérifier et de post-traiter de bout en bout le **Cas de Validation 6 (`06_cas_nappe_liquide_kelvin_helmholtz_cvi_markers`)** sur les moteurs précompilés **C++17** (`COVOF_2D-C++/bin/covof2d_cpp`) et **Fortran 2018** (`COVOF_2D-FORTRAN/bin/covof2d_fortran`).

### Interdictions absolues et règles de fidélité
1. **INTERDICTION DE MODIFIER LE CODE SOURCE OU LES BINAIRES** : Ne modifiez aucun fichier dans `COVOF_2D-C++/` ni `COVOF_2D-FORTRAN/` et ne recompilez pas les binaires.
2. **CONFORMITÉ DES CHEMINS (`output_dir <= 512` car.)** : Exécutez avec `cwd = VALIDATION_TOTALE/Run_CASES/Run_FC-WY-CVI-Marker/06_cas_nappe_liquide_kelvin_helmholtz_cvi_markers`, pré-créez `output/{cpp,fortran}/<cas>`, et utilisez des chemins courts relatifs `"output/cpp/<cas>"` et `"output/fortran/<cas>"`.
3. **PARAMÈTRES CODÉS EN DUR DANS LES BINAIRES POUR `LIQUID_SHEET`** :
   Dans `Simulation.cpp` et `module_simulation.f90`, le profil de vitesse initial ($A_{\mathrm{gaz}} = 30.0\,\mathrm{m/s}$, $A_{\mathrm{liq}} = 2.0\,\mathrm{m/s}$, amplitude `0.04`) est codé en dur pour `LIQUID_SHEET`, et la fraction volumique initiale utilise un profil lissé sur `eps = 1.5*max(dx,dy)`, ce qui ensemence les marqueurs sur ~3 mailles d'épaisseur par interface. Reportez fidèlement ces valeurs réelles dans le rapport.
4. **VÉRIFICATION PÉRIODIQUE ET DU LISSAGE SVOF** :
   - Vérifiez objectivement la conservation de masse $|\Delta V|/V_0$ entre `covof2d_cpp` et `covof2d_fortran` ainsi qu'entre `np=1`, `np=2` et `np=4` en condition `PERIODIC` (sans supposition a priori).
   - À $\sigma = 0$, `enable_diffusion_smoothing` (`C_tilde`) est sans effet ; seul **`enable_svof_smoothing`** ($C_s$) modifie $\rho(C_s)$ et $\mu(C_s)$.
5. **ZÉRO RÉSULTAT EN DUR, ZÉRO DISSIMULATION** : Toute divergence, asymétrie périodique ou écart C++/Fortran observé dans les `.csv` / `.vtr` / `.vtk` doit être rapporté tel quel.

---

## 2. Arborescence et Théorie de la Nappe Liquide Cisaillée

```text
VALIDATION_TOTALE/Run_CASES/Run_FC-WY-CVI-Marker/06_cas_nappe_liquide_kelvin_helmholtz_cvi_markers/
├── case/{cpp,fortran}/
├── output/{cpp,fortran}/
├── figures/
├── scripts/{create_cases.py,run_all_simulations.py,postprocess_analysis.py}
└── rapport/
```

Dans $\Omega = [0, L] \times [0, L]$ ($L = 3\times 10^{-3}\,\mathrm{m}$), périodique selon $x$ et `FREE_SLIP` selon $y$, une nappe liquide horizontale ($C = 1$, fluide 1 dans `COVOF_2D` pour `LIQUID_SHEET`, $\rho_1 = 1.0\,\mathrm{kg/m}^3$, $\mu_1 = 10^{-12}\,\mathrm{Pa\cdot s}$) d'épaisseur $\delta = 3\times 10^{-4}\,\mathrm{m}$ centrée en $y_c = 1.5\times 10^{-3}\,\mathrm{m}$ est plongée dans le fluide 2 ($\rho_2 \in \{100, 1000\}\,\mathrm{kg/m}^3$, $\mu_2 = 10^{-12}\,\mathrm{Pa\cdot s}$), avec $\sigma = 0$.
*(Remarque : vérifiez dans `Density` du `.vtr` à $t=0$ quel fluide porte $C=1$ et quel fluide porte $\rho_2$, afin que l'analyse physique soit cohérente avec le champ réel.)*

---

## 3. Matrice Complète des Paramètres Physiques et Numériques

Tous les essais utilisent `test_case: "CAS_10_LIQUID_SHEET"`, `MONOLITHIC` (`BICGSTAB2`, `BLOCK_TRIANGULAR`, `tol: 1e-10`, `max_iter: 1000`), `GEAR`, `ADAMS_BASHFORTH`, `CENTERED`, `WEYMOUTH_YUE`, `PLIC`, `tvd_limiter: "NONE"`, $\Delta t = 10^{-6}\,\mathrm{s}$, `t_final = 3e-4 s` (et `5e-4 s` pour étudier l'étirement maximal) :

| Axe d'étude | Grilles | Rapports $\rho_2/\rho_1$ | Lissage SVOF (`svof_filter_width`) | Marqueurs (`use_markers` / `scheme` / `marker_normals_plic`) |
| :--- | :---: | :---: | :--- | :--- |
| **Axe A : Conservation $E_k(t)/E_k(0)$ ($r_\rho = 100, 1000$)** | $64\times 64$, $128\times 128$ | $100,\, 1000$ | `false` vs `true` ($k \in \{0.5, 1.0\}$) | `false` (`WY+PLIC` seul) vs `true` (`"CVI_RK2"`, `marker_curvature: false`, `marker_normals_plic: true`) |
| **Axe B : CVI vs Bilinéaire sur cisaillement** | $64\times 64$ | $1000$ | `true` ($k=0.5$) | `"CVI_RK2"` vs `"BILINEAR_RK2"` (avec même `marker_normals_plic: false`) |
| **Axe C : Parité dual-code (`np=1`)** | $64\times 64$, $128\times 128$ | $100,\, 1000$ | `false` et `true` ($k=0.5$) | `covof2d_cpp` vs `covof2d_fortran` |

### Prototype `case/cpp/cas_06_sheet_cvi_svof_r1000_n64.json`
```json
{
  "mode": "NAVIER_STOKES",
  "test_case": "CAS_10_LIQUID_SHEET",
  "mesh": {
    "nx": 64,
    "ny": 64,
    "lx": 0.003,
    "ly": 0.003
  },
  "physics": {
    "enable_two_phase": true,
    "enable_surface_tension": false,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_svof_smoothing": true,
    "svof_filter_width": 0.5,
    "enable_diffusion_smoothing": false,
    "sigma": 0.0,
    "gx": 0.0,
    "gy": 0.0,
    "fluid1": { "rho": 1.0,    "mu": 1e-12 },
    "fluid2": { "rho": 1000.0, "mu": 1e-12 }
  },
  "initial_condition": {
    "type": "LIQUID_SHEET",
    "center_y": 0.0015,
    "thickness": 0.0003
  },
  "boundary_conditions": {
    "left":   { "type": "PERIODIC" },
    "right":  { "type": "PERIODIC" },
    "bottom": { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "top":    { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 }
  },
  "time": {
    "time_scheme": "GEAR",
    "dt": 1e-06,
    "t_final": 0.0003,
    "cfl_max": 0.7,
    "adaptive_dt": false
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.45,
    "use_markers": true,
    "markers_per_cell": 4,
    "marker_advection_scheme": "CVI_RK2",
    "marker_curvature": false,
    "marker_normals_plic": true,
    "marker_max_per_cell": 8,
    "write_markers_vtk": true
  },
  "solver": {
    "type": "MONOLITHIC",
    "linear_solver": "BICGSTAB2",
    "tol": 1e-10,
    "max_iter": 1000,
    "use_precond": true,
    "precond_type": "BLOCK_TRIANGULAR"
  },
  "output": {
    "output_dir": "output/cpp/sheet_cvi_svof_r1000_n64",
    "case_name": "cas_06_sheet_r1000_n64",
    "vtk_interval": 50,
    "log_interval": 25,
    "write_csv_history": true
  }
}
```

---

## 4. Livrables, Post-Traitement et Rédaction du Rapport

1. **Extraction fidèle** : Tracer $E_k(t)/E_k(0)$ (colonne `kinetic_energy` de `*_history.csv`), `max_div`, itérations Krylov depuis `.log`, nombre de marqueurs `n_markers(t)` et traversée périodique $x=0 \leftrightarrow x=L$, dérive volumique $|V(t)-V(0)|/V(0)$ depuis `.vtr`.
2. **Livrables** : `figures/conservation_energie_cinetique_nappe.{png,pdf}`, `figures/nappe_liquide_rideaux_marqueurs_cvi.{png,pdf}`, `figures/dual_code_parity_cas06.{png,pdf}`, `rapport/rapport_validation_06_nappe_liquide_kelvin_helmholtz_cvi_markers.{md,pdf}` avec table complète des paramètres JSON et compte-rendu fidèle des éventuels écarts C++/Fortran ou de périodicité.
