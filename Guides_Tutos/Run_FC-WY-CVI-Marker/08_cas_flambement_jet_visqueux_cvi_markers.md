# INSTRUCTIONS D'EXÉCUTION POUR AGENT / MODÈLE IA — CAS 8 (FC-WY-CVI-MARKER) : FLAMBEMENT D'UN JET LIQUIDE VISQUEUX 2D (*FLUIDS* 2021, CAS 3.3) (`08_cas_flambement_jet_visqueux_cvi_markers`)

## 1. Mission, Contraintes Critiques et Interdictions Absolues

Vous êtes un agent autonome chargé de préparer, d'exécuter, de vérifier et de post-traiter de bout en bout le **Cas de Validation 8 (`08_cas_flambement_jet_visqueux_cvi_markers`)** sur les moteurs précompilés **C++17** (`COVOF_2D-C++/bin/covof2d_cpp`) et **Fortran 2018** (`COVOF_2D-FORTRAN/bin/covof2d_fortran`).

### Interdictions absolues et règles de fidélité
1. **INTERDICTION DE MODIFIER LE CODE SOURCE OU LES BINAIRES** : Ne modifiez aucun fichier dans `COVOF_2D-C++/` ni `COVOF_2D-FORTRAN/` et ne recompilez pas les binaires.
2. **CONFORMITÉ DES CHEMINS (`output_dir <= 512` car.)** : Exécutez avec `cwd = VALIDATION_TOTALE/Run_CASES/Run_FC-WY-CVI-Marker/08_cas_flambement_jet_visqueux_cvi_markers`, pré-créez `output/{cpp,fortran}/<cas>`, et utilisez des chemins courts relatifs `"output/cpp/<cas>"` et `"output/fortran/<cas>"`.
3. **VÉRIFICATION DES MARQUEURS ET DU RÉ-ENSEMENCEMENT À L'INJECTEUR `JET_INLET`** :
   Dans les binaires C++ et Fortran actualisés, lorsque `vof.use_markers` est activé, la condition initiale `"VISCOUS_JET"` utilise un sous-échantillonnage $16\times 16$ par cellule autour du bouchon d'injection (`center_x`, `jet_width`) et `prune_and_reseed` ré-ensemence les flancs du jet entrant au bord supérieur `JET_INLET` ($y = L_y$).
   **Obligation absolue :** Le script de post-traitement doit vérifier et rapporter l'évolution réelle `n_markers(t)` depuis `*_markers_*.vtk` et vérifier que les marqueurs suivent bien le développement du jet entrant (ou signaler tout comportement réel observé).
4. **RÔLE DU LISSAGE SVOF ($\sigma = 0$)** : Puisque $\sigma = 0$, `enable_diffusion_smoothing` (`C_tilde`) n'a aucun effet dans les équations ; testez uniquement l'effet de **`enable_svof_smoothing`** (`false` vs `true` pour $k \in \{0.5, 1.0, 1.5\}$).
5. **ZÉRO RÉSULTAT EN DUR, ZÉRO TIMEOUT MASQUÉ** : Toute valeur de divergence, itérations Krylov depuis `.log`, loi de volume injecté $V_{\mathrm{jet}}(t)$ ou non-convergence doit être rapportée fidèlement.

---

## 2. Arborescence et Théorie du Flambement Visqueux

```text
VALIDATION_TOTALE/Run_CASES/Run_FC-WY-CVI-Marker/08_cas_flambement_jet_visqueux_cvi_markers/
├── case/{cpp,fortran}/
├── output/{cpp,fortran}/
├── figures/
├── scripts/{create_cases.py,run_all_simulations.py,postprocess_analysis.py}
└── rapport/
```

Dans une cavité carrée $\Omega = [0, 1] \times [0, 1]\,\mathrm{m}^2$ initialement remplie d'air ($\rho_1 = 1.1768\,\mathrm{kg/m}^3$, $\mu_1 = 1.85\times 10^{-3}\,\mathrm{Pa\cdot s}$, $C = 1$), un jet très visqueux ($\rho_2 = 1800.0\,\mathrm{kg/m}^3$, $\mu_2 \in \{300, 500\}\,\mathrm{Pa\cdot s}$, $C = 0$) de largeur $D$ est injecté en $(x_c = 0.5\,\mathrm{m}, y = 1.0\,\mathrm{m})$ à vitesse $v_{\mathrm{in}} = -V_0$ sous $g_y = -9.81\,\mathrm{m/s}^2$, avec $\sigma = 0$ et débit compensé exact sur `top` (`JET_INLET`, $v_{\mathrm{out}} = +V_0 \frac{D}{L_x - D}$).
- Critère de flambement : $Re = \frac{\rho_2 D V_0}{\mu_2} \le Re_c \approx 0.52$ et $H/D > 3\pi \approx 9.4248$.
- Volume théorique de fluide 2 injecté : $V_{\mathrm{jet}}^{\mathrm{th}}(t) = V_{\mathrm{jet}}(0) + D \, V_0 \, t$ (vérifié contre `bubble_volume` de `*_history.csv` et l'intégrale de `1-C` sur `.vtr`).

---

## 3. Matrice Complète des Paramètres Physiques et Numériques

Tous les essais utilisent `GEAR`, `ADAMS_BASHFORTH`, `CENTERED`, `WEYMOUTH_YUE`, `PLIC`, `tvd_limiter: "NONE"`, `VISCOUS_JET`, `JET_INLET`, $\Delta t = 5\times 10^{-4}\,\mathrm{s}$, grille $80\times 80$ et $160\times 160$ (privilégier `np=4` en MPI si le budget temps séquentiel l'exige pour `MONOLITHIC`/`AAL`) :

| Sous-cas | $D$ [m] | $\mu_2$ [Pa$\cdot$s] | $V_0$ [m/s] | $Re$ | $t_{\mathrm{final}}$ [s] | Configurations testées (`tol: 1e-8, max_iter: 1000`) |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **Sous-cas 1 (Flambement)** | $0.08$ | $500.0$ | $0.60$ | $0.1728$ | $2.0$ | `MONOLITHIC` + `SVOF` ($k \in \{0.5, 1.0, 1.5\}$ vs `false`) ; `use_markers: true` (`CVI_RK2`) vs `false` ; comparaison avec `AAL` (`k_aal=100, uzawa_iter=4`) et `STANDARD_PROJECTION` (seq) ; parité C++/Fortran |
| **Sous-cas 2 (Stable)** | $0.10$ | $300.0$ | $1.00$ | $0.6000$ | $2.0$ | `MONOLITHIC` + `SVOF (k=1.0)` + `WY + PLIC` (`use_markers: true` et `false`) |
| **Sous-cas 3 (Inertiel)** | $0.10$ | $300.0$ | $6.00$ | $3.6000$ | $0.9$ | `MONOLITHIC` + `SVOF (k=1.0)` + `WY + PLIC` (`use_markers: true` et `false`) |

### Prototype `case/cpp/cas_08_re01728_fc_svof_k10_160.json`
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
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_svof_smoothing": true,
    "svof_filter_width": 1.0,
    "enable_diffusion_smoothing": false
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
    "output_dir": "output/cpp/re01728_svof_k10_160",
    "case_name": "cas_08_re01728_160",
    "vtk_interval": 100,
    "log_interval": 20,
    "write_csv_history": true
  }
}
```

---

## 4. Livrables, Post-Traitement et Rédaction du Rapport

1. **Vérifications obligatoires** : `V_jet(t)` mesuré sur `.vtr` vs `V_jet(0) + D*V0*t`, `max_div(t)`, itérations Krylov depuis `.log`, `n_markers(t)` réel depuis `*_markers_*.vtk`, amplitude latérale de flambement $x_{\mathrm{jet}}(y,t)$ sur l'iso-contour $C=0.5$.
2. **Livrables** : `figures/sequence_flambement_re01728.{png,pdf}`, `figures/volume_injecte_et_divergence.{png,pdf}`, `figures/diagramme_stabilite_et_svof.{png,pdf}`, `figures/dual_code_parity_cas08.{png,pdf}`, `rapport/rapport_validation_08_flambement_jet_visqueux_cvi_markers.{md,pdf}` avec table complète des paramètres JSON et transparence totale sur `n_markers`, les solveurs et la parité.
