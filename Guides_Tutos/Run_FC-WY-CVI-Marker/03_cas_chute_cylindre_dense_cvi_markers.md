# INSTRUCTIONS D'EXÉCUTION POUR AGENT / MODÈLE IA — CAS 3 (FC-WY-CVI-MARKER) : CHUTE LIBRE D'UN CYLINDRE DENSE (*FLUIDS* 2021, CAS 3.1) (`03_cas_chute_cylindre_dense_cvi_markers`)

## 1. Mission, Contraintes Critiques et Interdictions Absolues

Vous êtes un agent autonome chargé de préparer, d'exécuter, de vérifier et de post-traiter de bout en bout le **Cas de Validation 3 (`03_cas_chute_cylindre_dense_cvi_markers`)** sur les moteurs précompilés **C++17** (`COVOF_2D-C++/bin/covof2d_cpp`) et **Fortran 2018** (`COVOF_2D-FORTRAN/bin/covof2d_fortran`).

### Interdictions absolues et règles de fidélité
1. **INTERDICTION DE MODIFIER LE CODE SOURCE OU LES BINAIRES** : Ne modifiez aucun fichier dans `COVOF_2D-C++/` ni `COVOF_2D-FORTRAN/` et ne recompilez pas les binaires.
2. **CONFORMITÉ DES CHEMINS (`output_dir <= 512` car.)** : Exécutez avec `cwd = VALIDATION_TOTALE/Run_CASES/Run_FC-WY-CVI-Marker/03_cas_chute_cylindre_dense_cvi_markers`, pré-créez `output/{cpp,fortran}/<cas>`, et utilisez des chemins courts relatifs `"output/cpp/<cas>"` et `"output/fortran/<cas>"`.
3. **ZÉRO TIMEOUT MASQUÉ, ZÉRO COMPARAISON À TEMPS DIFFÉRENTS** : Allouez un budget CPU réaliste (ou privilégiez MPI `np=4` pour les solveurs monolithiques si compatible, sans `STANDARD_PROJECTION` en MPI car son Poisson est séquentiel ; ou adaptez `dt` / `t_final` / grilles en le documentant explicitement dans les paramètres). Si une simulation n'atteint pas $t_{\mathrm{final}}$, il est **strictement interdit** de la comparer avec des simulations ayant atteint un autre temps dans les mêmes tableaux ou figures sans recalage temporel exact au pas commun maximal, et il faut explicitement la marquer `TIMEOUT / INCOMPLET`.
4. **ZÉRO CONSTANTE CODÉE EN DUR** : Le champ `Iter Krylov` doit être lu **réellement** depuis la colonne 4 des `.log` stdout (jamais de valeur par défaut `12.0` ou `10.5` !). Sauvegardez `stdout`/`stderr` en flux continu dans `<cas>.log` même si un processus expire.
5. **NE PAS INVENTER D'EFFET DU LISSAGE DIFFUSIF À $\sigma = 0$** : Dans le code, `enable_diffusion_smoothing` (`C_tilde`) n'intervient **que** dans le calcul capillaire CSF quand `enable_surface_tension: true`. Puisque $\sigma = 0$ dans ce cas, seul **`enable_svof_smoothing`** ($C_s$) modifie $\rho(C_s)$ et $\mu(C_s)$. Ne créez donc pas de cas factice `enable_diffusion_smoothing` à $\sigma = 0$ en prétendant tester son effet rhéologique.

---

## 2. Arborescence et Spécification Réelle des Sorties

```text
VALIDATION_TOTALE/Run_CASES/Run_FC-WY-CVI-Marker/03_cas_chute_cylindre_dense_cvi_markers/
├── case/{cpp,fortran}/
├── output/{cpp,fortran}/
├── figures/
├── scripts/{create_cases.py,run_all_simulations.py,postprocess_analysis.py}
└── rapport/
```

### Rappel des conventions d'indexation du code
- Pour `"type": "CYLINDER"`, le cylindre dense est le **fluide 2 (`C = 0`)** et le fluide ambiant léger est le **fluide 1 (`C = 1`)**, contrairement à `DROPLET`. Le script doit utiliser `fluid1` pour le milieu ambiant et `fluid2` pour le cylindre.
- Dans `*_history.csv`, `center_of_mass_y`, `rise_velocity` ($v_c < 0$ vers le bas) et `bubble_volume` ($\sum (1-C)\Delta x\Delta y$) suivent exactement le fluide 2 (`1 - C`).
- La conservation du volume $|\Delta V(t)| / V_0$ doit être calculée depuis `PhaseIndicator` dans les `.vtr` (pas depuis la colonne CSV tronquée à 6–7 chiffres).
- Pour le défaut de circularité $\Delta R_{\max}(t)$, calculez **symétriquement pour toutes les méthodes** l'écart sur l'isocontour $C = 0.5$ du `.vtr` autour du centroïde réel $(x_c, y_c)$ pondéré par `(1-C)`, **et** rapportez séparément la métrique sur `*_markers_*.vtk` autour du même $(x_c, y_c)$ (pour éviter tout biais de barycentre des marqueurs).

---

## 3. Théorie : Chute d'un Cylindre Dense en Canal Vertical Confiné

Dans $\Omega = [0, 0.05] \times [0, 0.15]\,\mathrm{m}^2$, un cylindre fluide de rayon $R_0 = 0.01\,\mathrm{m}$ ($2 R_0 / L_x = 0.40$, confinement latéral fort) initialement centré en $(0.025, 0.12)\,\mathrm{m}$ tombe depuis le repos sous $\mathbf{g} = (0, -9.81)\,\mathrm{m/s}^2$, avec $\sigma = 0$.
- Parois latérales et haute : `FREE_SLIP` ; paroi inférieure : `DIRICHLET` ($u=v=0$).
- *Note physique sur la vitesse de chute à court temps :* Sans viscosité ni confinement, un cylindre 2D subit une masse ajoutée égale au fluide déplacé ($C_m = 1$), donnant une accélération initiale non visqueuse $a_{\mathrm{pot}} = -\frac{r_\rho - 1}{r_\rho + 1} g$ (et encore réduite par le confinement du canal $2R_0/L_x=0.40$ ainsi que par la traînée visqueuse). Le modèle $v_c(t) = -(1 - r_\rho^{-1}) g t$ néglige la masse ajoutée ($C_m=1$) et constitue une borne asymptotique haute pour $r_\rho \gg 1$.

---

## 4. Matrice Complète des Paramètres Physiques et Numériques

Tous les essais utilisent `GEAR`, `ADAMS_BASHFORTH`, `CENTERED`, `WEYMOUTH_YUE`, `PLIC`, `tvd_limiter: "NONE"`, `CYLINDER`, $\Delta t = 2\times 10^{-4}\,\mathrm{s}$, `t_final = 0.10 s` ($500$ pas), `tol = 1e-8` (ou `1e-10`), `max_iter = 1000` (sans dupliquer inutilement la même configuration sous 4 noms distincts : si une même configuration sert plusieurs axes, référencez son identifiant unique).

| Axe d'étude | Grilles $N_x \times N_y$ | Rapports $(\rho_2/\rho_1,\, \mu_2/\mu_1)$ | Lissage SVOF (`svof_filter_width`) | Solveurs et Marqueurs comparés |
| :--- | :---: | :---: | :--- | :--- |
| **Axe 1 : Convergence spatiale** | $32\times 96$, $64\times 192$ (et $128\times 384$ si budget CPU suffisant) | $(1000,\, 1000)$ | `enable_svof_smoothing: true, k: 0.5` | `WY + PLIC` seul vs `WY + PLIC + CVI_RK2 (marker_normals_plic: true)` avec `MONOLITHIC` (`BT`) |
| **Axe 2 : Effet du lissage SVOF de Helmholtz** | $64\times 192$ | $(1000,\, 1000)$ | `false` vs `true` ($k \in \{0.25, 0.5, 1.0\}$) | `MONOLITHIC` + `CVI_RK2` (mesure itérations BiCGStab2 réelles, $\Delta R_{\max}$, divergence) |
| **Axe 3 : Contrastes rhéologiques** | $64\times 192$ | $(10,\, 10)$, $(1000,\, 1000)$, $(10^6,\, 1000)$ | `true, k: 0.5` | `MONOLITHIC` + `CVI_RK2` |
| **Axe 4 : Solveurs vitesse-pression** | $64\times 192$ | $(1000,\, 1000)$ | `true, k: 0.5` | `MONOLITHIC` (`BT`), `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (`k_aal=100, uzawa_iter=4`), `STANDARD_PROJECTION` (en séquentiel) |

### Prototype `case/cpp/cas_03_axe1_cvi_64x192.json`
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
    "enable_diffusion_smoothing": false,
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
    "time_scheme": "GEAR",
    "dt": 0.0002,
    "t_final": 0.10,
    "cfl_max": 0.3,
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
    "output_dir": "output/cpp/axe1_cvi_64x192",
    "case_name": "cas_03_cvi_64x192",
    "vtk_interval": 25,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 5. Exigences de Post-Traitement, Livrables et Rédaction du Rapport

1. **Extractions strictement fidèles** :
   - Parser les itérations Krylov réelles (colonne `Iter Krylov` du `.log`), la fréquence de saturation à `max_iter`, `max_div`, $y_c(t)$, $v_c(t)$, la dérive volumique sur `.vtr` et le défaut de circularité à l'instant exact indiqué sur chaque figure et tableau.
2. **Livrables attendus** :
   - `figures/cinematique_et_circularite_cylindre.{png,pdf}` ;
   - `figures/comparaison_axe2_lissage_svof.{png,pdf}` (avec les vraies itérations Krylov et divergence) ;
   - `figures/axe3_fort_contraste_densite.{png,pdf}` ;
   - `figures/champ_cylindre_cs_et_marqueurs_64x192.{png,pdf}` (titre affichant le **temps réel exact** du fichier `.vtr`/`.vtk` tracé) ;
   - `figures/dual_code_parity_cas03.{png,pdf}`.
3. **Rapport (`rapport/rapport_validation_03_chute_cylindre_dense_cvi_markers.{md,pdf}`)** :
   - Tableau complet des paramètres physiques/numériques de toutes les simulations, tous les cas (y compris échecs/timeouts s'il y en a), sans omission sélective de lignes C++/Fortran dans le tableau de parité.
