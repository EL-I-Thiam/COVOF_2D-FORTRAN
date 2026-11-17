# INSTRUCTIONS D'EXÉCUTION POUR AGENT / MODÈLE IA — CAS 2 (FC-WY-CVI-MARKER) : TOURBILLON CISAILLÉ RÉVERSIBLE DE RIDER-KOTHE (`02_cas_tourbillon_cisaille_rider_kothe_cvi_markers`)

## 1. Mission, Contraintes Critiques et Interdictions Absolues

Vous êtes un agent autonome chargé de préparer, d'exécuter, de vérifier et de post-traiter de bout en bout le **Cas de Validation 2 (`02_cas_tourbillon_cisaille_rider_kothe_cvi_markers`)** sur les moteurs précompilés **C++17** (`COVOF_2D-C++/bin/covof2d_cpp`) et **Fortran 2018** (`COVOF_2D-FORTRAN/bin/covof2d_fortran`).

### Interdictions absolues et règles de fidélité
1. **INTERDICTION DE MODIFIER LE CODE SOURCE OU LES BINAIRES** : Ne modifiez aucun fichier dans `COVOF_2D-C++/` ni `COVOF_2D-FORTRAN/` et ne recompilez pas les binaires.
2. **CONFORMITÉ DES CHEMINS (`output_dir <= 512` car.)** : Exécutez chaque simulation avec `cwd = VALIDATION_TOTALE/Run_CASES/Run_FC-WY-CVI-Marker/02_cas_tourbillon_cisaille_rider_kothe_cvi_markers`, pré-créez `output/{cpp,fortran}/<cas>`, et utilisez des chemins courts relatifs `"output/cpp/<cas>"` et `"output/fortran/<cas>"`. Purge totale requise avant lancement (aucun fichier doublon `* 2.*`).
3. **SPÉCIFICATION EXPLICITE DE LA GÉOMÉTRIE (`axis_a`, `axis_b`, `pert_amp`)** : En mode `VORTEX_ADVECTION`, spécifiez explicitement **`"axis_a": 0.125, "axis_b": 0.125, "pert_amp": 0.0`** dans `"initial_condition"` avec `"type": "DROPLET"` (pour comparer de manière strictement iso-IC à $t=0$). Vérifiez à $t=0$ dans les `.vtr` que le disque initial est identique pour toutes les configurations.
4. **COMPARAISON NON BIAISÉE SUR L'AXE B (`CVI_RK2` vs `BILINEAR_RK2`)** : Gardez strictement la même valeur de `marker_normals_plic` (testez d'abord `false` pour isoler purement l'effet CVI vs Bilinéaire sur les trajectoires, puis éventuellement `true`) entre `CVI_RK2` et `BILINEAR_RK2`.
5. **ZÉRO RÉSULTAT EN DUR, ZÉRO CHERRY-PICKING** : Ne supprimez jamais de lignes ni de colonnes défavorables des tableaux CSV dans le rapport final. Reportez systématiquement $E_{R,\mathrm{mean}}$, $E_{R,\max}$, la distance de Hausdorff exacte $d_H$, les écarts C++/Fortran et les écarts `NP=1` vs `NP=4` tels qu'ils sortent des simulations.

---

## 2. Arborescence et Spécification Réelle des Sorties en Mode `VORTEX_ADVECTION`

```text
VALIDATION_TOTALE/Run_CASES/Run_FC-WY-CVI-Marker/02_cas_tourbillon_cisaille_rider_kothe_cvi_markers/
├── case/{cpp,fortran}/
├── output/{cpp,fortran}/
├── figures/
├── scripts/{create_cases.py,run_all_simulations.py,postprocess_analysis.py}
└── rapport/
```

### Sorties réellement produites par les binaires en mode `VORTEX_ADVECTION`
- **`vortex_error.csv`** (`13` colonnes) :
  `nx,ny,dt,T,advection,reconstruction,limiter,err_l1,err_l2,err_linf,mass_init,mass_final,mass_conservation`
  avec $E_{L_1} = \frac{1}{N_x N_y}\sum_{i,j}|C_{i,j}(T)-C_{i,j}(0)|$ et `mass_conservation = mass_final / mass_init`.
- **`<case_name>_history.csv`** : en mode `VORTEX_ADVECTION`, les deux binaires créent l'en-tête uniquement (**0 ligne de données**).
- **`<case_name>_p<rank>_<step:06d>.vtr` & `<case_name>_<step:06d>.pvtr`** + **`<case_name>_markers_<step:06d>.vtk`** : écrits à `step = 0`, tous les `vtk_interval`, à $t \approx T/2 = 3.0\,\mathrm{s}$ et à $t \ge T = 6.0\,\mathrm{s}$. *(Vérifiez qu'en `np=4`, toutes les pièces `p0..p3` existent bien aux instants analysés.)*

---

## 3. Théorie et Spécificités Implémentées dans `COVOF_2D`

Dans $\Omega = [0, 1] \times [0, 1]\,\mathrm{m}^2$, un disque de fluide 1 ($C = 1$) de rayon $R_0 = 0.125\,\mathrm{m}$ est centré en $(x_c, y_c) = (0.5, 0.75)\,\mathrm{m}$ à $t = 0$. Le champ de vitesse analytique solénoïdal est imposé à chaque pas de temps avec $T = 6.0\,\mathrm{s}$ :
$$u(x, y, t) = -\sin^2(\pi x)\sin(2\pi y)\cos\left(\frac{\pi t}{T}\right), \qquad v(x, y, t) = +\sin(2\pi x)\sin^2(\pi y)\cos\left(\frac{\pi t}{T}\right).$$

### Comportement réel des binaires actuels à prendre en compte dans l'analyse
1. **Intégration temporelle au point milieu dans `advect_rk2`** : En `VORTEX_ADVECTION`, les vitesses analytiques `u, v` sont évaluées à $t^n$ puis réévaluées au point milieu $t^n + 0.5\Delta t$ pour l'étape 2 de `advect_rk2` (ordre 2 en temps).
2. **Cohérence angulaire, fusion $2\times 2$ & synchronisation MPI** :
   - `apply_marker_normals_to_segments` filtre les normales par rapport au gradient de Youngs `3x3` avec le seuil `dot_prod >= 0.3` (tandis que `compute_curvature_from_markers` utilise `0.5`).
   - En cas de surpopulation (`> marker_max_per_cell`), C++ et Fortran utilisent tous deux une fusion lagrangienne pondérée par sous-quadrants $2\times 2$ préservant la cohérence des normales.
   - En MPI, l'étape 2 de `advect_rk2` utilise une évaluation CVI en deux passes à travers les raccords de sous-domaines et synchronise les segments/normales des mailles fantômes pour `advect_weymouth_yue`. **Vérifiez objectivement** la conservation de masse (`mass_conservation`) et l'invariance `np=1` vs `np=4`.

---

## 4. Matrice Complète des Paramètres Physiques et Numériques

En `VORTEX_ADVECTION`, `physics` et `boundary_conditions` sont lus pour conformité JSON mais la cinématique est pilotée par `T = 6.0 s` et `dt`. Exécutez chaque configuration sur `covof2d_cpp` et `covof2d_fortran` :

| Axe d'étude | Grilles $N\times N$ | Pas de temps $\Delta t$ [s] (`t_final=6.0`) | `vtk_interval` | `use_markers` / `scheme` / `marker_normals_plic` | MPI `np` |
| :--- | :---: | :---: | :---: | :--- | :---: |
| **Axe A1 (Référence VOF)** | $64, 128, 256, 512$ | $0.005,\, 0.0025,\, 0.00125,\, 0.0005$ | $60,\, 120,\, 240,\, 600$ | `false` / — / `false` (`WY + PLIC`, `tvd_limiter: "NONE"`) | `1` |
| **Axe A2 (Hybride CVI)** | $64, 128, 256, 512$ | $0.005,\, 0.0025,\, 0.00125,\, 0.0005$ | $60,\, 120,\, 240,\, 600$ | `true` / `"CVI_RK2"` / `true` (`markers_per_cell: 4`, `max: 8`, `marker_curvature: false`) | `1` |
| **Axe B (CVI vs Bilin iso-config)** | $128, 256$ | $0.0025,\, 0.00125$ | $120,\, 240$ | Comparer `"CVI_RK2"` vs `"BILINEAR_RK2"` **avec le même `marker_normals_plic: false`** (puis `true` si testé) | `1` |
| **Axe C (Invariance MPI)** | $128, 256$ | $0.0025,\, 0.00125$ | $120,\, 240$ | `WY+PLIC` seul (`false`) et `CVI+Marqueurs` (`true`, `"CVI_RK2"`, `marker_normals_plic: true`) | `1` vs `4` (`--oversubscribe`) |

### Prototype `case/cpp/cas_02_axeA2_cvi_n128.json`
```json
{
  "mode": "VORTEX_ADVECTION",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 1.0,
    "ly": 1.0
  },
  "physics": {
    "enable_two_phase": true,
    "enable_surface_tension": false,
    "enable_convection": true,
    "sigma": 0.0,
    "fluid1": { "rho": 1000.0, "mu": 0.001 },
    "fluid2": { "rho": 1.0,    "mu": 1.0e-5 }
  },
  "initial_condition": {
    "type": "DROPLET",
    "center_x": 0.5,
    "center_y": 0.75,
    "radius": 0.125,
    "axis_a": 0.125,
    "axis_b": 0.125
  },
  "time": {
    "dt": 0.0025,
    "t_final": 6.0,
    "cfl_max": 0.32,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "top":    { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 }
  },
  "vof": {
    "advection_method": "WEYMOUTH_YUE",
    "reconstruction_method": "PLIC",
    "tvd_limiter": "NONE",
    "cfl_vof": 0.32,
    "use_markers": true,
    "markers_per_cell": 4,
    "marker_advection_scheme": "CVI_RK2",
    "marker_curvature": false,
    "marker_normals_plic": true,
    "marker_max_per_cell": 8,
    "write_markers_vtk": true
  },
  "output": {
    "output_dir": "output/cpp/axeA2_cvi_n128",
    "case_name": "cas_02_cvi_n128",
    "vtk_interval": 120,
    "log_interval": 120
  }
}
```

---

## 5. Exigences de Post-Traitement, Livrables et Rédaction du Rapport

1. **Métriques calculées sans approximation grossière ni fallback** :
   - Lire `err_l1`, `err_l2`, `err_linf`, `mass_init`, `mass_final`, `mass_conservation` depuis `vortex_error.csv`, et vérifier la cohérence par intégration des `.vtr` à $t = 0$ et $t = 6.0\,\mathrm{s}$.
   - Sur le dernier fichier `*_markers_*.vtk` exact (`t = 6.0 s`, vérifier le numéro de pas maximal sans mélanger de doublons), calculer $E_{R,\mathrm{mean}}(T)$, $E_{R,\mathrm{med}}(T)$, $E_{R,\max}(T)$ et la distance de Hausdorff symétrique exacte $d_H$ face à un cercle discrétisé à $\ge 1000$ points.
2. **Livrables attendus** :
   - `figures/spirale_t3s_et_retour_t6s_vof_vs_cvi.{png,pdf}` (superposition $C=0.5$ et marqueurs CVI à $t=3.0\,\mathrm{s}$ et $t=6.0\,\mathrm{s}$, signalant tout marqueur orphelin) ;
   - `figures/axeA_convergence_spatiale_l1.{png,pdf}` (avec ordres de convergence réels calculés entre grilles consécutives) ;
   - `figures/axeB_cvi_vs_bilinear_radial_error.{png,pdf}` ;
   - `figures/axeC_mpi_invariance_profils.{png,pdf}` ;
   - `figures/dual_code_parity_cas02.{png,pdf}`.
3. **Rapport (`rapport/rapport_validation_02_tourbillon_cisaille_rider_kothe_cvi_markers.{md,pdf}`)** :
   - Inclure la table complète des paramètres physiques/numériques de tous les `.json`, les tableaux intégraux (aucune ligne/colonne omise), les temps CPU réels, et une critique honnête des écarts observés (ex. présence de marqueurs orphelins, pertes de conservation MPI ou divergences `CVI` vs `BILINEAR` et C++ vs Fortran si elles apparaissent).
