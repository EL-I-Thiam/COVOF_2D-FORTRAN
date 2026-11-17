# INSTRUCTIONS D'EXÉCUTION POUR AGENT / MODÈLE IA — CAS 1 (FC-WY-CVI-MARKER) : GOUTTE STATIQUE DE YOUNG-LAPLACE ET RÉDUCTION DES COURANTS PARASITES (`01_cas_goutte_statique_young_laplace_cvi_markers`)

## 1. Mission, Contraintes Critiques et Interdictions Absolues

Vous êtes un agent autonome chargé de préparer, d'exécuter, de vérifier et de post-traiter de bout en bout le **Cas de Validation 1 (`01_cas_goutte_statique_young_laplace_cvi_markers`)** de la plateforme **COVOF 2D** sur les deux moteurs précompilés **C++17** (`COVOF_2D-C++/bin/covof2d_cpp`) et **Fortran 2018** (`COVOF_2D-FORTRAN/bin/covof2d_fortran`).

### Interdictions absolues et règles de fidélité
1. **INTERDICTION DE MODIFIER LE CODE SOURCE OU LES BINAIRES** : Ne modifiez aucun fichier dans `COVOF_2D-C++/` ni `COVOF_2D-FORTRAN/` et ne recompilez pas les binaires. Exécutez exclusivement les exécutables déjà compilés.
2. **CONFORMITÉ DES CHEMINS (`output_dir <= 512` car.)** : Dans `COVOF_2D-FORTRAN`, `cfg%output_dir` est dimensionné à `512` caractères (`case_name <= 128`) et `COVOF_2D-C++` crée récursivement `output_dir` (`std::filesystem::create_directories`). Par bonne pratique et lisibilité, exécutez chaque simulation avec `cwd = VALIDATION_TOTALE/Run_CASES/Run_FC-WY-CVI-Marker/01_cas_goutte_statique_young_laplace_cvi_markers`, créez au préalable `output/{cpp,fortran}/<cas>` via `os.makedirs(..., exist_ok=True)`, et utilisez des chemins courts relatifs `"output/cpp/<cas>"` et `"output/fortran/<cas>"`.
3. **SPÉCIFICATION EXPLICITE DE LA GÉOMÉTRIE (`axis_a`, `axis_b`, `pert_amp`)** : Spécifiez explicitement **`"axis_a": 0.25, "axis_b": 0.25, "pert_amp": 0.0`** dans `"initial_condition"` de tous les JSON afin de garantir de manière explicite et auto-documentée une goutte strictement circulaire de rayon $R_0 = 0.25\,\mathrm{m}$ (`is_circular_ic = true`).
4. **ZÉRO RÉSULTAT CODÉ EN DUR, ZÉRO AFFIRMATION NON JUSTIFIÉE** : Il est strictement interdit d'écrire une valeur numérique simulée, une annotation graphique, un ordre de convergence ou une conclusion en dur dans les scripts ou le rapport. Toutes les grandeurs doivent être calculées dynamiquement depuis les fichiers d'entrée `.json` et les sorties brutes réelles (`*_history.csv`, `*.vtr`, `*_markers_*.vtk`, `*.log`).
5. **FIDÉLITÉ AUX EXÉCUTIONS RÉELLES** : Vérifiez pour chaque run que `returncode == 0`, `t_end >= t_final - 0.5*dt`, et qu'aucun `NaN`/`Inf` n'est présent. En cas d'échec, de timeout, de non-convergence Krylov (`Iter Krylov == max_iter`) ou d'écart C++/Fortran, signalez-le explicitement sans fallback silencieux ni omission sélective. Comparez les méthodes **uniquement à des instants identiques**.

---

## 2. Arborescence Obligatoire et Spécification Réelle des Sorties

```text
VALIDATION_TOTALE/Run_CASES/Run_FC-WY-CVI-Marker/01_cas_goutte_statique_young_laplace_cvi_markers/
├── case/{cpp,fortran}/
├── output/{cpp,fortran}/
├── figures/
├── scripts/{create_cases.py,run_all_simulations.py,postprocess_analysis.py}
└── rapport/
```

### Sorties réellement produites par les binaires
- **`<case_name>_history.csv`** (écrit à chaque pas si `write_csv_history: true`) :
  - **C++ (`9` colonnes)** : `step,time,dt,center_of_mass_y,rise_velocity,bubble_volume,kinetic_energy,max_u,max_div`.
  - **Fortran (`11` colonnes)** : ajoute `avg_rho_gas,avg_p_gas`.
  - *Attention* : Les CSV étant formatés à 6–7 chiffres significatifs (`ES14.6`), l'erreur de conservation temporelle du volume $|V(t) - V(0)| / V(0)$ doit être calculée par sommation directe de `PhaseIndicator * dx * dy` dans les fichiers `.vtr`, en la distinguant de l'erreur géométrique initiale $|V(0) - \pi R_0^2|/(\pi R_0^2)$.
- **`<case_name>_p<rank>_<step:06d>.vtr` & `<case_name>_<step:06d>.pvtr`** : contiennent `PhaseIndicator`, `PhaseIndicator_Smoothed` (lorsque `C_s` est alloué), `Pressure`, `Density`, `Viscosity`, `Divergence`, `Temperature`, `Curvature`, `Velocity` (3 comp.).
- **`<case_name>_markers_<step:06d>.vtk`** (si `vof.use_markers && vof.write_markers_vtk`) : `POLYDATA` ASCII contenant `POINTS`, `VERTICES`, `VECTORS MarkerNormal double`, `SCALARS MarkerCurvature double 1`.
- **Journal `stdout` (`.log` à capturer intégralement même en cas de timeout)** : tableau toutes les `log_interval` itérations donnant `Pas | Temps | dt | Iter Krylov | Vitesse Max | Div Max` (+ `Nu Moyen` en Fortran). Aucun résidu Krylov `||r||` n'est imprimé en mode `NAVIER_STOKES`.

---

## 3. Théorie et Spécificités des Binaires sur ce Cas

Dans $\Omega = [0, 1] \times [0, 1]\,\mathrm{m}^2$ sans gravité ($\mathbf{g} = \mathbf{0}$), une goutte circulaire de rayon $R_0 = 0.25\,\mathrm{m}$ ($C = 1$) est centrée en $(0.5, 0.5)\,\mathrm{m}$ au repos ($\mathbf{u}(\mathbf{x}, 0) = \mathbf{0}$). À l'équilibre exact :
$$\mathbf{u}_{\mathrm{ex}} \equiv \mathbf{0}, \qquad \kappa_{\mathrm{ex}} = \frac{1}{R_0} = 4.0000\,\mathrm{m}^{-1}, \qquad \Delta P_{\mathrm{ex}} = P_{\mathrm{int}} - P_{\mathrm{ext}} = \frac{\sigma}{R_0} = 4.0000\,\mathrm{Pa} \quad (\sigma = 1.0\,\mathrm{N/m}).$$

### Spécificités critiques du code à respecter
1. **Contrainte de pas de temps capillaire explicite (Brackbill)** : Dans `COVOF_2D`, la force capillaire $\mathbf{F}_\sigma$ est traitée explicitement dans le second membre. Le pas de temps doit respecter :
   $$\Delta t < \Delta t_{\mathrm{cap}} = \sqrt{\frac{(\rho_1 + \rho_2)\, \Delta x^3}{4\pi \sigma}}.$$
   Pour $\rho_1 = \rho_2 = 1.0\,\mathrm{kg/m}^3$, $\sigma = 1.0\,\mathrm{N/m}$, $\Delta t_{\mathrm{cap}} \approx \{2.20\times 10^{-3},\, 7.78\times 10^{-4},\, 2.75\times 10^{-4},\, 9.72\times 10^{-5}\}\,\mathrm{s}$ pour $N \in \{32, 64, 128, 256\}$. Un $\Delta t$ supérieur à $\Delta t_{\mathrm{cap}}$ déstabilise le CSF eulérien avec le raffinement du maillage.
2. **Clés `vof.cfl_vof`, lissage diffusif `C_tilde` et condition initiale** :
   - Dans `COVOF_2D-C++`, `vof.cfl_vof` est ignoré (seul `time.cfl_max` est utilisé si `adaptive_dt: true`).
   - Lorsque `enable_diffusion_smoothing: true` (sans `enable_svof_smoothing`), `C_tilde` est lissé dès $t=0$ avec conditions de Neumann aux bords puis à chaque pas $n \ge 1$.
   - Lorsque `use_markers: false`, `DROPLET` utilise un profil sinusoïdal diffus sur $3h$, alors que `use_markers: true` (avec `axis_a == axis_b == 0.25`) utilise un sous-échantillonnage net $16\times 16$ : mentionnez explicitement cette différence de champ initial dans l'analyse.
3. **Condition limite `FREE_SLIP` sur les 4 bords** : La pression étant définie à une constante près (pinning sur cellule `(1,1)` en `MONOLITHIC`), relevez systématiquement le nombre d'itérations `Iter Krylov` depuis les `.log` (s'il sature à `max_iter`, rapportez-le).
4. **Lissage SVOF de Helmholtz** (`enable_svof_smoothing: true`, $C_s - k\Delta x\Delta y \nabla^2 C_s = C$) : En mode `use_markers: false`, $C_s$ remplace $C$ dans le calcul CSF de $\kappa$. En mode `use_markers: true, marker_curvature: true`, $C_s$ ne lisse que $\rho(C_s), \mu(C_s), \mu_{\mathrm{node}}(C_s)$, tandis que $\kappa^{\mathrm{marker}}$ et $\nabla C$ restent raides.

---

## 4. Matrice Complète des Paramètres Physiques et Numériques

Tous les essais sont exécutés sur `covof2d_cpp` et `covof2d_fortran` avec `ADAMS_BASHFORTH`, `GEAR`, `WEYMOUTH_YUE`, `PLIC`, `tvd_limiter: "NONE"`, `t_final = 0.02 s` (ou `0.05 s` avec $\Delta t \le 0.5 \Delta t_{\mathrm{cap}}$ afin que toutes les grilles dont $N=256$ terminent exactement à $t_{\mathrm{final}}$).

| Axe d'étude | Grilles $N\times N$ | Pas de temps $\Delta t$ [s] | Propriétés physiques | Configurations comparées |
| :--- | :---: | :---: | :--- | :--- |
| **Axe A : Convergence spatiale** | $32, 64, 128, 256$ | $10^{-3},\, 3\cdot 10^{-4},\, 10^{-4},\, 4\cdot 10^{-5}$ ($\le 0.5\,\Delta t_{\mathrm{cap}}$) | $\rho_1=\rho_2=1.0$, $\mu_1=\mu_2=1.0$, $\sigma=1.0$, $\mathbf{g}=\mathbf{0}$ | `CSF Brut` (`use_markers: false`, aucun lissage) vs `Hybride CVI-Marqueurs` (`CVI_RK2`, `marker_curvature: true`, `marker_normals_plic: false`, `markers_per_cell: 4`, `marker_max_per_cell: 8`) |
| **Axe B : Effet des lissages & Schéma marqueurs** | $64, 128$ | $3\cdot 10^{-4},\, 10^{-4}$ | Iso-propriétés ($\rho=1, \mu=1, \sigma=1$) | (1) `CSF Brut`, (2) `CSF Diffusif` (`diffusion_smoothing_steps: 3`), (3) `CSF SVOF` ($k \in \{0.5, 1.0\}$), (4) `Marqueurs BILINEAR_RK2`, (5) `Marqueurs CVI_RK2` sans SVOF, (6) `Marqueurs CVI_RK2` + `SVOF` ($k=0.5$) |
| **Axe C : Fort contraste eau/air** | $64, 128$ | $3\cdot 10^{-4},\, 10^{-4}$ (respectant $\Delta t_{\mathrm{cap,gaz}}$ si $\rho_2=1$) | $\rho_1=1000, \rho_2=1.0$, $\mu_1=0.1, \mu_2=0.001$, $\sigma=1.0$ | `CSF Brut`, `CSF SVOF (k=0.5)`, `Marqueurs CVI` sans SVOF, `Marqueurs CVI + SVOF (k=0.5)` |
| **Axe D : Solveurs vitesse-pression** | $64$ | $3\cdot 10^{-4}$ | Iso-propriétés ($\rho=1, \mu=1, \sigma=1$) | `MONOLITHIC` (`BT`), `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (`k_aal=100, uzawa_iter=4, uzawa_tol=1e-6`), `STANDARD_PROJECTION` (avec Marqueurs CVI) |

### Prototype `case/cpp/cas_01_axeA_cvi_n128.json`
*(Pour Fortran : remplacer `"output/cpp/axeA_cvi_n128"` par `"output/fortran/axeA_cvi_n128"`)*
```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 1.0,
    "ly": 1.0
  },
  "physics": {
    "enable_two_phase": true,
    "enable_surface_tension": true,
    "enable_convection": false,
    "enable_energy": false,
    "enable_svof_smoothing": false,
    "svof_filter_width": 0.5,
    "enable_diffusion_smoothing": false,
    "diffusion_smoothing_steps": 0,
    "sigma": 1.0,
    "gx": 0.0,
    "gy": 0.0,
    "fluid1": { "rho": 1.0, "mu": 1.0 },
    "fluid2": { "rho": 1.0, "mu": 1.0 }
  },
  "initial_condition": {
    "type": "DROPLET",
    "center_x": 0.5,
    "center_y": 0.5,
    "radius": 0.25,
    "axis_a": 0.25,
    "axis_b": 0.25
  },
  "time": {
    "time_scheme": "GEAR",
    "dt": 0.0001,
    "t_final": 0.02,
    "cfl_max": 0.3,
    "adaptive_dt": false
  },
  "boundary_conditions": {
    "left":   { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "right":  { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
    "bottom": { "type": "FREE_SLIP", "u": 0.0, "v": 0.0 },
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
    "marker_curvature": true,
    "marker_normals_plic": false,
    "marker_max_per_cell": 8,
    "write_markers_vtk": true
  },
  "output": {
    "output_dir": "output/cpp/axeA_cvi_n128",
    "case_name": "cas_01_cvi_n128",
    "vtk_interval": 20,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 5. Exigences de Post-Traitement, Livrables et Rédaction du Rapport

1. **Vérifications automatiques obligatoires** (`scripts/postprocess_analysis.py`) :
   - Rejeter ou marquer explicitement `INCOMPLET`/`ECHEC` tout cas n'ayant pas atteint `t_final` ou contenant `NaN`/`Inf`.
   - Extraire toutes les métriques au **même instant exact** (`t_final`) : $\langle\kappa\rangle$, $E_2(\kappa)$, $E_\infty(\kappa)$ sur cellules mixtes (`0.05 < C < 0.95`), $\Delta P = \langle p \rangle_{r<0.15} - \langle p \rangle_{r>0.38}$, $\|\mathbf{u}\|_\infty$, $Ca_{\max} = \mu_1 \|\mathbf{u}\|_\infty / \sigma$, $\|B_u \mathbf{u}\|_\infty$, $|V(t_{\mathrm{final}}) - V(0)| / V(0)$ depuis les `.vtr`, `n_markers(t_final)`, `Iter Krylov` moyen et taux de saturation à `max_iter` depuis `.log`, temps CPU mesuré réel.
   - Ne jamais utiliser de valeur factice (`1.0 s`, constante codée en dur) ; en cas d'absence de donnée, afficher `N/A`.
2. **Figures (Matplotlib / ParaView fidèles aux sorties, sans légendes superposées)** :
   - `figures/convergence_courbure_et_courants_parasites.{png,pdf}` (Axe A + Axe B à $t_{\mathrm{final}}$ identique) ;
   - `figures/comparaison_axeB_lissages.{png,pdf}` ;
   - `figures/champ_pression_cs_et_marqueurs_n128.{png,pdf}` (valeurs de $\Delta P$ calculées dynamiquement depuis `.vtr`) ;
   - `figures/axeC_fort_contraste_eau_air.{png,pdf}` ;
   - `figures/dual_code_parity_cas01.{png,pdf}`.
3. **Rapport (`rapport/rapport_validation_01_goutte_statique_young_laplace_cvi_markers.{md,pdf}`)** :
   - Inclure une **table exhaustive des paramètres physiques et numériques** lus dans les `.json` exécutés ;
   - Ne pas utiliser de balises `<kappa>` non échappées dans les tableaux Markdown (utiliser `$\langle \kappa \rangle$`) et formater les tableaux LaTeX/Pandoc sans chevauchement de colonnes ;
   - Analyser objectivement les résultats obtenus (y compris si le SVOF seul dégrade la courbure CSF ou si $\|\mathbf{u}\|_\infty$ croît encore à $t_{\mathrm{final}}$).
