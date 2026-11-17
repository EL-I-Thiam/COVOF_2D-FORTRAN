# INSTRUCTIONS D'EXÉCUTION POUR AGENT / MODÈLE IA — CAS 4 (FC-WY-CVI-MARKER) : OSCILLATIONS CAPILLAIRES D'UNE BULLE DE LAMB À GRAND RAPPORT DE DENSITÉ (`04_cas_oscillation_bulle_lamb_cvi_markers`)

## 1. Mission, Contraintes Critiques et Interdictions Absolues

Vous êtes un agent autonome chargé de préparer, d'exécuter, de vérifier et de post-traiter de bout en bout le **Cas de Validation 4 (`04_cas_oscillation_bulle_lamb_cvi_markers`)** sur les moteurs précompilés **C++17** (`COVOF_2D-C++/bin/covof2d_cpp`) et **Fortran 2018** (`COVOF_2D-FORTRAN/bin/covof2d_fortran`).

### Interdictions absolues et règles de fidélité
1. **INTERDICTION DE MODIFIER LE CODE SOURCE OU LES BINAIRES** : Ne modifiez aucun fichier dans `COVOF_2D-C++/` ni `COVOF_2D-FORTRAN/` et ne recompilez pas les binaires.
2. **CONFORMITÉ DES CHEMINS (`output_dir <= 512` car.)** : Exécutez avec `cwd = VALIDATION_TOTALE/Run_CASES/Run_FC-WY-CVI-Marker/04_cas_oscillation_bulle_lamb_cvi_markers`, pré-créez `output/{cpp,fortran}/<cas>`, et utilisez des chemins courts relatifs `"output/cpp/<cas>"` et `"output/fortran/<cas>"`.
3. **INITIALISATION ELLIPTIQUE EXPLICITE PAR `axis_a` ET `axis_b`** : Spécifiez explicitement `"axis_a"`, `"axis_b"` et `"pert_amp": 0.0` dans `"initial_condition"` :
   - Pour $\varepsilon = 0.05$ ($R_0 = 0.02\,\mathrm{m}$) : `"axis_a": 0.021`, `"axis_b": 0.01904761904761905` ;
   - Pour $\varepsilon = 0.20$ (régime non linéaire) : `"axis_a": 0.024`, `"axis_b": 0.016666666666666666`.
   Vérifiez dans le `.vtr` à $t=0$ que $a_0/R_0 \approx 1+\varepsilon$ et $b_0/R_0 \approx 1/(1+\varepsilon)$.
4. **ZÉRO FALLBACK ANALYTIQUE, ZÉRO DISSIMULATION D'INSTABILITÉ** : En cas de crash (`SIGSEGV`), `NaN`, perte de masse $> 1\,\%$, ou dispersion chaotique des marqueurs dans le domaine, **il est strictement interdit** de renvoyer $f_{\mathrm{th}}$ ou un cosinus synthétique, ou de prétendre dans le rapport que les simulations ont convergé. Tous les comportements réels (y compris divergences ou écarts entre codes) doivent être fidèlement exposés.

---

## 2. Théorie Exacte : Dispersion Capillaire de Lamb (Mode $n = 2$) et Contrainte de Stabilité

Dans une cavité carrée $\Omega = [0, 0.1] \times [0, 0.1]\,\mathrm{m}^2$ sans gravité ($\mathbf{g} = \mathbf{0}$) à parois `FREE_SLIP`, une bulle de gaz ($\rho_2 = 1.17683\,\mathrm{kg/m}^3$, $\mu_2 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$, $C = 0$) de rayon équivalent $R_0 = 0.02\,\mathrm{m}$ est immergée dans un liquide lourd ($\rho_1 = 7000.0\,\mathrm{kg/m}^3$, $\mu_1 = 4.0\times 10^{-3}\,\mathrm{Pa\cdot s}$, $C = 1$) avec $\sigma = 1.8\,\mathrm{N/m}$.

### 2.1 Fréquence théorique exacte de Lamb ($n=2$, milieu infini)
La pulsation propre théorique vaut exactement :
$$\omega_2 = \sqrt{\frac{(2^3 - 2)\,\sigma}{(\rho_1 + \rho_2)\,R_0^3}} = \sqrt{\frac{6 \times 1.8}{(7000.0 + 1.17683)\times (0.02)^3}} = \sqrt{\frac{10.8}{0.056009415}} = \mathbf{13.88613\,\mathrm{rad/s}},$$
soit une fréquence et une période théoriques exactes de :
$$\mathbf{f_{\mathrm{th}} = \frac{\omega_2}{2\pi} = 2.21005\,\mathrm{Hz}}, \qquad \mathbf{T_{\mathrm{th}} = 0.45248\,\mathrm{s}}.$$
*(Attention : l'horizon de simulation doit donc être $t_{\mathrm{final}} \ge 1.35\,\mathrm{s}$ pour observer au moins $3$ périodes complètes ; un calcul limité à $0.60\,\mathrm{s}$ ne couvre que $1.33$ période, et $0.35\,\mathrm{s}$ moins d'un cycle.)*

Pour le taux d'amortissement visqueux de Lamb en 2D ($\nu_1 = \mu_1/\rho_1 = 5.7143\times 10^{-7}\,\mathrm{m}^2/\mathrm{s}$), notez le modèle retenu lors de la comparaison (ex. formule asymptotique 2D $2n(n-1)\nu_1/R_0^2 = 5.714\times 10^{-3}\,\mathrm{s}^{-1}$ ou $2(n^2-1)\nu_1/R_0^2 = 8.571\times 10^{-3}\,\mathrm{s}^{-1}$).

### 2.2 Contrainte de pas de temps capillaire côté gaz ($\rho_2 = 1.17683$)
Dans une formulation à force capillaire explicite et rapport de densité $\rho_1/\rho_2 \approx 5948$, le critère local de stabilité capillaire dans le gaz léger vaut :
$$\Delta t_{\sigma,\mathrm{gaz}} = \sqrt{\frac{2 \rho_2 \, \Delta x^3}{4\pi \sigma}} = \begin{cases} 6.31\times 10^{-4}\,\mathrm{s} & (64^2), \\ 2.23\times 10^{-4}\,\mathrm{s} & (128^2), \\ 7.87\times 10^{-5}\,\mathrm{s} & (256^2). \end{cases}$$
L'utilisation d'un pas $\Delta t = 2.5\times 10^{-4}\,\mathrm{s}$ viole donc le critère côté gaz dès $128^2$ et $256^2$. Choisissez des $\Delta t \le 0.5\,\Delta t_{\sigma,\mathrm{gaz}}$ (ou testez comparativement $\Delta t = 2.5\times 10^{-4}\,\mathrm{s}$ vs $\Delta t$ résolu) et utilisez un critère linéaire strict (`tol: 1e-10`, `max_iter: 1000`). Par ailleurs, si le couplage $\rho(C_s)$ lissé + $\nabla C$ raide + $\kappa^{\mathrm{marker}}$ reste instable sur ce binaire à $\rho_1/\rho_2 = 5950$, documentez-le objectivement !

---

## 3. Matrice Complète des Paramètres Physiques et Numériques

| Configuration | Description | `use_markers` | `marker_advection_scheme` | `enable_svof_smoothing` ($k$) | `enable_diffusion_smoothing` (`steps`) |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **Config 1a** | VOF + CSF diffusif | `false` | — | `false` | `true` (`3`) |
| **Config 1b** | VOF + CSF SVOF Helmholtz | `false` | — | `true` (`0.5` et `1.0`) | `false` (`0`) |
| **Config 2a** | Hybride CVI sans SVOF | `true` | `"CVI_RK2"` | `false` | `false` (`0`) |
| **Config 2b** | Hybride CVI + SVOF Helmholtz | `true` | `"CVI_RK2"` | `true` (`0.5`) | `false` (`0`) |
| **Config 3** | Hybride Bilinéaire + SVOF | `true` | `"BILINEAR_RK2"` | `true` (`0.5`) | `false` (`0`) |

- **Grilles et pas de temps** : $N \in \{64, 128\}$ (et $256$ si stable dans le budget CPU) avec $\Delta t \in \{1.0\times 10^{-4},\, 5.0\times 10^{-5},\, 2.5\times 10^{-5}\}\,\mathrm{s}$ ($\le 0.5\,\Delta t_{\sigma,\mathrm{gaz}}$) jusqu'à $t_{\mathrm{final}} = 1.40\,\mathrm{s}$ ($\approx 3.1$ périodes de Lamb à $2.210\,\mathrm{Hz}$).
- **Amplitudes initiales** : $\varepsilon = 0.05$ (`axis_a = 0.021, axis_b = 0.01904761904761905`) et $\varepsilon = 0.20$ (`axis_a = 0.024, axis_b = 0.016666666666666666`).
- **Parité dual-code** : Exécuter à minima `Config 1b` et `Config 2b` sur $N \in \{64, 128\}$ à la fois sur `covof2d_cpp` et `covof2d_fortran`.

### Prototype `case/cpp/cas_04_cfg2b_cvi_svof_n128.json`
```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 128,
    "lx": 0.1,
    "ly": 0.1
  },
  "physics": {
    "enable_two_phase": true,
    "compressible": false,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_surface_tension": true,
    "enable_svof_smoothing": true,
    "svof_filter_width": 0.5,
    "enable_diffusion_smoothing": false,
    "diffusion_smoothing_steps": 0,
    "sigma": 1.8,
    "gx": 0.0,
    "gy": 0.0,
    "fluid1": { "rho": 7000.0,  "mu": 4.0e-3 },
    "fluid2": { "rho": 1.17683, "mu": 1.85e-5 }
  },
  "initial_condition": {
    "type": "BUBBLE",
    "center_x": 0.05,
    "center_y": 0.05,
    "radius": 0.02,
    "axis_a": 0.021,
    "axis_b": 0.01904761904761905
  },
  "time": {
    "dt": 0.00005,
    "t_final": 1.4,
    "time_scheme": "GEAR",
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
    "marker_normals_plic": true,
    "marker_max_per_cell": 8,
    "write_markers_vtk": true
  },
  "output": {
    "output_dir": "output/cpp/cfg2b_cvi_svof_n128",
    "case_name": "cas_04_lamb_cvi_n128",
    "vtk_interval": 200,
    "log_interval": 50,
    "write_csv_history": true
  }
}
```

---

## 4. Exigences de Post-Traitement, Livrables et Rédaction du Rapport

1. **Extraction géométrique homogène** : Calculer les demi-axes $a(t), b(t)$ depuis le champ $C$ des fichiers `.vtr` (isocontour $C=0.5$ ou moments d'inertie de `1 - C`), et **jamais** par moyenne spatiale brute des marqueurs ni par génération analytique de secours.
2. **Analyse spectrale (FFT + passages par zéro / extrema)** : Calculer la fréquence $f_{\mathrm{num}}$ uniquement si le signal est resté physiquement oscillatoire et borné ($|\Delta V|/V_0 < 1\,\%$, $a(t)/R_0 < 1.5$, aucun `NaN`). Dans le cas contraire, renseigner `INSTABLE / DIVERGE` sans jamais substituer $f_{\mathrm{th}} = 2.2100\,\mathrm{Hz}$.
3. **Livrables attendus** :
   - `figures/comparaison_demi_axes_et_frequence_lamb.{png,pdf}` (référence à $f_{\mathrm{th}} = 2.2100\,\mathrm{Hz}$, $T_{\mathrm{th}} = 0.4525\,\mathrm{s}$) ;
   - `figures/conservation_volume_bulle.{png,pdf}` (calculée depuis `.vtr`, avec seuils réels et sans forcer une ligne plane à $10^{-15}$) ;
   - `figures/comparaison_cvi_vs_bilineaire_poles.{png,pdf}` ;
   - `figures/dual_code_parity_cas04.{png,pdf}`.
4. **Rapport (`rapport/rapport_validation_04_oscillation_bulle_lamb_cvi_markers.{md,pdf}`)** :
   - Reproduire intégralement les paramètres physiques/numériques de chaque `.json`, la fréquence théorique $2.2100\,\mathrm{Hz}$, et exactement les données des CSV de synthèse générés par le script de post-traitement (y compris échecs/divergences s'il y en a).
