# INSTRUCTIONS D'EXÉCUTION POUR AGENT / MODÈLE IA — CAS 4 (FC-WY-CVI-MARKER) : OSCILLATIONS CAPILLAIRES AMORTIES D'UNE BULLE DE LAMB À GRAND RAPPORT DE DENSITÉ (`04_cas_oscillation_bulle_lamb_cvi_markers`)

## 1. Mission et Rôle de l'Agent

Vous êtes un agent autonome chargé d'exécuter, de vérifier et de post-traiter de bout en bout le **Cas de Validation 4 de la méthode hybride FC-WY-CVI-Marker (`04_cas_oscillation_bulle_lamb_cvi_markers`)** sur les moteurs **C++17** (`covof2d_cpp`) et **Fortran 2018** (`covof2d_fortran`) de la plateforme **COVOF 2D**.

Ce cas évalue la capacité du couplage **Solveur Monolithique (FC) + VOF-PLIC (Weymouth & Yue) + Interpolation CVI + Marqueurs Non Connectés** à restituer sans amortissement numérique parasite la fréquence propre et le taux d'amortissement visqueux des oscillations capillaires d'une bulle de gaz dans un liquide très dense ($\rho_1/\rho_2 \approx 5950$).

**Règles impératives :**
1. Isolez toute la campagne dans `Validation_Totale/Run_CASES/Run_FC-WY-CVI-Marker/04_cas_oscillation_bulle_lamb_cvi_markers/` avec l'arborescence `case/`, `output/`, `figures/`, `scripts/`, `rapport/`.
2. Comparez systématiquement sur les grilles $N \in \{64, 128, 256\}$ :
   - **Configuration 1a (Référence VOF + CSF à lissage diffusif Laplacien $\widetilde{C}$)** : `use_markers: false`, `enable_diffusion_smoothing: true` (`diffusion_smoothing_steps: 3`), `enable_svof_smoothing: false`.
   - **Configuration 1b (VOF + CSF à lissage elliptique SVOF de Helmholtz $C_s$)** : `use_markers: false`, `enable_svof_smoothing: true` (`svof_filter_width` $k \in \{0.5, 1.0\}$), `enable_diffusion_smoothing: false` (où $C_s - k\Delta x\Delta y\nabla^2 C_s = C$ sert à la fois au lissage rhéologique $\rho(C_s), \mu(C_s)$ et au calcul de $\kappa = -\nabla \cdot (\nabla C_s / \|\nabla C_s\|)$).
   - **Configuration 2a (Hybride FC + WY + CVI + Marqueurs sans lissage SVOF)** : `use_markers: true`, `marker_curvature: true`, `marker_advection_scheme: "CVI_RK2"`, `enable_svof_smoothing: false`, `enable_diffusion_smoothing: false`.
   - **Configuration 2b (Synergie optimale : Hybride FC + WY + CVI + Marqueurs + Lissage rhéologique SVOF de Helmholtz)** : `use_markers: true`, `marker_curvature: true`, `marker_advection_scheme: "CVI_RK2"`, `enable_svof_smoothing: true` (`svof_filter_width: 0.5`), `enable_diffusion_smoothing: false` (régularisation elliptique de $\rho(C_s), \mu(C_s)$ à $\rho_1/\rho_2 = 5950$ tout en conservant la courbure parabolique sous-maille $\kappa^{\mathrm{marker}}$ et le saut capillaire raide $\nabla C$).
   - **Configuration 3 (Hybride avec interpolation bilinéaire non conservative)** : `use_markers: true`, `marker_advection_scheme: "BILINEAR_RK2"`, `enable_svof_smoothing: true` (`svof_filter_width: 0.5`), afin de démontrer le gain de la condition de Liouville ($\nabla \cdot \mathbf{u}_{\mathrm{CVI}} = 0$) sur la distribution des marqueurs aux pôles de l'ellipse après 3 périodes d'oscillation.

---

## 2. Théorie : Dispersion Capillaire de Lamb et Couplage CVI-Marqueurs

Dans une cavité carrée $\Omega = [0, 0.1] \times [0, 0.1]\,\mathrm{m}^2$ sans gravité ($\mathbf{g} = \mathbf{0}$) avec parois `FREE_SLIP`, une bulle de gaz ($\rho_2 = 1.17683\,\mathrm{kg/m}^3$, $\mu_2 = 1.85\times 10^{-5}\,\mathrm{Pa\cdot s}$) est entourée d'un liquide lourd ($\rho_1 = 7000.0\,\mathrm{kg/m}^3$, $\mu_1 = 4.0\times 10^{-3}\,\mathrm{Pa\cdot s}$) avec $\sigma = 1.8\,\mathrm{N/m}$. À $t=0$, la bulle de rayon équivalent $R_0 = 0.02\,\mathrm{m}$ est déformée sur le mode quadrupolaire $n=2$ ($a_0 = R_0(1+\varepsilon)$, $b_0 = R_0/(1+\varepsilon)$, avec $\varepsilon = 0.05$ en régime linéaire et $\varepsilon = 0.20$ en régime non linéaire).

### 2.1 Pulsation propre et amortissement analytique de Lamb (mode $n=2$)
En milieu infini, la pulsation propre théorique $\omega_2$ et la fréquence $f_{\mathrm{th}}$ valent :
$$\omega_2 = \sqrt{\frac{(2^3 - 2)\,\sigma}{(\rho_1 + \rho_2)\,R_0^3}} = \sqrt{\frac{6\,\sigma}{(\rho_1 + \rho_2)\,R_0^3}} = 35.8539\,\mathrm{rad/s}, \qquad f_{\mathrm{th}} = \frac{\omega_2}{2\pi} = 5.70632\,\mathrm{Hz} \quad (T_{\mathrm{th}} = 0.17524\,\mathrm{s}).$$

### 2.2 Pourquoi la méthode FC + WY + CVI + Marqueurs lève le verrou du Cas 13
Dans le cas VOF eulérien standard, le lissage laplacien $\widetilde{C}$ requis pour stabiliser les dérivées secondes de $\kappa$ étale le saut de pression capillaire sur 3 à 4 mailles, ce qui réduit le rappel capillaire effectif aux sommets de l'ellipse (abaissant la fréquence numérique autour de $4.95\,\mathrm{Hz}$ sur grille $128^2$) et dissipe prématurément l'énergie cinétique.
Grâce à l'estimation parabolique locale de $\kappa$ sur le nuage de marqueurs advecté par **CVI**, la courbure aux pôles ($\kappa_{\max} = a_0 / b_0^2$) et à l'équateur ($\kappa_{\min} = b_0 / a_0^2$) est capturée sans étaler la fraction volumique $C$, tandis que le solveur **Monolithique FC** avec lissage rhéologique SVOF ($k = 0.5$) absorbe le rapport de densité $\rho_1/\rho_2 = 5950$ avec une divergence discrète $\|B_u \mathbf{u}^{n+1}\|_\infty \le 10^{-10}$.

---

## 3. Paramètres Physiques, Numériques et Prototype `.json`

| Catégorie | Paramètre | Valeur / Spécification |
| :--- | :--- | :--- |
| Domaine et Bulle initiale | $[0, 0.1]^2\,\mathrm{m}^2$, $(x_c, y_c)$, $R_0$, $\varepsilon$ | $(0.05, 0.05)\,\mathrm{m}$, $R_0 = 0.02\,\mathrm{m}$, $\varepsilon \in \{0.05, 0.20\}$ |
| Contrastes rhéologiques | $\rho_1/\rho_2$, $\mu_1/\mu_2$, $\sigma$ | $7000 / 1.17683 \approx 5948$, $4\times 10^{-3} / 1.85\times 10^{-5} \approx 216.2$, $\sigma = 1.8\,\mathrm{N/m}$ |
| Grilles et Temps | $N \in \{64, 128, 256\}$, $\Delta t = 2.5\times 10^{-4}\,\mathrm{s}$ | $t_{\mathrm{final}} = 0.60\,\mathrm{s}$ ($2400$ pas de temps, $>3.4$ périodes) |
| Solveur Monolithique | `MONOLITHIC` + `BICGSTAB2` + `BLOCK_TRIANGULAR` | `enable_svof_smoothing: true`, `svof_filter_width: 0.5`, `tol: 1e-10` |

### Prototype `case/cpp/cas_04_fc_wy_cvi_markers_n128.json`
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
    "fluid1": { "rho": 7000.0, "mu": 4.0e-3 },
    "fluid2": { "rho": 1.17683, "mu": 1.85e-5 }
  },
  "initial_condition": {
    "type": "BUBBLE",
    "center_x": 0.05,
    "center_y": 0.05,
    "radius": 0.02
  },
  "time": {
    "dt": 0.00025,
    "t_final": 0.6,
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
    "output_dir": "Validation_Totale/Run_CASES/Run_FC-WY-CVI-Marker/04_cas_oscillation_bulle_lamb_cvi_markers/output/cpp/cvi_n128",
    "case_name": "cas_04_lamb_cvi_n128",
    "vtk_interval": 20,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 4. Livrables et Vérifications Attendues

1. **Analyse spectrale et temporelle (`figures/comparaison_demi_axes_et_frequence_lamb.pdf`)** : Tracer l'évolution temporelle du demi-grand axe $a(t)/R_0$ et de l'énergie cinétique $E_k(t)$ avec et sans marqueurs CVI, et tabuler l'erreur relative $|f_{\mathrm{num}} - f_{\mathrm{th}}| / f_{\mathrm{th}}$.
2. **Conservation exacte du volume (`figures/conservation_volume_bulle.pdf`)** : Vérifier que $|\Delta V(t)| / V_0 \le 10^{-12}$ sur les 2400 pas de temps.
3. **Rapport PDF complet (`rapport/rapport_validation_04_oscillation_bulle_lamb_cvi_markers.pdf`)**.
