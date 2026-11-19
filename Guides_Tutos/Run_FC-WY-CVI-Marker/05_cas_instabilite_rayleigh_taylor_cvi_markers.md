# INSTRUCTIONS D'EXÉCUTION POUR AGENT / MODÈLE IA — CAS 5 (FC-WY-CVI-MARKER) : INSTABILITÉ GRAVITAIRE DE RAYLEIGH-TAYLOR DIPHASIQUE (`05_cas_instabilite_rayleigh_taylor_cvi_markers`)

## 1. Mission, Contraintes Critiques et Interdictions Absolues

Vous êtes un agent autonome chargé de préparer, d'exécuter, de vérifier et de post-traiter de bout en bout le **Cas de Validation 5 (`05_cas_instabilite_rayleigh_taylor_cvi_markers`)** sur les moteurs précompilés **C++17** (`COVOF_2D-C++/bin/covof2d_cpp`) et **Fortran 2018** (`COVOF_2D-FORTRAN/bin/covof2d_fortran`).

### Interdictions absolues et règles de fidélité
1. **INTERDICTION DE MODIFIER LE CODE SOURCE OU LES BINAIRES** : Ne modifiez aucun fichier dans `COVOF_2D-C++/` ni `COVOF_2D-FORTRAN/` et ne recompilez pas les binaires.
2. **CONFORMITÉ DES CHEMINS (`output_dir <= 512` car.)** : Exécutez avec `cwd = VALIDATION_TOTALE/Run_CASES/RUN_FC-WY-CVI-Marker/05_cas_instabilite_rayleigh_taylor_cvi_markers`, pré-créez `output/{cpp,fortran}/<cas>`, et utilisez des chemins courts relatifs `"output/cpp/<cas>"` et `"output/fortran/<cas>"`.
3. **VÉRIFICATION OBLIGATOIRE DES MARQUEURS SUR L'INITIALISATION `RAYLEIGH_TAYLOR`** :
   Dans les binaires C++ et Fortran actualisés, lorsque `vof.use_markers` est activé, la condition initiale `"RAYLEIGH_TAYLOR"` utilise un sous-échantillonnage $16\times 16$ par cellule afin de générer des cellules mixtes (`0.01 < C < 0.99`) et d'ensemencer des marqueurs dès $t = 0$.
   **Obligation absolue :** Le script de post-traitement doit compter et rapporter le nombre réel `n_markers(t)` dans `*_markers_*.vtk` pour chaque simulation (`n_markers(0) > 0` attendu lorsque `use_markers: true`). Si `n_markers(t) == 0` sur une exécution donnée, signalez-le explicitement sans attribuer aux marqueurs CVI un résultat issu de `WY + PLIC (+ SVOF)` seul.
4. **RÔLE EXACT DU LISSAGE DIFFUSIF VS SVOF** : À $\sigma = 0$ (`enable_surface_tension: false`), `enable_diffusion_smoothing` (`C_tilde`) n'a aucun effet dans les équations ; seul `enable_svof_smoothing` ($C_s$) modifie $\rho(C_s)$ et $\mu(C_s)$. `enable_diffusion_smoothing` n'est pertinent qu'à $\sigma > 0$ lorsque `use_markers: false` (Axe B/C).
5. **ZÉRO RÉSULTAT EN DUR, ZÉRO FALLBACK** : Toutes les valeurs (taux de croissance $\gamma_{\mathrm{num}}$, itérations Krylov depuis `.log`, `max_div`, conservation `.vtr`) doivent provenir exclusivement des sorties brutes.

---

## 2. Arborescence et Théorie de l'Instabilité de Rayleigh-Taylor

```text
VALIDATION_TOTALE/Run_CASES/RUN_FC-WY-CVI-Marker/05_cas_instabilite_rayleigh_taylor_cvi_markers/
├── case/{cpp,fortran}/
├── output/{cpp,fortran}/
├── figures/
├── scripts/{create_cases.py,run_all_simulations.py,postprocess_analysis.py}
└── rapport/
```

Dans $\Omega = [0, L_x] \times [0, L_y]$ ($L_x = 0.5\,\mathrm{m}, L_y = 1.5\,\mathrm{m}$), un fluide lourd ($\rho_1 = 1.225\,\mathrm{kg/m}^3$, en haut, $C=1$) surmonte un fluide léger ($\rho_2 = 0.1694\,\mathrm{kg/m}^3$, en bas, $C=0$) sous $g_y = -9.81\,\mathrm{m/s}^2$, soit $At = \frac{\rho_1 - \rho_2}{\rho_1 + \rho_2} = 0.75705$, avec $\mu_1 = \mu_2 = 3.13\times 10^{-3}\,\mathrm{Pa\cdot s}$.
Pour une interface perturbée par $y_{\mathrm{int}}(x) = y_c - \eta_0 \cos\left(\frac{2\pi\, \texttt{pert\_k}\, x}{L_x}\right)$ ($\eta_0 = 0.01\,\mathrm{m}$, `pert_k = 1.0`, $k = 2\pi/L_x = 4\pi\,\mathrm{m}^{-1}$), le taux de croissance linéaire théorique non visqueux vaut :
$$\gamma_{\mathrm{th}} = \sqrt{At \, g \, k - \frac{\sigma \, k^3}{\rho_1 + \rho_2}} \approx \begin{cases} 9.660\,\mathrm{s}^{-1} & (\sigma = 0), \\ 9.289\,\mathrm{s}^{-1} & (\sigma = 0.005\,\mathrm{N/m}), \\ 8.075\,\mathrm{s}^{-1} & (\sigma = 0.02\,\mathrm{N/m}), \end{cases}$$
*(en incluant la viscosité de Bellman-Pennington $\gamma_{\mathrm{visc}}$, on pourra aussi comparer la racine positive de $\gamma^2 + 2\nu k^2 \gamma - \left(At g k - \frac{\sigma k^3}{\rho_1+\rho_2}\right) = 0$).*

---

## 3. Matrice Complète des Paramètres Physiques et Numériques

Toutes les simulations utilisent `MONOLITHIC` (`BICGSTAB2`, `BLOCK_TRIANGULAR`, `tol: 1e-10`, `max_iter: 1000`), `GEAR`, `ADAMS_BASHFORTH`, `CENTERED`, `WEYMOUTH_YUE`, `PLIC`, `tvd_limiter: "NONE"`, `t_final = 0.9 s`, $\Delta t = 0.001\,\mathrm{s}$ (vérifier $\mathrm{CFL}$ et critère capillaire pour $\sigma > 0$), `FREE_SLIP` à gauche/droite, `DIRICHLET` ($u=v=0$) en bas/haut :

| Axe d'étude | Grilles $N_x \times N_y$ | Tension superficielle $\sigma$ | Lissage / Schéma |
| :--- | :---: | :---: | :--- |
| **Axe A : Croissance linéaire ($\sigma = 0$)** | $64\times 192$, $128\times 384$ | $\sigma = 0.0\,\mathrm{N/m}$ | `WY + PLIC` (`use_markers: false`) vs `use_markers: true` (`CVI_RK2`, vérifier `n_markers` réel !) avec `SVOF k=0.5` et sans SVOF |
| **Axe B : Effet capillaire ($\sigma > 0$)** | $64\times 192$, $128\times 384$ | $\sigma \in \{0.005, 0.02\}\,\mathrm{N/m}$ | `CSF Diffusif` (`3 pas`) vs `CSF SVOF` ($k \in \{0.5, 1.0\}$) vs `use_markers: true` (vérifier `n_markers` réel !) |
| **Axe C : Parité dual-code** | $64\times 192$ | $\sigma \in \{0.0, 0.005\}\,\mathrm{N/m}$ | `covof2d_cpp` vs `covof2d_fortran` |

### Prototype `case/cpp/cas_05_axeB_rt_s0005_128x384.json`
```json
{
  "mode": "NAVIER_STOKES",
  "mesh": {
    "nx": 128,
    "ny": 384,
    "lx": 0.5,
    "ly": 1.5
  },
  "physics": {
    "enable_two_phase": true,
    "enable_surface_tension": true,
    "enable_convection": true,
    "convective_scheme": "CENTERED",
    "convective_linearization": "ADAMS_BASHFORTH",
    "enable_svof_smoothing": true,
    "svof_filter_width": 0.5,
    "enable_diffusion_smoothing": false,
    "diffusion_smoothing_steps": 0,
    "sigma": 0.005,
    "gx": 0.0,
    "gy": -9.81,
    "fluid1": { "rho": 1.225,  "mu": 3.13e-3 },
    "fluid2": { "rho": 0.1694, "mu": 3.13e-3 }
  },
  "initial_condition": {
    "type": "RAYLEIGH_TAYLOR",
    "center_y": 0.75,
    "pert_amp": 0.01,
    "pert_k": 1.0
  },
  "time": {
    "dt": 0.001,
    "t_final": 0.9,
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
    "output_dir": "output/cpp/axeB_rt_s0005_128",
    "case_name": "cas_05_rt_s0005_128",
    "vtk_interval": 25,
    "log_interval": 10,
    "write_csv_history": true
  }
}
```

---

## 4. Livrables, Post-Traitement et Rédaction du Rapport

1. **Grandeurs extraites depuis `.vtr` et `.log`** :
   - Position de la pointe du panache descendant ($y_{\mathrm{spike}}(t)$ sur colonne centrale) et de la bulle montante ($y_{\mathrm{bubble}}(t)$ sur bords latéraux) le long des iso-contours $C = 0.5$ ; ajustement du taux de croissance linéaire $\gamma_{\mathrm{num}}$ sur $\ln(\eta(t)/\eta_0)$ dans la phase linéaire ($t \in [0.05, 0.25]\,\mathrm{s}$).
   - Relevé exact de `n_markers(0)` et `n_markers(t_final)` depuis `*_markers_*.vtk`, conservation volumique depuis `.vtr`, itérations Krylov et `max_div` depuis `.log` / `.csv`.
2. **Livrables** : `figures/croissance_lineaire_rt.{png,pdf}`, `figures/volutes_champignon_t06_t08_t09.{png,pdf}`, `figures/dual_code_parity_cas05.{png,pdf}`, `rapport/rapport_validation_05_instabilite_rayleigh_taylor_cvi_markers.{md,pdf}` contenant la table complète des paramètres JSON et les résultats réels sans falsification.
