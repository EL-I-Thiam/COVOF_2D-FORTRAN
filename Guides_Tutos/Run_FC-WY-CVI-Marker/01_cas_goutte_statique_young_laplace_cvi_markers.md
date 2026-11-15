# INSTRUCTIONS D'EXÉCUTION POUR AGENT / MODÈLE IA — CAS 1 (FC-WY-CVI-MARKER) : GOUTTE STATIQUE DE YOUNG-LAPLACE ET RÉDUCTION DES COURANTS PARASITES (`01_cas_goutte_statique_young_laplace_cvi_markers`)

## 1. Mission et Rôle de l'Agent

Vous êtes un agent autonome chargé d'exécuter, de vérifier et de post-traiter de bout en bout le **Cas de Validation 1 de la méthode hybride FC-WY-CVI-Marker (`01_cas_goutte_statique_young_laplace_cvi_markers`)** de la plateforme CFD diphasique incompressible **COVOF 2D** dans ses deux implémentations **C++17** (`covof2d_cpp`) et **Fortran 2018** (`covof2d_fortran`).

Ce premier cas fondamental quantifie le gain apporté par la nouvelle architecture hybride **Solveur Monolithique (FC) + VOF-PLIC (Weymouth & Yue) + Interpolation Conservative de Vitesse (CVI) + Marqueurs Lagrangiens Non Connectés** sur la précision géométrique de la courbure $\kappa$ et la suppression des courants parasites capillaires (*spurious currents*).

**Règles impératives :**
1. Créez et isolez l'intégralité de la campagne dans `Validation_Totale/Run_CASES/Run_FC-WY-CVI-Marker/01_cas_goutte_statique_young_laplace_cvi_markers/` avec l'arborescence stricte `case/`, `output/`, `figures/`, `scripts/`, `rapport/`.
2. Exécutez les trois axes d'étude :
   - **Axe A (Convergence de la courbure et des courants parasites avec Marqueurs CVI)** : grilles $N \in \{32, 64, 128, 256\}$ avec `use_markers: true`, `marker_curvature: true`, `marker_advection_scheme: "CVI_RK2"` **sans lissage diffusif artificiel** (`enable_diffusion_smoothing: false`).
   - **Axe B (Comparaison directe VOF seul vs Hybride CVI-Marqueurs)** : comparaison à $N = 64$ et $N = 128$ entre (1) VOF brut sans lissage, (2) VOF avec lissage diffusif $\widetilde{C}$ (`diffusion_smoothing_steps: 3`), (3) Marqueurs avec interpolation bilinéaire (`BILINEAR_RK2`), et (4) Marqueurs avec interpolation conservative `CVI_RK2`.
   - **Axe C (Synergie avec le Solveur Monolithique FC vs Projection / Lagrangien Augmenté)** : comparaison de `MONOLITHIC` (FC), `ADAPTATIVE_AUGMENTED_LAGRANGIAN` (AAL) et `STANDARD_PROJECTION` couplés aux marqueurs CVI.

---

## 2. Compilation des Moteurs et Arborescence du Dossier de Simulation

```bash
mkdir -p COVOF_2D-C++/build && cd COVOF_2D-C++/build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j4 && cd ../..
mkdir -p COVOF_2D-FORTRAN/build && cd COVOF_2D-FORTRAN/build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j4 && cd ../..
```

Arborescence obligatoire à créer :
```text
Validation_Totale/Run_CASES/Run_FC-WY-CVI-Marker/01_cas_goutte_statique_young_laplace_cvi_markers/
├── case/{cpp,fortran}/
├── output/{cpp,fortran}/
├── figures/
├── scripts/{create_cases.py,run_all_simulations.py,postprocess_analysis.py,paraview_render.py}
└── rapport/
```

---

## 3. Théorie : Équations Physiques, Courbure Lagrangienne et Loi de Young-Laplace

Dans un domaine carré fermé $\Omega = [0, 1] \times [0, 1]\,\mathrm{m}^2$ sans gravité ($\mathbf{g} = \mathbf{0}$) avec parois à glissement libre (`FREE_SLIP`), une goutte circulaire de rayon $R_0 = 0.25\,\mathrm{m}$ est centrée en $(x_c, y_c) = (0.5, 0.5)\,\mathrm{m}$ au repos initial $\mathbf{u}(\mathbf{x}, 0) = \mathbf{0}$, avec $\rho_1 = \rho_2 = 1.0\,\mathrm{kg/m}^3$, $\mu_1 = \mu_2 = 1.0\,\mathrm{Pa\cdot s}$ et $\sigma = 1.0\,\mathrm{N/m}$.

### 3.1 Solution exacte de Young-Laplace
À l'équilibre exact, la vitesse est identiquement nulle ($\mathbf{u}_{\mathrm{ex}}(\mathbf{x}, t) \equiv \mathbf{0}$), la courbure sur le cercle vaut $\kappa_{\mathrm{ex}} = 1/R_0 = 4.0000\,\mathrm{m}^{-1}$, et le saut de pression capillaire vaut :
$$\Delta P_{\mathrm{ex}} = P_{\mathrm{int}} - P_{\mathrm{ext}} = \frac{\sigma}{R_0} = 4.0000\,\mathrm{Pa}.$$

### 3.2 Verrou du CSF eulérien standard et levée du verrou par les Marqueurs CVI
Dans l'approche CSF eulérienne standard, la courbure $\kappa = -\nabla \cdot (\nabla \widetilde{C}/\|\nabla \widetilde{C}\|)$ est évaluée par différences finies sur la fraction volumique, ce qui engendre des erreurs d'ordre $O(1)$ sur $\kappa$ et génère des tourbillons parasites de nombre capillaire :
$$Ca_{\max}(t) = \frac{\mu \, \|\mathbf{u}(\cdot, t)\|_\infty}{\sigma}.$$
Dans la nouvelle méthode hybride :
1. Les marqueurs non connectés $(x_p, y_p)$ distribués sur l'interface sont advectés par le schéma RK2-Heun utilisant l'interpolateur **CVI** ($\nabla \cdot \mathbf{u}_{\mathrm{CVI}} \equiv 0$).
2. Dans chaque cellule interfaciale $(i,j)$, les marqueurs du voisinage $3\times 3$ filtrés par cohérence angulaire ($\mathbf{n}_p \cdot \mathbf{n}_{i,j}^{\mathrm{ref}} \ge 0.5$) sont projetés dans le repère local tangent-normal $(\tau, \eta)$ où une parabole pondérée $\eta(\tau) = a\tau^2 + b\tau + c$ fournit la courbure d'ordre 2 :
   $$\kappa_{i,j}^{\mathrm{marker}} = \frac{2a}{(1 + b^2)^{3/2}}.$$
3. La force capillaire *Balanced-Force* aux faces MAC s'écrit alors directement sans diffusion artificielle de $C$ :
   $$f_{\sigma, x, i,j} = \sigma \, \frac{\kappa_{i-1,j}^{\mathrm{marker}} + \kappa_{i,j}^{\mathrm{marker}}}{2} \, \frac{C_{i,j} - C_{i-1,j}}{\Delta x}.$$

---

## 4. Paramètres Physiques, Numériques et Matrice des Jeux de Paramètres

| Catégorie | Paramètre | Valeur / Spécification |
| :--- | :--- | :--- |
| Domaine et Goutte | $\Omega$, $(x_c, y_c)$, $R_0$ | $[0,1]\times[0,1]\,\mathrm{m}^2$, $(0.5, 0.5)\,\mathrm{m}$, $R_0 = 0.25\,\mathrm{m}$ |
| Propriétés physiques | $\rho_1=\rho_2, \mu_1=\mu_2, \sigma, \mathbf{g}$ | $\rho = 1.0\,\mathrm{kg/m}^3$, $\mu = 1.0\,\mathrm{Pa\cdot s}$, $\sigma = 1.0\,\mathrm{N/m}$, $\mathbf{g} = \mathbf{0}$ |
| Séquence de grilles et $\Delta t$ | $N \in \{32, 64, 128, 256\}$, $t_{\mathrm{final}}=0.05\,\mathrm{s}$ | $\Delta t \in \{0.005, 0.005, 0.0025, 0.00125\}\,\mathrm{s}$ |
| Paramètres Marqueurs & CVI | `use_markers`, `markers_per_cell`, `marker_curvature` | `true`, $N_{\mathrm{target}} = 4$, `marker_advection_scheme = "CVI_RK2"` |
| Métriques à extraire | Courbure, Pression, Courants parasites, Volume | $E_\infty(\kappa) = \|\kappa^{\mathrm{marker}} - 4.0\|_\infty / 4.0$, $|\Delta P - 4.0|/4.0$, $\|\mathbf{u}\|_\infty$, $Ca_{\max}$, $|\Delta V|/V_0$ |

---

## 5. Prototype de Fichier `.json` (C++17 et Fortran 2018)

### Prototype `case/cpp/cas_01_axeA_fc_wy_cvi_markers_n128.json`
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
    "radius": 0.25
  },
  "time": {
    "dt": 0.0025,
    "t_final": 0.05,
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
    "output_dir": "Validation_Totale/Run_CASES/Run_FC-WY-CVI-Marker/01_cas_goutte_statique_young_laplace_cvi_markers/output/cpp/axeA_n128",
    "case_name": "cas_01_fc_wy_cvi_n128",
    "vtk_interval": 5,
    "log_interval": 1,
    "write_csv_history": true
  }
}
```

---

## 6. Post-Traitement et Livrables Attendus

1. **Graphique de convergence log-log (`figures/convergence_courbure_et_courants_parasites.pdf`)** : Tracer $\|\mathbf{u}\|_\infty(t_{\mathrm{final}})$ et l'erreur relative sur $\Delta P$ en fonction de $N \in \{32, 64, 128, 256\}$ en comparant le CSF standard et la courbure par marqueurs CVI.
2. **Visualisation ParaView (`figures/paraview_champ_pression_et_marqueurs_n128.png`)** : Superposer le champ de pression $p(x,y)$ sur la grille MAC et le nuage de marqueurs `*_markers_*.vtk` coloré par la courbure locale $\kappa_p$.
3. **Rapport PDF complet (`rapport/rapport_validation_01_goutte_statique_young_laplace_cvi_markers.pdf`)** compilé sans débordement de marge.
