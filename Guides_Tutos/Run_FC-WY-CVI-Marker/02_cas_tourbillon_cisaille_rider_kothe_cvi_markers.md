# INSTRUCTIONS D'EXÉCUTION POUR AGENT / MODÈLE IA — CAS 2 (FC-WY-CVI-MARKER) : TOURBILLON CISAILLÉ RÉVERSIBLE DE RIDER-KOTHE ET PRÉSERVATION SOUS-MAILLE PAR MARQUEURS CVI (`02_cas_tourbillon_cisaille_rider_kothe_cvi_markers`)

## 1. Mission et Rôle de l'Agent

Vous êtes un agent autonome chargé d'exécuter, de vérifier et de post-traiter de bout en bout le **Cas de Validation 2 de la méthode hybride FC-WY-CVI-Marker (`02_cas_tourbillon_cisaille_rider_kothe_cvi_markers`)** sur les moteurs **C++17** (`covof2d_cpp`) et **Fortran 2018** (`covof2d_fortran`) de la plateforme **COVOF 2D**.

Ce cas cinématique de référence teste la capacité du couplage **VOF-PLIC (Weymouth & Yue) + Interpolation Conservative de Vitesse (CVI) + Marqueurs Non Connectés** à :
1. Satisfaire le **théorème de Liouville discret** ($\nabla \cdot \mathbf{u}_{\mathrm{CVI}} \equiv 0$) sur une longue intégration temporelle ($T = 6.0\,\mathrm{s}$) ;
2. Conserver la trace continue sous-maille de la queue de la spirale à l'étirement maximal ($t = T/2 = 3.0\,\mathrm{s}$) là où le VOF-PLIC eulérien seul subit un pincement numérique (*flotsam*) sur maillages $64^2$ et $128^2$ ;
3. Réduire l'erreur géométrique $L_1(T)$ et la distance de Hausdorff $d_H(T)$ des marqueurs après l'inversion complète du tourbillon à $t = T = 6.0\,\mathrm{s}$.

---

## 2. Théorie : Champ Solénoïdal de Rider-Kothe et Réversibilité Lagrangienne CVI

Dans la cavité unitaire $\Omega = [0, 1] \times [0, 1]\,\mathrm{m}^2$, un disque de fluide 1 ($C = 1$) de rayon $R_0 = 0.125\,\mathrm{m}$ est centré en $(x_c, y_c) = (0.5, 0.75)\,\mathrm{m}$ à $t = 0$. Le champ de vitesse analytique solénoïdal dérive de la fonction de courant :
$$\Psi(x, y, t) = \frac{1}{\pi} \sin^2(\pi x) \sin^2(\pi y) \cos\left(\frac{\pi t}{T}\right), \qquad T = 6.0\,\mathrm{s},$$
soit aux faces décalées de la grille MAC :
$$u(x, y, t) = -\sin^2(\pi x)\sin(2\pi y)\cos\left(\frac{\pi t}{T}\right), \qquad v(x, y, t) = +\sin(2\pi x)\sin^2(\pi y)\cos\left(\frac{\pi t}{T}\right).$$

### 2.1 Comparaison théorique : Interpolation Bilinéaire vs Interpolation CVI
- Avec une **interpolation bilinéaire standard** (`BILINEAR_RK2`), la divergence interne $\nabla \cdot \mathbf{u}_{\mathrm{bilin}}(\xi, \eta) \neq 0$ dévie les marqueurs transversalement aux lignes de courant de $\Psi$, provoquant un regroupement artificiel (*clustering*) et une erreur irréversible à $t = T$.
- Avec l'**interpolation conservative CVI** (`CVI_RK2`) couplée au schéma de Heun d'ordre 2, $\nabla \cdot \mathbf{u}_{\mathrm{CVI}}(\xi, \eta) = \frac{\Delta u}{\Delta x} + \frac{\Delta v}{\Delta y} \equiv 0$. Les marqueurs suivent fidèlement les iso-valeurs de $\Psi$ à l'aller ($0 \to T/2$) comme au retour ($T/2 \to T$).
- De plus, l'utilisation des normales des marqueurs (`marker_normals_plic: true`) avec **filtrage angulaire** ($\mathbf{n}_p \cdot \mathbf{n}_{i,j}^{\mathrm{ref}} \ge 0.5$) oriente les segments PLIC dans la queue fine de la spirale sans confondre les deux bords opposés du ligament.

---

## 3. Matrice des Jeux de Paramètres et Axes d'Étude

| Axe d'étude | Grilles $N\times N$ | Pas de temps $\Delta t$ | Configurations comparées |
| :--- | :---: | :---: | :--- |
| **Axe A : Convergence spatiale** | $64^2, 128^2, 256^2, 512^2$ | $5\cdot 10^{-3}, 2.5\cdot 10^{-3}, 1.25\cdot 10^{-3}, 5\cdot 10^{-4}\,\mathrm{s}$ | `WEYMOUTH_YUE + PLIC` seul vs `WEYMOUTH_YUE + PLIC + CVI_RK2 (marker_normals_plic: true)` |
| **Axe B : Qualité Lagrangienne CVI vs Bilinéaire** | $128^2, 256^2$ | $2.5\cdot 10^{-3}, 1.25\cdot 10^{-3}\,\mathrm{s}$ | Erreur radiale de retour $E_R(T) = \frac{1}{M}\sum_p \left|\|\mathbf{x}_p(T) - \mathbf{x}_c\| - R_0\right|$ entre `CVI_RK2` et `BILINEAR_RK2` |
| **Axe C : Invariance MPI ($NP=1$ vs $NP=4$)** | $128^2$ | $2.5\cdot 10^{-3}\,\mathrm{s}$ | Vérification du routage MPI en 2 passes (Passe X puis Passe Y) lorsque la spirale franchit les 4 sous-domaines |

### Prototype `case/cpp/cas_02_wy_plic_cvi_markers_n128.json`
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
    "fluid2": { "rho": 1.0,    "mu": 1.0e-05 }
  },
  "initial_condition": {
    "type": "DROPLET",
    "center_x": 0.5,
    "center_y": 0.75,
    "radius": 0.125
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
    "output_dir": "Validation_Totale/Run_CASES/Run_FC-WY-CVI-Marker/02_cas_tourbillon_cisaille_rider_kothe_cvi_markers/output/cpp/cvi_n128",
    "case_name": "cas_02_vortex_cvi_n128",
    "vtk_interval": 120,
    "log_interval": 120
  }
}
```

---

## 4. Livrables et Post-Traitement Attendus

1. **Superposition de l'interface VOF ($C = 0.5$) et du nuage de marqueurs CVI à $t = 3.0\,\mathrm{s}$ et $t = 6.0\,\mathrm{s}$** (`figures/spirale_t3s_et_retour_t6s_vof_vs_cvi.pdf`) illustrant la continuité sous-maille de la queue spiralée par les marqueurs CVI.
2. **Tableau et courbes de convergence log-log** de l'erreur géométrique $E_{L_1}(T) = \sum_{i,j} |C_{i,j}(T) - C_{i,j}(0)|\Delta x \Delta y$, de l'erreur lagrangienne $E_R(T)$ (`CVI_RK2` vs `BILINEAR_RK2`) et de la conservation de masse $|\Delta V(T)|/V_0$.
3. **Rapport PDF complet (`rapport/rapport_validation_02_tourbillon_cisaille_rider_kothe_cvi_markers.pdf`)**.
