# Correspondance Théorie $\leftrightarrow$ Discrétisation $\leftrightarrow$ Code Source

Ce chapitre constitue la **section pivot** du Guide de Développement. Il établit le lien explicite entre les équations continues de la mécanique des fluides, leur discrétisation spatiotemporelle sur grille MAC et leur implémentation exacte dans `COVOF_2D-C++` (C++17) et `COVOF_2D-FORTRAN` (Fortran 2018).

Chaque opérateur est présenté selon le triptyque rigoureux :
1. **Équation discrète exacte** (indicée selon les conventions du code).
2. **Localisation exacte** (fichiers, classes/modules, méthodes/subroutines et plages de lignes).
3. **Analyse de l'implémentation** (boucles critiques, régularisations et invariants).

---

## 4.1. Diffusion Visqueuse et Moyenne Harmonique aux Nœuds ($\mu_{\text{node}}$)

### 4.1.1. Équation Discrète Exacte

Sur la grille MAC, la composante $u_{i,j}$ (face Ouest de la cellule $(i,j)$) intègre deux types de flux visqueux sur le volume de contrôle $\Delta x \Delta y$ :
- Les contraintes normales d'élongation $2\mu \partial_x u$ évaluées aux centres $(i,j)$ et $(i-1,j)$ :
  $$\text{Diff}_{\text{elong}} = \left[ 2 \mu_{i,j} \frac{u_{i+1,j} - u_{i,j}}{\Delta x} - 2 \mu_{i-1,j} \frac{u_{i,j} - u_{i-1,j}}{\Delta x} \right] \Delta y$$
- Les contraintes tangentielles de cisaillement $\mu_{\text{node}}(\partial_y u + \partial_x v)$ évaluées aux nœuds Nord $(i, j+1)$ et Sud $(i, j)$ :
  $$\text{Shear} = \left[ \mu_{\text{node}, i, j+1} \left( \frac{u_{i,j+1} - u_{i,j}}{\Delta y} + \frac{v_{i,j+1} - v_{i-1,j+1}}{\Delta x} \right) - \mu_{\text{node}, i, j} \left( \frac{u_{i,j} - u_{i,j-1}}{\Delta y} + \frac{v_{i,j} - v_{i-1,j}}{\Delta x} \right) \right] \Delta x$$

Pour assurer la continuité des contraintes tangentielles à travers une interface diphasique à fort saut $\mu_1/\mu_2$, $\mu_{\text{node}}$ est calculée par **moyenne harmonique** (ou interpolation SVOF harmonique de $C_{s,\text{node}}$ si `enable_svof_smoothing = true`) :

$$\mu_{\text{node}, i, j} = \frac{4}{\frac{1}{\mu_{i-1, j-1}} + \frac{1}{\mu_{i, j-1}} + \frac{1}{\mu_{i-1, j}} + \frac{1}{\mu_{i, j}}} \quad \text{ou} \quad \frac{\mu_1 \mu_2}{\mu_2 C_{s,\text{node}} + \mu_1 (1 - C_{s,\text{node}})}$$

### 4.1.2. Localisation dans le Code Source

- **Calcul de $\mu_{\text{node}}$ :**
  - C++17 : `COVOF_2D-C++/src/InterfaceTracker.cpp` $\rightarrow$ `InterfaceTracker::update_physical_properties` (lignes 558–582).
  - Fortran 2018 : `COVOF_2D-FORTRAN/src/module_interface_vof.f90` $\rightarrow$ `update_physical_properties` (lignes 479–626).
- **Application dans l'opérateur linéaire :**
  - C++17 : `COVOF_2D-C++/src/NavierStokesSolver.cpp` $\rightarrow$ `NavierStokesSolver::apply_system_matrix` (lignes 334–358 pour $u$, lignes 480–510 pour $v$).
  - Fortran 2018 : `COVOF_2D-FORTRAN/src/module_navier_stokes.f90` $\rightarrow$ `ns_apply_system_matrix` (lignes 380–809).

### 4.1.3. Analyse de l'Implémentation C++ (`NavierStokesSolver.cpp`, lignes 336–358)

```cpp
// (a) Viscosité d'élongation 2*mu d2u/dx2
double mu_cell_E = fields_old.mu[mesh.cell_idx(i, j)];
double mu_cell_W = fields_old.mu[mesh.cell_idx(i - 1, j)];
double diff_elong_E = 2.0 * mu_cell_E * (tf.u[mesh.u_idx(i + 1, j)] - tf.u[mesh.u_idx(i, j)]) / dx * dy;
double diff_elong_W = 2.0 * mu_cell_W * (tf.u[mesh.u_idx(i, j)] - tf.u[mesh.u_idx(i - 1, j)]) / dx * dy;
val -= (diff_elong_E - diff_elong_W);

// (b) Viscosité de cisaillement et rotation mu (d2u/dy2 + d2v/dxdy)
double mu_node_N = fields_old.mu_node[mesh.node_idx(i, j + 1)];
double mu_node_S = fields_old.mu_node[mesh.node_idx(i, j)];
double shear_N = mu_node_N * ((tf.u[mesh.u_idx(i, j + 1)] - tf.u[mesh.u_idx(i, j)]) / dy +
                              (tf.v[mesh.v_idx(i, j + 1)] - tf.v[mesh.v_idx(i - 1, j + 1)]) / dx) * dx;
double shear_S = mu_node_S * ((tf.u[mesh.u_idx(i, j)] - tf.u[mesh.u_idx(i, j - 1)]) / dy +
                              (tf.v[mesh.v_idx(i, j)] - tf.v[mesh.v_idx(i - 1, j)]) / dx) * dx;
val -= (shear_N - shear_S);
```

Notice : L'opérateur complet étant multiplié par le volume cellulaire $\Delta x \Delta y$, les flux surfaciques sont directement multipliés par $\Delta y$ (faces $X$) et $\Delta x$ (faces $Y$), garantissant la conservation symétrique.

---

## 4.2. Advection Convective et Linéarisations (Picard, AB2, Newton)

### 4.2.1. Équation Discrète Exacte

Dans l'opérateur implicite `apply_system_matrix`, la forme semi-implicite de la convection sur la composante $u_{i,j}$ avec vitesse advectante de référence $(\mathbf{u}_{\text{ref}}, \mathbf{v}_{\text{ref}})$ s'écrit avec flux aux faces $F_E = \frac{1}{2}(u_{\text{ref}, i,j} + u_{\text{ref}, i+1,j})\Delta y$, $F_W = \frac{1}{2}(u_{\text{ref}, i-1,j} + u_{\text{ref}, i,j})\Delta y$, $F_N = \bar{v}_{\text{ref}, N}\Delta x$, $F_S = \bar{v}_{\text{ref}, S}\Delta x$ :

$$\mathcal{C}_u(\delta \mathbf{u}) = \frac{\rho_{u, i, j}}{2} \left[ F_E (\delta u_{i+1,j} + \delta u_{i,j}) - F_W (\delta u_{i-1,j} + \delta u_{i,j}) + F_N (\delta u_{i,j+1} + \delta u_{i,j}) - F_S (\delta u_{i,j-1} + \delta u_{i,j}) \right]$$

- **Mode `"PICARD"` :** $\mathbf{u}_{\text{ref}} = \mathbf{u}^n$.
- **Mode `"ADAMS_BASHFORTH"` / `"AB2"` :** lorsque `has_nm1 == true`, $\mathbf{u}_{\text{ref}} = 2 \mathbf{u}^n - \mathbf{u}^{n-1}$ (extrapolation à $t^{n+1}$ depuis $n$ et $n-1$).
- **Mode `"NEWTON"` :** ajoute le jacobien croisé $\rho_u \left( \delta u_{i,j} \frac{u^n_{i+1,j} - u^n_{i-1,j}}{2\Delta x} + \bar{\delta v}_{i,j} \frac{u^n_{i,j+1} - u^n_{i,j-1}}{2\Delta y} \right) \Delta x \Delta y$, et compense le terme quadratique dans `assemble_rhs`.

### 4.2.2. Localisation dans le Code Source

- **C++17 :** `COVOF_2D-C++/src/NavierStokesSolver.cpp` $\rightarrow$ `apply_system_matrix` (lignes 230–332) et `assemble_rhs` (lignes 553–839).
- **Fortran 2018 :** `COVOF_2D-FORTRAN/src/module_navier_stokes.f90` $\rightarrow$ `ns_apply_system_matrix` (lignes 380–809) et `ns_assemble_rhs` (lignes 810–1124).

---

## 4.3. Lagrangien Augmenté Adaptatif (AL, AAL, C-AAL) et Mise à Jour d'Uzawa

### 4.3.1. Équation Discrète Exacte

1. **Pénalité spatio-temporelle locale AAL (Vincent et al. 2004, Éq. 14) :**
   Dans chaque cellule $(i, j)$ :
   $$r_{i,j} = K_{\text{aal}} \max\left( \alpha_{\text{time}} \rho_{i,j} \frac{L_0^2}{t_0}, \; \rho_{i,j} u_0 L_0, \; \rho_{i,j} \frac{L_0^2}{u_0}\|\mathbf{g}\|, \; \frac{p_0 L_0}{u_0}, \; \mu_{i,j}, \; \mu_{i,j} \frac{\Delta x \Delta y}{\mathcal{K}_{i,j}}, \; \frac{\sigma}{u_0} \right)$$
   En mode compressible (`C-AAL` ou `compressible && !is_aal_incompressible`), régularisation acoustique :
   $$r_{\text{eff}, i,j} = \frac{r_{i,j}}{1 + r_{i,j} \frac{\chi_{T, i,j}}{\Delta t}}$$

2. **Terme de divergence pénalisée dans $A \mathbf{x}$ (sur face $u_{i,j}$) :**
   $$-\left( r_{\text{eff}, i,j} D_{i,j} - r_{\text{eff}, i-1,j} D_{i-1,j} \right) \Delta y$$
   où $D_{i,j} = \frac{u_{i+1,j} - u_{i,j}}{\Delta x} + \frac{v_{i,j+1} - v_{i,j}}{\Delta y}$ (en mode `C-AAL` ou `MONOLITHIC` compressible, la ligne de pression inclut également $+\frac{\chi_T}{\Delta t} p^{n+1} \Delta x \Delta y$).

3. **Actualisation duale d'Uzawa :**
   $$p_{i,j}^{m+1} = p_{i,j}^m - r_{\text{eff}, i,j} \left[ D_{i,j}^{*, m+1} + \frac{\chi_{T,i,j}}{\Delta t}(p_{i,j}^m - p_{i,j}^n) \right]$$

### 4.3.2. Localisation dans le Code Source

- **C++17 :** `COVOF_2D-C++/src/NavierStokesSolver.cpp` $\rightarrow$ `compute_adaptative_r_al` (lignes 841–927), `update_pressure_al` (lignes 929–951), `apply_system_matrix` (lignes 361–365).
- **Fortran 2018 :** `COVOF_2D-FORTRAN/src/module_navier_stokes.f90` $\rightarrow$ `ns_compute_adaptative_r_al` (lignes 1125–1241), `ns_update_pressure_al` (lignes 1242–1286).

---

## 4.4. Pénalisation d'Obstacles Immergés (Darcy-Brinkman)

### 4.4.1. Équation Discrète Exacte

Pour toute face $u_{i,j}$ (ou $v_{i,j}$) détectée dans un obstacle (`BOX` ou `CYLINDER`) de perméabilité $\mathcal{K}$ :
- Si $\mathcal{K} \le 10^{-12}$ (solide imperméable) : la ligne de la matrice est remplacée par la condition de Dirichlet homogène $\left(\frac{\alpha_{\text{time}} \rho_u}{\Delta t} \Delta x \Delta y\right) \delta u_{i,j} = 0$.
- Si $\mathcal{K} > 10^{-12}$ (milieu poreux) : ajout implicite du terme de Darcy $+\frac{\mu_u}{\mathcal{K}} \Delta x \Delta y \, \delta u_{i,j}$ à l'opérateur et à la diagonale `diag_u` / `u_aP`.

### 4.4.2. Localisation dans le Code Source

- **C++17 :** `COVOF_2D-C++/include/covof/NavierStokesSolver.hpp` (`is_point_in_obstacle`, `is_u_face_in_obstacle`, `get_u_obstacle_perm`, lignes 101–202) et `src/NavierStokesSolver.cpp` (lignes 278–287 et 376–382).
- **Fortran 2018 :** `COVOF_2D-FORTRAN/src/module_mesh.f90` (`is_in_obstacle`, lignes 203–222) et `module_navier_stokes.f90` (`is_u_face_in_obstacle`, `is_v_face_in_obstacle`, lignes 82–131).

---

## 4.5. Suivi d'Interface VOF-PLIC et Advection de Weymouth-Yue

### 4.5.1. Équation Discrète Exacte

1. **Reconstruction affine par morceaux (`InterfaceSegment` : $n_x x + n_y y + \alpha = 0$ orienté vers $C = 0$) :**
   - `reconstruct_plic` : normale de Youngs $3\times 3$ pondérée $[1, 2, 1]$.
   - `reconstruct_elvira` : sélection du candidat minimisant l'erreur volumique $L_2$ sur le bloc $3\times 3$ parmi 6 fonctions de hauteur horizontales et verticales.
   - `compute_alpha` : bisection exacte sur l'intervalle $[\min_{4\text{ coins}} (-\mathbf{n}\cdot\mathbf{x}), \max_{4\text{ coins}} (-\mathbf{n}\cdot\mathbf{x})]$ utilisant l'aire polygonale exacte `cut_box_area` (clipping de Sutherland-Hodgman + formule de Gauss/Shoelace).

2. **Advection conservative de Weymouth & Yue (2010) avec splitting alterné de Strang :**
   - En début de pas : $C_{d,i,j} = 1$ si $C_{i,j}^n \ge 0.5$, sinon $0$.
   - Passe X : $C_{i,j}^* = C_{i,j}^n - \frac{f_{x, i+1, j} - f_{x, i, j}}{\Delta x \Delta y} + C_{d,i,j} \frac{\Delta t}{\Delta x}(u_{i+1,j} - u_{i,j})$ (sans écrêtage intermédiaire : `clamp_out = false`).
   - Passe Y : $C_{i,j}^{n+1} = C_{i,j}^* - \frac{f_{y, i, j+1} - f_{y, i, j}}{\Delta x \Delta y} + C_{d,i,j} \frac{\Delta t}{\Delta y}(v_{i,j+1} - v_{i,j})$ (avec écrêtage final `[0, 1]`).

### 4.5.2. Localisation dans le Code Source

- **C++17 :**
  - `COVOF_2D-C++/src/InterfaceReconstruction.cpp` : `cut_box_area` (l. 37), `compute_alpha` (l. 106), `clip_segment_in_cell` (l. 147), `reconstruct_plic` (l. 252), `reconstruct_elvira` (l. 321), `reconstruct_lvira` (l. 436).
  - `COVOF_2D-C++/src/InterfaceTracker.cpp` : `advect` (l. 849–894), `advect_weymouth_yue` (l. 1148–1327).
- **Fortran 2018 :**
  - `COVOF_2D-FORTRAN/src/module_interface_reconstruction.f90` : `cut_box_area` (l. 33), `compute_alpha` (l. 105), `clip_segment_in_cell` (l. 152), `reconstruct_elvira` (l. 346).
  - `COVOF_2D-FORTRAN/src/module_interface_vof.f90` : `advect_vof` (l. 792), `advect_vof_weymouth_yue` (l. 1162).

---

## 4.6. Méthode Hybride FC-WY-CVI-Marker (Variante B Exacte sur Boîtes MAC)

### 4.6.1. Formulation Mathématique Discrète

Pour tout point matériel $(x_p, y_p)$ localisé dans la cellule $(i, j)$, la **Variante B exacte** (Wang et al., 2015 ; Pusok et al., 2017) définit deux boîtes décalées distinctes restant strictement dans le halo `(0:nx+1, 0:ny+1)` :
- Boîte horizontale $u$ : $j_u = (y_p < y_j^c) \;?\; j-1 : j$, $\xi_u = \frac{x_p - x_i^f}{\Delta x}$, $\eta_u = \frac{y_p - y_{j_u}^c}{\Delta y}$.
- Boîte verticale $v$ : $i_v = (x_p < x_i^c) \;?\; i-1 : i$, $\xi_v = \frac{x_p - x_{i_v}^c}{\Delta x}$, $\eta_v = \frac{y_p - y_j^f}{\Delta y}$.

Avec $\Delta u_S = u_{i+1, j_u} - u_{i, j_u}$, $\Delta u_N = u_{i+1, j_u+1} - u_{i, j_u+1}$, $\Delta v_W = v_{i_v, j+1} - v_{i_v, j}$ et $\Delta v_E = v_{i_v+1, j+1} - v_{i_v+1, j}$ :

$$u_{\text{CVI}} = u_{\text{bilin}}(\xi_u, \eta_u) + \frac{\Delta x}{2 \Delta y} \xi_u (1 - \xi_u)(\Delta v_E - \Delta v_W)$$
$$v_{\text{CVI}} = v_{\text{bilin}}(\xi_v, \eta_v) + \frac{\Delta y}{2 \Delta x} \eta_v (1 - \eta_v)(\Delta u_N - \Delta u_S)$$

L'advection temporelle utilise le schéma RK2 de Heun (prédicteur $\mathbf{x}_p^*$ avec $\mathbf{u}^n$ puis correcteur avec $\mathbf{u}^{n+1}$ ou le champ courant en mode cinématique), suivie de la migration MPI 2 passes, de l'élagage par sous-octants et de l'ajustement parabolique pondéré dans le repère tourné $(\tau, \eta)$ :

$$\min_{a,b,c} \sum_{p=1}^{M} w_p \left(a \tau_p^2 + b \tau_p + c - \eta_p\right)^2 \quad \text{sous } \mathbf{n}_p \cdot \mathbf{n}^{\text{ref}} \ge 0.5 \implies \kappa_{i,j}^{\text{marker}} = \frac{2a}{(1 + b^2)^{3/2}}$$

### 4.6.2. Localisation dans le Code Source

- **C++17 :** `COVOF_2D-C++/src/InterfaceMarkers.cpp` :
  - `init_from_vof` (l. 16–109)
  - `interpolate_velocity_cvi` (l. 123–185)
  - `advect_rk2` (l. 239–495, incluant routage MPI aller-retour à 4 passes)
  - `migrate_mpi` (l. 497–679, migration X puis Y)
  - `prune_and_reseed` (l. 681–935)
  - `compute_curvature_from_markers` (l. 937–1384, incluant halo $2.5 h$ et lissage tangentiel)
  - `write_vtk` (l. 1390–1459)
- **Fortran 2018 :** `COVOF_2D-FORTRAN/src/module_markers.f90` :
  - `init_markers_from_vof` (l. 110), `interpolate_velocity_cvi` (l. 232), `advect_markers_rk2` (l. 327), `migrate_markers_mpi` (l. 763), `prune_and_reseed_markers` (l. 995), `compute_curvature_from_markers` (l. 1277), `write_markers_vtk` (l. 1785).
- **Couplage dans VOF et CSF :**
  - C++17 : `InterfaceTracker::advect` (`InterfaceTracker.cpp`, l. 887–893), `compute_curvature_and_CSF` (l. 658–733), `apply_marker_normals_to_segments` (l. 736–847).
  - Fortran 2018 : `module_interface_vof.f90` $\rightarrow$ `advect_vof` (l. 792+), `compute_curvature_and_csf` (l. 627+).

### 4.6.3. Analyse de l'Implémentation CVI (`InterfaceMarkers.cpp`, lignes 145–185)

```cpp
int j_u = (yp < mesh.yc(j)) ? (j - 1) : j;
int i_v = (xp < mesh.xc(i)) ? (i - 1) : i;

double xi_u  = std::clamp((xp - mesh.xu(i))   / dx, 0.0, 1.0);
double eta_u = std::clamp((yp - mesh.yc(j_u)) / dy, 0.0, 1.0);
double xi_v  = std::clamp((xp - mesh.xc(i_v)) / dx, 0.0, 1.0);
double eta_v = std::clamp((yp - mesh.yv(j))   / dy, 0.0, 1.0);
// ... lecture des 4 sommets MAC u_{i..i+1, j_u..j_u+1} et v_{i_v..i_v+1, j..j+1} ...
double corr_u = 0.5 * (dx / dy) * xi_u  * (1.0 - xi_u)  * (dv_e - dv_w);
double corr_v = 0.5 * (dy / dx) * eta_v * (1.0 - eta_v) * (du_n - du_s);
up = u_bilin + corr_u;
vp = v_bilin + corr_v;
```

Notice sur les conventions de signes : dans `InterfaceReconstruction`, les segments PLIC $n_x^{\text{PLIC}} x + n_y^{\text{PLIC}} y + \alpha = 0$ ont leur normale orientée vers $C = 0$ ($-\nabla C / \|\nabla C\|$), tandis que les marqueurs et le modèle CSF orientent $\mathbf{n}_p$ vers $C = 1$ ($+\nabla C / \|\nabla C\|$). À chaque ensemencement et lors de l'injection `apply_marker_normals_to_segments`, le signe est rigoureusement inversé (`mk.nx = -seg.nx`, `snx = -fields.marker_nx[idx]`).

---

## 4.7. Lissage SVOF Helmholtz et Préconditionneur Bloc Triangulaire BT-SBP

### 4.7.1. Lissage SVOF Helmholtz (`smooth_svof_helmholtz`)

- **Équation discrète :** $(1 + a_E + a_W + a_N + a_S) C_{s,i,j} - \sum_{\text{nb}} a_{\text{nb}} C_{s,\text{nb}} = C_{i,j}$ avec $a_{E,W} = k \Delta y / \Delta x$, $a_{N,S} = k \Delta x / \Delta y$ (annulés sur les bords physiques de Neumann homogène) et 25 balayages Gauss-Seidel avant/arrière synchronisés en MPI.
- **Localisation :** `InterfaceTracker.cpp` (l. 288–380) et `module_interface_vof.f90` (l. 387–478).

### 4.7.2. Préconditionneur Triangulaire par Blocs BT-SBP

- **Équation discrète :** Pour le système monolithique $\begin{pmatrix} F_u & 0 & G_x \\ 0 & F_v & G_y \\ D_x & D_y & -C_p \end{pmatrix}$, `BlockTriangularPreconditioner` résout d'abord le complément de Schur $S \mathbf{r}_p = \mathbf{z}_p$ (8 balayages de Gauss-Seidel symétrique sur les coefficients $a_W, a_E = \Delta y^2 / u_{aP}$, $a_S, a_N = \Delta x^2 / v_{aP}$, $a_P = \sum a_{\text{nb}} + \frac{\chi_T}{\Delta t}\Delta x\Delta y + 10^{-12}\Delta x\Delta y$, avec retrait du mode constant en incompressible), puis calcule $\text{rhs}_u = \mathbf{z}_u - G_x \mathbf{r}_p$, $\text{rhs}_v = \mathbf{z}_v - G_y \mathbf{r}_p$ et résout $F_u \mathbf{r}_u = \text{rhs}_u$ et $F_v \mathbf{r}_v = \text{rhs}_v$ par 4 balayages SSOR.
- **Localisation :**
  - C++17 : `COVOF_2D-C++/src/Preconditioner.cpp` $\rightarrow$ `BlockTriangularPreconditioner::update` (l. 210–382) et `apply` (l. 384–579).
  - Fortran 2018 : `COVOF_2D-FORTRAN/src/module_navier_stokes.f90` $\rightarrow$ `ns_assemble_preconditioner_blocks` (l. 1287–1585) et `module_linear_solvers.f90` $\rightarrow$ `apply_preconditioner` (l. 47–489).
