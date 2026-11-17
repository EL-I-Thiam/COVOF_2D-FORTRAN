# Méthodes Numériques

## Discrétisation Spatiale sur Grille Cartésienne Décalée (MAC)

Le solveur **COVOF 2D** emploie la discrétisation par volumes finis sur grille cartésienne décalée dite *Marker-and-Cell* (MAC, Harlow & Welch 1965). Ce positionnement entrelacé des inconnues garantit la condition inf-sup (Ladyzhenskaya-Babuška-Brezzi discrète) et élimine intrinsèquement tout mode parasite en damier (*checkerboard*) sur la pression sans recourir à l'interpolation de Rhie-Chow ni à une stabilisation artificielle.

### Localisation des Variables Primordiales

Sur une grille uniforme de pas $\Delta x = L_x / N_x$ et $\Delta y = L_y / N_y$, les grandeurs sont stockées aux positions géométriques suivantes (indices locaux $i \in [1, N_{x,\text{loc}}]$, $j \in [1, N_{y,\text{loc}}]$ entourés d'une couche fantôme $0$ et $N+1$) :

1. **Centres des cellules $(i, j)$ :**
   Pression $p_{i,j}$, fraction volumique VOF $C_{i,j}$, versions régularisées $C_{s,i,j}$ et $\widetilde{C}_{i,j}$, masse volumique $\rho_{i,j}$, viscosité $\mu_{i,j}$, température $T_{i,j}$, énergie totale $e_{i,j}$ et $\rho e_{i,j}$, compressibilité $\chi_{T,i,j}$, paramètre AAL $r_{\text{aal},i,j}$, courbure $\kappa_{i,j}$ et normale $\mathbf{n}_{i,j}$.
2. **Faces verticales $(i - \frac{1}{2}, j)$ (notées `u(i, j)` dans le code, face Ouest de la cellule $(i,j)$) :**
   Vitesse horizontale $u_{i,j}$, masse volumique de face $\rho_{u,i,j}$, force capillaire horizontale $F_{\sigma, x, i, j}$, quantité de mouvement intermédiaire $(\rho u)^*_{i,j}$.
3. **Faces horizontales $(i, j - \frac{1}{2})$ (notées `v(i, j)` dans le code, face Sud de la cellule $(i,j)$) :**
   Vitesse verticale $v_{i,j}$, masse volumique de face $\rho_{v,i,j}$, force capillaire verticale $F_{\sigma, y, i, j}$, flottabilité de Boussinesq et quantité de mouvement intermédiaire $(\rho v)^*_{i,j}$.
4. **Sommets / Nœuds $(i - \frac{1}{2}, j - \frac{1}{2})$ (notés `mu_node(i, j)` dans le code) :**
   Viscosité harmonique nodale $\mu_{\text{node}, i, j}$ utilisée pour discrétiser les dérivées croisées de cisaillement $\mu(\partial_y u + \partial_x v)$.

```text
                       v(i, j+1)  [face Nord de (i,j)]
                 +---------^---------+
                 |                   |
                 |  p(i,j), C(i,j)   |
      u(i, j)  ---->  rho, mu, T   ---->  u(i+1, j)
   [face Ouest]  |  kappa, chi_T     |  [face Est]
                 |                   |
                 +---------^---------+
                       v(i, j)    [face Sud de (i,j)]
     mu_node(i, j) au coin Sud-Ouest (i-1/2, j-1/2)
```

### Dualité Divergence-Gradient sur Grille MAC

Avec cette convention d'indiçage exacte du code (`u(i,j)` à gauche et `u(i+1,j)` à droite de la cellule `(i,j)`), la divergence discrète $(\nabla_h \cdot \mathbf{u})_{i,j}$ et le gradient de pression s'écrivent :

$$(\nabla_h \cdot \mathbf{u})_{i,j} = \frac{u_{i+1,j} - u_{i,j}}{\Delta x} + \frac{v_{i,j+1} - v_{i,j}}{\Delta y}$$

$$\left(\frac{\partial p}{\partial x}\right)_{i,j} = \frac{p_{i,j} - p_{i-1,j}}{\Delta x} \quad (\text{sur face } u_{i,j}), \qquad \left(\frac{\partial p}{\partial y}\right)_{i,j} = \frac{p_{i,j} - p_{i,j-1}}{\Delta y} \quad (\text{sur face } v_{i,j})$$

Cette dualité adjointe ($\mathbf{D} = -\mathbf{G}^T$ à un facteur volumique $\Delta x \Delta y$ près) assure la symétrie parfaite de la matrice de Stokes et l'annulation exacte de la divergence dans les solveurs monolithiques et AAL.

---

## Schémas Convectifs et Linéarisations de l'Inertie

Le traitement du terme non linéaire $(\mathbf{u} \cdot \nabla)\mathbf{u}$ est configuré par deux clés complémentaires dans `physics` (ou au niveau racine) :
- `physics.convective_scheme` : `"UPWIND"`, `"CENTERED"`, `"TVD_SUPERBEE"`.
- `physics.convective_linearization` : `"PICARD"`, `"ADAMS_BASHFORTH"` (alias `"AB2"`), `"NEWTON"`.

### Schémas Spatiaux (`convective_scheme`)

1. **Décentré Amont d'Ordre 1 (`"UPWIND"`) :**
   Inconditionnellement stable et monotone même à grand nombre de Reynolds de maille $Re_\Delta > 2$, au prix d'une diffusion numérique $\mathcal{O}(\Delta x)$.
2. **Centré d'Ordre 2 (`"CENTERED"`) :**
   Évalue les vitesses aux faces par moyenne arithmétique symétrique, atteignant l'ordre $\mathcal{O}(\Delta x^2)$ sans dissipation artificielle, idéal pour la DNS lorsque $Re_\Delta \le 2$.
3. **Haute Résolution TVD (`"TVD_SUPERBEE"`) :**
   Combine la précision d'ordre 2 dans les zones régulières et le limiteur non linéaire Superbee près des gradients raides.

### Linéarisations Temporelles du Terme Convectif (`convective_linearization`)

1. **Linéarisation de Picard (`"PICARD"`, défaut) :**
   La vitesse advectante est figée au pas connu $\mathbf{u}^n$ : $(\mathbf{u}^n \cdot \nabla)\mathbf{u}^{n+1}$. Simple et robuste, précise à l'ordre 1 en temps sur l'advection.
2. **Extrapolation d'Adams-Bashforth d'Ordre 2 (`"ADAMS_BASHFORTH"` / `"AB2"`) :**
   Lorsque l'historique $\mathbf{u}^{n-1}$ est disponible (`fields.has_nm1 == true`), la vitesse advectrice est extrapolée à l'ordre 2 selon l'Équation 10 de *Fluids* 2021 :
   $$\mathbf{u}_{\text{adv}} = 2 \mathbf{u}^n - \mathbf{u}^{n-1} \quad \text{ou} \quad \text{flux extrapolés } \frac{3}{2}\mathcal{C}(\mathbf{u}^n) - \frac{1}{2}\mathcal{C}(\mathbf{u}^{n-1})$$
   préservant la linéarité du système implicite tout en atteignant l'ordre 2 temporel.
3. **Linéarisation Complète de Newton (`"NEWTON"`) :**
   Inclut le terme jacobien croisé $(\delta \mathbf{u} \cdot \nabla)\mathbf{u}^n$ en plus de $(\mathbf{u}^n \cdot \nabla)\mathbf{u}^{n+1}$ dans l'opérateur matriciel `apply_system_matrix` (et le terme quadratique $-(\mathbf{u}^n \cdot \nabla)\mathbf{u}^n$ au second membre `assemble_rhs`), accélérant la convergence sur les écoulements stationnaires ou à grands pas de temps.

---

## Discrétisation Temporelle et Pas de Temps Adaptatif

### Schémas Temporels Implémentés (`time.time_scheme`)

1. **Euler Implicite d'Ordre 1 (`"EULER"`, défaut) :**
   $$\frac{\mathbf{u}^{n+1} - \mathbf{u}^n}{\Delta t}, \qquad \alpha_{\text{time}} = 1.0$$
2. **Schéma de Gear / BDF2 d'Ordre 2 (`"GEAR"` ou `"BDF2"`) :**
   Dès le pas 2 (`fields.has_nm1 == true`), le terme instationnaire utilise la formule rétrograde d'ordre 2 à 3 niveaux temporels :
   $$\frac{3 \mathbf{u}^{n+1} - 4 \mathbf{u}^n + \mathbf{u}^{n-1}}{2 \Delta t}, \qquad \alpha_{\text{time}} = 1.5$$
   (avec amorçage automatique par Euler implicite $\alpha_{\text{time}} = 1.0$ au pas 1).
3. **Transport Inertiel Conservatif `SSP-RK3` (`physics.enable_ssp_rk3_momentum = true`) :**
   Activé automatiquement lorsque `compressible = true` (ou manuellement), un schéma de Runge-Kutta TVD à 3 étages transporte la quantité de mouvement $(\rho u, \rho v)$ de façon conservative avant l'étape de couplage vitesse-pression (étape 1 de *Fluids* 2024).

### Pas de Temps Adaptatif (`time.adaptive_dt = true`)

À chaque itération, `Simulation::compute_time_step()` (`sim_compute_time_step` en Fortran) évalue la vitesse maximale globale $U_{\max} = \max_{i,j} \sqrt{|u_{i,j}|^2 + |v_{i,j}|^2}$ via `MPI_Allreduce(..., MPI_MAX)` et borne le pas de temps par :

$$\Delta t_{\text{CFL}} = \text{CFL}_{\max} \frac{\min(\Delta x, \Delta y)}{\max(U_{\max}, 10^{-6})}$$

$$\Delta t_{\sigma} = \begin{cases} \sqrt{\frac{(\rho_1 + \rho_2) \, \min(\Delta x, \Delta y)^3}{4 \pi \sigma}} & \text{si } \text{enable\_surface\_tension} \text{ et } \sigma > 0 \\ +\infty & \text{sinon} \end{cases}$$

$$\Delta t = \min\left(\Delta t_{\text{init}}, \, \Delta t_{\text{CFL}}, \, \Delta t_{\sigma}\right)$$

où $\Delta t_{\text{init}}$ est `time.dt` et $\text{CFL}_{\max}$ est `time.cfl_max` (le pas est également tronqué sur le dernier pas pour atteindre exactement `time.t_final`).

---

## Suivi d'Interface VOF et Méthode Hybride FC-WY-CVI-Marker

### Reconstruction Géométrique et Advection VOF (`vof.*`)

Le module `InterfaceTracker` (`module_interface_vof` en Fortran) met à disposition deux familles d'advection de la fraction volumique $C$ :

1. **Schémas géométriques (`"WEYMOUTH_YUE"` par défaut, et `"CIAM"`) :**
   - Couplés aux algorithmes de reconstruction affine $n_x x + n_y y + \alpha = 0$ : `"PLIC"` (Youngs), `"ELVIRA"` (6 candidats de fonctions de hauteur, ordre 2 exact sur droites), `"LVIRA"` (moindres carrés non linéaires sur $3\times 3$) ou `"SLIC"`.
   - **Schéma conservatif de Weymouth & Yue (2010)** avec splitting alterné de Strang (X puis Y aux pas pairs, Y puis X aux pas impairs) : fixe l'indicatrice binaire de dilatation $C_{d,i,j} = \mathbf{1}_{C_{i,j}^n \ge 0.5}$ au début du pas et transporte :
     $$C_{i,j}^* = C_{i,j}^n - \frac{\Delta V_{i+1,j}^x - \Delta V_{i,j}^x}{\Delta x \Delta y} + C_{d,i,j} \frac{\Delta t}{\Delta x}(u_{i+1,j}^{n+1} - u_{i,j}^{n+1})$$
     $$C_{i,j}^{n+1} = C_{i,j}^* - \frac{\Delta V_{i,j+1}^y - \Delta V_{i,j}^y}{\Delta x \Delta y} + C_{d,i,j} \frac{\Delta t}{\Delta y}(v_{i,j+1}^{n+1} - v_{i,j}^{n+1})$$
     où $\Delta V$ est calculé exactement par découpe géométrique Sutherland-Hodgman + Shoelace (`cut_box_area`). Sous $\text{CFL}_{\text{vof}} < 0.5$ et $\nabla_h \cdot \mathbf{u}^{n+1} = 0$, la somme des flux télescope et **conserve le volume à la précision du solveur linéaire**.
2. **Schéma algébrique (`"VOF_TVD"` / `"TVD"`) :**
   Sans reconstruction géométrique (`reconstruction_method` forcée à `"NONE"`), piloté par les limiteurs `"SUPERBEE"`, `"VAN_LEER"`, `"MUSCL"` ou `"HRIC"`.

### Méthode Hybride FC-WY-CVI-Marker (Variante B Exacte sur Boîtes MAC Décalées)

Lorsqu'on active `vof.use_markers = true` (ou `markers.enabled = true`), COVOF 2D associe au transport VOF-WY un **nuage de marqueurs lagrangiens non connectés** (`InterfaceMarkers` / `module_markers`) sans arêtes ni maillage surfacique, totalement insensible aux ruptures et coalescences.

#### 1. Interpolation Conservative de la Vitesse (CVI — Variante B exacte)
Selon le théorème de Liouville, l'interpolation bilinéaire standard sur grille MAC produit une divergence continue non nulle d'ordre $\mathcal{O}(h)$ qui crée des amas artificiels de marqueurs. COVOF 2D implémente la **Variante B exacte de Wang et al. (2015) / Pusok et al. (2017)** (`interpolate_velocity_cvi`) en interpolant chaque composante sur **sa propre boîte décalée MAC** :
- Boîte pour $u$ : $\Omega_u(i, j_u) = [x_i^f, x_{i+1}^f] \times [y_{j_u}^c, y_{j_u+1}^c]$ avec $j_u = j - 1$ si $y_p < y_j^c$ sinon $j$, et coordonnées réduites $(\xi_u, \eta_u) \in [0, 1]^2$.
- Boîte pour $v$ : $\Omega_v(i_v, j) = [x_{i_v}^c, x_{i_v+1}^c] \times [y_j^f, y_{j+1}^f]$ avec $i_v = i - 1$ si $x_p < x_i^c$ sinon $i$, et coordonnées réduites $(\xi_v, \eta_v) \in [0, 1]^2$.
- Formule CVI exacte avec corrections quadratiques solénoïdales ($\Delta v_E = v_{i_v+1,j+1} - v_{i_v+1,j}$, $\Delta v_W = v_{i_v,j+1} - v_{i_v,j}$, $\Delta u_N = u_{i+1,j_u+1} - u_{i,j_u+1}$, $\Delta u_S = u_{i+1,j_u} - u_{i,j_u}$) :
  $$u_{\text{CVI}}(x_p, y_p) = u_{\text{bilin}}(\xi_u, \eta_u) + \frac{\Delta x}{2 \Delta y} \, \xi_u (1 - \xi_u) (\Delta v_E - \Delta v_W)$$
  $$v_{\text{CVI}}(x_p, y_p) = v_{\text{bilin}}(\xi_v, \eta_v) + \frac{\Delta y}{2 \Delta x} \, \eta_v (1 - \eta_v) (\Delta u_N - \Delta u_S)$$
En dérivant analytiquement dans chacun des 4 sous-quadrants de $(i, j)$, tous les termes linéaires croisés s'annulent identiquement :
$$\frac{\partial u_{\text{CVI}}}{\partial x} + \frac{\partial v_{\text{CVI}}}{\partial y} \equiv \frac{u_{i+1,j} - u_{i,j}}{\Delta x} + \frac{v_{i,j+1} - v_{i,j}}{\Delta y} = D_{i,j}$$
tout en restant strictement dans la couche unique de halo `(0:nx+1, 0:ny+1)`. *(L'option `marker_advection_scheme = "BILINEAR_RK2"` permet de désactiver la correction quadratique à des fins d'étude comparative).*

#### 2. Advection RK2-Heun, Migration MPI, Élagage et Courbure Parabolique
À l'intérieur de `InterfaceTracker::advect` (`advect_vof`), juste après l'advection VOF et le lissage SVOF/diffusif et juste avant `compute_curvature_and_CSF` :
1. **Advection RK2-Heun (`advect_rk2` / `advect_markers_rk2`) :** Prédit $\mathbf{x}_p^* = \mathbf{x}_p^n + \Delta t \, \mathbf{u}_{\text{CVI}}(\mathbf{x}_p^n; \mathbf{u}^n)$, effectue si nécessaire un routage MPI aller-retour à 4 passes pour évaluer $\mathbf{u}_{\text{CVI}}(\mathbf{x}_p^*; \mathbf{u}^{n+1})$ sur le rang propriétaire, puis corrige $\mathbf{x}_p^{n+1} = \mathbf{x}_p^n + \frac{\Delta t}{2}[\mathbf{u}_1 + \mathbf{u}_2]$.
2. **Migration MPI en 2 passes (`migrate_mpi` / `migrate_markers_mpi`) :** Passe X (Ouest $\leftrightarrow$ Est avec repliement périodique éventuel) puis Passe Y (Sud $\leftrightarrow$ Nord) sans nécessiter de voisins diagonaux.
3. **Élagage et ré-ensemencement sélectif (`prune_and_reseed`) :**
   - Supprime les marqueurs tombés en phase pure ($C < 10^{-6}$ ou $C > 1 - 10^{-6}$ dans le voisinage $3\times 3$) et élague par sous-octants si une cellule dépasse `marker_max_per_cell` (défaut `8`).
   - Ensemence uniquement les cellules interfaciales dépeuplées ($N_p = 0$) avec `markers_per_cell` marqueurs (défaut `4`) sans jamais écraser ni projeter les marqueurs actifs existants.
4. **Courbure et normale dans le repère tourné $(\tau, \eta)$ (`compute_curvature_from_markers`) :**
   - Échange les marqueurs fantômes dans une bande $2.5 h$ entre rangs MPI.
   - Filtre les marqueurs voisins par proximité et cohérence angulaire $\mathbf{n}_p \cdot \mathbf{n}_{i,j}^{\text{ref}} \ge 0.5$ (protégeant les filaments et jupes minces contre le mélange des faces opposées).
   - Ajuste une parabole pondérée $\eta(\tau) = a \tau^2 + b \tau + c$ par moindres carrés $3\times 3$ :
     $$\mathbf{n}_{i,j}^{\text{marker}} = \frac{\mathbf{e}_\eta - b \mathbf{e}_\tau}{\sqrt{1 + b^2}}, \qquad \kappa_{i,j}^{\text{marker}} = \frac{2a}{(1 + b^2)^{3/2}}$$
     avec écrêtage Nyquist $|\kappa| \le 2/\min(\Delta x, \Delta y)$, passe conservative de lissage tangentiel et repli sur la courbure CSF si $M < 4$.
   - Si `marker_normals_plic = true`, injecte $\mathbf{n}^{\text{marker}}$ (sous réserve de cohérence angulaire $\ge 0.5$ avec le gradient VOF) dans les segments PLIC pour l'étape suivante (`apply_marker_normals_to_segments`).

---

## Formulations du Couplage Vitesse-Pression

COVOF 2D propose cinq formulations pilotées par `solver.type` :

### 1. Lagrangien Augmenté Standard (`"AUGMENTED_LAGRANGIAN"`, défaut)
Ajoute le terme de pénalité $-r_{\text{al}} \nabla(\nabla \cdot \mathbf{u}^{k+1})$ à l'équation de quantité de mouvement avec un scalaire uniforme `solver.r_al`, puis actualise la pression par boucle d'Uzawa (`solver.uzawa_iter`) :

$$p^{k+1} = p^k - r_{\text{al}} (\nabla_h \cdot \mathbf{u}^{k+1})$$

avec arrêt anticipé dès que $\|\nabla_h \cdot \mathbf{u}^{k+1}\|_\infty \le \text{uzawa\_tol}$.

### 2. Lagrangien Augmenté Adaptatif Incompressible (`"AAL"` / `"ADAPTATIVE_AUGMENTED_LAGRANGIAN"`)
Détermine cellule par cellule le paramètre optimal $r(t, M)$ (Vincent et al., *Computers & Fluids* 33, 2004, Éq. 14) en pondérant les sept échelles physiques locales par `solver.k_aal` ($K_{\text{aal}} \in [10, 1000]$, défaut `100.0`) :

$$r_{i,j} = K_{\text{aal}} \max\left( \alpha_{\text{time}}\rho \frac{L_0^2}{t_0}, \; \rho u_0 L_0, \; \rho \frac{L_0^2}{u_0}\|\mathbf{g}\|, \; \frac{p_0 L_0}{u_0}, \; \mu, \; \mu \frac{\Delta x \Delta y}{\mathcal{K}}, \; \frac{\sigma}{u_0} \right)$$

Les échelles de référence (`l0_ref`, `t0_ref`, `u0_ref`, `p0_ref`) sont calculées automatiquement si laissées à `0.0`.

### 3. Lagrangien Augmenté Compressible Bas-Mach (`"C-AAL"` / `"CAAL"`)
Active automatiquement `compressible = true` et `adaptative_al = true`. Régularise la pénalité dans les zones compressibles ($\chi_T > 0$) et actualise la pression selon la contrainte barotrope :

$$r_{\text{eff}, i,j} = \frac{r_{i,j}}{1 + r_{i,j} \frac{\chi_{T,i,j}}{\Delta t}}, \qquad p_{i,j}^{k+1} = p_{i,j}^k - r_{\text{eff}, i,j} \left[ (\nabla_h \cdot \mathbf{u}^{k+1})_{i,j} + \frac{\chi_{T,i,j}}{\Delta t}(p_{i,j}^k - p_{i,j}^n) \right]$$

### 4. Résolution Monolithique Tout-Couplé (`"MONOLITHIC"` / `"FULLY_COUPLED"` / `"FC"`)
Assemble simultanément la vitesse et la pression dans le système point-selle sans pénalité ($r = 0$) :

$$\begin{pmatrix} \mathbf{A}_u & \mathbf{A}_{uv} & \mathbf{G}_x \\ \mathbf{A}_{vu} & \mathbf{A}_v & \mathbf{G}_y \\ \mathbf{D}_x & \mathbf{D}_y & -\mathbf{C}_p \end{pmatrix} \begin{pmatrix} \mathbf{u}^{n+1} \\ \mathbf{v}^{n+1} \\ \mathbf{p}^{n+1} \end{pmatrix} = \begin{pmatrix} \mathbf{f}_u \\ \mathbf{f}_v \\ \mathbf{g}_p \end{pmatrix}$$

où $\mathbf{C}_p = \text{diag}(\chi_T \Delta x \Delta y / \Delta t)$ en compressible ($\mathbf{0}$ en incompressible). La ligne de continuité n'étant pas divisée par $\rho$, ce système garantit une divergence discrète à la tolérance linéaire `solver.tol` (typiquement $10^{-10}$) quel que soit le saut $\rho_1 / \rho_2$.

### 5. Méthode de Projection Standard (`"STANDARD_PROJECTION"` / `"PROJECTION"`)
Résout successivement la prédiction visqueuse $\mathbf{u}^*$, l'équation de Poisson en pression $\nabla \cdot (\frac{1}{\rho}\nabla \phi) = \frac{1}{\Delta t}\nabla \cdot \mathbf{u}^*$, puis corrige $\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho}\nabla \phi$.

---

## Solveurs Linéaires de Krylov et Préconditionneurs

COVOF 2D implémente une architecture **Matrix-Free** (l'opérateur complet $A \mathbf{x}$ est appliqué par parcours des pochoirs MAC avec synchronisation MPI des halos via `apply_system_matrix`) complétée par des blocs locaux pour le préconditionnement (`Fuu`, `Fvv`, `SchurLaplace` en CSR côté C++, stencils 5-points côté Fortran).

### Solveurs Linéaires (`solver.linear_solver`)

| Mot-Clé | Algorithme | Domaine d'Emploi Recommandé |
| :--- | :--- | :--- |
| `"BICGSTAB2"` (défaut) | BiCGSTAB(2) (Sleijpen & van der Vorst 1993) | Opérateurs asymétriques convectifs, AAL, C-AAL et Monolithique |
| `"BICGSTAB"` | BiCGSTAB standard (van der Vorst 1992) | Systèmes asymétriques modérément raides |
| `"CG"` / `"CONJUGATE_GRADIENT"` | Gradient Conjugué Préconditionné | Opérateurs symétriques définis positifs (Stokes pur, Poisson, diffusion) |
| `"GMRES"` | GMRES(30) redémarré (Arnoldi + Givens) | Systèmes fortement non normaux |
| `"MINRES"` | Minimal Residual (Paige & Saunders 1975) | Systèmes point-selle symétriques indéfinis |
| `"DIRECT"` / `"MUMPS"` | Factorisation LU LAPACK (`dgetrf`/`dgetrs`) | Petits maillages de référence ou validation exacte (avec mise en cache si $\Delta t$ constant) |

*(Remarque : C++ dispose également des passerelles `PetscSolverBridge` (`"PETSC"`) et `HypreSolverBridge` (`"HYPRE"`) qui basculent automatiquement sur `BiCGSTAB2Solver` si PETSc/HYPRE ne sont pas liés).*

### Préconditionneurs (`solver.precond_type` / `solver.preconditioner`)

| Mot-Clé (et Alias) | Formulation | Recommandation |
| :--- | :--- | :--- |
| `"JACOBI"` / `"DIAGONAL"` (défaut) | Inversion de la diagonale locale enrichie (`diag_u`, `diag_v`, `diag_p`) | Optimal et 100 % parallèle pour `AAL` et `C-AAL` |
| `"BLOCK_TRIANGULAR"` / `"BT_SBP"` / `"SCHUR"` | Factorisation par blocs triangulaires $P = \begin{pmatrix} D_u & G \\ 0 & S \end{pmatrix}$ avec 8 balayages SSOR sur le complément de Schur $S$ et 4 balayages sur $F_u, F_v$ | Indispensable pour le mode `"MONOLITHIC"` |
| `"SSOR"` / `"GAUSS_SEIDEL"` | Lissage avant/arrière symétrique ($\omega = 1.0$) | Écoulements modérément visqueux |
| `"ILU0"` / `"ILU"` | Factorisation incomplète ILU(0) locale par sous-domaine | Systèmes raides sur sous-domaines compacts |
| `"CUSTOM"` | Hook utilisateur (`CustomPreconditioner::set_custom_operator` / `custom_precond_hook`) | Prototypage de nouveaux préconditionneurs |
| `"NONE"` | Identité ($P = I$) | Diagnostic de conditionnement |
