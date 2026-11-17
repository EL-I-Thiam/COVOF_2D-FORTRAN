# Formulation Mathématique et Modèles Physiques

## Système d'Équations Gouvernantes

Le solveur **COVOF 2D** résout les équations fondamentales de la mécanique des milieux continus régissant les écoulements 2D de fluides visqueux, newtoniens, monophasiques ou diphasiques, incompressibles ou faiblement compressibles (*All-Mach*), avec ou sans transfert thermique et en présence éventuelle d'obstacles solides ou poreux immergés.

### Équations de Navier-Stokes sous Forme Continue (*One-Fluid*)

Dans un domaine rectangulaire ouvert et borné $\Omega = [0, L_x] \times [0, L_y] \subset \mathbb{R}^2$ de frontière $\partial\Omega$ et pour $t \in [0, t_{\text{final}}]$, la formulation *One-Fluid* regroupe les deux phases en un milieu continu unique de masse volumique $\rho(\mathbf{x}, t)$ et de viscosité dynamique $\mu(\mathbf{x}, t)$ :

$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0$$

$$\frac{\partial (\rho \mathbf{u})}{\partial t} + \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u}) = -\nabla p + \nabla \cdot \boldsymbol{\tau} + \rho \mathbf{g} + \mathbf{f}_\sigma + \mathbf{f}_{\text{brinkman}} + \mathbf{f}_{\text{boussinesq}}$$

où :
- $\rho = \rho(\mathbf{x}, t)$ est la masse volumique locale du mélange $[\text{kg}\cdot\text{m}^{-3}]$.
- $\mathbf{u} = (u, v)^T$ est le vecteur vitesse eulérien $[\text{m}\cdot\text{s}^{-1}]$.
- $p = p(\mathbf{x}, t)$ est la pression hydrodynamique ou thermodynamique $[\text{Pa}]$.
- $\boldsymbol{\tau} = \mu \left( \nabla \mathbf{u} + (\nabla \mathbf{u})^T \right)$ est le tenseur des contraintes visqueuses newtoniennes en écoulement incompressible (avec contribution de dilatation $\lambda (\nabla \cdot \mathbf{u}) \mathbf{I}$ négligée en régime bas-Mach isotherme ou prise en compte dans les flux d'élongation).
- $\mathbf{g} = (g_x, g_y)^T$ est le vecteur accélération de la pesanteur $[\text{m}\cdot\text{s}^{-2}]$ (`physics.gx`, `physics.gy`).
- $\mathbf{f}_\sigma$ est la force volumique de tension superficielle interfaciale $[\text{N}\cdot\text{m}^{-3}]$.
- $\mathbf{f}_{\text{brinkman}}$ est la force de traînée de pénalisation de Darcy-Brinkman dans les obstacles ou milieux poreux $[\text{N}\cdot\text{m}^{-3}]$.
- $\mathbf{f}_{\text{boussinesq}}$ est la force de flottabilité thermique sous approximation de Boussinesq $[\text{N}\cdot\text{m}^{-3}]$.

### Formulation Advective Utilisée sur la Grille Décalée MAC

En exploitant la conservation de la masse, la forme advective non conservative résolue par le noyau implicite de COVOF 2D s'écrit, projetée selon $x$ et $y$ :

$$\rho \left( \frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} \right) = -\frac{\partial p}{\partial x} + \frac{\partial}{\partial x}\left( 2\mu \frac{\partial u}{\partial x} \right) + \frac{\partial}{\partial y}\left[ \mu \left( \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \right) \right] + \rho g_x + f_{\sigma, x} + f_{\text{brinkman}, x}$$

$$\rho \left( \frac{\partial v}{\partial t} + u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} \right) = -\frac{\partial p}{\partial y} + \frac{\partial}{\partial x}\left[ \mu \left( \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \right) \right] + \frac{\partial}{\partial y}\left( 2\mu \frac{\partial v}{\partial y} \right) + \rho g_y + f_{\sigma, y} + f_{\text{brinkman}, y} + f_{\text{boussinesq}, y}$$

Lorsque l'option `physics.enable_ssp_rk3_momentum = true` est activée (automatique par défaut en régime `compressible = true`), une étape préliminaire explicite **SSP-RK3** assure le transport inertiel conservatif de $(\rho u, \rho v)$ avant le couplage implicite vitesse-pression.

---

## Modèle Diphasique VOF, Lissage SVOF et Marqueurs Hybrides

### Fonction Indicatrice de Phase (Fraction Volumique $C$)

Dans la formulation *One-Fluid*, les deux fluides immiscibles — **Fluide 1** (`physics.fluid1`, généralement le fluide lourd/porteur ou ambiant) et **Fluide 2** (`physics.fluid2`, généralement la goutte, bulle, cylindre dense ou gaz selon le cas initial) — sont repérés par la fonction indicatrice de couleur $C(\mathbf{x}, t) \in [0, 1]$ (`fields.C` en C++, `fields%C` en Fortran) :

$$C(\mathbf{x}, t) = \begin{cases} 1 & \text{dans la Phase 1 (`fluid1`)} \\ 0 & \text{dans la Phase 2 (`fluid2`)} \\ ]0, 1[ & \text{dans les mailles mixtes coupées par l'interface } \Gamma(t) \end{cases}$$

Sans changement de phase, le transport lagrangien de $C$ obéit à :

$$\frac{\partial C}{\partial t} + \mathbf{u} \cdot \nabla C = 0 \quad \Longleftrightarrow \quad \frac{\partial C}{\partial t} + \nabla \cdot (C \mathbf{u}) = C (\nabla \cdot \mathbf{u})$$

### Régularisation Rhéologique SVOF par Filtre de Helmholtz

Conformément à El Ouafa, Vincent & Le Chenadec (*Fluids* 2021), lorsque le contraste de masse volumique ($\rho_1 / \rho_2 \gg 1$) ou de viscosité ($\mu_1 / \mu_2 \gg 1$) est élevé, COVOF 2D distingue rigoureusement :

1. **La fraction volumique advectée $C$** (`fields.C` / `fields%C`), qui demeure **parfaitement raide** (épaisseur d'une maille) afin de préserver la netteté géométrique de l'interface.
2. **La fonction couleur lissée rhéologique $C_s$** (`fields.C_s` / `fields%C_s`), activée par `physics.enable_svof_smoothing = true` et calculée de manière purement diagnostique après l'advection en résolvant l'équation elliptique de Helmholtz avec conditions aux limites de Neumann homogènes ($\partial_n C_s = 0$ sur $\partial\Omega$) :
   $$C_s - D \nabla^2 C_s = C, \qquad D = k \, \Delta x \Delta y$$
   où $k$ est le paramètre `physics.svof_filter_width` (défaut `1.0`). Ce système diagonalement dominant est résolu par 25 itérations de Gauss-Seidel symétrique parallèle avec bornage strict $C_s \in [0, 1]$.
3. **La fonction couleur filtrée capillaire $\widetilde{C}$** (`fields.C_tilde` / `fields%C_tilde`), activée par `physics.enable_diffusion_smoothing = true` (`diffusion_smoothing_steps`, défaut `3`), obtenue par $N_{\text{diff}}$ pas explicites de diffusion laplacienne ($\Delta \tau = 0.2 \min(\Delta x^2, \Delta y^2)$) pour régulariser le calcul classique de la courbure CSF lorsque les marqueurs ne sont pas activés.

### Lois de Mélange des Propriétés Matérielles sur Grille MAC

En notant $C_{\text{src}} = C_s$ si `enable_svof_smoothing = true` (sinon $C_{\text{src}} = C$) :

1. **Centres des cellules de pression $(i, j)$** — moyenne arithmétique :
   $$\rho_{i,j} = C_{\text{src}, i,j} \rho_1 + (1 - C_{\text{src}, i,j}) \rho_2, \qquad \mu_{i,j} = C_{\text{src}, i,j} \mu_1 + (1 - C_{\text{src}, i,j}) \mu_2$$
2. **Faces décalées de vitesse $(i - 1/2, j)$ et $(i, j - 1/2)$** — moyenne arithmétique des cellules adjacentes :
   $$\rho_{u, i, j} = \frac{1}{2}\left(\rho_{i-1, j} + \rho_{i, j}\right), \qquad \rho_{v, i, j} = \frac{1}{2}\left(\rho_{i, j-1} + \rho_{i, j}\right)$$
3. **Sommets / Nœuds de cisaillement $(i - 1/2, j - 1/2)$** — **moyenne harmonique** (`fields.mu_node` / `fields%mu_node`), indispensable pour assurer la continuité de la contrainte tangentielle à travers les interfaces à grand saut de viscosité :
   - En mode VOF standard (`enable_svof_smoothing = false`) :
     $$\mu_{\text{node}, i, j} = \frac{4}{\frac{1}{\mu_{i-1, j-1}} + \frac{1}{\mu_{i, j-1}} + \frac{1}{\mu_{i-1, j}} + \frac{1}{\mu_{i, j}}}$$
   - En mode SVOF (`enable_svof_smoothing = true`, Éqs. 8–9 *Fluids* 2021) :
     $$\mu_{\text{node}, i, j} = \frac{\mu_1 \mu_2}{\mu_2 \, C_{s,\text{node}} + \mu_1 (1 - C_{s,\text{node}})}, \qquad C_{s,\text{node}} = \frac{1}{4}\left( C_{s,i-1,j-1} + C_{s,i,j-1} + C_{s,i-1,j} + C_{s,i,j} \right)$$

En régime **monophasique** (`enable_two_phase = false`), COVOF 2D propose en outre les modèles de viscosité spatiale `physics.viscosity_model` :
- `"CONSTANT"` : $\mu = \mu_1$.
- `"EXPONENTIAL"` : $\mu(y) = \mu_0 \exp\left( a \left( \frac{y - e}{e} \right)^2 \right)$ (avec `visc_mu0`, `visc_a`, `visc_e`).
- `"ASYMMETRIC_EXPONENTIAL"` (en C++) : $\mu(y) = \mu_0 \exp\left( -a \frac{y - e}{e} \right)$.

### Modèle de Tension Superficielle CSF et Couplage Marqueurs Lagrangiens CVI

La tension superficielle est modélisée par la formulation volumique *Balanced-Force CSF* (Brackbill et al., 1992) discrétisée directement aux faces décalées MAC :

$$F_{\sigma, x, i-1/2, j} = \sigma \, \bar{\kappa}_{i-1/2, j} \, \frac{C_{i, j} - C_{i-1, j}}{\Delta x}, \qquad F_{\sigma, y, i, j-1/2} = \sigma \, \bar{\kappa}_{i, j-1/2} \, \frac{C_{i, j} - C_{i, j-1}}{\Delta y}$$

où $\sigma$ est la tension superficielle (`physics.sigma` $[\text{N}\cdot\text{m}^{-1}]$) et $\kappa$ la courbure algébrique plane $[\text{m}^{-1}]$. Deux modes de calcul de $(\mathbf{n}, \kappa)$ sont disponibles :

1. **Mode VOF eulérien standard :**
   À partir du champ couleur $C_{\text{curv}}$ ($\widetilde{C}$ si `enable_diffusion_smoothing`, sinon $C_s$ si `enable_svof_smoothing`, sinon $C$) :
   $$\mathbf{n}_{i,j} = \frac{\nabla_h C_{\text{curv}, i,j}}{\|\nabla_h C_{\text{curv}, i,j}\|}, \qquad \kappa_{i,j} = -\nabla_h \cdot \mathbf{n}_{i,j}$$
2. **Mode hybride FC-WY-CVI-Marker (`vof.use_markers = true` et `vof.marker_curvature = true`) :**
   Un nuage de marqueurs lagrangiens non connectés $(x_p, y_p, n_{x,p}, n_{y,p}, \kappa_p)$ est advecté à l'ordre 2 par le schéma RK2-Heun à l'aide de l'**Interpolation Conservative de la Vitesse (CVI — Variante B exacte sur boîtes décalées MAC)** garantissant $\nabla \cdot \mathbf{u}_{\text{CVI}} \equiv D_{i,j}$ en tout point. Dans chaque cellule mixte, la normale et la courbure $\kappa_{i,j}^{\text{marker}}$ sont obtenues par régression parabolique pondérée dans le repère local tangent-normal $(\tau, \eta)$ avec filtrage angulaire $\mathbf{n}_p \cdot \mathbf{n}^{\text{ref}} \ge 0.5$ pour protéger les ligaments minces. Cette approche élimine la majeure partie des courants parasites et préserve la courbure avec une convergence d'ordre 2.

---

## Modèle Faiblement Compressible, Régime All-Mach et Équations d'État (EOS)

Pour les écoulements diphasiques compressibles ou acoustiques à bas nombre de Mach (compression d'une bulle sous saut de pression, onde de choc acoustique faible, piston isotherme), COVOF 2D implémente la formulation **All-Mach C-AAL** (El Ouafa et al., *Fluids* 2024).

### Équation de Continuité Barotrope et Compressibilité Isotherme

En introduisant la compressibilité isotherme $\chi_T = \frac{1}{\rho}\left(\frac{\partial \rho}{\partial p}\right)_T = \frac{1}{\rho c_s^2}$ (`fields.chi_T` $[\text{Pa}^{-1}]$), la contrainte de continuité discrète s'écrit :

$$\chi_{T, i,j} \frac{p_{i,j}^{n+1} - p_{i,j}^n}{\Delta t} + (\nabla_h \cdot \mathbf{u}^{n+1})_{i,j} = 0$$

La masse volumique cellulaire est actualisée en fonction de la pression et de la température via l'équation d'état propre à chaque phase (`physics.fluid1.eos`, `physics.fluid2.eos`), puis moyennée par $C_{\text{src}}$ : $\rho = C_{\text{src}} \rho_1(p, T) + (1 - C_{\text{src}}) \rho_2(p, T)$ et $\chi_T = C_{\text{src}} \chi_{T,1} + (1 - C_{\text{src}}) \chi_{T,2}$.

### Catalogue des Quatre Lois d'État Implémentées (`EOS.hpp` / `module_eos.f90`)

| Mot-Clé `eos` (et Alias) | Loi de Masse Volumique $\rho(p, T)$ | Compressibilité Isotherme $\chi_T$ | Célérité du Son $c_s$ | Paramètres JSON Associés |
| :--- | :--- | :--- | :--- | :--- |
| `"INCOMPRESSIBLE"` | $\rho = \rho_0$ | $\chi_T = 0$ | $10^5\,\text{m}\cdot\text{s}^{-1}$ | `rho` |
| `"LINEAR"` / `"ACOUSTIC"` | $\rho = \rho_0 \left[ 1 + \chi_{T,0}(p - p_{\text{ref}}) \right]$ | $\chi_T = \chi_{T,0}$ | $1 / \sqrt{\rho \chi_{T,0}}$ | `rho`, `chi_T` ($p_{\text{ref}} = 1.01325\times 10^5\,\text{Pa}$) |
| `"IDEAL_GAS"` / `"PERFECT_GAS"` | $\rho = \frac{p}{r_{\text{gas}} T}$ | $\chi_T = \frac{1}{\max(1, p)}$ | $\sqrt{\gamma \, r_{\text{gas}} T}$ | `r_gas`, `gamma`, `cv` |
| `"NASG"` / `"STIFFENED_GAS"` | $\rho = \frac{p + \gamma p_\infty}{(\gamma - 1) c_v T + b_{\text{nasg}}(p + \gamma p_\infty)}$ | $\chi_T = \frac{1 - \rho b_{\text{nasg}}}{\rho \max(1, p + p_\infty)}$ | $1500\,\text{m}\cdot\text{s}^{-1}$ (repli liquide) | `gamma`, `cv`, `p_inf`, `b_nasg`, `q_nasg` |

En régime compressible non isotherme (`compressible = true` et `enable_energy = true`), la routine `solve_total_energy` résout l'équation de l'énergie totale spécifique $e = u_{\text{int}} + \frac{1}{2}\|\mathbf{u}\|^2$ (`fields.energy`, `fields.rho_e`), puis reconstruit la température via $T = \frac{\gamma - 1}{r_{\text{gas}}}\left(e - \frac{1}{2}\|\mathbf{u}\|^2\right)$.

---

## Modèle de Milieu Poreux et Pénalisation de Darcy-Brinkman

Pour modéliser des obstacles solides ou des matrices poreuses sans déformer la grille cartésienne MAC, COVOF 2D implémente la pénalisation volumique de Darcy-Brinkman (`obstacles` dans le fichier JSON) :

$$\mathbf{f}_{\text{brinkman}} = -\frac{\mu}{\mathcal{K}(\mathbf{x})} \mathbf{u}$$

Deux comportements distincts sont automatiquement sélectionnés selon la valeur de `permeability` ($\mathcal{K}$) de chaque obstacle (`BOX` ou `CYLINDER`) :
- **Solide rigide imperméable ($\mathcal{K} \le 10^{-12}\,\text{m}^2$, défaut $10^{-14}$ en C++ / $10^{-40}$ en Fortran) :** Les faces de vitesse $u$ et $v$ situées dans l'obstacle reçoivent directement la condition de Dirichlet nulle ($u = 0, v = 0$) dans l'opérateur linéaire et le second membre.
- **Milieu poreux perméable ($\mathcal{K} > 10^{-12}\,\text{m}^2$) :** Le terme de traînée implicite $+\frac{\mu}{\mathcal{K}} \Delta x \Delta y \, u$ est ajouté à la diagonale de la matrice de quantité de mouvement (et dans l'échelle $\tau_{\text{brinkman}} = \mu \frac{\Delta x \Delta y}{\mathcal{K}}$ du mode AAL).

---

## Transfert Thermique et Approximation de Boussinesq

En régime incompressible avec `physics.enable_energy = true` (et `compressible = false`), COVOF 2D résout l'équation de transport-diffusion de la température $T$ (`solve_energy` / `ns_solve_energy`) :

$$\rho C_p \left( \frac{\partial T}{\partial t} + \mathbf{u} \cdot \nabla T \right) = \nabla \cdot (k_{\text{th}} \nabla T)$$

Le couplage avec l'équation de quantité de mouvement verticale s'effectue par le terme de flottabilité de Boussinesq utilisant le nombre de Rayleigh `physics.Ra`, le nombre de Prandtl `physics.Pr`, le coefficient de dilatation `physics.beta` et la température de référence `physics.T_ref` :

$$\mathbf{f}_{\text{boussinesq}} = -\rho \, \beta \, (T - T_{\text{ref}}) \, \mathbf{g}$$

À chaque pas de temps, le solveur évalue et affiche le **nombre de Nusselt moyen sur la paroi chaude ($x = 0$)** via `compute_mean_nusselt_hot_wall`.

---

## Nombres Adimensionnels et Table Générale de Nomenclature

| Nombre Adimensionnel | Définition | Signification Physique |
| :--- | :--- | :--- |
| **Reynolds ($Re$)** | $Re = \frac{\rho U L}{\mu}$ | Rapport forces d'inertie / forces visqueuses |
| **Froude ($Fr$)** | $Fr = \frac{U}{\sqrt{g L}}$ | Rapport forces d'inertie / forces gravitationnelles |
| **Weber ($We$)** | $We = \frac{\rho U^2 L}{\sigma}$ | Rapport forces d'inertie / forces capillaires |
| **Bond / Eötvös ($Bo$)** | $Bo = \frac{\Delta \rho g L^2}{\sigma}$ | Rapport flottabilité gravitaire / tension de surface |
| **Capillaire ($Ca$)** | $Ca = \frac{\mu U}{\sigma}$ | Rapport forces visqueuses / tension de surface |
| **Ohnesorge ($Oh$)** | $Oh = \frac{\mu}{\sqrt{\rho \sigma L}}$ | Amortissement visqueux des ondes capillaires |
| **Prandtl ($Pr$)** | $Pr = \frac{\mu C_p}{k_{\text{th}}}$ | Rapport diffusivité cinématique / diffusivité thermique |
| **Rayleigh ($Ra$)** | $Ra = \frac{g \beta \Delta T L^3}{\nu \alpha}$ | Intensité motrice de la convection naturelle |
| **Darcy ($Da$)** | $Da = \frac{\mathcal{K}}{L^2}$ | Perméabilité adimensionnelle du milieu poreux |

### Table Générale de Nomenclature (Symboles $\leftrightarrow$ Variables C++ & Fortran)

| Symbole | Signification Physique | Unités SI | Identifiant C++ (`COVOF_2D-C++`) | Identifiant Fortran (`COVOF_2D-FORTRAN`) |
| :--- | :--- | :--- | :--- | :--- |
| $\mathbf{u} = (u,v)$ | Vitesse sur faces décalées MAC | $\text{m}\cdot\text{s}^{-1}$ | `fields.u`, `fields.v` | `fields%u`, `fields%v` |
| $p, p^n$ | Pression cellulaire et état précédent | $\text{Pa}$ | `fields.p`, `fields.p_old` | `fields%p`, `fields%p_old` |
| $C$ | Fraction volumique VOF raide | $-$ | `fields.C` | `fields%C` |
| $C_s$ | Fraction volumique lissée SVOF (Helmholtz) | $-$ | `fields.C_s` | `fields%C_s` |
| $\widetilde{C}$ | Fraction volumique diffusée pour CSF | $-$ | `fields.C_tilde` | `fields%C_tilde` |
| $\rho, \rho_u, \rho_v$ | Masse volumique (centres et faces) | $\text{kg}\cdot\text{m}^{-3}$ | `fields.rho`, `fields.rho_u`, `fields.rho_v` | `fields%rho`, `fields%rho_u`, `fields%rho_v` |
| $\mu, \mu_{\text{node}}$ | Viscosité dynamique (centres et nœuds) | $\text{Pa}\cdot\text{s}$ | `fields.mu`, `fields.mu_node` | `fields%mu`, `fields%mu_node` |
| $\mathbf{n}, \kappa$ | Normale unitaire et courbure d'interface | $-, \text{m}^{-1}$ | `fields.norm_x`, `fields.norm_y`, `fields.kappa` | `fields%curvature` (et normales locales) |
| $\mathbf{f}_\sigma$ | Force capillaire CSF aux faces | $\text{N}\cdot\text{m}^{-3}$ | `fields.F_sigma_x`, `fields.F_sigma_y` | `fields%fx_sigma`, `fields%fy_sigma` |
| $(x_p, y_p, \mathbf{n}_p, \kappa_p)$ | Marqueur lagrangien non connecté | $\text{m}, -, \text{m}^{-1}$ | `tracker.markers.cloud[m]` (`Marker`) | `global_marker_cloud%list(m)` (`marker_t`) |
| $\chi_T$ | Compressibilité isotherme | $\text{Pa}^{-1}$ | `fields.chi_T` | `fields%chi_T` |
| $r_{\text{al}}, r(t, M)$ | Pénalité AL constante et champ local AAL | $\text{Pa}\cdot\text{s}$ | `config.solver.r_al`, `fields.r_aal` | `cfg%r_al`, `fields%r_aal` |
| $T, C_p, k_{\text{th}}, \beta$ | Température et paramètres thermiques | $\text{K}, \text{J}/(\text{kg}\cdot\text{K}), \text{W}/(\text{m}\cdot\text{K}), \text{K}^{-1}$ | `fields.T`, `config.fluid1.Cp`, `config.fluid1.k_thermal`, `config.beta` | `fields%T`, `cfg%Cp`, `cfg%k_thermal`, `cfg%beta` |
| $e, \rho e$ | Énergie totale spécifique et volumique | $\text{J}\cdot\text{kg}^{-1}, \text{J}\cdot\text{m}^{-3}$ | `fields.energy`, `fields.rho_e` | `fields%energy`, `fields%rho_e` |
| $\mathcal{K}$ | Perméabilité de Darcy-Brinkman | $\text{m}^2$ | `config.obstacles[k].permeability` | `cfg%obs_perm` |
| $\sigma, \mathbf{g}$ | Tension de surface et gravité | $\text{N}\cdot\text{m}^{-1}, \text{m}\cdot\text{s}^{-2}$ | `config.sigma`, `config.gx`, `config.gy` | `cfg%sigma`, `cfg%gx`, `cfg%gy` |
| $\Delta x, \Delta y, \Delta t$ | Pas d'espace et pas de temps | $\text{m}, \text{s}$ | `mesh.dx`, `mesh.dy`, `dt` | `mesh%dx`, `mesh%dy`, `dt` |
