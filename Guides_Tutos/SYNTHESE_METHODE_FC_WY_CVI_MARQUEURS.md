
---

# 1. Positionnement Scientifique, État de l'Art et Originalité de l'Approche

## 1.1 Le Défi Fondamental des Écoulements Diphasiques Incompressibles à Forts Contrastes

La simulation numérique directe (DNS) d'écoulements diphasiques incompressibles mettant en jeu de grands rapports de masse volumique ($\rho_1 / \rho_2 \in [10^2, 10^4]$), de grands sauts de viscosité ($\mu_1 / \mu_2 \in [10^2, 10^6]$) et des effets capillaires intenses ($\sigma > 0$) se heurte historiquement à une triple difficulté :

1. **La conservation exacte de la masse de chaque phase** : Les méthodes de suivi d'interface implicites comme le *Level-Set* souffrent d'une perte ou d'un gain systématique de masse lors de l'étape de réinitialisation par équation eikonale $|\nabla \phi| = 1$. Les méthodes *Volume-of-Fluid* (VOF) à reconstruction géométrique affine par morceaux (*Piecewise Linear Interface Calculation* — PLIC) conservent la masse à condition que le schéma d'advection directionnel traite exactement les termes de dilatation croisés et que le champ de vitesse discret soit rigoureusement à divergence nulle ($\nabla_h \cdot \mathbf{u} = 0$).
2. **La résolution sous-maille des structures filamentaires et des fortes courbures** : Dans une méthode VOF-PLIC eulérienne pure, l'interface est représentée par un unique segment de droite par maille :
   $$n_x x + n_y y + \alpha = 0.$$
   Dès que l'épaisseur caractéristique $\delta$ d'une nappe liquide, d'une jupe gazeuse ou d'une vrille tourbillonnaire devient inférieure au pas de grille ($\delta < \Delta x$), la reconstruction eulérienne brise artificiellement le ligament en gouttelettes numériques (*numerical pinch-off* ou *flotsam*). De surcroît, l'évaluation de la courbure $\kappa = -\nabla \cdot \mathbf{n}$ par différentiation finie de la fraction volumique $C$ (ou de sa régularisée $\widetilde{C}$) introduit des erreurs d'ordre $O(1)$ à l'échelle de la maille, générant des **courants parasites** (*spurious currents*) qui polluent l'équilibre capillaire.

3. **Le couplage vitesse-pression à grands contrastes rhéologiques** : Dans les solveurs découplés à projection (Chorin-Temam, Bell-Colella-Glaz), la contrainte d'incompressibilité est ramenée à une équation de Poisson à coefficients discontinus $\nabla \cdot \left(\frac{1}{\rho}\nabla p\right) = \frac{1}{\Delta t}\nabla \cdot \mathbf{u}^*$. À grand rapport de densité ($\rho_1/\rho_2 \ge 1000$), le mauvais conditionnement de cet opérateur dégrade la divergence discrète finale $\nabla_h \cdot \mathbf{u}^{n+1}$, ce qui compromet à la fois la stabilité cinétique et la conservation du volume VOF.

## 1.2 Architecture Hybride Proposée pour la Plateforme `COVOF_2D`

Pour lever simultanément ces trois verrous au sein des moteurs **Fortran 2018** (`COVOF_2D-FORTRAN`) et **C++17** (`COVOF_2D-C++`), la présente synthèse établit les fondements théoriques, l'analyse critique et le cadre algorithmique d'un couplage à quatre piliers :

1. **Le solveur Navier-Stokes monolithique entièrement couplé (*Fully Coupled* — FC)** avec pénalisation visqueuse harmonique aux nœuds et lissage rhéologique SVOF par équation de Helmholtz ([El Ouafa, Vincent & Le Chenadec, *Fluids* 2021](https://doi.org/10.3390/fluids6010023)) ;
2. **Le schéma d'advection VOF-PLIC à conservation exacte de Weymouth & Yue (2010)** avec séparation directionnelle de Strang alternée ([Weymouth & Yue, 2010](https://doi.org/10.1016/j.jcp.2009.12.018) ; [Arrufat et al., 2018](https://arxiv.org/html/1811.12327v2)) ;
3. **Un nuage de marqueurs lagrangiens non connectés** (*Unconnected Interface Marker Cloud*) distribués uniquement au voisinage immédiat de l'interface (portant position $(x_p, y_p)$, normale orientée $(n_{x,p}, n_{y,p})$ et courbure locale $\kappa_p$) ;
4. **L'Interpolation Conservative de Vitesse (*Conservative Velocity Interpolation* — CVI)** sur grille décalée MAC ([Wang, Agrusta & van Hunen, 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC5089062/) ; [Pusok, Kaus & Popov, 2017](https://www.researchgate.net/publication/311165935_On_the_Quality_of_Velocity_Interpolation_Schemes_for_Marker-in-Cell_Method_and_Staggered_Grids)), garantissant analytiquement $\nabla \cdot \mathbf{u}_{\text{CVI}}(\mathbf{x}) \equiv 0$ en tout point continu $(x,y)$ de chaque maille.

## 1.3 Comparaison Critique avec l'État de l'Art International

Le tableau ci-dessous positionne l'architecture proposée vis-à-vis des grandes familles de méthodes hybrides VOF/Particules/Front-Tracking de la littérature :

| Méthode de la littérature | Représentation de l'interface | Conservation du volume | Interpolation de vitesse | Couplage Navier-Stokes | Complexité topologique & MPI |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Mixed Markers & VOF** ([Aulisa et al., 2003–2025](https://pubs.aip.org/aip/acp/article/3315/1/400022/3363173/Marker-reconstruction-by-using-best-fit-quadric)) | Marqueurs connectés par chaînes + quadriques locales | Très bonne via correction VOF | Bilinéaire / quadratique locale | Découplé (Projection / SIMPLE) | Reconnexion explicite lourde lors des ruptures |
| **Hybrid Particle-VOF** ([Karnakov et al., 2020](https://arxiv.org/html/1906.00314v2)) | VOF-WY + particules non connectées | Exacte (WY) | Bilinéaire sur quad/octree | Projection (Bell-Colella-Glaz) | Divergence continue non nulle sous interpolation bilinéaire |
| **Intersection Marker (ISM)** ([Ho et al., 2018](https://www.emerald.com/hff/article/28/1/118/193408/Bubble-flow-simulations-using-the-intersection)) | Marqueurs aux intersections grille-interface | Locale par marqueurs | Linéaire sur faces | Découplé | Reconnexion topologique explicite requise |
| **EBIT (Edge-Based Tracking)** ([Chirco & Zaleski, 2022–2025](https://arxiv.org/html/2509.26361v2)) | Marqueurs contraints sur les arêtes | Quasi-conservative (sans VOF) | Linéaire 1D sur arêtes | Projection | Topologie implicite par règles de coloriage |
| **One-Layer Particle Level-Set** ([OPLS, 2019](https://www.emerald.com/insight/content/doi/10.1108/hff-07-2019-0547/full/html)) | Level-Set + 1 couche de particules | Corrigée mais non exacte ($10^{-4}$) | Bilinéaire + RK3 | Projection | Dérive de masse résiduelle du Level-Set |
| **Sharp Front-Tracking + Jump** ([Gorges et al., 2022](https://arxiv.org/abs/2202.10223) ; [Evrard et al., 2025](https://arxiv.org/html/2505.06932v1)) | Maillage surfacique connecté | Exacte sous correction globale | Interpolation à divergence nulle | Projection / IIM | Chirurgie de maillage surfacique complexe en MPI |
| **Méthode `COVOF_2D` (FC + WY + CVI + Marqueurs non connectés)** | **VOF-PLIC (WY) + Nuage de marqueurs non connectés** | **Exacte ($10^{-12}$)** via WY + $B_u\mathbf{u}^{n+1}=0$ monolithique | **CVI analytique** ($\nabla\cdot\mathbf{u}_{\text{CVI}}\equiv 0$) | **Monolithique FC** (Vitesse-Pression couplées) | **100 % locale, sans connectivité, routage MPI X-Y à 4 voisins** |

## 1.4 Les Trois Avantages Décisifs de l'Architecture

1. **Synergie mathématique entre le solveur Monolithique (FC), Weymouth-Yue et la CVI** :
   Dans `COVOF_2D`, le solveur monolithique résout le système selle-point complet où la ligne de continuité $B_u \mathbf{u}^{n+1} = 0$ est purement géométrique (non divisée par $\rho$). Le résidu de divergence discret $D_{i,j}$ est donc nul à la précision du solveur linéaire ($\|B_u \mathbf{u}^{n+1}\|_\infty \le 10^{-10}$), même lorsque $\rho_1/\rho_2 = 10^4$ et $\mu_1/\mu_2 = 10^5$. Cette divergence discrète nulle verrouille **simultanément** la conservation exacte du volume dans Weymouth & Yue et l'annulation analytique ponctuelle $\nabla \cdot \mathbf{u}_{\text{CVI}}(x,y) \equiv 0$ dans l'interpolation CVI.

2. **Absence totale de connectivité lagrangienne** :
   Les marqueurs forment un simple nuage de points $(x_p, y_p, n_{x,p}, n_{y,p}, \kappa_p)$ sans arêtes, sans listes chaînées et sans triangulation. Lors d'une coalescence de bulles ou du pincement d'un ligament, aucune chirurgie topologique n'est nécessaire : la topologie globale est portée par la fraction volumique $C$ (VOF-PLIC), tandis que les marqueurs fournissent la géométrie fine sous-maille (normale et courbure).

3. **Parallélisme MPI natif sur halo de largeur 1** :
   Toute l'interpolation CVI et la migration des marqueurs s'appuient exclusivement sur la couche unique de cellules fantômes `(0:nx+1, 0:ny+1)` déjà présente dans `COVOF_2D` et sur les 4 voisins cardinaux MPI existants via un routage dimensionnel en deux passes (Passe X puis Passe Y).

---

# 2. Fondements Théoriques et Mathématiques Complets

## 2.1 Modèle « One-Fluid » et Conditions de Saut Interfaciales

Soit un domaine borné $\Omega = \Omega_1(t) \cup \Omega_2(t) \subset \mathbb{R}^2$ séparé par une interface mobile $\Gamma(t) = \partial\Omega_1(t) \cap \partial\Omega_2(t)$ de normale unitaire $\mathbf{n}$ orientée du fluide 1 vers le fluide 2. La fonction caractéristique de phase $C(\mathbf{x}, t) \in \{0, 1\}$ (dont la moyenne par cellule est la fraction volumique VOF $C_{i,j} \in [0, 1]$) satisfait l'équation de transport lagrangienne :

$$\frac{\partial C}{\partial t} + \mathbf{u} \cdot \nabla C = 0.$$

Sous l'hypothèse d'incompressibilité de chaque phase pure, le système de Navier-Stokes *One-Fluid* s'écrit au sens des distributions sur $\Omega$ :

$$\nabla \cdot \mathbf{u} = 0,$$

$$\rho(C)\left( \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} \right) = -\nabla p + \rho(C)\mathbf{g} + \nabla \cdot \left[ \mu(C)\left( \nabla \mathbf{u} + \nabla \mathbf{u}^T \right) \right] + \mathbf{F}_\sigma,$$

où la force capillaire singulière est modélisée par $\mathbf{F}_\sigma = \sigma \kappa \mathbf{n} \delta_\Gamma$, avec $\sigma$ le coefficient de tension superficielle constant, $\kappa = -\nabla_s \cdot \mathbf{n}$ la courbure algébrique plane et $\delta_\Gamma$ la mesure de Dirac surfacique portée par $\Gamma(t)$.

À la traversée de l'interface $\Gamma(t)$ sans changement de phase, les conditions de saut physiques imposent :

- **Continuité cinématique de la vitesse** :
  $$[\![ \mathbf{u} ]\!]_\Gamma = \mathbf{0} \implies [\![ \nabla_s \mathbf{u} ]\!]_\Gamma = \mathbf{0}.$$

- **Continuité de la contrainte visqueuse tangentielle** :
  $$\left[\!\left[ \mu \left( \frac{\partial u_\tau}{\partial n} + \frac{\partial u_n}{\partial \tau} \right) \right]\!\right]_\Gamma = 0 \implies \left[\!\left[ \mu \frac{\partial u_\tau}{\partial n} \right]\!\right]_\Gamma = -[\![ \mu ]\!]_\Gamma \frac{\partial u_n}{\partial \tau}.$$
  Lorsque $\mu_1 / \mu_2 \gg 1$, la dérivée normale de la vitesse tangentielle $\frac{\partial u_\tau}{\partial n}$ subit donc une discontinuité franche (un « coude » du profil de vitesse à l'interface).

- **Saut de contrainte normale (Loi de Young-Laplace dynamique)** :
  $$[\![ p ]\!]_\Gamma - 2 \left[\!\left[ \mu \frac{\partial u_n}{\partial n} \right]\!\right]_\Gamma = \sigma \kappa.$$

## 2.2 Régularisation Rhéologique SVOF (Helmholtz) et Pénalisation Visqueuse dans `COVOF_2D`

Conformément à [El Ouafa, Vincent & Le Chenadec (*Fluids* 2021)](https://doi.org/10.3390/fluids6010023), pour stabiliser le système monolithique en présence de forts contrastes $\rho_1/\rho_2$ et $\mu_1/\mu_2$, `COVOF_2D` distingue rigoureusement :

1. **La fraction volumique advectée $C$** (`fields%C`), qui demeure **parfaitement raide** (discrétisée sur une seule cellule d'épaisseur par VOF-PLIC) et qui n'est **jamais diffusée** ;
2. **La fonction couleur lissée rhéologique $C_s$** (`fields%C_s`), recalculée de manière purement diagnostique à chaque pas de temps par la résolution de l'équation elliptique de Helmholtz (`smooth_svof_helmholtz`) :
   $$C_s - D \nabla^2 C_s = C, \qquad D = k \, \Delta x \Delta y, \quad k \in [0.2, 1.5],$$
   avec conditions aux limites de Neumann homogènes $\frac{\partial C_s}{\partial n} = 0$ sur $\partial\Omega$.

3. **La fonction couleur filtrée capillaire $\widetilde{C}$** (`fields%C_tilde`), obtenue par $N_{\text{diff}} \in \{2, 3\}$ pas explicites de diffusion laplacienne (`vof_smooth_diffusion`) pour l'évaluation classique du modèle *Continuum Surface Force* (CSF) lorsque les marqueurs ne sont pas activés.

Les propriétés physiques sont ensuite injectées sur la grille décalée MAC selon la **méthode de pénalisation visqueuse** :

- **Centres des cellules de pression $(i,j)$** (où sont évaluées la divergence et les dérivées diagonales $\mu \partial_x u$, $\mu \partial_y v$) — moyenne arithmétique :
  $$\rho_{i,j} = C_{s,i,j}\rho_1 + (1 - C_{s,i,j})\rho_2, \qquad \mu_{i,j} = C_{s,i,j}\mu_1 + (1 - C_{s,i,j})\mu_2.$$

- **Faces de vitesse $(i-1/2, j)$ et $(i, j-1/2)$** — moyenne arithmétique de la densité :
  $$\rho_{u, i,j} = \frac{1}{2}(\rho_{i-1,j} + \rho_{i,j}), \qquad \rho_{v, i,j} = \frac{1}{2}(\rho_{i,j-1} + \rho_{i,j}).$$

- **Sommets / Nœuds de cisaillement $(i-1/2, j-1/2)$** (où est discrétisé le flux visqueux croisé $\mu(\partial_y u + \partial_x v)$) — **moyenne harmonique** (`fields%mu_node`) assurant la continuité du flux de cisaillement tangentiel :
  $$\mu_{\text{node}, i,j} = \left( \frac{C_{s,\text{node}}}{\mu_1} + \frac{1 - C_{s,\text{node}}}{\mu_2} \right)^{-1}, \qquad C_{s,\text{node}} = \frac{1}{4}\left( C_{s,i-1,j-1} + C_{s,i,j-1} + C_{s,i-1,j} + C_{s,i,j} \right).$$

## 2.3 Système Monolithique Entièrement Couplé (FC) et Invariance de la Divergence

Avec le schéma temporel d'Euler implicite ou de Gear d'ordre 2 (BDF2 : $\frac{\alpha_0 \mathbf{u}^{n+1} + \alpha_1 \mathbf{u}^n + \alpha_2 \mathbf{u}^{n-1}}{\Delta t}$) et la linéarisation convective d'Adams-Bashforth, la discrétisation spatiale par volumes finis sur grille MAC conduit au système algébrique selle-point monolithique (`SOLVER_MONOLITHIC = 2`) :

$$\begin{pmatrix} \mathbf{A}_u(\rho^n, \mu^n, \mathbf{u}^n) & \mathbf{A}_{uv}(\mu^n) & \mathbf{G}_x \\ \mathbf{A}_{vu}(\mu^n) & \mathbf{A}_v(\rho^n, \mu^n, \mathbf{u}^n) & \mathbf{G}_y \\ \mathbf{D}_x & \mathbf{D}_y & \mathbf{0} \end{pmatrix} \begin{pmatrix} \mathbf{u}^{n+1} \\ \mathbf{v}^{n+1} \\ \mathbf{p}^{n+1} \end{pmatrix} = \begin{pmatrix} \mathbf{f}_u^n + \mathbf{f}_{\sigma, x}^n \\ \mathbf{f}_v^n + \mathbf{f}_{\sigma, y}^n \\ \mathbf{0} \end{pmatrix}.$$

En notant $\mathbf{B}_u = (\mathbf{D}_x, \mathbf{D}_y)$ l'opérateur de divergence discret MAC, la troisième ligne du système impose directement dans chaque cellule $(i,j) \in [1, N_x] \times [1, N_y]$ :

$$D_{i,j}^{n+1} = \frac{u_{i+1,j}^{n+1} - u_{i,j}^{n+1}}{\Delta x} + \frac{v_{i,j+1}^{n+1} - v_{i,j}^{n+1}}{\Delta y} = 0,$$

*(selon la convention d'indiçage de `COVOF_2D` où `u(i,j)` est la face Ouest $x_{i-1/2}$ et `u(i+1,j)` la face Est $x_{i+1/2}$ de la cellule $(i,j)$).*
Résolu par **BiCGStab(2)** préconditionné par le préconditionneur bloc-triangulaire à complément de Schur approché (**`BLOCK_TRIANGULAR` / `BT_SBP`**), ce système garantit $\| \mathbf{D}^{n+1} \|_\infty \le \epsilon_{\text{tol}}$ (typiquement $10^{-10}$ à $10^{-12}$) indépendamment du rapport $\rho_1/\rho_2$.

## 2.4 Schéma VOF-PLIC Conservatif de Weymouth & Yue (2010)

Pour transporter $C$ sans perte ni gain de masse, Weymouth & Yue (2010) intègrent l'équation d'advection sous forme conservative corrigée de la dilatation directionnelle en fixant l'indicatrice binaire de dilatation au début du pas de temps $t^n$ (`c_d(i,j)` dans `advect_vof_weymouth_yue`) :

$$C_{d, i,j} = \begin{cases} 1 & \text{si } C_{i,j}^n \ge 0.5, \\ 0 & \text{si } C_{i,j}^n < 0.5. \end{cases}$$

Pour un pas pair du splitting de Strang (Passe X puis Passe Y) :

1. **Passe X** (avec reconstruction PLIC $n_x x + n_y y + \alpha = 0$ sur $C^n$) :
   $$C_{i,j}^* = C_{i,j}^n - \frac{\Delta V_{i+1,j}^x - \Delta V_{i,j}^x}{\Delta x \Delta y} + C_{d, i,j} \frac{\Delta t}{\Delta x}\left( u_{i+1,j}^{n+1} - u_{i,j}^{n+1} \right),$$
   où $\Delta V_{i,j}^x$ est l'aire géométrique exacte de fluide 1 coupée par le demi-plan PLIC dans la boîte de flux amont de largeur $|u_{i,j}^{n+1}|\Delta t$ (`cut_box_area`).

2. **Passe Y** (avec nouvelle reconstruction PLIC sur $C^*$, mais en conservant le **même** $C_{d, i,j}$) :
   $$C_{i,j}^{n+1} = C_{i,j}^* - \frac{\Delta V_{i,j+1}^y - \Delta V_{i,j}^y}{\Delta x \Delta y} + C_{d, i,j} \frac{\Delta t}{\Delta y}\left( v_{i,j+1}^{n+1} - v_{i,j}^{n+1} \right).$$

En additionnant les deux passes et en sommant sur l'ensemble des cellules $(i,j)$, les flux géométriques de faces $\Delta V^x$ et $\Delta V^y$ se télescopent exactement, laissant :

$$\sum_{i,j} C_{i,j}^{n+1} - \sum_{i,j} C_{i,j}^n = \Delta t \sum_{i,j} C_{d, i,j} \underbrace{\left( \frac{u_{i+1,j}^{n+1} - u_{i,j}^{n+1}}{\Delta x} + \frac{v_{i,j+1}^{n+1} - v_{i,j}^{n+1}}{\Delta y} \right)}_{=\, D_{i,j}^{n+1}}.$$

Puisque le solveur monolithique FC assure $|D_{i,j}^{n+1}| \le 10^{-10}$ et que $C_{d,i,j} \in \{0, 1\}$ évite tout dépassement ($C_{i,j}^{n+1} \in [0, 1]$ sans écrêtage destructif sous $\text{CFL} = \max(|u|\Delta t/\Delta x, |v|\Delta t/\Delta y) < 0.5$), **la conservation du volume est exacte à la précision du solveur**.

## 2.5 Théorème de Liouville et Dérivation Rigoureuse de l'Interpolation CVI sur Grille MAC

### A. Condition de Liouville pour les marqueurs lagrangiens
Considérons un marqueur matériel $\mathbf{x}_p(t)$ advecté par le champ de vitesse continu interpolé $\mathbf{u}_h(\mathbf{x}, t) = (u_h(x,y,t), v_h(x,y,t))^T$ :

$$\frac{d\mathbf{x}_p}{dt} = \mathbf{u}_h(\mathbf{x}_p(t), t).$$

Le déterminant jacobien $J_p(t) = \det\left(\frac{\partial \mathbf{x}_p(t)}{\partial \mathbf{x}_p(0)}\right)$ mesurant la dilatation d'un élément de surface matériel obéit à l'identité d'Euler-Liouville :

$$\frac{dJ_p}{dt} = \left( \frac{\partial u_h}{\partial x} + \frac{\partial v_h}{\partial y} \right)_{\mathbf{x}=\mathbf{x}_p} J_p(t).$$

Pour que l'advection lagrangienne ne crée ni amas artificiels (*clustering*), ni raréfaction dans les zones de cisaillement ou de point d'arrêt, il est nécessaire et suffisant que la divergence analytique du champ interpolé soit identiquement nulle en tout point $(x,y)$ à l'intérieur de chaque cellule :

$$\nabla \cdot \mathbf{u}_h(x, y) = \frac{\partial u_h(x,y)}{\partial x} + \frac{\partial v_h(x,y)}{\partial y} \equiv 0, \qquad \forall (x,y) \in \Omega_{i,j}.$$

### B. Pourquoi l'interpolation bilinéaire standard viole le théorème de Liouville
Sur la cellule $(i,j) = [x_i^f, x_{i+1}^f] \times [y_j^f, y_{j+1}^f]$ de coordonnées réduites :

$$\xi = \frac{x - x_i^f}{\Delta x} \in [0, 1], \qquad \eta = \frac{y - y_j^f}{\Delta y} \in [0, 1],$$

notons $u_{SW}, u_{SE}, u_{NW}, u_{NE}$ les vitesses horizontales aux 4 coins (ou nœuds de la boîte d'interpolation) et $v_{SW}, v_{SE}, v_{NW}, v_{NE}$ les vitesses verticales correspondantes. L'interpolation bilinéaire standard s'écrit :

$$u_{\text{bilin}}(\xi, \eta) = (1 - \eta)\left[ (1 - \xi)u_{SW} + \xi u_{SE} \right] + \eta\left[ (1 - \xi)u_{NW} + \xi u_{NE} \right],$$

$$v_{\text{bilin}}(\xi, \eta) = (1 - \xi)\left[ (1 - \eta)v_{SW} + \eta v_{NW} \right] + \xi\left[ (1 - \eta)v_{SE} + \eta v_{NE} \right].$$

Calculons ses dérivées partielles exactes ($\partial_x = \frac{1}{\Delta x}\partial_\xi$ et $\partial_y = \frac{1}{\Delta y}\partial_\eta$) :

$$\frac{\partial u_{\text{bilin}}}{\partial x}(\xi, \eta) = \frac{(1 - \eta)\Delta u_S + \eta \Delta u_N}{\Delta x} = \frac{\Delta u_S + \Delta u_N}{2 \Delta x} + \left(\eta - \frac{1}{2}\right)\frac{\Delta u_N - \Delta u_S}{\Delta x},$$

$$\frac{\partial v_{\text{bilin}}}{\partial y}(\xi, \eta) = \frac{(1 - \xi)\Delta v_W + \xi \Delta v_E}{\Delta y} = \frac{\Delta v_W + \Delta v_E}{2 \Delta y} + \left(\xi - \frac{1}{2}\right)\frac{\Delta v_E - \Delta v_W}{\Delta y},$$

où $\Delta u_S = u_{SE} - u_{SW}$, $\Delta u_N = u_{NE} - u_{NW}$, $\Delta v_W = v_{NW} - v_{SW}$ et $\Delta v_E = v_{NE} - v_{SE}$.
En sommant ces deux dérivées, même si la divergence moyenne est nulle, il subsiste un **résidu de divergence non nul dépendant linéairement de $\xi$ et $\eta$** :

$$\nabla \cdot \mathbf{u}_{\text{bilin}}(\xi, \eta) = \left(\eta - \frac{1}{2}\right)\frac{\Delta u_N - \Delta u_S}{\Delta x} + \left(\xi - \frac{1}{2}\right)\frac{\Delta v_E - \Delta v_W}{\Delta y} \neq 0.$$

### C. Formulation Mathématique Exacte de l'Interpolation CVI pour `COVOF_2D`
Pour annuler identiquement ce résidu en tout point $(\xi, \eta) \in [0, 1]^2$, on remarque que :

$$\frac{\partial}{\partial x}\left[ \frac{1}{2}\xi(1 - \xi) \right] = \frac{1 - 2\xi}{2\Delta x} = -\frac{\xi - \frac{1}{2}}{\Delta x}, \qquad \frac{\partial}{\partial y}\left[ \frac{1}{2}\eta(1 - \eta) \right] = \frac{1 - 2\eta}{2\Delta y} = -\frac{\eta - \frac{1}{2}}{\Delta y}.$$

L'**Interpolation Conservative de Vitesse (CVI)** ([Wang et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC5089062/) ; [Pusok et al., 2017](https://www.researchgate.net/publication/311165935_On_the_Quality_of_Velocity_Interpolation_Schemes_for_Marker-in-Cell_Method_and_Staggered_Grids)) consiste donc à ajouter à $u_{\text{bilin}}$ le terme quadratique en $\xi(1-\xi)$ qui compense la variation en $\xi$ de $\partial_y v_{\text{bilin}}$, et à ajouter à $v_{\text{bilin}}$ le terme quadratique en $\eta(1-\eta)$ qui compense la variation en $\eta$ de $\partial_x u_{\text{bilin}}$ :

$$u_{\text{CVI}}(\xi, \eta) = u_{\text{bilin}}(\xi, \eta) + \frac{\Delta x}{2 \Delta y} \, \xi(1 - \xi) \left( \Delta v_E - \Delta v_W \right),$$

$$v_{\text{CVI}}(\xi, \eta) = v_{\text{bilin}}(\xi, \eta) + \frac{\Delta y}{2 \Delta x} \, \eta(1 - \eta) \left( \Delta u_N - \Delta u_S \right).$$

**Preuve directe de divergence analytique nulle :**
En dérivant analytiquement $u_{\text{CVI}}$ par rapport à $x$ et $v_{\text{CVI}}$ par rapport à $y$ :

$$\frac{\partial u_{\text{CVI}}}{\partial x}(\xi, \eta) = \frac{\Delta u_S + \Delta u_N}{2 \Delta x} + \left(\eta - \frac{1}{2}\right)\frac{\Delta u_N - \Delta u_S}{\Delta x} - \left(\xi - \frac{1}{2}\right)\frac{\Delta v_E - \Delta v_W}{\Delta y},$$

$$\frac{\partial v_{\text{CVI}}}{\partial y}(\xi, \eta) = \frac{\Delta v_W + \Delta v_E}{2 \Delta y} + \left(\xi - \frac{1}{2}\right)\frac{\Delta v_E - \Delta v_W}{\Delta y} - \left(\eta - \frac{1}{2}\right)\frac{\Delta u_N - \Delta u_S}{\Delta x}.$$

Par addition terme à terme, les variations linéaires en $(\xi - 1/2)$ et $(\eta - 1/2)$ s'éliminent **exactement** pour tout $(\xi, \eta) \in [0, 1]^2$, laissant uniquement la divergence discrète :

$$\nabla \cdot \mathbf{u}_{\text{CVI}}(\xi, \eta) \equiv \frac{\Delta u_S + \Delta u_N}{2 \Delta x} + \frac{\Delta v_W + \Delta v_E}{2 \Delta y} = 0.$$

**Adaptation exacte aux cellules de pression $(i,j)$ de `COVOF_2D` (Halo de largeur 1) :**
Dans `COVOF_2D`, où la cellule $(i,j)$ possède un halo unique `(0:nx+1, 0:ny+1)`, on peut construire la CVI :

- **Variante A (Par cellule de pression $(i,j)$, $C^0$ aux sommets)** : Les valeurs aux 4 coins $(SW, SE, NW, NE)$ de la cellule $(i,j)$ sont obtenues par moyenne sur les faces adjacentes :
  $$u_{SW} = \frac{u_{i,j-1} + u_{i,j}}{2}, \quad u_{SE} = \frac{u_{i+1,j-1} + u_{i+1,j}}{2}, \quad u_{NW} = \frac{u_{i,j} + u_{i,j+1}}{2}, \quad u_{NE} = \frac{u_{i+1,j} + u_{i+1,j+1}}{2},$$
  $$v_{SW} = \frac{v_{i-1,j} + v_{i,j}}{2}, \quad v_{SE} = \frac{v_{i,j} + v_{i+1,j}}{2}, \quad v_{NW} = \frac{v_{i-1,j+1} + v_{i,j+1}}{2}, \quad v_{NE} = \frac{v_{i,j+1} + v_{i+1,j+1}}{2}.$$
  Alors $\frac{\Delta u_S + \Delta u_N}{2\Delta x} + \frac{\Delta v_W + \Delta v_E}{2\Delta y} = \frac{1}{4}\left(D_{i,j-1} + 2 D_{i,j} + D_{i,j+1}\right) + \dots = 0$ dès que $D_{k,l} = 0$ dans la cellule et ses 4 voisines ! De plus, tous les indices $(i\pm 1, j\pm 1)$ pour $(i,j) \in [1, N_x] \times [1, N_y]$ appartiennent **strictement** à `[0:nx+1, 0:ny+1]` sans jamais déborder en `nx+2` !

- **Variante B (Par boîtes décalées de Pusok et al., 2017)** : On interpole $u$ dans la cellule décalée en $y$ encadrant le marqueur ($y \in [y_j^c, y_{j+1}^c]$ si $\eta \ge 0.5$, ou $[y_{j-1}^c, y_j^c]$ si $\eta < 0.5$) avec $\Delta u_S = u_{i+1,j} - u_{i,j}$ et $\Delta u_N = u_{i+1,j+1} - u_{i,j+1}$, et symétriquement pour $v$ dans la cellule décalée en $x$ ($x \in [x_i^c, x_{i+1}^c]$ ou $[x_{i-1}^c, x_i^c]$) :
  $$u_{\text{CVI}}(\xi_u, \eta_u) = (1 - \eta_u)\left[(1-\xi_u)u_{i,j_u} + \xi_u u_{i+1,j_u}\right] + \eta_u\left[(1-\xi_u)u_{i,j_u+1} + \xi_u u_{i+1,j_u+1}\right],$$
  corrigé par le terme quadratique en $\xi_u(1-\xi_u)$ issu de la divergence nulle de $D_{i,j_u}$ et $D_{i,j_u+1}$. Les deux variantes utilisent uniquement les indices `0` à `nx+1` et garantissent $\nabla \cdot \mathbf{u}_{\text{CVI}} \equiv 0$.

---

# 3. Analyse Critique d'Expert : Verrous de `COVOF_2D`, Pièges Numériques et Solutions

L'audit croisé de la théorie et du code source de `PROJET_COVOF-2D` permet d'identifier **six pièges numériques critiques** et de formuler les solutions exactes qui les lèvent :

## 3.1 Verrou n°1 : Exactitude polynomiale du terme correctif CVI

- **Problème identifié** : Dans certaines retranscriptions de la littérature, un facteur supplémentaire $(\eta - 1/2)$ est ajouté par erreur devant le terme quadratique $\xi(1-\xi)$ de $u_{\text{CVI}}$. Comme démontré à la Section 2.5.C, cela introduit une dérivée parasite en $(1-2\xi)(\eta - 1/2)$ qui détruit la propriété $\nabla \cdot \mathbf{u}_{\text{CVI}} = 0$.
- **Solution retenue dans `COVOF_2D`** : Coder la paire symétrique exacte où la correction de $u_{\text{CVI}}$ est purement quadratique en $\xi(1-\xi)$ et la correction de $v_{\text{CVI}}$ est purement quadratique en $\eta(1-\eta)$, assortie d'un test unitaire automatique vérifiant par différences finies centrées à $\delta = 10^{-7}$ que $\left|\frac{\partial u_{\text{CVI}}}{\partial x} + \frac{\partial v_{\text{CVI}}}{\partial y}\right| \le 10^{-11}$ en tout point aléatoire $(x,y)$ du domaine.

## 3.2 Verrou n°2 : Le Paradoxe du *Reseeding* PLIC (« Effet Esclave »)

- **Problème identifié** : Si l'on ré-ensemence les marqueurs à chaque pas de temps sur les segments PLIC ou si l'on projette orthogonalement les marqueurs existants sur la droite PLIC reconstruite à partir de $C^{n+1}$, les marqueurs perdent toute mémoire lagrangienne sous-maille et redeviennent une simple copie discrète de la reconstruction VOF-PLIC.
- **Solution retenue dans `COVOF_2D`** :
  1. **Aucune projection sur PLIC** pour les marqueurs actifs situés dans une cellule interfaciale ($\epsilon_{\text{vof}} < C_{i,j} < 1 - \epsilon_{\text{vof}}$) : leurs positions évoluent librement selon $\mathbf{u}_{\text{CVI}}$.
  2. **Ensemencement sélectif uniquement dans les cellules interfaciales dépeuplées ($N_p = 0$)** : Lorsqu'une interface pénètre dans une nouvelle cellule ou qu'un étirement vide une cellule ($N_p = 0$), on y sème $N_{\text{target}}$ marqueurs :
     - Soit par **interpolation parabolique locale** appuyée sur les marqueurs des cellules voisines si au moins 3 marqueurs cohérents existent dans le voisinage $3\times 3$ ;
     - Soit, en repli (*fallback*), uniformément le long du segment PLIC local $(n_x, n_y, \alpha)$.

## 3.3 Verrou n°3 : Précision Temporelle d'Ordre 2 du Schéma RK2 (Heun)

- **Problème identifié** : Advecter les marqueurs entre $t^n$ et $t^{n+1}$ par un schéma RK2 utilisant uniquement le champ figé $\mathbf{u}^{n+1}$ n'est que d'ordre 1 en temps pour un écoulement instationnaire.
- **Solution retenue dans `COVOF_2D`** : Exploiter les tableaux `fields%u_old` ($\mathbf{u}^n$) et `fields%v_old` ($\mathbf{v}^n$) déjà présents dans `module_fields.f90` (lignes 60–61) et `Simulation.hpp`. Puisque $\nabla_h \cdot \mathbf{u}^n = 0$ et $\nabla_h \cdot \mathbf{u}^{n+1} = 0$, le schéma prédicteur-correcteur de Heun :
  $$\mathbf{x}_p^* = \mathbf{x}_p^n + \Delta t \, \mathbf{u}_{\text{CVI}}\left(\mathbf{x}_p^n; \mathbf{u}^n\right),$$
  $$\mathbf{x}_p^{n+1} = \mathbf{x}_p^n + \frac{\Delta t}{2} \left[ \mathbf{u}_{\text{CVI}}\left(\mathbf{x}_p^n; \mathbf{u}^n\right) + \mathbf{u}_{\text{CVI}}\left(\mathbf{x}_p^*; \mathbf{u}^{n+1}\right) \right],$$
  est **exactement d'ordre 2 en temps et en espace** tout en utilisant à chaque sous-étape un champ de vitesse à divergence analytique nulle.
  En outre, la normale matérielle $\mathbf{n}_p$ portée par chaque marqueur peut être soit réévaluée par l'ajustement local, soit tournée cinématiquement par la formule de transport d'un élément de surface matériel :
  $$\frac{d\mathbf{n}_p}{dt} = -\left(\nabla \mathbf{u}_{\text{CVI}}^T \cdot \mathbf{n}_p\right) + \left[\mathbf{n}_p \cdot \left(\nabla \mathbf{u}_{\text{CVI}}^T \cdot \mathbf{n}_p\right)\right]\mathbf{n}_p.$$

## 3.4 Verrou n°4 : Estimation de la Courbure dans un Repère Local Tourné et Protection des Filaments Minces

- **Problème identifié** :
  1. Un ajustement cartésien $y = a x^2 + b x + c$ devient singulier lorsque l'interface est quasi-verticale ($|n_x| \to 1$).
  2. Dans un filament mince d'épaisseur $\delta < 2\Delta x$ (ex. Cas 10 *Liquid Sheet* ou Cas 7 *Vortex Shear*), le voisinage $3\times 3$ d'une cellule contient des marqueurs appartenant aux **deux interfaces opposées** du ligament ($\mathbf{n}_p \approx -\mathbf{n}_q$). Un ajustement global sur tous les marqueurs du bloc $3\times 3$ échoue alors catastrophiquement.
- **Solution mathématique complète à implémenter** :
  Pour toute cellule interfaciale $(i,j)$ :
  1. **Normale de référence VOF** : Calculer l'orientation de référence unitaire $\mathbf{n}_{i,j}^{\text{ref}} = \frac{\nabla \widetilde{C}_{i,j}}{\|\nabla \widetilde{C}_{i,j}\|}$ (pointant vers $C=1$, conformément à la convention de `compute_curvature_and_csf` dans `module_interface_vof.f90`, lignes 594–610).
  2. **Collecte et filtrage angulaire sur le voisinage $3\times 3$** : Rassembler les marqueurs $p \in \mathcal{V}_{3\times 3}(i,j)$ situés à une distance $r_p = \|\mathbf{x}_p - \mathbf{x}_{i,j}^c\| \le 1.8 \max(\Delta x, \Delta y)$ **et vérifiant la condition de cohérence directionnelle** :
     $$\mathbf{n}_p \cdot \mathbf{n}_{i,j}^{\text{ref}} \ge \cos(60^\circ) = 0.5.$$
     Ce filtre élimine instantanément les marqueurs situés sur la face opposée d'un filament mince !
  3. **Repère local tangent-normal $(\tau, \eta)$** : Si le nombre de marqueurs filtrés est $M \ge 4$, définir le centre local $\bar{\mathbf{x}}$ (moyenne pondérée des marqueurs de la cellule $(i,j)$), le vecteur normal $\mathbf{e}_\eta = \mathbf{n}_{i,j}^{\text{ref}} = (n_x^{\text{ref}}, n_y^{\text{ref}})$ et le vecteur tangent direct $\mathbf{e}_\tau = (-n_y^{\text{ref}}, n_x^{\text{ref}})$. Chaque marqueur $p$ a pour coordonnées locales :
     $$\tau_p = (\mathbf{x}_p - \bar{\mathbf{x}}) \cdot \mathbf{e}_\tau, \qquad \eta_p = (\mathbf{x}_p - \bar{\mathbf{x}}) \cdot \mathbf{e}_\eta.$$
  4. **Moindres carrés pondérés $3\times 3$** : Minimiser sur $(a, b, c) \in \mathbb{R}^3$ la fonctionnelle quadratique pondérée :
     $$\mathcal{J}(a, b, c) = \sum_{p=1}^M w_p \left( a \tau_p^2 + b \tau_p + c - \eta_p \right)^2, \qquad w_p = \exp\left( -\frac{\tau_p^2 + \eta_p^2}{(1.2 \Delta h)^2} \right),$$
     ce qui revient à inverser analytiquement un système linéaire symétrique $3\times 3$ par la formule de Cramer ou pivot de Gauss.
  5. **Extraction de la normale corrigée et de la courbure signée** :
     En $\tau = 0$, la pente tangentielle est $\eta'(0) = b$ et la dérivée seconde est $\eta''(0) = 2a$. La normale unitaire corrigée (pointant vers $C=1$) et la courbure algébrique (avec la convention $\kappa = -\nabla \cdot \mathbf{n}$ de `compute_curvature_and_csf`) valent :
     $$\mathbf{n}_{i,j}^{\text{marker}} = \frac{\mathbf{e}_\eta - b \, \mathbf{e}_\tau}{\sqrt{1 + b^2}}, \qquad \kappa_{i,j}^{\text{marker}} = \frac{2 a}{(1 + b^2)^{3/2}}.$$
     *(Vérification de signe : pour une goutte circulaire de fluide 1 ($C=1$ au centre) de rayon $R$, $\mathbf{e}_\eta = \nabla C / \|\nabla C\|$ pointe vers l'intérieur du cercle ; la courbe s'incurve vers l'intérieur ($\eta_p = +a \tau_p^2$ avec $a = +1/(2R) > 0$), donc $\kappa_{i,j}^{\text{marker}} = +1/R > 0$, et $F_{\sigma,x} = \sigma \kappa \partial_x C$ pointe vers l'intérieur dans la moitié droite où $\partial_x C < 0$, ce qui donne bien $\Delta P = P_{\text{int}} - P_{\text{ext}} = +\sigma / R > 0$ !)*
  6. **Condition de repli robuste (*Fallback*)** : Si $M < 4$ ou si le déterminant du système $3\times 3$ est mal conditionné ($\det < 10^{-14}$), ou si $|\kappa_{i,j}^{\text{marker}}| > \kappa_{\max} = \frac{2}{\min(\Delta x, \Delta y)}$ (limite de Nyquist de la grille), la cellule conserve la courbure CSF standard issue de $\widetilde{C}$.

## 3.5 Verrou n°5 : Point d'Insertion Exact dans `advect_vof` et Conventions de Signes PLIC vs CSF

- **Audit du code existant (`module_interface_vof.f90` et `InterfaceTracker.cpp`)** :
  1. Dans `module_interface_vof.f90` (lignes 669–696), la routine `advect_vof` exécute successivement : `advect_vof_weymouth_yue`, `exchange_halos_scalar`, `apply_vof_bc`, `smooth_svof_helmholtz`, `update_physical_properties`, `vof_smooth_diffusion`, **puis appelle directement `compute_curvature_and_csf` à la ligne 695**.
  2. Si l'advection des marqueurs était placée dans `sim_run` *après* `advect_vof`, alors `compute_curvature_and_csf` calculerait la force capillaire au temps $t^{n+1}$ avec les marqueurs encore situés au temps $t^n$ !
  3. **Solution** : L'étape d'advection CVI, de migration MPI et de maintenance des marqueurs doit impérativement être appelée **à l'intérieur de `advect_vof` (et `InterfaceTracker::advect`), juste après `vof_smooth_diffusion` et juste avant `compute_curvature_and_csf`**.
  4. **Attention à l'opposition de signes entre `module_interface_reconstruction.f90` et `compute_curvature_and_csf`** :
     - Dans `module_interface_reconstruction.f90` (`reconstruct_plic`, lignes 180–210), le segment PLIC $n_x^{\text{PLIC}} x + n_y^{\text{PLIC}} y + \alpha = 0$ utilise $(n_x^{\text{PLIC}}, n_y^{\text{PLIC}}) = -\frac{\nabla C}{\|\nabla C\|}$ (pointant vers $C = 0$).
     - Dans `compute_curvature_and_csf` (lignes 594–610), la normale utilise $+\frac{\nabla C}{\|\nabla C\|}$ (pointant vers $C = 1$).
     - Lorsqu'on ensemence un marqueur sur un segment PLIC $(n_x^{\text{PLIC}}, n_y^{\text{PLIC}}, \alpha)$, il faut initialiser sa normale orientée vers $C=1$ en posant $\mathbf{n}_p = (-n_x^{\text{PLIC}}, -n_y^{\text{PLIC}})$ !

## 3.6 Verrou n°6 : Migration MPI par Routage Dimensionnel en 2 Passes (4 Voisins au lieu de 8)

- **Audit du code MPI (`module_mpi_comm.f90` et `MpiComm.cpp`)** :
  La structure `mesh_2d_t` stocke uniquement les 4 voisins cardinaux (`neighbor_left`, `neighbor_right`, `neighbor_bottom`, `neighbor_top`). Un marqueur situé dans un coin de sous-domaine peut, sous une condition $\text{CFL} \le 0.45$, franchir simultanément une frontière verticale et une frontière horizontale au cours d'un pas $\Delta t$.

- **Solution HPC sans ajout de voisins diagonaux** :
  Effectuer la migration en **deux passes séquentielles** (exactement comme pour `exchange_halos_scalar`) :
  - **Passe 1 (Axe X)** : Tout marqueur avec $x_p < x_{\min,\text{local}}$ est envoyé à `neighbor_left` ; tout marqueur avec $x_p \ge x_{\max,\text{local}}$ est envoyé à `neighbor_right` (avec décalage périodique $\pm L_x$ si `bc_left == BC_PERIODIC`). Les marqueurs reçus sont immédiatement ajoutés au tableau local.
  - **Passe 2 (Axe Y)** : On balaie le tableau local mis à jour (incluant les marqueurs reçus lors de la Passe 1 !). Tout marqueur avec $y_p < y_{\min,\text{local}}$ est envoyé à `neighbor_bottom` ; tout marqueur avec $y_p \ge y_{\max,\text{local}}$ est envoyé à `neighbor_top` (avec décalage périodique $\pm L_y$ si `bc_bottom == BC_PERIODIC`).
  Ainsi, un marqueur diagonal transite par le voisin Est/Ouest en Passe 1, puis est automatiquement relayé vers le voisin Nord-Est/Sud-Ouest en Passe 2, avec **zéro modification de la topologie MPI** !

---

# 4. Algorithme Complet du Pas de Temps Hybride dans `COVOF_2D`

Au pas de temps $n \to n+1$ ($t^{n+1} = t^n + \Delta t$) :

1. **Sauvegarde de l'état antérieur** :
   $\mathbf{u}_{\text{old}} \leftarrow \mathbf{u}^n$, $\mathbf{v}_{\text{old}} \leftarrow \mathbf{v}^n$ (tous deux à divergence discrète nulle : $\mathbf{B}_u \mathbf{u}^n = \mathbf{0}$).

2. **Résolution du système monolithique Navier-Stokes (FC)** :
   - Assemblage des matrices de quantité de mouvement avec $\rho(C_s^n)$, $\mu(C_s^n)$, $\mu_{\text{node}}(C_s^n)$ et du second membre incluant la force capillaire $\mathbf{F}_\sigma^n$.
   - Résolution monolithique par BiCGStab(2) + préconditionneur `BT_SBP` jusqu'à $\|\mathbf{B}_u \mathbf{u}^{n+1}\|_\infty \le 10^{-10}$.
   - Synchronisation MPI des halos de vitesse (`sync_u`, `sync_v`) et application des conditions aux limites de vitesse (`apply_velocity_bc`).
3. **Advection VOF-PLIC conservative (Weymouth & Yue)** :
   - Reconstruction des segments PLIC $(n_x^{\text{PLIC}}, n_y^{\text{PLIC}}, \alpha)$ (guidée optionnellement par les normales des marqueurs $\mathbf{n}^{\text{marker}, n}$).
   - Advection directionnelle de Strang (`advect_vof_weymouth_yue`) avec $\mathbf{u}^{n+1}$ donnant $C^{n+1}$.
   - Échange MPI des halos de $C^{n+1}$ (`exchange_halos_scalar`) et conditions aux limites (`apply_vof_bc`).
4. **Advection Lagrangienne CVI et Maintenance des Marqueurs** *(si `use_markers == .true.`)* :
   - **Étape 4a (Advection RK2-CVI)** : Pour chaque marqueur actif $p$, calculer $\mathbf{u}_{\text{CVI}}(\mathbf{x}_p^n; \mathbf{u}^n)$, prédire $\mathbf{x}_p^* = \mathbf{x}_p^n + \Delta t \, \mathbf{u}_{\text{CVI}}(\mathbf{x}_p^n; \mathbf{u}^n)$, appliquer les conditions aux limites géométriques sur $\mathbf{x}_p^*$, puis corriger :
     $$\mathbf{x}_p^{n+1} = \mathbf{x}_p^n + \frac{\Delta t}{2}\left[ \mathbf{u}_{\text{CVI}}(\mathbf{x}_p^n; \mathbf{u}^n) + \mathbf{u}_{\text{CVI}}(\mathbf{x}_p^*; \mathbf{u}^{n+1}) \right].$$
   - **Étape 4b (Conditions aux limites physiques et Migration MPI X-Y)** :
     - Réflexion symétrique sur parois `FREE_SLIP` / `SYMMETRY`, suppression sur bord `OUTFLOW`, injection de nouveaux marqueurs sur bord `JET_INLET` aux abscisses $x_{\text{jet}} \pm D/2$.
     - Exécution de `migrate_markers_mpi` en 2 passes (Axe X puis Axe Y).
   - **Étape 4c (Élagage / *Pruning*)** :
     - Mettre à jour les indices de cellule $(i_p, j_p) = \left(\lfloor (x_p - x_{\min})/\Delta x \rfloor + 1, \lfloor (y_p - y_{\min})/\Delta y \rfloor + 1\right)$.
     - Supprimer tout marqueur situé dans une cellule pure ($C_{i,j}^{n+1} < 10^{-6}$ ou $C_{i,j}^{n+1} > 1 - 10^{-6}$) dont les 8 voisines sont également pures, ou si la densité locale dépasse $N_{\max} = 8$ marqueurs par cellule (en conservant les marqueurs les mieux espacés).
   - **Étape 4d (Ensemencement sélectif / *Seeding*)** :
     - Pour toute cellule interfaciale ($10^{-4} \le C_{i,j}^{n+1} \le 1 - 10^{-4}$) ayant $N_p = 0$ marqueur, reconstruire le segment PLIC local $(n_x, n_y, \alpha)$ de $C_{i,j}^{n+1}$, calculer ses deux extrémités d'intersection $\mathbf{P}_1, \mathbf{P}_2$ avec le bord de la maille $(i,j)$, et insérer $N_{\text{target}}$ marqueurs aux abscisses curvilignes $s_k = \frac{k - 0.5}{N_{\text{target}}} \in (0, 1)$ ($k = 1, \dots, N_{\text{target}}$) avec normale $\mathbf{n}_p = (-n_x, -n_y)$.
5. **Régularisations et Calcul de la Force Capillaire Hybride** :
   - Calcul de $C_s^{n+1}$ (`smooth_svof_helmholtz`), mise à jour de $\rho, \mu, \mu_{\text{node}}$ (`update_physical_properties`), et calcul de $\widetilde{C}^{n+1}$ (`vof_smooth_diffusion`).
   - Dans `compute_curvature_and_csf` : si `marker_curvature == .true.`, évaluer $\kappa_{i,j}$ par ajustement parabolique dans le repère local $(\tau, \eta)$ avec filtrage angulaire sur les marqueurs du voisinage $3\times 3$ (Section 3.4), puis calculer les forces capillaires *Balanced-Force* aux faces :
     $$f_{\sigma, x, i,j} = \sigma \, \frac{\kappa_{i-1,j} + \kappa_{i,j}}{2} \, \frac{C_{i,j}^{n+1} - C_{i-1,j}^{n+1}}{\Delta x}, \qquad f_{\sigma, y, i,j} = \sigma \, \frac{\kappa_{i,j-1} + \kappa_{i,j}}{2} \, \frac{C_{i,j}^{n+1} - C_{i,j-1}^{n+1}}{\Delta y}.$$

---

# 5. Synthèse des Cas de Validation de `COVOF_2D` Améliorés par la Méthode

Le déploiement de cette méthode dans `COVOF_2D` cible huit campagnes majeures répertoriées dans `Validation_Totale/Run_CASES/` :

| Réf. Cas dans `Run_CASES` | Nom du Cas Test | Phénomène Physique Dominant | Verrou Actuel du VOF Seul | Gain Apporté par `FC + WY + CVI + Marqueurs` |
| :--- | :--- | :--- | :--- | :--- |
| **Cas 10** | `validation_cas_10_liquid_sheet` | Cisaillement Kelvin-Helmholtz d'une nappe mince ($r_\rho = 1000$) | Rupture numérique du ligament dès que $\delta < \Delta x$ | Suivi continu sous-maille des deux interfaces de la nappe + conservation d'énergie cinétique FC |
| **Cas 21** | `validation_cas_21_ascension_bulle_2d` | Ascension de bulle 2D ($Bo=10$ ellipsoïdal et $Bo=125$ à jupes) | Fragmentation des jupes fines à $Bo=125$ et bruit sur $\kappa$ | Courbure précise au nez de la bulle et maintien de l'intégrité des jupes latérales |
| **Cas 24** | `validation_cas_24_viscous_jet_buckling` | Flambement d'un jet visqueux ($\mu_2/\mu_1 = 2.7\times 10^5$, $Re=0.17$) | Érosion des replis serrés de l'accordéon visqueux | Normales sous-maille fidèles dans les lobes de flambement + injection continue `JET_INLET` |
| **Cas 25** | `validation_cas_25_laplace_droplet_cvi_markers` | Équilibre capillaire de Young-Laplace ($\Delta P = \sigma / R_0$) | Courants parasites $\|\mathbf{u}\|_\infty$ dus aux différences finies sur $\widetilde{C}$ | Convergence d'ordre 2 de $\kappa^{\text{marker}}$ et chute de 1 à 2 ordres de grandeur des courants parasites |
| **Cas 26** | `validation_cas_26_oscillation_bulle_lamb_cvi_markers` | Oscillations capillaires amorties de Lamb ($r_\rho = 5950$) | Amortissement numérique excessif et décalage de fréquence $f_2$ | Restitution fidèle de la fréquence de Lamb et du décrément logarithmique visqueux |
| **Cas 27** | `validation_cas_27_vortex_shear_cvi_markers` | Tourbillon cisaillé réversible de Rider-Kothe ($T = 6\,\text{s}$) | Fragmentation de la queue spiralée à $t = T/2 = 3\,\text{s}$ | Réversibilité lagrangienne CVI ($\nabla\cdot\mathbf{u}_{\text{CVI}}=0$) et réduction de l'erreur géométrique $L_1(T)$ |
| **Cas 28** | `validation_cas_28_dense_cylinder_cvi_markers` | Chute libre d'un cylindre dense (Cas 3.1 *Fluids* 2021) | Facettage angulaire du cylindre sur maillage modéré | Maintien de la circularité géométrique et conservation exacte de l'aire ($\Delta V/V_0 < 10^{-12}$) |
| **Cas 29** | `validation_cas_29_rayleigh_taylor_cvi_markers` | Instabilité gravitaire de Rayleigh-Taylor diphasique | Diffusion des enroulements secondaires du panache | Résolution fine des volutes de Kelvin-Helmholtz le long des flancs du champignon |

---

# 6. Références Bibliographiques

1. **El Ouafa, S., Vincent, S., & Le Chenadec, V. (2021).** *Monolithic Solvers for Incompressible Two-Phase Flows at Large Density and Viscosity Ratios*. Fluids, 6(1), 23. [https://doi.org/10.3390/fluids6010023](https://doi.org/10.3390/fluids6010023)
2. **Weymouth, G. D., & Yue, D. K.-P. (2010).** *Conservative Volume-of-Fluid method for free-surface simulations on Cartesian-grids*. Journal of Computational Physics, 229(8), 2853–2865. [https://doi.org/10.1016/j.jcp.2009.12.018](https://doi.org/10.1016/j.jcp.2009.12.018)
3. **Wang, H., Agrusta, R., & van Hunen, J. (2015).** *Advantages of a conservative velocity interpolation (CVI) scheme for particle-in-cell methods in geodynamic modeling*. Geochemistry, Geophysics, Geosystems, 16(6), 2015–2023. [https://pmc.ncbi.nlm.nih.gov/articles/PMC5089062/](https://pmc.ncbi.nlm.nih.gov/articles/PMC5089062/)
4. **Pusok, A. E., Kaus, B. J. P., & Popov, A. A. (2017).** *On the Quality of Velocity Interpolation Schemes for Marker-in-Cell Method and Staggered Grids*. Pure and Applied Geophysics, 174, 1071–1089. [https://www.researchgate.net/publication/311165935](https://www.researchgate.net/publication/311165935_On_the_Quality_of_Velocity_Interpolation_Schemes_for_Marker-in-Cell_Method_and_Staggered_Grids)
5. **Gerya, T. (2019).** *Introduction to Numerical Geodynamic Modelling* (Chapitre 8 : *Advection Equation and Marker-in-Cell Method*). Cambridge University Press.
6. **Arrufat, T., Crialesi-Esposito, M., Fuster, D., Ling, Y., Malan, L., Scardovelli, R., Tryggvason, G., & Zaleski, S. (2021).** *A mass-momentum consistent, Volume-of-Fluid method for incompressible flow on staggered grids*. Computers & Fluids, 215, 104785. [https://arxiv.org/html/1811.12327v2](https://arxiv.org/html/1811.12327v2)
7. **Karnakov, P., Litvinov, S., & Koumoutsakos, P. (2020).** *A hybrid particle volume-of-fluid method for curvature estimation in multiphase flows*. International Journal of Multiphase Flow, 125, 103209. [https://arxiv.org/html/1906.00314v2](https://arxiv.org/html/1906.00314v2)
8. **Aulisa, E., Manservisi, S., & Scardovelli, R. (2003, 2025).** *A mixed markers and volume-of-fluid method for the reconstruction and advection of interfaces in two-phase and free-boundary flows*. Journal of Computational Physics, 188(2), 611–639 ; & AIP Conf. Proc. 3315, 040022 (2025).
9. **Chirco, L., & Zaleski, S. (2022, 2025).** *An Edge-Based Interface Tracking (EBIT) method for multiphase-flow simulations*. Int. J. Numer. Meth. Fluids, 94, 1468–1492 ; & *EBIT 3D*, [arXiv:2509.26361](https://arxiv.org/html/2509.26361v2).
10. **Gorges, C., Evrard, F., van Wachem, B., & Denner, F. (2022).** *Divergence-preserving velocity interpolation in Front-Tracking*. [https://arxiv.org/abs/2202.10223](https://arxiv.org/abs/2202.10223)
11. **Evrard, F., Chirco, L., Zaleski, S., & van Wachem, B. (2025).** *Sharp Front-Tracking with geometric interface reconstruction*. [https://arxiv.org/html/2505.06932v1](https://arxiv.org/html/2505.06932v1)
