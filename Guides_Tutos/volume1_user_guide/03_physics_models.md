# Formulation Mathématique et Modèles Physiques

## Système d'Équations Gouvernantes

Le solveur **COVOF 2D** résout les équations fondamentales de la mécanique des milieux continus régissant les écoulements de fluides visqueux, newtoniens, monophasiques ou diphasiques, sous diverses approximations physiques (incompressible, faiblement compressible à bas nombre de Mach, avec ou sans transfert thermique, et en présence d'obstacles solides ou poreux).

### Équations de Navier-Stokes sous Forme Continue

Dans un domaine spatial ouvert et borné $\Omega \subset \mathbb{R}^2$ de frontière $\partial \Omega$ et pour un intervalle de temps $t \in [0, t_{\text{end}}]$, les équations de conservation de la masse et de la quantité de mouvement s'écrivent :

$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0$$

$$\frac{\partial (\rho \mathbf{u})}{\partial t} + \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u}) = -\nabla p + \nabla \cdot \boldsymbol{\tau} + \rho \mathbf{g} + \mathbf{f}_\sigma + \mathbf{f}_{\text{brinkman}} + \mathbf{f}_{\text{boussinesq}}$$

où :
- $\rho = \rho(\mathbf{x}, t)$ est la masse volumique locale du milieu $[\text{kg}\cdot\text{m}^{-3}]$.
- $\mathbf{u} = (u, v)^T$ est le vecteur vitesse eulérien $[\text{m}\cdot\text{s}^{-1}]$.
- $p = p(\mathbf{x}, t)$ est le champ de pression hydrodynamique $[\text{Pa}]$.
- $\boldsymbol{\tau}$ est le tenseur des contraintes visqueuses de déviation $[\text{Pa}]$.
- $\mathbf{g} = (g_x, g_y)^T$ est le vecteur accélération de la pesanteur $[\text{m}\cdot\text{s}^{-2}]$.
- $\mathbf{f}_\sigma$ est la force volumique de tension superficielle concentrée à l'interface $[\text{N}\cdot\text{m}^{-3}]$.
- $\mathbf{f}_{\text{brinkman}}$ est la force de traînée volumique de pénalisation poreuse ou solide de Brinkman $[\text{N}\cdot\text{m}^{-3}]$.
- $\mathbf{f}_{\text{boussinesq}}$ est la force volumique de flottabilité thermique sous approximation de Boussinesq $[\text{N}\cdot\text{m}^{-3}]$.

### Formulation Advective Non Conservative

En exploitant l'équation de continuité, la dérivée matérielle de la vitesse s'exprime sous forme advective (forme utilisée pour la discrétisation par volumes finis décalés dans COVOF 2D) :

$$\rho \left( \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} \right) = -\nabla p + \nabla \cdot \boldsymbol{\tau} + \rho \mathbf{g} + \mathbf{f}_\sigma + \mathbf{f}_{\text{brinkman}} + \mathbf{f}_{\text{boussinesq}}$$

Soit, projetée selon les composantes cartésiennes horizontales $x$ et verticales $y$ :

$$\rho \left( \frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} \right) = -\frac{\partial p}{\partial x} + \frac{\partial \tau_{xx}}{\partial x} + \frac{\partial \tau_{xy}}{\partial y} + \rho g_x + f_{\sigma, x} + f_{\text{brinkman}, x}$$

$$\rho \left( \frac{\partial v}{\partial t} + u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} \right) = -\frac{\partial p}{\partial y} + \frac{\partial \tau_{yx}}{\partial x} + \frac{\partial \tau_{yy}}{\partial y} + \rho g_y + f_{\sigma, y} + f_{\text{brinkman}, y} + f_{\text{boussinesq}, y}$$

---

## Modèle 1-Fluide Diphasique et Suivi d'Interface VOF

Dans la formulation 1-fluide de Kataoka (1986), le milieu diphasique est modélisé comme un fluide unique continu équivalent dont les propriétés thermodynamiques et de transport varient de manière discontinue à la traversée de l'interface séparant la phase 1 (fluide porteur ou liquide) de la phase 2 (goutte, bulle ou gaz).

### Fonction Indicatrice de Couleur (Fraction Volumique VOF)

La distribution spatiale des deux phases est représentée par la fonction indicatrice de couleur $C(\mathbf{x}, t) \in [0, 1]$ :

$$C(\mathbf{x}, t) = \begin{cases} 1 & \text{si } \mathbf{x} \in \text{Phase 1} \\ 0 & \text{si } \mathbf{x} \in \text{Phase 2} \\ ]0, 1[ & \text{dans les mailles traversées par l'interface} \end{cases}$$

L'équation cinématique de transport de l'interface sans changement de phase s'énonce :

$$\frac{\partial C}{\partial t} + \mathbf{u} \cdot \nabla C = 0$$

Sous forme conservative pour un écoulement à divergence nulle :

$$\frac{\partial C}{\partial t} + \nabla \cdot (C \mathbf{u}) = 0$$

### Lois de Mélange des Propriétés Matérielles

Dans chaque cellule de calcul, la masse volumique effective $\rho$ et la viscosité dynamique $\mu$ sont reconstruites à partir de la fraction volumique locale $C$ :

1. **Masse volumique cellulaire (Moyenne arithmétique) :**
   $$\rho(C) = C \rho_1 + (1 - C) \rho_2$$
   où $\rho_1$ est la masse volumique de la phase 1 (`cfg%rho1` / `config.rho1`) et $\rho_2$ celle de la phase 2 (`cfg%rho2` / `config.rho2`).

2. **Viscosité dynamique cellulaire et nodale :**
   Dans le centre de la maille :
   $$\mu(C) = C \mu_1 + (1 - C) \mu_2$$
   Pour le calcul des termes de cisaillement aux sommets du maillage décalé MAC ($\mu_{xy}$ aux nœuds $(i+\frac{1}{2}, j+\frac{1}{2})$), COVOF 2D utilise la **moyenne harmonique** des quatre mailles adjacentes, qui est la seule formulation rigoureusement conservative et stable pour les rapports de viscosité élevés ($\mu_1/\mu_2 \gg 1$) :
   $$\mu_{i+1/2, j+1/2} = \left[ \frac{1}{4} \left( \frac{1}{\mu_{i, j}} + \frac{1}{\mu_{i+1, j}} + \frac{1}{\mu_{i, j+1}} + \frac{1}{\mu_{i+1, j+1}} \right) \right]^{-1}$$

### Modèle de Tension Superficielle CSF (Continuum Surface Force)

La tension superficielle agissant sur l'interface est convertie en un terme source volumique équivalent selon la méthode CSF de Brackbill, Kothe et Zemach (1992) :

$$\mathbf{f}_\sigma = \sigma \kappa \nabla C$$

où :
- $\sigma$ est le coefficient de tension superficielle interfaciale $[\text{N}\cdot\text{m}^{-1}]$ (`cfg%sigma` / `config.sigma`).
- $\kappa$ est la courbure moyenne locale de l'interface $[\text{m}^{-1}]$.
- $\nabla C$ est le gradient spatial de la fraction volumique, orienté selon la normale unitaire à l'interface.

La normale locale unitaire $\mathbf{n}$ et la courbure $\kappa$ sont calculées par dérivation numérique du champ de fraction volumique régularisé :

$$\mathbf{n} = \frac{\nabla C}{\|\nabla C\| + \epsilon}, \quad \kappa = -\nabla \cdot \mathbf{n} = -\left( \frac{\partial n_x}{\partial x} + \frac{\partial n_y}{\partial y} \right)$$

où $\epsilon = 10^{-30}$ est un seuil de régularisation infinitésimal évitant les divisions par zéro dans les phases pures uniformes.

---

## Modèle Faiblement Compressible et Régime All-Mach

Pour les écoulements diphasiques présentant une compressibilité acoustique ou thermodynamique (notamment la compression d'inclusions gazeuses ou la propagation d'ondes de pression à bas nombre de Mach), COVOF 2D intègre la formulation **All-Mach C-AAL** (El Ouafa et al., *Fluids* 2024).

### Équation de Continuité Généralisée

L'équation de conservation de la masse est réécrite en introduisant le coefficient de compressibilité isotherme $\chi_T$ :

$$\chi_T = \frac{1}{\rho}\left( \frac{\partial \rho}{\partial p} \right)_T = \frac{1}{\rho c_s^2}$$

où $c_s$ désigne la célérité du son dans le fluide considéré $[\text{m}\cdot\text{s}^{-1}]$. L'équation de continuité prend alors la forme :

$$\chi_T \left( \frac{\partial p}{\partial t} + \mathbf{u} \cdot \nabla p \right) + \nabla \cdot \mathbf{u} = 0$$

À bas nombre de Mach ($Ma = \|\mathbf{u}\|/c_s \ll 1$), le terme convectif $\mathbf{u} \cdot \nabla p$ est d'ordre $\mathcal{O}(Ma^2)$ devant la dérivée instationnaire et est négligé, ce qui conduit à la relation de contrainte barotrope discrète :

$$\chi_T \frac{\partial p}{\partial t} + \nabla \cdot \mathbf{u} = 0 \quad \iff \quad \nabla \cdot \mathbf{u}^{n+1} + \frac{\chi_T}{\Delta t} (p^{n+1} - p^n) = 0$$

Cette relation modifie directement l'étape de mise à jour d'Uzawa dans la méthode du Lagrangien augmenté adaptatif compressible (C-AAL), régularisant le paramètre de pénalité local :

$$r_{\text{eff}}(\mathbf{x}) = \frac{r_{\text{AAL}}(\mathbf{x})}{1 + r_{\text{AAL}}(\mathbf{x}) \frac{\chi_T(\mathbf{x})}{\Delta t}}$$

---

## Modèle de Milieu Poreux et Pénalisation de Brinkman

Pour représenter des obstacles rigides immergés (pylônes, cylindres, marches, vannes) sans faire appel à un remaillage non structuré conforme, COVOF 2D implémente la méthode de pénalisation volumique de Brinkman (Angot, Bruneau & Fabrie, 1999).

### Équation de Brinkman-Navier-Stokes

Un terme de traînée de Darcy-Brinkman est ajouté au membre de droite de l'équation de quantité de mouvement :

$$\mathbf{f}_{\text{brinkman}} = -\frac{\mu}{\mathcal{K}(\mathbf{x})} \mathbf{u}$$

où $\mathcal{K}(\mathbf{x})$ est le champ spatial de perméabilité intrinsèque de Darcy $[\text{m}^2]$ :
- **Dans le fluide libre ($\mathbf{x} \in \Omega_f$) :** $\mathcal{K} \to +\infty \implies \mathbf{f}_{\text{brinkman}} = \mathbf{0}$, les équations de Navier-Stokes standard sont recouvrées exactement.
- **Dans un milieu poreux perméable ($\mathbf{x} \in \Omega_p$) :** $\mathcal{K} = K_{\text{poreux}} \in [10^{-6}, 10^{-2}]\,\text{m}^2$, la loi de Darcy avec cisaillement de Brinkman régit la filtration.
- **Dans un obstacle solide rigide imperméable ($\mathbf{x} \in \Omega_s$) :** $\mathcal{K} = \eta \to 0$ ($\eta \le 10^{-12}\,\text{m}^2$, paramètre `obs_perm` / `PERMEABILITY_SOLID_DEFAULT`). Par pénalisation asymptotique, l'équation impose mathématiquement :
  $$\mathbf{u} = \mathbf{0} + \mathcal{O}(\sqrt{\eta}) \quad \text{dans } \Omega_s$$
  garantissant le respect de la condition d'adhérence sans coût géométrique supplémentaire.

---

## Transfert Thermique et Approximation de Boussinesq

Pour les écoulements thermoconvectifs (cavité différentiellement chauffée, instabilités de Rayleigh-Bénard), le solveur couple l'équation de quantité de mouvement à une équation de transport-diffusion de l'énergie thermique.

### Équation de l'Énergie en Température

En négligeant la dissipation visqueuse et les variations de pression :

$$\rho C_p \left( \frac{\partial T}{\partial t} + \mathbf{u} \cdot \nabla T \right) = \nabla \cdot (k_{\text{th}} \nabla T) + S_T$$

où :
- $T(\mathbf{x}, t)$ est le champ de température absolue $[\text{K}]$ (`fields.T` / `fields%T`).
- $C_p$ est la capacité thermique isobare $[\text{J}\cdot\text{kg}^{-1}\cdot\text{K}^{-1}]$ (`cfg%cp1` / `config.cp1`).
- $k_{\text{th}}$ est la conductivité thermique du fluide $[\text{W}\cdot\text{m}^{-1}\cdot\text{K}^{-1}]$ (`cfg%th_cond1` / `config.th_cond1`).
- $S_T$ est un terme source volumique de chaleur $[\text{W}\cdot\text{m}^{-3}]$.

### Couplage de Flottabilité de Boussinesq

La variation de la masse volumique avec la température est linéarisée autour d'un état de référence $(T_{\text{ref}}, \rho_0)$ :

$$\rho(T) \approx \rho_0 \left( 1 - \beta (T - T_{\text{ref}}) \right)$$

où $\beta$ est le coefficient de dilatation volumique thermique $[\text{K}^{-1}]$ (`cfg%beta_th` / `config.beta_th`). Le terme source de gravité prend alors la forme :

$$\mathbf{f}_{\text{boussinesq}} = -\rho_0 \beta (T - T_{\text{ref}}) \mathbf{g}$$

Ce terme est orienté verticalement selon l'axe de la pesanteur, induisant la mise en mouvement naturelle du fluide chaud plus léger vers le haut et du fluide froid plus dense vers le bas.

---

## Nombres Adimensionnels et Régimes d'Écoulement

La physique des écoulements simulés dans COVOF 2D est caractérisée par les nombres sans dimension classiques, définis à partir d'une longueur de référence $L$, d'une vitesse de référence $U$, d'un écart de température $\Delta T = T_{\text{hot}} - T_{\text{cold}}$ et de la viscosité cinématique $\nu = \mu/\rho$ :

| Nombre Sans Dimension | Expression Mathématique | Signification Physique | Modules COVOF 2D Concernés |
| :--- | :--- | :--- | :--- |
| **Reynolds ($Re$)** | $Re = \frac{\rho U L}{\mu} = \frac{U L}{\nu}$ | Rapport forces d'inertie / forces visqueuses | Tous (Cavité, Poiseuille, Cylindre) |
| **Froude ($Fr$)** | $Fr = \frac{U}{\sqrt{g L}}$ | Rapport forces d'inertie / forces gravitationnelles | Surface libre, Rupture de barrage |
| **Weber ($We$)** | $We = \frac{\rho U^2 L}{\sigma}$ | Rapport forces d'inertie / forces capillaires | Rupture de nappe, Goutte, Jet |
| **Capillaire ($Ca$)** | $Ca = \frac{\mu U}{\sigma} = \frac{We}{Re}$ | Rapport forces visqueuses / forces capillaires | Microfluidique, Instabilité de jet |
| **Bond / Eötvös ($Bo$)** | $Bo = \frac{\Delta \rho g L^2}{\sigma}$ | Rapport forces gravitationnelles / capillaires | Ascension de bulle, Sédimentation |
| **Ohnesorge ($Oh$)** | $Oh = \frac{\mu}{\sqrt{\rho \sigma L}} = \frac{\sqrt{We}}{Re}$ | Dissipation visqueuse relative aux ondes capillaires | Oscillation de bulle, Atomisation |
| **Prandtl ($Pr$)** | $Pr = \frac{\nu}{\alpha} = \frac{\mu C_p}{k_{\text{th}}}$ | Rapport diffusivité visqueuse / diffusivité thermique | Rayleigh-Bénard, Convection naturelle |
| **Rayleigh ($Ra$)** | $Ra = \frac{g \beta \Delta T L^3}{\nu \alpha} = Gr \cdot Pr$ | Intensité motrice de la convection thermique | Cavité différentiellement chauffée |
| **Darcy ($Da$)** | $Da = \frac{K}{L^2}$ | Perméabilité adimensionnelle du milieu poreux | Milieux poreux, Pénalisation Brinkman |

---

## Table Générale de Nomenclature

Le tableau ci-dessous recense l'intégralité des grandeurs scalaires, vectorielles et tensorielles utilisées dans la formulation continue et leur dénomination exacte dans les codes sources C++ et Fortran :

| Symbole | Signification Physique | Unités SI | Variable C++ (`COVOF_2D-C++`) | Variable Fortran (`COVOF_2D-FORTRAN`) |
| :--- | :--- | :--- | :--- | :--- |
| $\mathbf{u}$ | Vecteur vitesse eulérien | $\text{m}\cdot\text{s}^{-1}$ | `fields.u`, `fields.v` | `fields%u`, `fields%v` |
| $p$ | Champ de pression scalaire | $\text{Pa} = \text{N}\cdot\text{m}^{-2}$ | `fields.p` | `fields%p` |
| $p^n$ | Pression au pas de temps antérieur | $\text{Pa}$ | `fields.p_old` | `fields%p_old` |
| $\rho$ | Masse volumique du fluide | $\text{kg}\cdot\text{m}^{-3}$ | `fields.rho` | `fields%rho` |
| $\mu$ | Viscosité dynamique | $\text{Pa}\cdot\text{s} = \text{kg}\cdot\text{m}^{-1}\cdot\text{s}^{-1}$ | `fields.mu` | `fields%mu` |
| $C$ | Fraction volumique de couleur (VOF) | Sans dimension | `fields.c` | `fields%c` |
| $\mathbf{n}$ | Vecteur normal unitaire à l'interface | Sans dimension | `fields.nx`, `fields.ny` | `fields%nx`, `fields%ny` |
| $\kappa$ | Courbure locale de l'interface | $\text{m}^{-1}$ | `fields.kappa` | `fields%kappa` |
| $\sigma$ | Tension superficielle | $\text{N}\cdot\text{m}^{-1}$ | `config.sigma` | `cfg%sigma` |
| $\mathbf{f}_\sigma$ | Force de tension superficielle CSF | $\text{N}\cdot\text{m}^{-3}$ | `fields.fx_surf`, `fields.fy_surf` | `fields%fx_surf`, `fields%fy_surf` |
| $\chi_T$ | Compressibilité isotherme | $\text{Pa}^{-1} = \text{m}\cdot\text{s}^2\cdot\text{kg}^{-1}$ | `fields.chi_T` | `fields%chi_T` |
| $r_{\text{al}}$ | Paramètre de pénalité Lagrangien | $\text{Pa}\cdot\text{s}\cdot\text{m}^{-1}$ | `config.r_al`, `fields.r_aal` | `cfg%r_al`, `fields%r_aal` |
| $K$ | Perméabilité de Darcy | $\text{m}^2$ | `config.obs_perm` | `cfg%obs_perm` |
| $T$ | Température absolue | $\text{K}$ | `fields.T` | `fields%T` |
| $k_{\text{th}}$ | Conductivité thermique | $\text{W}\cdot\text{m}^{-1}\cdot\text{K}^{-1}$ | `fields.k_therm` | `fields%k_therm` |
| $C_p$ | Capacité calorifique isobare | $\text{J}\cdot\text{kg}^{-1}\cdot\text{K}^{-1}$ | `fields.cp_therm` | `fields%cp_therm` |
| $\beta$ | Dilatation thermique volumique | $\text{K}^{-1}$ | `config.beta_th` | `cfg%beta_th` |
| $\mathbf{g}$ | Accélération de la pesanteur | $\text{m}\cdot\text{s}^{-2}$ | `config.gx`, `config.gy` | `cfg%gx`, `cfg%gy` |
| $\Delta t$ | Pas de temps d'intégration | $\text{s}$ | `dt` | `dt` |
| $\Delta x, \Delta y$ | Pas d'espace de la grille MAC | $\text{m}$ | `mesh.dx`, `mesh.dy` | `mesh%dx`, `mesh%dy` |
| $N_x, N_y$ | Nombre de mailles de calcul locales | Entier | `mesh.nx_local`, `mesh.ny_local` | `mesh%nx_local`, `mesh%ny_local` |
