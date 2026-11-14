# Méthodes Numériques

## Discrétisation Spatiale sur Grille Cartésienne Décalée (MAC)

Le solveur **COVOF 2D** utilise la discrétisation spatiale sur grille décalée dite *Marker-and-Cell* (MAC), introduite par Harlow et Welch (1965). Ce choix topologique fondamental est universellement reconnu en mécanique des fluides incompressibles pour ses propriétés mathématiques intrinsèques d'inf-sup stabilité (condition de Ladyzhenskaya-Babuška-Brezzi ou LBB discrète) et l'élimination naturelle des oscillations parasites en damier (*checkerboard oscillations*) de pression, sans requérir de dissipation artificielle ni de correction de Rhie-Chow.

### Positionnement Staggered des Variables Primordiales

Le domaine de calcul bidimensionnel est pavé par un maillage régulier de pas $\Delta x$ et $\Delta y$. Les différentes grandeurs physiques sont localisées comme suit :

1. **Centres des cellules $(i, j)$ :**
   - Pression hydrodynamique $p_{i, j}$
   - Fraction volumique VOF $C_{i, j}$
   - Masse volumique $\rho_{i, j}$ et viscosité dynamique cellulaire $\mu_{i, j}$
   - Température scalaire $T_{i, j}$
   - Paramètre adaptatif de pénalité de Lagrangien augmenté $r_{i, j}$
   - Compressibilité isotherme $\chi_{T, i, j}$

2. **Faces verticales (Normales à $x$) $(i + \frac{1}{2}, j)$ :**
   - Composante horizontale de la vitesse $u_{i+1/2, j}$
   - Masse volumique transportée aux faces $\rho_{u, i+1/2, j}$
   - Force volumique superficielle horizontale $f_{\sigma, x, i+1/2, j}$

3. **Faces horizontales (Normales à $y$) $(i, j + \frac{1}{2})$ :**
   - Composante verticale de la vitesse $v_{i, j+1/2}$
   - Masse volumique transportée aux faces $\rho_{v, i, j+1/2}$
   - Force volumique superficielle verticale $f_{\sigma, y, i, j+1/2}$
   - Terme de flottabilité thermique de Boussinesq $f_{\text{boussinesq}, y, i, j+1/2}$

4. **Sommets / Nœuds de coin $(i + \frac{1}{2}, j + \frac{1}{2})$ :**
   - Viscosité dynamique harmonique nodale $\mu_{i+1/2, j+1/2}$
   - Taux de déformation croisé $\left( \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \right)$
   - Contrainte de cisaillement visqueuse $\tau_{xy, i+1/2, j+1/2}$

```
               (i, j+1)
                  |
        +---------v---------+
        |                   |
        |      p(i, j)      |
   u    |      C(i, j)      |    u
(i-1/2,j)----> T(i, j)     ----> (i+1/2, j)
        |     rho(i, j)     |
        |                   |
        +---------^---------+
                  |
                (i, j)
              v(i, j-1/2)
```

### Gradients Discrets et Conservation de la Divergence

Sur la grille MAC, l'opérateur de divergence discrète de la vitesse $(\nabla \cdot \mathbf{u})_{i, j}$ au centre de la cellule $(i, j)$ s'écrit de manière compacte et exacte au second ordre :

$$(\nabla \cdot \mathbf{u})_{i, j} = \frac{u_{i+1/2, j} - u_{i-1/2, j}}{\Delta x} + \frac{v_{i, j+1/2} - v_{i, j-1/2}}{\Delta y}$$

Réciproquement, les composantes du gradient de pression agissant respectivement sur les équations de quantité de mouvement en $u$ et en $v$ s'évaluent sans aucune interpolation :

$$\left( \frac{\partial p}{\partial x} \right)_{i+1/2, j} = \frac{p_{i+1, j} - p_{i, j}}{\Delta x}, \quad \left( \frac{\partial p}{\partial y} \right)_{i, j+1/2} = \frac{p_{i, j+1} - p_{i, j}}{\Delta y}$$

Cette dualité adjointe parfaite assure que le gradient discret $B^T$ est l'opposé exact de l'adjoint de l'opérateur de divergence discrète $B$, conférant au système de Stokes une structure point-selle symétrique parfaite.

---

## Schémas d'Advection et Linéarisations Convectives

Le terme d'advection non linéaire $(\mathbf{u} \cdot \nabla) \mathbf{u}$ dans les équations de quantité de mouvement représente le couplage cinématique d'inertie. COVOF 2D met à disposition quatre formulations sélectionnables par la clé `convective_linearization` dans le fichier JSON.

### Schéma Décentré Amont d'Ordre 1 (`UPWIND`)

Le schéma décentré amont utilise le sens de propagation de l'écoulement pour évaluer la vitesse sur les faces de contrôle :

$$\left( u \frac{\partial u}{\partial x} \right)_{i+1/2, j} \approx \begin{cases} u_{i+1/2, j} \frac{u_{i+1/2, j} - u_{i-1/2, j}}{\Delta x} & \text{si } u_{i+1/2, j} \ge 0 \\ u_{i+1/2, j} \frac{u_{i+3/2, j} - u_{i+1/2, j}}{\Delta x} & \text{si } u_{i+1/2, j} < 0 \end{cases}$$

- **Avantages :** Stabilité inconditionnelle, positivité stricte, absence d'oscillations numériques parasites même à très haut nombre de Reynolds de maille $Re_{\Delta} = \frac{\rho U \Delta x}{\mu} > 2$.
- **Inconvénients :** Diffusion numérique artificielle d'ordre 1 proportionnelle à $\frac{u \Delta x}{2}$, tronquant les tourbillons fins.

### Schéma Centré d'Ordre 2 (`CENTERED`)

Le schéma centré calcule la moyenne arithmétique de la vitesse aux interfaces de chaque volume de contrôle :

$$\left( u \frac{\partial u}{\partial x} \right)_{i+1/2, j} \approx \frac{1}{\Delta x} \left[ \left(\frac{u_{i+1/2, j} + u_{i+3/2, j}}{2}\right)^2 - \left(\frac{u_{i-1/2, j} + u_{i+1/2, j}}{2}\right)^2 \right]$$

- **Avantages :** Précision spatiale rigoureusement d'ordre 2 $\mathcal{O}(\Delta x^2, \Delta y^2)$, aucune diffusion numérique artificielle, conservation stricte de l'énergie cinétique globale.
- **Précautions :** Peut générer des oscillations de dispersion lorsque $Re_{\Delta} > 2$, nécessitant un maillage suffisamment résolu.

### Linéarisation Complète de Newton (`NEWTON`)

Pour les simulations monolithiques stationnaires ou à grands pas de temps, COVOF 2D dispose d'une linéarisation complète de Newton-Raphson de l'opérateur convectif :

$$\mathbf{u}^{n+1} = \mathbf{u}^* + \delta \mathbf{u}$$

$$(\mathbf{u}^{n+1} \cdot \nabla) \mathbf{u}^{n+1} \approx (\mathbf{u}^* \cdot \nabla) \delta \mathbf{u} + (\delta \mathbf{u} \cdot \nabla) \mathbf{u}^* + (\mathbf{u}^* \cdot \nabla) \mathbf{u}^*$$

Cette formulation introduit les dérivées croisées $\frac{\partial u}{\partial x}$, $\frac{\partial u}{\partial y}$, $\frac{\partial v}{\partial x}$, $\frac{\partial v}{\partial y}$ dans la matrice jacobienne du système linéaire, garantissant un taux de convergence asymptotique **quadratique** indispensable pour la convergence rapide des solveurs de Krylov couplés.

### Schéma d'Adams-Bashforth d'Ordre 2 (`ADAMS_BASHFORTH` / `AB2`)

Dans ce schéma, la vitesse advectante est extrapolée explicitement à partir des deux pas de temps antérieurs $(n)$ et $(n-1)$ :

$$\mathbf{u}_{\text{adv}} = \frac{3}{2}\mathbf{u}^n - \frac{1}{2}\mathbf{u}^{n-1}$$

L'équation en vitesse devient alors un problème linéaire à coefficients constants durant le pas de temps, découplant l'advection du problème de diffusion-continuité.

---

## Discrétisation Temporelle et Stabilité

Le choix de l'algorithme d'avancement en temps est piloté par la clé `time.scheme` dans le fichier de configuration.

### Schémas Temporels Implémentés

1. **Euler Explicite d'Ordre 1 (`EULER_EXPLICIT`) :**
   $$\frac{\mathbf{u}^{n+1} - \mathbf{u}^n}{\Delta t} = \mathbf{F}(\mathbf{u}^n)$$
   Utilisé principalement pour la vérification numérique ou les régimes purement dominés par la convection avec CFL strict.

2. **Schéma de Gear / BDF2 d'Ordre 2 (`GEAR` / `BDF2`) :**
   Schéma à trois pas de temps implicite, inconditionnellement stable pour les termes linéaires :
   $$\frac{3 \mathbf{u}^{n+1} - 4 \mathbf{u}^n + \mathbf{u}^{n-1}}{2 \Delta t} = \mathbf{F}(\mathbf{u}^{n+1})$$
   Dans COVOF 2D, le coefficient $\alpha_{\text{time}}$ est automatiquement basculé à $1.5$ dès que l'historique temporel $\mathbf{u}^{n-1}$ (`fields.has_nm1`) est disponible, avec repli sur Euler implicite ($\alpha_{\text{time}} = 1.0$) lors du tout premier pas de temps.

3. **Runge-Kutta à Préservation de Stabilité Forte (`SSP_RK3`) :**
   Schéma d'ordre 3 en trois sous-étapes pour l'advection conservative conjointe de la masse et de la quantité de mouvement (modèle Consistent Mass-Momentum d'El Ouafa et al., *Fluids* 2021) :
   $$\mathbf{u}^{(1)} = \mathbf{u}^n + \Delta t \mathcal{L}(\mathbf{u}^n)$$
   $$\mathbf{u}^{(2)} = \frac{3}{4}\mathbf{u}^n + \frac{1}{4}\mathbf{u}^{(1)} + \frac{1}{4}\Delta t \mathcal{L}(\mathbf{u}^{(1)})$$
   $$\mathbf{u}^{n+1} = \frac{1}{3}\mathbf{u}^n + \frac{2}{3}\mathbf{u}^{(2)} + \frac{2}{3}\Delta t \mathcal{L}(\mathbf{u}^{(2)})$$

### Critères de Stabilité et Pas de Temps Adaptatif (CFL)

Lorsque `time.adaptative_dt` est activé, COVOF 2D recalcule à chaque itération le pas de temps maximal admissible $\Delta t$ sous la triple contrainte de Courant-Friedrichs-Lewy (CFL) :

1. **Condition CFL Convective :**
   $$\Delta t_{\text{conv}} = \text{CFL} \cdot \left[ \max_{i, j} \left( \frac{|u_{i+1/2, j}|}{\Delta x} + \frac{|v_{i, j+1/2}|}{\Delta y} \right) \right]^{-1}$$

2. **Condition CFL Diffusive Visqueuse (pour intégration explicite) :**
   $$\Delta t_{\text{diff}} = \frac{1}{2} \min_{i, j} \left( \frac{\rho_{i, j}}{\mu_{i, j}} \right) \left[ \frac{1}{\Delta x^2} + \frac{1}{\Delta y^2} \right]^{-1}$$

3. **Condition Capillaire d'Ondes de Surface (Brackbill et al. 1992) :**
   $$\Delta t_{\text{cap}} = \sqrt{\frac{(\rho_1 + \rho_2) \min(\Delta x, \Delta y)^3}{4 \pi \sigma}}$$

Le pas de temps effectif retenu est :
$$\Delta t^{n+1} = \min \left( \Delta t_{\text{max}}, \Delta t_{\text{conv}}, \Delta t_{\text{diff}}, \Delta t_{\text{cap}} \right)$$

---

## Formulations du Couplage Vitesse-Pression

Le traitement de la contrainte d'incompressibilité $\nabla \cdot \mathbf{u} = 0$ constitue le cœur algorithmique de tout solveur Navier-Stokes. COVOF 2D intègre cinq formulations majeures.

### 1. Lagrangien Augmenté Standard (`AUGMENTED_LAGRANGIAN`)

La méthode du Lagrangien Augmenté (Fortin & Glowinski, 1982 ; Caltagirone & Breil, 1999) pénalise l'équation de quantité de mouvement par la divergence de la vitesse tout en conservant le gradient de pression :

$$\rho \left( \frac{\alpha_{\text{time}} \mathbf{u}^{k+1}}{\Delta t} + (\mathbf{u}^n \cdot \nabla)\mathbf{u}^{k+1} \right) - \nabla \cdot \boldsymbol{\tau}^{k+1} - r \nabla(\nabla \cdot \mathbf{u}^{k+1}) = -\nabla p^k + \mathbf{f}$$

où $r$ est un paramètre scalaire uniforme de pénalité de Lagrangien augmenté ($r \sim 10^2 - 10^4\,\text{Pa}\cdot\text{s}\cdot\text{m}^{-1}$, configuré via `solver.r_al`). La pression est ensuite réactualisée par une étape d'Uzawa :

$$p^{k+1} = p^k - r \nabla \cdot \mathbf{u}^{k+1}$$

Si `solver.uzawa_iter = 1`, la méthode se comporte comme un pseudo-pas fractionnaire très rapide. Si `solver.uzawa_iter > 1`, une boucle itérative d'Uzawa est exécutée jusqu'à ce que $\|\nabla \cdot \mathbf{u}^{k+1}\|_\infty \le \epsilon_{\text{uzawa}}$ (`solver.uzawa_tol`), assurant une divergence rigoureusement nulle au pas de temps $n+1$.

### 2. Lagrangien Augmenté Adaptatif Incompressible (`AAL`)

Dans les écoulements diphasiques présentant de forts contrastes de masse volumique ($\rho_1/\rho_2 \ge 10^3$) et de viscosité ($\mu_1/\mu_2 \ge 10^2$), une valeur scalaire uniforme de $r$ dégrade considérablement le conditionnement du système linéaire. La méthode AAL (*Adaptative Augmented Lagrangian*, Vincent et al., *Computers & Fluids* 2004) détermine automatiquement un champ local de pénalité $r(t, \mathbf{x})$ cellule par cellule :

$$r_{i, j} = K_{\text{aal}} \cdot \max \left( \tau_{\text{instat}}, \tau_{\text{inert}}, \tau_{\text{grav}}, \tau_{\text{pres}}, \tau_{\text{visc}}, \tau_{\text{cap}}, \tau_{\text{brink}} \right)$$

où les temps caractéristiques physiques sont construits à partir des échelles globales du problème :
- Terme instationnaire : $\tau_{\text{instat}} = \alpha_{\text{time}} \rho_{i, j} \frac{L_0^2}{t_0}$
- Terme inertiel : $\tau_{\text{inert}} = \rho_{i, j} u_0 L_0$
- Terme gravitationnel : $\tau_{\text{grav}} = \rho_{i, j} \frac{L_0^2}{u_0} \|\mathbf{g}\|$
- Terme de pression : $\tau_{\text{pres}} = \frac{p_0 L_0}{u_0}$
- Terme visqueux : $\tau_{\text{visc}} = \mu_{i, j}$
- Terme capillaire : $\tau_{\text{cap}} = \frac{\sigma}{u_0}$
- Terme de Brinkman : $\tau_{\text{brink}} = \mu_{i, j} \frac{\Delta x \Delta y}{\mathcal{K}_{i, j}}$

Ce préconditionnement physique optimal assure une convergence quasi insensible aux sauts de propriétés interfaciales.

### 3. Lagrangien Augmenté Compressible Bas-Mach (`C-AAL`)

Pour les régimes compressibles bas-Mach (`solver.type = "C-AAL"`), le terme de compressibilité isotherme $\chi_T = \frac{1}{\rho c_s^2}$ modifie la contrainte de fermeture. Le paramètre adaptatif de pénalité est régularisé selon :

$$r_{\text{eff}, i, j} = \frac{r_{\text{AAL}, i, j}}{1 + r_{\text{AAL}, i, j} \frac{\chi_{T, i, j}}{\Delta t}}$$

et la mise à jour d'Uzawa devient :

$$p^{k+1}_{i, j} = p^k_{i, j} - r_{\text{eff}, i, j} \left[ (\nabla \cdot \mathbf{u}^{k+1})_{i, j} + \frac{\chi_{T, i, j}}{\Delta t} (p^k_{i, j} - p^n_{i, j}) \right]$$

### 4. Résolution Monolithique Tout-Couplé (`MONOLITHIC` / `FULLY_COUPLED`)

Dans la formulation monolithique, la vitesse et la pression sont assemblées simultanément au sein d'un système point-selle global sans paramètre de pénalité ($r = 0$) :

$$\begin{pmatrix} A_{\mathbf{u}} & B^T \\ B & -C_p \end{pmatrix} \begin{pmatrix} \mathbf{u}^{n+1} \\ p^{n+1} \end{pmatrix} = \begin{pmatrix} \mathbf{f}_{\mathbf{u}} \\ g_p \end{pmatrix}$$

où :
- $A_{\mathbf{u}}$ est le bloc de quantité de mouvement vectoriel discret.
- $B^T$ est la matrice du gradient de pression discret.
- $B$ est la matrice de divergence discrète ($B = - (B^T)^T$).
- $C_p = \text{diag}\left( \frac{\chi_T \Delta x \Delta y}{\Delta t} \right)$ est la matrice de compressibilité barotrope ($C_p = 0$ en incompressible strict).

Ce système est résolu de manière unifiée à chaque pas de temps via les solveurs de Krylov couplés au préconditionneur par blocs triangulaires.

### 5. Méthode de Projection Standard (`STANDARD_PROJECTION`)

La méthode de projection de Chorin-Temam sépare la résolution en deux sous-étapes découplées :
1. **Étape prédictive :** Calcul d'un champ de vitesse intermédiaire $\mathbf{u}^*$ sans gradient de pression.
2. **Étape de projection :** Résolution d'une équation de Poisson pour le pseudo-potentiel de pression $\phi$ :
   $$\nabla \cdot \left( \frac{1}{\rho} \nabla \phi \right) = \frac{\nabla \cdot \mathbf{u}^*}{\Delta t}$$
3. **Étape de correction solénoïdale :**
   $$\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho} \nabla \phi, \quad p^{n+1} = p^n + \phi$$

---

## Solveurs Linéaires de Krylov Matrix-Free

COVOF 2D met en œuvre une infrastructure de solveurs de Krylov **sans matrice explicite (*Matrix-Free*)**. L'opérateur complet $A \mathbf{x}$ n'est jamais assemblé sous forme de matrice creuse globale en mémoire, mais calculé à la volée par passage direct des pochoirs volumes finis décalés sur les tableaux locaux avec synchronisation MPI des mailles fantômes.

### Algorithmes de Résolution Linéaire

1. **BiCGSTAB(2) (`LINEAR_BICGSTAB2`) :**
   Algorithme de Sleijpen et van der Vorst (1993). Il combine deux étapes de BiCG avec une étape de minimisation polynomiale locale d'ordre 2, offrant une robustesse exceptionnelle contre les stagnations sur les systèmes asymétriques à fort gradient convectif.
2. **BiCGSTAB (`LINEAR_BICGSTAB`) :**
   Algorithme standard de van der Vorst (1992), rapide et efficace pour les matrices asymétriques modérément raides.
3. **Gradient Conjugué Préconditionné (`LINEAR_CG`) :**
   Algorithme classique de Hestenes et Stiefel (1952). Strictement dédié aux opérateurs symétriques définis positifs (système de Poisson en pression ou étape de diffusion pure).
4. **GMRES(m) (`LINEAR_GMRES`) :**
   *Generalized Minimal Residual* redémarré de Saad et Schultz (1986). Minimise la norme euclidienne du résidu dans le sous-espace de Krylov orthogonalisé par Modified Gram-Schmidt.
5. **MINRES (`LINEAR_MINRES`) :**
   Algorithme de Paige et Saunders (1975). Dédié aux systèmes linéaires **symétriques indéfinis** tels que les systèmes point-selle Navier-Stokes de Stokes ou monolithiques. Repose sur la tridiagonalisation de Lanczos à trois termes et les rotations de Givens, garantissant une minimisation monotone de la norme du résidu avec un coût mémoire minimal.
6. **Solveur Direct Exact (`LINEAR_DIRECT` / `MUMPS`) :**
   Factorisation LU complète par blocs via les routines LAPACK standards `dgetrf` (décomposition $P A = L U$) et `dgetrs` (descente et remontée triangulaires). Pour les régimes stationnaires linéaires, la factorisation factorisée est mise en cache pour éliminer tout surcoût de factorisation aux pas de temps ultérieurs.

---

## Préconditionneurs Linéaires Parallèles

Le choix du préconditionneur est essentiel pour accélérer la convergence des méthodes de Krylov. COVOF 2D propose cinq préconditionneurs natifs :

1. **Jacobi Diagonal (`PRECOND_JACOBI` / `DIAGONAL`) :**
   Inversion de la diagonale enrichie de l'opérateur local. Cette diagonale intègre la contribution instationnaire, la diffusion visqueuse $\mathcal{O}(\mu/\Delta x^2)$, la pénalité de Lagrangien augmenté adaptatif $\mathcal{O}(r/\Delta x^2)$, la résistance de Brinkman $\mu/\mathcal{K}$ et les termes diagonaux de Newton.
2. **SSOR / Gauss-Seidel Symétrique (`PRECOND_SSOR`) :**
   Balayage avant (*forward sweep*) puis balayage arrière (*backward sweep*) sur la grille locale avec coefficient de relaxation optimal $\omega = 1.0$, offrant un lissage haute fréquence efficace.
3. **ILU(0) (`PRECOND_ILU0`) :**
   Factorisation incomplète sans remplissage sur les blocs de vitesse et pression.
4. **Blocs Triangulaires avec Complément de Schur (`PRECOND_BLOCK_TRIANGULAR` / `BT_SBP`) :**
   Préconditionneur spécifique dédié au système monolithique tout-couplé (El Ouafa et al., 2021/2023). Il s'appuie sur une factorisation triangulaire supérieure par blocs :
   $$P = \begin{pmatrix} D_{\mathbf{u}} & B^T \\ 0 & S_p \end{pmatrix}$$
   où $D_{\mathbf{u}}$ est l'approximation diagonale du bloc de quantité de mouvement et $S_p$ est le complément de Schur approché $S_p \approx B D_{\mathbf{u}}^{-1} B^T - C_p$, résolu par balayages de Gauss-Seidel symétriques.
5. **Préconditionneur Personnalisé (`PRECOND_CUSTOM`) :**
   Point d'entrée permettant à l'utilisateur d'injecter une fonction lambda C++ ou un pointeur de procédure Fortran (`custom_precond_hook`) pour tester de nouvelles heuristiques spectrales ou des couplages multigrilles externes.
