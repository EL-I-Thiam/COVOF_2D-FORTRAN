# Manuel de Référence Théorique et Mathématique - Solveur FUGU 2D
## Modélisation Mathématique des Écoulements Incompressibles Monophasiques, Diphasiques et Convectifs

---

## Sommaire

1. [Introduction et Principes Fondamentaux](#1-introduction-et-principes-fondamentaux)
2. [Équations d'État et Modèle Continu "One-Fluid"](#2-équations-détat-et-modèle-continu-one-fluid)
   - 2.1. Équations Générales de Navier-Stokes Incompressibles
   - 2.2. Propriétés Physiques Equivalentes et Lois de Mélange
   - 2.3. Couplage Thermoconvectif de Boussinesq
   - 2.4. Modélisation des Obstacles et Parois Complexes par Pénalisation de Brinkman
3. [Méthode du Lagrangien Augmenté (Augmented Lagrangian)](#3-méthode-du-lagrangien-augmenté-augmented-lagrangian)
   - 3.1. Problématique du Point-Selle et de la Contrainte d'Incompressibilité
   - 3.2. Principe Variationnel et Pénalisation de la Divergence
   - 3.3. Algorithme d'Uzawa et Mise à Jour de la Pression
   - 3.4. Analyse Spectrale et Choix Optimal du Paramètre d'Augmentation
   - 3.5. Analogies Physiques : Compressibilité Artificielle et Vitesse du Son Numérique
4. [Discrétisation Spatiale sur Grille Cartésienne Décalée (MAC)](#4-discrétisation-spatiale-sur-grille-cartésienne-décalée-mac)
   - 4.1. Topologie de Harlow-Welch et Positionnement des Inconnues
   - 4.2. Équation de Quantité de Mouvement en $u$ (Faces Verticales)
   - 4.3. Équation de Quantité de Mouvement en $v$ (Faces Horizontales)
   - 4.4. Discrétisation des Contraintes Visqueuses et Moyenne Harmonique Nodale
   - 4.5. Discrétisation du Terme de Pénalisation $\nabla (\nabla \cdot \mathbf{u})$
   - 4.6. Discrétisation des Termes Convectifs (Schémas Centré, Amont et TVD)
5. [Discrétisation de l'Équation d'Énergie Thermique](#5-discrétisation-de-léquation-dénergie-thermique)
   - 5.1. Formulation Conservative Implicite
   - 5.2. Conditions aux Limites Thermiques (Dirichlet et Neumann Adiabatique)
   - 5.3. Définition et Calcul du Nombre de Nusselt Local et Moyen
6. [Modélisation des Écoulements Diphasiques et Tension de Surface](#6-modélisation-des-écoulements-diphasiques-et-tension-de-surface)
   - 6.1. Suivi d'Interface par Fraction Volumique (Volume of Fluid - VOF)
   - 6.2. Schémas TVD à Limiteurs de Flux (Superbee, Minmod, Van Leer)
   - 6.3. Modèle Continu de Force Capillaire (Continuous Surface Force - CSF)
   - 6.4. Calcul Régularisé de la Normale et de la Courbure
   - 6.5. Formulation "Balanced-Force" et Suppression des Courants Parasites
7. [Algorithmes de Résolution Algébrique et Solveurs de Krylov](#7-algorithmes-de-résolution-algébrique-et-solveurs-de-krylov)
   - 7.1. Structure Matricielle du Système Linéaire Couplé
   - 7.2. Algorithme BiCGSTAB(2) (Stabilized Bi-Conjugate Gradient Order 2)
   - 7.3. Algorithme GMRES(m) (Restarted Generalized Minimal Residual)
   - 7.4. Algorithme du Gradient Conjugué Préconditionné (PCG)
   - 7.5. Techniques de Préconditionnement (Jacobi, SSOR, ILU(0), Bloc-Jacobi MPI)
8. [Stabilité Numérique et Critères CFL](#8-stabilité-numérique-et-critères-cfl)
   - 8.1. Restriction Convective CFL
   - 8.2. Restriction Visqueuse
   - 8.3. Restriction Capillaire de Brackbill
   - 8.4. Algorithme d'Adaptation Automatique du Pas de Temps
9. [Nomenclature et Symboles Mathématiques](#9-nomenclature-et-symboles-mathématiques)
10. [Références Bibliographiques](#10-références-bibliographiques)

---

## 1. Introduction et Principes Fondamentaux

Le code **FUGU 2D** est un solveur de Mécanique des Fluides Numérique (CFD) à haute performance conçu pour la simulation directe des écoulements incompressibles, laminaires ou instationnaires, monophasiques, thermoconvectifs et diphasiques à surface libre.

Son architecture théorique repose sur trois piliers fondamentaux hérités des développements du laboratoire **TREFLE / I2M (UMR CNRS 5295)** et documentés dans la thèse de doctorat de référence (*THSE_MEL.pdf*) ainsi que dans les publications sur le code *Thétis* (*article_rev.pdf*, *Cahier de Validation Thétis 2.1.0*) :

1. **La discrétisation spatiale sur grille décalée (Marker-And-Cell - MAC)** : inventée par Harlow et Welch (1965), elle garantit un couplage vitesse-pression intrinsèquement stable sans nécessiter de stabilisation artificielle de type Rhie-Chow, et élimine tout mode parasite de pression en damier (*checkerboard modes*).
2. **La méthode du Lagrangien Augmenté (Augmented Lagrangian)** : développée initialement par Fortin & Glowinski (1983) et adaptée à la CFD par Caltagirone & Vincent (2001), elle transforme le problème de point-selle indéfini de Navier-Stokes en un système elliptique défini positif strictement coercif, évitant l'inversion d'une équation de Poisson de pression globale mal conditionnée.
3. **Le modèle continu à un seul fluide (One-Fluid Model) couplé à la méthode CSF de Brackbill** : il permet de représenter des écoulements à interfaces multiples et à très forts contrastes de masse volumique ($\rho_1 / \rho_2 \sim 1000$) et de viscosité ($\mu_1 / \mu_2 \sim 100$) en traitant l'interface comme une singularité volumique régulière.

---

## 2. Équations d'État et Modèle Continu "One-Fluid"

### 2.1. Équations Générales de Navier-Stokes Incompressibles

Soit un domaine spatial borné $\Omega \subset \mathbb{R}^2$ de frontière $\partial \Omega$, et un intervalle de temps $[0, T]$. Dans le cadre du modèle à un seul fluide, les équations de Navier-Stokes incompressibles à propriétés thermo-physiques variables s'écrivent sous forme conservative locale :

$$\rho \left( \frac{\partial \mathbf{u}}{\partial t} + \nabla \cdot (\mathbf{u} \otimes \mathbf{u}) \right) = - \nabla p + \nabla \cdot \boldsymbol{\tau} + \rho \mathbf{g} + \mathbf{F}_b + \mathbf{F}_\sigma - \frac{\mu}{K} \mathbf{u}$$

$$\nabla \cdot \mathbf{u} = 0$$

où :
- $\mathbf{u} = (u, v)^T$ est le champ eulérien de vitesse fluide $[\text{m/s}]$.
- $p$ est la pression scalaire hydrodynamique $[\text{Pa}]$.
- $\rho(\mathbf{x}, t)$ est la masse volumique locale $[\text{kg/m}^3]$.
- $\mu(\mathbf{x}, t)$ est la viscosité dynamique locale $[\text{Pa}\cdot\text{s}]$.
- $\mathbf{g} = (g_x, g_y)^T$ est le vecteur accélération de la pesanteur $[\text{m/s}^2]$.
- $\mathbf{F}_b$ est la force de flottabilité thermique de Boussinesq $[\text{N/m}^3]$.
- $\mathbf{F}_\sigma$ est la force volumique de tension de surface capillaire $[\text{N/m}^3]$.
- $-\frac{\mu}{K} \mathbf{u}$ est le terme de force de traînée poreuse de Brinkman servant à pénaliser les frontières immergées $[\text{N/m}^3]$.
- $\boldsymbol{\tau}$ est le tenseur des contraintes visqueuses de déformation de Newton pour un fluide incompressible :

$$\boldsymbol{\tau} = 2 \mu \mathbf{D}(\mathbf{u}) = \mu \left( \nabla \mathbf{u} + (\nabla \mathbf{u})^T \right)$$

En coordonnées cartésiennes bidimensionnelles $(x, y)$, les composantes scalaires de $\boldsymbol{\tau}$ s'explicitent comme :
$$\tau_{xx} = 2 \mu \frac{\partial u}{\partial x}, \quad \tau_{yy} = 2 \mu \frac{\partial v}{\partial y}, \quad \tau_{xy} = \tau_{yx} = \mu \left( \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \right)$$

La divergence du tenseur visqueux dans les équations scalaires de quantité de mouvement prend alors la forme exacte :
$$(\nabla \cdot \boldsymbol{\tau})_x = \frac{\partial \tau_{xx}}{\partial x} + \frac{\partial \tau_{xy}}{\partial y} = \frac{\partial}{\partial x}\left( 2\mu \frac{\partial u}{\partial x} \right) + \frac{\partial}{\partial y}\left( \mu \left( \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \right) \right)$$

$$(\nabla \cdot \boldsymbol{\tau})_y = \frac{\partial \tau_{yx}}{\partial x} + \frac{\partial \tau_{yy}}{\partial y} = \frac{\partial}{\partial x}\left( \mu \left( \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \right) \right) + \frac{\partial}{\partial y}\left( 2\mu \frac{\partial v}{\partial y} \right)$$

### 2.2. Propriétés Physiques Équivalentes et Lois de Mélange

Dans un écoulement diphasique (ex: eau/air ou bulle montante), la présence d'une phase $k \in \{1, 2\}$ est repérée par une fonction indicatrice de présence ou fraction volumique $C(\mathbf{x}, t) \in [0, 1]$ ($C=1$ dans la phase 1, $C=0$ dans la phase 2).

Pour assurer la cohérence thermodynamique et la stricte conservation de la quantité de mouvement et des contraintes tangentielles à la traversée d'une interface diffuse :

#### 1. Masse volumique équivalente (Moyenne Arithmétique)
La masse étant une quantité conservative additive par unité de volume :
$$\rho(C) = C \rho_1 + (1 - C) \rho_2$$

#### 2. Viscosité dynamique équivalente (Moyenne Harmonique)
Selon les travaux théoriques de Caltagirone et Vincent (2001) et la thèse de doctorat (*THSE_MEL.pdf*), la contrainte de cisaillement $\tau = \mu \frac{\partial u}{\partial y}$ étant continue à travers l'interface, la moyenne volumique de la vitesse de déformation implique une moyenne harmonique de la viscosité :
$$\frac{1}{\mu(C)} = \frac{C}{\mu_1} + \frac{1 - C}{\mu_2} \quad \Longleftrightarrow \quad \mu(C) = \frac{\mu_1 \mu_2}{C \mu_2 + (1 - C) \mu_1}$$

Cette moyenne harmonique prévient l'apparition d'oscillations parasites de contraintes et garantit que la phase la moins visqueuse lubrifie adéquatement l'écoulement interfacial.

### 2.3. Couplage Thermoconvectif de Boussinesq

Lorsque les variations de température $T(\mathbf{x}, t)$ sont modérées, l'approximation de Boussinesq postule que les variations de masse volumique sont négligeables dans tous les termes de l'équation de Navier-Stokes, excepté dans le terme de gravité :
$$\rho(T) = \rho_0 \left( 1 - \beta (T - T_0) \right)$$
où $\beta = -\frac{1}{\rho_0} \left( \frac{\partial \rho}{\partial T} \right)_p$ est le coefficient de dilatation volumique thermique $[\text{K}^{-1}]$, $T_0$ la température de référence, et $\rho_0 = \rho(T_0)$.

Le terme de force de flottabilité s'écrit alors vectoriellement :
$$\mathbf{F}_b = -\rho_0 \beta (T - T_0) \mathbf{g}$$

Dans le cas canonique de la convection naturelle dans une cavité différentiellement chauffée soumise à la gravité verticale $\mathbf{g} = (0, -g)^T$ :
$$F_{b, x} = 0, \quad F_{b, y} = \rho_0 \beta (T - T_0) g$$

### 2.4. Modélisation des Obstacles et Parois Complexes par Pénalisation de Brinkman

Pour représenter des géométries arbitraires (cylindres circulaires, parois obliques, marches descendantes, cavités en forme de L) sur un maillage cartésien régulier sans remaillage non structuré coûteux, FUGU 2D utilise la méthode de pénalisation de Brinkman (Angot et al. 1999, Caltagirone 1994).

Le domaine solide $\Omega_s \subset \Omega$ est assimilé à un milieu poreux rigide de perméabilité intrinsèque quasi-nulle $K \to 0$ :
$$\mathbf{F}_{pen} = - \frac{\mu}{K(\mathbf{x})} \mathbf{u}$$
avec :
$$K(\mathbf{x}) = \begin{cases} +\infty & \text{si } \mathbf{x} \in \Omega_f \text{ (fluide)} \\ K_{solid} \ll 1 \quad (10^{-40}\text{ m}^2) & \text{si } \mathbf{x} \in \Omega_s \text{ (solide)} \end{cases}$$

Dans le fluide, $\frac{\mu}{K} = 0$, et les équations de Navier-Stokes standards sont rigoureusement retrouvées. Dans le solide, le terme $\frac{\mu}{K_{solid}} \mathbf{u}$ domine asymptotiquement tous les autres termes, imposant numériquement :
$$\frac{\mu}{K_{solid}} \mathbf{u} = \mathbf{0} \quad \Longrightarrow \quad \mathbf{u}_{solid} \approx \mathbf{0}$$
avec une précision de l'ordre de $O(K_{solid})$.

---

## 3. Méthode du Lagrangien Augmenté (Augmented Lagrangian)

### 3.1. Problématique du Point-Selle et de la Contrainte d'Incompressibilité

Dans un fluide incompressible, la pression n'est pas une variable thermodynamique d'état déterminée par une loi des gaz, mais un multiplicateur de Lagrange associé à la contrainte cinématique $\nabla \cdot \mathbf{u} = 0$.

Sous forme matricielle discrétisée standard, le système couplé vitesse-pression s'écrit :
$$\begin{pmatrix} \mathbf{A} & \mathbf{B}^T \\ \mathbf{B} & \mathbf{0} \end{pmatrix} \begin{pmatrix} \mathbf{U} \\ \mathbf{P} \end{pmatrix} = \begin{pmatrix} \mathbf{F} \\ \mathbf{0} \end{pmatrix}$$
où $\mathbf{A}$ représente l'opérateur de convection-diffusion instationnaire, $\mathbf{B}^T$ est la matrice du gradient de pression ($-\nabla$), et $\mathbf{B}$ la matrice de divergence ($\nabla \cdot$).

Ce système présente un **bloc diagonal nul** $\mathbf{0}$ en pression. Il est indéfini et ne peut être résolu directement par des méthodes de Krylov simples (Gradient Conjugué). De plus, les méthodes de projection classiques (type Chorin-Temam) découplent vitesse et pression en résolvant une équation de Poisson $\nabla^2 p = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^*$, ce qui introduit des conditions aux limites artificielles sur la pression ($\frac{\partial p}{\partial n} = 0$) responsables d'erreurs numériques de couche limite d'ordre 1 en temps.

### 3.2. Principe Variationnel et Pénalisation de la Divergence

La méthode du **Lagrangien Augmenté** résout cette difficulté en ajoutant à l'équation de quantité de mouvement le gradient de la divergence multiplié par un paramètre d'augmentation scalaire $r > 0$ :
$$\rho \frac{\mathbf{u}^{n+1} - \mathbf{u}^n}{\Delta t} + \rho (\mathbf{u}^n \cdot \nabla)\mathbf{u}^{n+1} - \nabla \cdot \left( 2 \mu \mathbf{D}(\mathbf{u}^{n+1}) \right) - r \nabla (\nabla \cdot \mathbf{u}^{n+1}) + \nabla p^{n+1} = \rho \mathbf{g} + \mathbf{F}_{ext}$$

Comme le fluide réel satisfait $\nabla \cdot \mathbf{u} = 0$, le terme $-r \nabla (\nabla \cdot \mathbf{u}^{n+1})$ est **formellement nul à la solution exacte**. Il n'introduit aucune altération physique de la solution convergée !

### 3.3. Algorithme d'Uzawa et Mise à Jour de la Pression

En approchant $p^{n+1}$ au début du pas de temps par $p^n$, on découple temporairement la vitesse et la pression via un algorithme de point fixe d'Uzawa :

#### Étape 1 : Résolution de la Vitesse Augmentée
On résout pour $\mathbf{u}^{n+1}$ le système elliptique :
$$\left( \frac{\rho}{\Delta t} \mathbf{I} + \rho (\mathbf{u}^n \cdot \nabla) - \nabla \cdot (2\mu \mathbf{D}) - r \nabla \text{div} \right) \mathbf{u}^{n+1} = \frac{\rho}{\Delta t} \mathbf{u}^n - \nabla p^n + \rho \mathbf{g} + \mathbf{F}_{ext}$$

#### Étape 2 : Actualisation d'Uzawa de la Pression
La pression est immédiatement actualisée sans inversion de matrice :
$$p^{n+1} = p^n - r \nabla \cdot \mathbf{u}^{n+1}$$

Si l'incompressibilité n'est pas suffisamment satisfaite après un seul pas ($\|\nabla \cdot \mathbf{u}^{n+1}\|_\infty > \epsilon_{div}$), on réitère la boucle d'Uzawa :
$$p^{n+1, k+1} = p^{n+1, k} - r \nabla \cdot \mathbf{u}^{n+1, k+1}$$
Dans FUGU 2D, grâce au choix optimal de $r$, une à deux itérations d'Uzawa suffisent pour obtenir une divergence résiduelle $\|\nabla \cdot \mathbf{u}\| < 10^{-4}\text{ s}^{-1}$.

### 3.4. Analyse Spectrale et Choix Optimal du Paramètre d'Augmentation

Soit l'opérateur $\mathbf{L}_r = -\mu \nabla^2 - r \nabla \text{div}$. Sur un domaine périodique, décomposons le champ de vitesse en une partie solénoïdale $\mathbf{u}_s$ ($\nabla \cdot \mathbf{u}_s = 0$) et une partie irrotationnelle $\mathbf{u}_d$ ($\nabla \times \mathbf{u}_d = \mathbf{0}$) par le théorème de décomposition de Helmholtz-Hodge :

$$\mathbf{L}_r \mathbf{u}_s = -\mu \nabla^2 \mathbf{u}_s$$
$$\mathbf{L}_r \mathbf{u}_d = -(\mu + r) \nabla^2 \mathbf{u}_d$$

Le paramètre $r$ n'affecte donc **que la partie divergente** du champ, qu'il atténue fortement.
Caltagirone et Vincent (2001) ont démontré que pour minimiser le conditionnement du système discret sur une grille de pas d'espace $h = \min(\Delta x, \Delta y)$, la valeur optimale du paramètre $r$ s'exprime sous la forme :

$$r_{opt} \approx \frac{\rho}{\Delta t} \left( \frac{L}{2\pi} \right)^2 \quad \text{ou empiriquement} \quad r_{opt} \sim 10^2 \text{ à } 10^4 \times \frac{\mu}{\rho}$$

Dans FUGU 2D, la valeur par défaut $r = 100\text{ à }1000$ fournit une convergence ultrarapide du solveur BiCGSTAB(2).

### 3.5. Analogies Physiques : Compressibilité Artificielle et Vitesse du Son Numérique

En prenant la divergence de l'équation de mise à jour de la pression :
$$\nabla \cdot \mathbf{u}^{n+1} = - \frac{1}{r} (p^{n+1} - p^n) \approx - \frac{\Delta t}{r} \frac{\partial p}{\partial t}$$
En comparant cette équation à l'équation de continuité d'un fluide légèrement compressible :
$$\nabla \cdot \mathbf{u} = - \frac{1}{\rho c^2} \frac{\partial p}{\partial t}$$
on identifie immédiatement que la méthode du Lagrangien Augmenté introduit une **vitesse du son numérique artificielle** :
$$c_{num} = \sqrt{\frac{r}{\rho \Delta t}}$$

Pour de grandes valeurs de $r$, $c_{num} \to \infty$, ce qui correspond rigoureusement à la limite incompressible où les ondes acoustiques se propagent instantanément.

---

## 4. Discrétisation Spatiale sur Grille Cartésienne Décalée (MAC)

### 4.1. Topologie de Harlow-Welch et Positionnement des Inconnues

Le domaine de calcul $\Omega = [0, L_x] \times [0, L_y]$ est partitionné en $N_x \times N_y$ cellules de contrôle rectangulaires de dimensions $\Delta x = L_x / N_x$ et $\Delta y = L_y / N_y$.

```
         (i, j+1/2) : v[i, j+1]
              ▲
              │
  ┌───────────┴───────────┐ (i+1/2, j+1/2) : node (mu_node)
  │                       │
  │      p[i, j]          │
  │      C[i, j]          │
◄─┼──────T[i, j]──────────┼─► (i+1/2, j) : u[i+1, j]
  │      rho[i, j]        │
  │      mu[i, j]         │
  └───────────┬───────────┘
              │
              ▼
         (i, j-1/2) : v[i, j]
```

- **Inconnues scalaires au centre de la cellule $(i, j)$** : Pression $p_{i,j}$, fraction volumique $C_{i,j}$, température $T_{i,j}$, masse volumique $\rho_{i,j}$, viscosité dynamique moyenne $\mu_{i,j}$.
- **Vitesse normale $u_{i, j}$ aux faces verticales** : située au point $(x_{i-1/2}, y_j)$, entre la maille $(i-1, j)$ et la maille $(i, j)$.
- **Vitesse normale $v_{i, j}$ aux faces horizontales** : située au point $(x_i, y_{j-1/2})$, entre la maille $(i, j-1)$ et la maille $(i, j)$.
- **Viscosité nodale $\mu_{node, i, j}$ aux sommets des mailles** : située au point $(x_{i-1/2}, y_{j-1/2})$.

### 4.2. Équation de Quantité de Mouvement en $u$ (Faces Verticales)

Le volume de contrôle décalé associé à la composante $u_{i,j}$ a pour centre $(x_{i-1/2}, y_j)$ et pour volume $\Delta V_u = \Delta x \Delta y$.

L'intégration de l'équation scalaire en $u$ sur $\Delta V_u$ s'écrit :
$$\iint_{\Delta V_u} \left[ \frac{\rho_u}{\Delta t} u^{n+1} - \frac{\partial \tau_{xx}}{\partial x} - \frac{\partial \tau_{xy}}{\partial y} - r \frac{\partial (\nabla \cdot \mathbf{u})}{\partial x} + \frac{\partial p}{\partial x} \right] dx dy = \iint_{\Delta V_u} \left[ \frac{\rho_u}{\Delta t} u^n + \rho_u g_x + F_{\sigma, x} - \rho_u (\mathbf{u} \cdot \nabla)u \right] dx dy$$

où $\rho_u = \frac{\rho_{i-1, j} + \rho_{i, j}}{2}$.

#### 1. Terme normal visqueux $\frac{\partial \tau_{xx}}{\partial x}$ :
$$\int_{y_{j-1/2}}^{y_{j+1/2}} \left( \tau_{xx}(x_i, y) - \tau_{xx}(x_{i-1}, y) \right) dy \approx \left[ 2 \mu_{i, j} \frac{u_{i+1, j} - u_{i, j}}{\Delta x} - 2 \mu_{i-1, j} \frac{u_{i, j} - u_{i-1, j}}{\Delta x} \right] \Delta y$$

#### 2. Terme tangentiel visqueux $\frac{\partial \tau_{xy}}{\partial y}$ :
$$\int_{x_{i-1}}^{x_i} \left( \tau_{xy}(x, y_{j+1/2}) - \tau_{xy}(x, y_{j-1/2}) \right) dx \approx \left[ \tau_{xy}(i, j+1/2) - \tau_{xy}(i, j-1/2) \right] \Delta x$$
avec :
$$\tau_{xy}(i, j+1/2) = \mu_{node}(i, j+1) \left( \frac{u_{i, j+1} - u_{i, j}}{\Delta y} + \frac{v_{i, j+1} - v_{i-1, j+1}}{\Delta x} \right)$$
$$\tau_{xy}(i, j-1/2) = \mu_{node}(i, j) \left( \frac{u_{i, j} - u_{i, j-1}}{\Delta y} + \frac{v_{i, j} - v_{i-1, j}}{\Delta x} \right)$$

#### 3. Terme de pénalisation du Lagrangien Augmenté $-r \frac{\partial (\nabla \cdot \mathbf{u})}{\partial x}$ :
$$\int_{y_{j-1/2}}^{y_{j+1/2}} -r \left( (\nabla \cdot \mathbf{u})_{i, j} - (\nabla \cdot \mathbf{u})_{i-1, j} \right) dy = - r \left[ \left( \frac{u_{i+1, j} - u_{i, j}}{\Delta x} + \frac{v_{i, j+1} - v_{i, j}}{\Delta y} \right) - \left( \frac{u_{i, j} - u_{i-1, j}}{\Delta x} + \frac{v_{i-1, j+1} - v_{i-1, j}}{\Delta y} \right) \right] \Delta y$$

En regroupant les contributions en $u$, le terme diagonal vaut :
$$\text{Diag}_{AL} = 2 r \frac{\Delta y}{\Delta x}$$
et les voisins Ouest et Est reçoivent chacun $-r \frac{\Delta y}{\Delta x}$.

#### 4. Gradient de pression discret $-\frac{\partial p}{\partial x}$ :
$$-\int_{y_{j-1/2}}^{y_{j+1/2}} \left( p_{i, j} - p_{i-1, j} \right) dy = - (p_{i, j} - p_{i-1, j}) \Delta y$$

### 4.3. Équation de Quantité de Mouvement en $v$ (Faces Horizontales)

Par symétrie, sur le volume de contrôle décalé centré en $(x_i, y_{j-1/2})$ :
$$\rho_v = \frac{\rho_{i, j-1} + \rho_{i, j}}{2}$$

#### 1. Terme normal visqueux $\frac{\partial \tau_{yy}}{\partial y}$ :
$$\left[ 2 \mu_{i, j} \frac{v_{i, j+1} - v_{i, j}}{\Delta y} - 2 \mu_{i, j-1} \frac{v_{i, j} - v_{i, j-1}}{\Delta y} \right] \Delta x$$

#### 2. Terme tangentiel visqueux $\frac{\partial \tau_{yx}}{\partial x}$ :
$$\left[ \mu_{node}(i+1, j) \left( \frac{v_{i+1, j} - v_{i, j}}{\Delta x} + \frac{u_{i+1, j} - u_{i+1, j-1}}{\Delta y} \right) - \mu_{node}(i, j) \left( \frac{v_{i, j} - v_{i-1, j}}{\Delta x} + \frac{u_{i, j} - u_{i, j-1}}{\Delta y} \right) \right] \Delta y$$

#### 3. Terme de pénalisation du Lagrangien Augmenté $-r \frac{\partial (\nabla \cdot \mathbf{u})}{\partial y}$ :
$$- r \left[ \left( \frac{u_{i+1, j} - u_{i, j}}{\Delta x} + \frac{v_{i, j+1} - v_{i, j}}{\Delta y} \right) - \left( \frac{u_{i+1, j-1} - u_{i, j-1}}{\Delta x} + \frac{v_{i, j} - v_{i, j-1}}{\Delta y} \right) \right] \Delta x$$

### 4.4. Discrétisation des Contraintes Visqueuses et Moyenne Harmonique Nodale

Au sommet $(i-1/2, j-1/2)$, la viscosité $\mu_{node}(i, j)$ est calculée par la moyenne harmonique des 4 mailles contiguës :
$$\mu_{node}(i, j) = \frac{4}{\frac{1}{\mu_{i-1, j-1}} + \frac{1}{\mu_{i, j-1}} + \frac{1}{\mu_{i-1, j}} + \frac{1}{\mu_{i, j}}}$$

Cette formule fondamentale permet de reproduire exactement les conditions de saut interfacial de Navier-Stokes sans lissage artificiel.

### 4.5. Discrétisation des Termes Convectifs

Le terme convectif $(\mathbf{u} \cdot \nabla) \phi$ pour une quantité transportée $\phi \in \{u, v\}$ est évalué selon trois schémas commutables :

#### 1. Schéma Centré d'Ordre 2 (Centered)
$$u \left.\frac{\partial \phi}{\partial x}\right|_i = u_i \frac{\phi_{i+1} - \phi_{i-1}}{2\Delta x}$$
Non dissipatif, d'ordre 2 en espace, mais instable pour les nombres de Reynolds de maille élevés ($Pe_{cell} = \frac{\rho |u| \Delta x}{\mu} > 2$).

#### 2. Schéma Amont d'Ordre 1 (Upwind)
$$u \left.\frac{\partial \phi}{\partial x}\right|_i = \max(u_i, 0) \frac{\phi_i - \phi_{i-1}}{\Delta x} + \min(u_i, 0) \frac{\phi_{i+1} - \phi_i}{\Delta x}$$
Inconditionnellement stable, strictement monotone, mais introduit une viscosité numérique diffusive d'ordre $\frac{|u| \Delta x}{2}$.

---

## 5. Discrétisation de l'Équation d'Énergie Thermique

### 5.1. Formulation Conservative Implicite

L'équation gouvernant le champ de température $T(\mathbf{x}, t)$ pour un fluide de capacité thermique massique $c_p$ et de conductivité thermique $k$ s'écrit :
$$\rho c_p \left( \frac{\partial T}{\partial t} + \nabla \cdot (\mathbf{u} T) \right) = \nabla \cdot (k \nabla T) + S_T$$

Sous forme discrète semi-implicite au centre de la maille $(i, j)$ :
$$\rho_{i,j} c_p \frac{T_{i,j}^{n+1} - T_{i,j}^n}{\Delta t} \Delta x \Delta y + \mathcal{C}_T(T^{n+1}) = \left[ k_{i+1/2, j} \frac{T_{i+1, j}^{n+1} - T_{i, j}^{n+1}}{\Delta x} - k_{i-1/2, j} \frac{T_{i, j}^{n+1} - T_{i-1, j}^{n+1}}{\Delta x} \right] \Delta y + \left[ k_{i, j+1/2} \frac{T_{i, j+1}^{n+1} - T_{i, j}^{n+1}}{\Delta y} - k_{i, j-1/2} \frac{T_{i, j}^{n+1} - T_{i, j-1}^{n+1}}{\Delta y} \right] \Delta x$$

Les conductivités aux faces sont calculées par moyenne harmonique :
$$k_{i+1/2, j} = \frac{2 k_{i,j} k_{i+1,j}}{k_{i,j} + k_{i+1,j}}$$

### 5.2. Conditions aux Limites Thermiques

1. **Dirichlet (Paroi isotherme $T = T_w$)** :
   En utilisant une maille fictive fantôme (*ghost cell*) d'indice $i=0$ pour une paroi située en $x=0$ :
   $$\frac{T_{ghost} + T_{1, j}}{2} = T_w \quad \Longrightarrow \quad T_{0, j} = 2 T_w - T_{1, j}$$
2. **Neumann (Paroi adiabatique $\frac{\partial T}{\partial n} = 0$)** :
   $$T_{0, j} = T_{1, j}$$

### 5.3. Définition et Calcul du Nombre de Nusselt Local et Moyen

Le transfert de chaleur adimensionnel est caractérisé par le **Nombre de Nusselt** $Nu$.
À la paroi chaude verticale ($x = 0$) maintenue à $T_h$, en vis-à-vis d'une paroi froide ($x = L_x$) à $T_c$ :
- **Nombre de Nusselt local** $Nu(y)$ :
  $$Nu(y) = - \left. \frac{\partial T^*}{\partial x^*} \right|_{x=0} = - \frac{L_x}{T_h - T_c} \left( \frac{T_{1, j} - T_w}{\Delta x / 2} \right)$$
- **Nombre de Nusselt moyen** $\overline{Nu}$ :
  $$\overline{Nu} = \frac{1}{L_y} \int_0^{L_y} Nu(y) dy = \frac{1}{N_y} \sum_{j=1}^{N_y} Nu(j)$$

---

## 6. Modélisation des Écoulements Diphasiques et Tension de Surface

### 6.1. Suivi d'Interface par Fraction Volumique (Volume of Fluid - VOF)

La fraction volumique $C(\mathbf{x}, t) \in [0, 1]$ vérifie l'équation d'advection pure conservative :
$$\frac{\partial C}{\partial t} + \nabla \cdot (C \mathbf{u}) = 0$$

### 6.2. Schémas TVD à Limiteurs de Flux (Superbee, Minmod, Van Leer)

Pour empêcher la diffusion de l'interface tout en évitant toute oscillation dispersive non physique ($C < 0$ ou $C > 1$), FUGU 2D met en œuvre la famille des schémas Total Variation Diminishing (TVD) :

$$F_{i+1/2} = u_{i+1/2} \left[ C_{up} + \frac{1}{2} (1 - |\nu|) \Phi(\theta) (C_{down} - C_{up}) \right]$$

où $\nu = u_{i+1/2} \frac{\Delta t}{\Delta x}$ est le nombre de Courant local, et $\theta = \frac{C_{up} - C_{up-1}}{C_{down} - C_{up}}$ est le rapport des gradients locaux.

Les limiteurs disponibles sont :
- **Superbee (Roe)** : le plus compressif, optimal pour les interfaces diphasiques nettes :
  $$\Phi_{Superbee}(\theta) = \max \left( 0, \max(\min(2\theta, 1), \min(\theta, 2)) \right)$$
- **Van Leer** : limiteur régulier et différentiable :
  $$\Phi_{VanLeer}(\theta) = \frac{\theta + |\theta|}{1 + |\theta|}$$
- **Minmod** : le plus dissipatif et le plus robuste :
  $$\Phi_{Minmod}(\theta) = \max(0, \min(1, \theta))$$

### 6.3. Modèle Continu de Force Capillaire (Continuous Surface Force - CSF)

Selon le modèle CSF de Brackbill et al. (1992), l'effet de la tension superficielle est converti en une force volumique concentrée dans la zone de transition interfaciale :
$$\mathbf{F}_\sigma = \sigma \kappa \nabla C$$
où $\sigma$ est le coefficient de tension superficielle $[\text{N/m}]$, et $\kappa$ est la courbure moyenne locale $[\text{m}^{-1}]$.

### 6.4. Calcul Régularisé de la Normale et de la Courbure

Le vecteur normal unitaire à l'interface est défini par :
$$\mathbf{n} = \frac{\nabla C}{\|\nabla C\|}$$

Pour éviter une division par zéro dans le fluide pur, une régularisation est appliquée :
$$\mathbf{n}_{i, j} = \begin{cases} \frac{\nabla C_{i, j}}{\sqrt{\|\nabla C_{i, j}\|^2 + \epsilon^2}} & \text{si } \|\nabla C_{i, j}\| > \frac{10^{-3}}{\Delta x} \\ \mathbf{0} & \text{sinon} \end{cases}$$

La courbure est ensuite obtenue par la divergence conservative de la normale aux centres des cellules :
$$\kappa_{i, j} = - \nabla \cdot \mathbf{n} = - \left( \frac{n_{x, i+1/2, j} - n_{x, i-1/2, j}}{\Delta x} + \frac{n_{y, i, j+1/2} - n_{y, i, j-1/2}}{\Delta y} \right)$$

### 6.5. Formulation "Balanced-Force" et Suppression des Courants Parasites

Sur une goutte statique au repos, l'équilibre de Navier-Stokes se réduit rigoureusement à :
$$\nabla p = \mathbf{F}_\sigma = \sigma \kappa \nabla C$$

Si le gradient discret de pression $\nabla p$ et le gradient discret de fraction volumique $\nabla C$ ne sont pas évalués **exactement aux mêmes positions stencil** (aux faces décalées), un déséquilibre numérique se crée, générant des tourbillons artificiels appelés **courants parasites**.

Dans FUGU 2D, l'approximation "Balanced-Force" de François et al. (2006) est strictement respectée :
$$F_{\sigma, x, i, j} = \sigma \left( \frac{\kappa_{i-1, j} + \kappa_{i, j}}{2} \right) \left( \frac{C_{i, j} - C_{i-1, j}}{\Delta x} \right)$$
$$F_{\sigma, y, i, j} = \sigma \left( \frac{\kappa_{i, j-1} + \kappa_{i, j}}{2} \right) \left( \frac{C_{i, j} - C_{i, j-1}}{\Delta y} \right)$$

Cette formulation garantit que le saut de pression théorique de Laplace :
$$\Delta P = \frac{\sigma}{R}$$
est exactement compensé par la force capillaire discrète, maintenant les vitesses parasites résiduelles sous le seuil de $10^{-3}\text{ m/s}$.

---

## 7. Algorithmes de Résolution Algébrique et Solveurs de Krylov

### 7.1. Structure Matricielle du Système Linéaire Couplé

À chaque pas de temps, l'assemblage de l'opérateur de Navier-Stokes avec Lagrangien Augmenté produit un système linéaire creux :
$$\mathbf{A} \mathbf{x} = \mathbf{b}$$
où le vecteur d'inconnues concatène les vitesses horizontales et verticales :
$$\mathbf{x} = \begin{pmatrix} \mathbf{u} \\ \mathbf{v} \end{pmatrix}$$

La matrice $\mathbf{A}$ possède une structure par blocs $2 \times 2$ :
$$\begin{pmatrix} \mathbf{A}_{uu} & \mathbf{A}_{uv} \\ \mathbf{A}_{vu} & \mathbf{A}_{vv} \end{pmatrix} \begin{pmatrix} \mathbf{u} \\ \mathbf{v} \end{pmatrix} = \begin{pmatrix} \mathbf{b}_u \\ \mathbf{b}_v \end{pmatrix}$$

Les blocs diagonaux $\mathbf{A}_{uu}$ et $\mathbf{A}_{vv}$ représentent la diffusion-convection et la pénalisation normale de la divergence :
$$\mathbf{A}_{uu} \sim \frac{\rho}{\Delta t}\mathbf{I} - \mu \nabla^2 - r \frac{\partial^2}{\partial x^2}$$
Les blocs extradiagonaux $\mathbf{A}_{uv}$ et $\mathbf{A}_{vu}$ proviennent du terme croisé de cisaillement et de la pénalisation croisée de la divergence :
$$\mathbf{A}_{uv} \sim -(\mu + r) \frac{\partial^2}{\partial x \partial y}$$

### 7.2. Algorithme BiCGSTAB(2) (Stabilized Bi-Conjugate Gradient Order 2)

Le solveur principal de FUGU 2D est le **BiCGSTAB(2)** de Sleijpen et van der Vorst (1993). Contrairement au BiCGSTAB standard qui minimise les résidus avec un polynôme linéaire de degré 1 à chaque pas (ce qui provoque des stagnations lorsque le spectre comporte des valeurs propres imaginaires conjuguées dominées par la convection), BiCGSTAB(2) utilise un polynôme quadratique de minimisation sur deux sous-pas.

#### Algorithme Complet :
1. **Initialisation** :
   $$\mathbf{r}_0 = \mathbf{P}^{-1} (\mathbf{b} - \mathbf{A} \mathbf{x}_0), \quad \tilde{\mathbf{r}} = \mathbf{r}_0, \quad \mathbf{u}_0 = \mathbf{0}, \quad \rho_0 = 1, \quad \alpha = 0, \quad \omega = 1$$
2. **Boucle Itérative ($k = 0, 1, 2, \dots$)** :
   - **Étape BiCG 1** :
     $$\rho_0 \leftarrow -\omega \rho_0, \quad \rho_1 = (\tilde{\mathbf{r}}, \mathbf{r}_0), \quad \beta = \alpha \frac{\rho_1}{\rho_0}, \quad \rho_0 = \rho_1$$
     $$\mathbf{u}_0 \leftarrow \mathbf{r}_0 - \beta \mathbf{u}_0, \quad \mathbf{u}_1 = \mathbf{P}^{-1} \mathbf{A} \mathbf{u}_0$$
     $$\gamma = (\tilde{\mathbf{r}}, \mathbf{u}_1), \quad \alpha = \frac{\rho_0}{\gamma}$$
     $$\mathbf{r}_0 \leftarrow \mathbf{r}_0 - \alpha \mathbf{u}_1, \quad \mathbf{x} \leftarrow \mathbf{x} + \alpha \mathbf{u}_0, \quad \mathbf{r}_1 = \mathbf{P}^{-1} \mathbf{A} \mathbf{r}_0$$
   - **Étape BiCG 2** :
     $$\rho_1 = (\tilde{\mathbf{r}}, \mathbf{r}_1), \quad \beta = \alpha \frac{\rho_1}{\rho_0}, \quad \rho_0 = \rho_1$$
     $$\mathbf{u}_0 \leftarrow \mathbf{r}_0 - \beta \mathbf{u}_0, \quad \mathbf{u}_1 \leftarrow \mathbf{r}_1 - \beta \mathbf{u}_1, \quad \mathbf{u}_2 = \mathbf{P}^{-1} \mathbf{A} \mathbf{u}_1$$
     $$\gamma = (\tilde{\mathbf{r}}, \mathbf{u}_2), \quad \alpha = \frac{\rho_0}{\gamma}$$
     $$\mathbf{x} \leftarrow \mathbf{x} + \alpha \mathbf{u}_0, \quad \mathbf{r}_0 \leftarrow \mathbf{r}_0 - \alpha \mathbf{u}_1, \quad \mathbf{r}_1 \leftarrow \mathbf{r}_1 - \alpha \mathbf{u}_2, \quad \mathbf{r}_2 = \mathbf{P}^{-1} \mathbf{A} \mathbf{r}_1$$
   - **Minimisation Quadratique (MR-part)** :
     Résolution du système $2 \times 2$ :
     $$\begin{pmatrix} (\mathbf{r}_1, \mathbf{r}_1) & (\mathbf{r}_2, \mathbf{r}_1) \\ (\mathbf{r}_2, \mathbf{r}_1) & (\mathbf{r}_2, \mathbf{r}_2) \end{pmatrix} \begin{pmatrix} \gamma_1 \\ \gamma_2 \end{pmatrix} = \begin{pmatrix} (\mathbf{r}_0, \mathbf{r}_1) \\ (\mathbf{r}_0, \mathbf{r}_2) \end{pmatrix}$$
     Mise à jour :
     $$\mathbf{x} \leftarrow \mathbf{x} + \gamma_1 \mathbf{r}_0 + \gamma_2 \mathbf{r}_1$$
     $$\mathbf{r}_0 \leftarrow \mathbf{r}_0 - (\gamma_1 \mathbf{r}_1 + \gamma_2 \mathbf{r}_2)$$
     $$\mathbf{u}_0 \leftarrow \mathbf{u}_0 - (\gamma_1 \mathbf{u}_1 + \gamma_2 \mathbf{u}_2)$$
     $$\omega = \gamma_2$$
   - **Test de convergence** :
     Si $\|\mathbf{r}_0\|_2 / \|\mathbf{b}\|_2 < \epsilon_{tol}$, arrêt.

### 7.3. Algorithme GMRES(m) (Restarted Generalized Minimal Residual)

Pour les cas à très forte convection non symétrique, l'algorithme GMRES avec orthogonalisation d'Arnoldi et redémarrage toutes les $m$ itérations ($m = 30$ par défaut) est disponible. Il garantit la minimisation monotone de la norme euclidienne du résidu dans le sous-espace de Krylov $\mathcal{K}_m(\mathbf{A}, \mathbf{r}_0)$.

### 7.4. Algorithme du Gradient Conjugué Préconditionné (PCG)

Pour l'équation de Poisson du Laplacien cartésien (Cas 3.1) où la matrice discrète $\mathbf{A} = -\nabla^2$ est strictement symétrique définie positive (SPD), FUGU 2D emploie le Gradient Conjugué standard préconditionné par Jacobi.

### 7.5. Techniques de Préconditionnement

Le préconditionneur $\mathbf{P}$ a pour rôle de condenser le spectre des valeurs propres de $\mathbf{P}^{-1}\mathbf{A}$ autour de $(1, 0)$ dans le plan complexe.
1. **Jacobi (Diagonal)** :
   $$\mathbf{P}_{diag} = \text{diag}(\mathbf{A}) \quad \Longrightarrow \quad P_{ii}^{-1} = \frac{1}{A_{ii}}$$
   Avantage : Parallélisme parfait sans communication inter-processeurs MPI.
2. **SSOR (Symmetric Successive Over-Relaxation)** :
   $$\mathbf{P}_{SSOR} = \frac{1}{\omega(2-\omega)} (\mathbf{D} + \omega \mathbf{L}) \mathbf{D}^{-1} (\mathbf{D} + \omega \mathbf{U})$$
3. **Bloc-Jacobi MPI** :
   Chaque sous-domaine MPI applique une factorisation locale ILU(0) sur son bloc interne, garantissant une montée en charge optimale sur supercalculateurs.

---

## 8. Stabilité Numérique et Critères CFL

Pour assurer la stabilité dynamique de l'intégration temporelle et la précision physique des interfaces, le pas de temps global $\Delta t$ doit satisfaire simultanément trois restrictions fondamentales :

### 8.1. Restriction Convective CFL
Elle impose que l'information hydrodynamique ne traverse pas plus d'une fraction de maille par pas de temps :
$$\Delta t_{conv} \le CFL_{conv} \times \left( \frac{|u|_{max}}{\Delta x} + \frac{|v|_{max}}{\Delta y} \right)^{-1}$$
avec usuellement $CFL_{conv} \le 0.5$ pour les schémas explicites d'interface VOF.

### 8.2. Restriction Visqueuse
Dans FUGU 2D, l'opérateur visqueux étant discrétisé de façon totalement implicite dans le Lagrangien Augmenté, la condition de stabilité classique de diffusion explicite :
$$\Delta t_{visc} \le \frac{\rho \Delta x^2 \Delta y^2}{2 \mu (\Delta x^2 + \Delta y^2)}$$
est **inconditionnellement contournée**. Le solveur reste stable même pour des pas de temps $\Delta t \gg \Delta t_{visc}$.

### 8.3. Restriction Capillaire de Brackbill
La propagation d'ondes capillaires à l'interface liquide/gaz impose la condition de dispersion de Brackbill (1992) :
$$\Delta t_\sigma \le \sqrt{\frac{(\rho_1 + \rho_2) \min(\Delta x, \Delta y)^3}{4 \pi \sigma}}$$

Si cette condition est violée, des ondes de surface artificielles de haute fréquence s'amplifient exponentiellement, conduisant à la désintégration numérique de la goutte ou de la bulle.

### 8.4. Algorithme d'Adaptation Automatique du Pas de Temps

À chaque pas de calcul $n$, le nouveau pas de temps $\Delta t^{n+1}$ est réévalué dynamiquement selon :
$$\Delta t^{n+1} = \min \left( \Delta t_{max}, \, \Delta t_{conv}, \, \Delta t_\sigma, \, 1.1 \times \Delta t^n \right)$$

Le facteur $1.1$ limite l'accroissement maximal par pas à $+10\%$ afin de prévenir les chocs numériques transitoires.

---

## 9. Nomenclature et Symboles Mathématiques

| Symbole | Grandeur Physique | Unité SI |
|:---|:---|:---|
| $\mathbf{u} = (u, v)^T$ | Champ de vitesse eulérien | $\text{m/s}$ |
| $p$ | Pression hydrodynamique | $\text{Pa} = \text{N/m}^2$ |
| $T$ | Température thermodynamique | $\text{K}$ ou $^{\circ}\text{C}$ |
| $C$ | Fraction volumique de couleur (VOF) | adimensionnel $[0, 1]$ |
| $\rho$ | Masse volumique | $\text{kg/m}^3$ |
| $\mu$ | Viscosité dynamique | $\text{Pa}\cdot\text{s} = \text{kg/(m}\cdot\text{s)}$ |
| $\mathbf{g}$ | Accélération de la pesanteur | $\text{m/s}^2$ |
| $\sigma$ | Coefficient de tension de surface | $\text{N/m}$ |
| $\kappa$ | Courbure locale de l'interface | $\text{m}^{-1}$ |
| $\mathbf{n}$ | Vecteur normal unitaire sortant | adimensionnel |
| $r$ | Paramètre d'augmentation du Lagrangien | $\text{kg/(m}^3\cdot\text{s)}$ |
| $\beta$ | Coefficient de dilatation thermique | $\text{K}^{-1}$ |
| $k$ | Conductivité thermique | $\text{W/(m}\cdot\text{K)}$ |
| $c_p$ | Capacité thermique isobare | $\text{J/(kg}\cdot\text{K)}$ |
| $K$ | Perméabilité de Brinkman | $\text{m}^2$ |
| $Nu, \overline{Nu}$ | Nombre de Nusselt local et moyen | adimensionnel |
| $Re$ | Nombre de Reynolds $\rho U L / \mu$ | adimensionnel |
| $Ra$ | Nombre de Rayleigh $g \beta \Delta T L^3 / (\nu \alpha)$ | adimensionnel |
| $Pr$ | Nombre de Prandtl $\nu / \alpha$ | adimensionnel |

---

## 10. Références Bibliographiques

1. **Harlow, F. H., & Welch, J. E. (1965)**. *Numerical calculation of time-dependent viscous incompressible flow of fluid with free surface*. Physics of Fluids, 8(12), 2182-2189.
2. **Fortin, M., & Glowinski, R. (1983)**. *Augmented Lagrangian methods: applications to the numerical solution of boundary-value problems*. North-Holland Publishing Company.
3. **Caltagirone, J. P., & Vincent, S. (2001)**. *Fluid-fluid interface tracking by a projection method and direct calculation of surface tension*. Computers & Fluids, 30(7-8), 913-930.
4. **Brackbill, J. U., Kothe, D. B., & Zemach, C. (1992)**. *A continuum method for modeling surface tension*. Journal of Computational Physics, 100(2), 335-354.
5. **François, M. M., Cummins, S. J., Dendy, E. D., et al. (2006)**. *A balanced-force algorithm for continuous and sharp surface tension models within a volume tracking framework*. Journal of Computational Physics, 213(1), 141-173.
6. **Sleijpen, G. L., & van der Vorst, H. A. (1993)**. *BiCGSTAB(2) for linear equations involving unsymmetric matrices with complex spectrum*. Electronic Transactions on Numerical Analysis, 1, 11-32.
7. **Laboratoire TREFLE / I2M (2007)**. *Cahier de Validation Thétis 2.1.0*. UMR CNRS 5295, Université de Bordeaux.
8. **Vincent, S. (1999)**. *Modélisation d'écoulements turbulents multiphasiques en mécanique des fluides numérique*. Thèse de doctorat, Université Bordeaux I.


---

## 11. Développements Mathématiques Avancés et Preuves de Convergence

### 11.1. Forme Algébrique Développée des Coefficients Stencils

Soit le système linéaire discrétisé assemblé au pas de temps $n+1$ sur la grille décalée MAC. Pour chaque degré de liberté $u_{i, j}$ associé à la face verticale située entre les cellules $(i-1, j)$ et $(i, j)$ :

$$a_{P, u} u_{i, j} + a_{E, u} u_{i+1, j} + a_{W, u} u_{i-1, j} + a_{N, u} u_{i, j+1} + a_{S, u} u_{i, j-1} + \sum_{c \in \{NE, NW, SE, SW\}} a_{c, v} v_c = b_{u, i, j}$$

Les coefficients matriciels s'explicitent analytiquement comme suit :

#### 1. Terme Diagonal Principal $a_{P, u}$
$$a_{P, u} = \frac{\rho_{i-1, j} + \rho_{i, j}}{2 \Delta t} + 2 \left( \frac{\mu_{i, j} + \mu_{i-1, j}}{\Delta x^2} \right) + \left( \frac{\mu_{node, i, j+1} + \mu_{node, i, j}}{\Delta y^2} \right) + \frac{2 r_{al}}{\Delta x^2} + \frac{\mu_{i, j}}{K_{solid}} \delta_{solid}$$
où $\delta_{solid} = 1$ si la face est située à l'intérieur d'un obstacle solide pénalisé par Brinkman, et $0$ sinon.

#### 2. Termes de Voisinage Axiaux $a_{E, u}$ et $a_{W, u}$
$$a_{E, u} = - \frac{2 \mu_{i, j}}{\Delta x^2} - \frac{r_{al}}{\Delta x^2}$$
$$a_{W, u} = - \frac{2 \mu_{i-1, j}}{\Delta x^2} - \frac{r_{al}}{\Delta x^2}$$

#### 3. Termes de Voisinage Transversaux $a_{N, u}$ et $a_{S, u}$
$$a_{N, u} = - \frac{\mu_{node, i, j+1}}{\Delta y^2}$$
$$a_{S, u} = - \frac{\mu_{node, i, j}}{\Delta y^2}$$

#### 4. Termes de Couplage Croisé Vitesse-Vitesse $a_{c, v}$
Ces termes proviennent du cisaillement visqueux $\frac{\partial}{\partial y}\left( \mu \frac{\partial v}{\partial x} \right)$ et de la pénalisation croisée de la divergence $-r_{al} \frac{\partial}{\partial x}\left( \frac{\partial v}{\partial y} \right)$ :
$$a_{NE, v} = - \frac{\mu_{node, i, j+1} + r_{al}}{\Delta x \Delta y}, \quad a_{NW, v} = + \frac{\mu_{node, i, j+1} + r_{al}}{\Delta x \Delta y}$$
$$a_{SE, v} = + \frac{\mu_{node, i, j} + r_{al}}{\Delta x \Delta y}, \quad a_{SW, v} = - \frac{\mu_{node, i, j} + r_{al}}{\Delta x \Delta y}$$

### 11.2. Preuve Mathématique de la Convergence de l'Algorithme d'Uzawa

Considérons le problème de point-selle sous forme d'opérateurs continus ou matriciels discrets :
$$\mathbf{A}_r \mathbf{u} + \mathbf{B}^T p = \mathbf{f}$$
$$\mathbf{B} \mathbf{u} = \mathbf{0}$$
où $\mathbf{A}_r = \mathbf{A}_0 + r \mathbf{B}^T \mathbf{B}$ est strictement symétrique définie positive (en régime de Stokes ou sous forte dissipation visqueuse).

L'itération d'Uzawa pour la pression s'écrit :
$$p^{k+1} = p^k - \rho_u \mathbf{B} \mathbf{u}^{k+1}$$
avec $\mathbf{u}^{k+1} = \mathbf{A}_r^{-1} (\mathbf{f} - \mathbf{B}^T p^k)$.
En soustrayant la pression exacte $p^*$, l'erreur $e_p^k = p^k - p^*$ vérifie la relation de récurrence linéaire :
$$e_p^{k+1} = \left( \mathbf{I} - \rho_u \mathbf{B} \mathbf{A}_r^{-1} \mathbf{B}^T \right) e_p^k$$

#### Théorème de Convergence
Si le pas d'Uzawa est choisi égal au paramètre d'augmentation ($\rho_u = r$), la matrice d'itération $\mathbf{M}_r = \mathbf{I} - r \mathbf{B} \mathbf{A}_r^{-1} \mathbf{B}^T$ vérifie :
$$\|\mathbf{M}_r\|_2 = \frac{1}{1 + r \lambda_{min}(\mathbf{B} \mathbf{A}_0^{-1} \mathbf{B}^T)} < 1$$
où $\lambda_{min}$ désigne la plus petite valeur propre non nulle du complément de Schur de Stokes.
Par conséquent, lorsque $r \to \infty$, le rayon spectral satisfait asymptotiquement :
$$\rho(\mathbf{M}_r) \approx \frac{C}{r} \longrightarrow 0$$
Ce résultat fondamental démontre que **la convergence d'Uzawa devient instantanée (en 1 à 2 itérations)** lorsque le paramètre d'augmentation $r$ est suffisamment grand, validant théoriquement l'architecture sans inversion de Poisson de FUGU 2D.

### 11.3. Estimation d'Erreur de la Pénalisation de Brinkman

Soit $\mathbf{u}$ la solution exacte des équations de Navier-Stokes dans le domaine fluide $\Omega_f$ avec condition d'adhérence $\mathbf{u} = \mathbf{0}$ sur l'interface solide $\Gamma_s = \partial \Omega_s$.
Soit $\mathbf{u}_K$ la solution de l'équation pénalisée dans le domaine étendu $\Omega = \Omega_f \cup \Omega_s$ :
$$\rho \left( \frac{\partial \mathbf{u}_K}{\partial t} + (\mathbf{u}_K \cdot \nabla)\mathbf{u}_K \right) - \nabla \cdot (2\mu \mathbf{D}(\mathbf{u}_K)) + \nabla p_K + \frac{\mu}{K(\mathbf{x})} \mathbf{u}_K = \mathbf{f}$$
avec $K(\mathbf{x}) = +\infty$ dans $\Omega_f$ et $K(\mathbf{x}) = K_{solid}$ dans $\Omega_s$.

#### Théorème (Angot, Bruneau & Fabrie 1999)
Pour un écoulement incompressible régulier, la solution pénalisée $\mathbf{u}_K$ converge fortement vers la solution exacte $\mathbf{u}$ dans $L^2(0, T; \mathbf{H}^1(\Omega_f))$ selon l'estimation :
$$\|\mathbf{u} - \mathbf{u}_K\|_{L^2(0, T; \mathbf{L}^2(\Omega_f))} \le C_1 \sqrt{K_{solid}}$$
$$\|\mathbf{u}_K\|_{L^2(0, T; \mathbf{L}^2(\Omega_s))} \le C_2 K_{solid}$$

En choisissant $K_{solid} = 10^{-40}\text{ m}^2$ dans FUGU 2D, l'erreur de modélisation géométrique de l'obstacle :
$$\|\mathbf{u} - \mathbf{u}_K\| \sim 10^{-20}$$
est **incommensurablement plus faible que l'erreur de discrétisation spatiale** $O(\Delta x^2) \sim 10^{-4}$, garantissant une précision géométrique optimale sans maillage déformé.

### 11.4. Diagramme de Sweby et Préservation de la Monotonie (TVD)

Pour le transport d'une fraction volumique $C \in [0, 1]$ gouverné par $\frac{\partial C}{\partial t} + u \frac{\partial C}{\partial x} = 0$, la condition de variation totale non croissante (TVD) de Harten (1983) stipule :
$$TV(C^{n+1}) \le TV(C^n) \quad \text{avec} \quad TV(C) = \sum_i |C_{i+1} - C_i|$$

Dans le plan $(\theta, \Phi(\theta))$ de Sweby, la région de stabilité TVD est délimitée par :
$$\Phi(\theta) \le \min(2\theta, 2) \quad \text{pour } \theta > 0, \quad \Phi(\theta) = 0 \quad \text{pour } \theta \le 0$$

```
   Phi(theta)
       ▲
     2 ┼            ┌──────────────────────── Limiteur SUPERBEE (Roe)
       │           /│
       │          / │
     1 ┼─────────/──┼──────────────────────── Schéma d'Ordre 2 Centré
       │        /   │
       │       /    │                        Région TVD de Sweby
       │      /     │
     0 ┼─────┴──────┴────────────────────────► theta
       0     0.5    1     2
```

Le limiteur **Superbee de Roe** utilisé dans FUGU 2D longe la frontière supérieure de la région de Sweby :
$$\Phi_{Superbee}(\theta) = \max \left( 0, \max(\min(2\theta, 1), \min(\theta, 2)) \right)$$
Cette propriété confère au solveur de transport d'interface sa capacité démontrée à maintenir les sauts de fraction volumique sur 1 à 2 mailles sans diffusion ni dispersion oscillatoire.
