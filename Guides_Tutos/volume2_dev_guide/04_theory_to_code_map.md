# Correspondance Théorie $\leftrightarrow$ Discrétisation $\leftrightarrow$ Code Source

Ce chapitre constitue la **section pivot** du Guide de Développement. Il établit la passerelle rigoureuse et exhaustive entre les équations aux dérivées partielles continues de la mécanique des fluides, leur discrétisation spatiotemporelle sur grille décalée MAC (*Marker-And-Cell*), et leur traduction concrète ligne par ligne dans les codes sources C++17 (`COVOF_2D-C++`) et Fortran 2018 (`COVOF_2D-FORTRAN`).

Chaque opérateur numérique majeur est disséqué selon un triptyque méthodologique strict :
1. **Équation discrète exacte :** Formulation mathématique complète indicée en espace et en temps.
2. **Localisation exacte dans le code :** Fichiers, classes, méthodes C++ et modules/routines Fortran.
3. **Analyse ligne par ligne de l'implémentation :** Pochoirs mémoire, régularisations physiques et optimisations SIMD.

---

## 4.1. Diffusion Visqueuse et Moyenne Harmonique aux Nœuds ($\mu_{xy}$)

### 4.1.1. Équation Discrète Exacte

Dans un écoulement à viscosité variable $\mu(\mathbf{x})$, le tenseur des contraintes visqueuses s'écrit $\boldsymbol{\tau} = \mu \left( \nabla \mathbf{u} + (\nabla \mathbf{u})^T \right)$. Sur la grille décalée MAC, la divergence du tenseur visqueux pour la composante horizontale de vitesse $u_{i+1/2, j}$ fait intervenir :
- La contrainte normale $\tau_{xx} = 2 \mu \frac{\partial u}{\partial x}$ évaluée aux centres de mailles $(i, j)$ et $(i+1, j)$.
- La contrainte tangentielle $\tau_{xy} = \mu \left( \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \right)$ évaluée aux sommets (nœuds) $(i+1/2, j+1/2)$ et $(i+1/2, j-1/2)$.

L'opérateur discret de diffusion horizontale $\mathcal{D}_u(u)_{i+1/2, j}$ s'écrit :
$$
\begin{aligned}
\mathcal{D}_u(u)_{i+1/2, j} = & \frac{2 \mu_{i+1, j} \frac{u_{i+3/2, j} - u_{i+1/2, j}}{\Delta x} - 2 \mu_{i, j} \frac{u_{i+1/2, j} - u_{i-1/2, j}}{\Delta x}}{\Delta x} \\
& + \frac{\mu_{xy, i+1/2, j+1/2} \left( \frac{u_{i+1/2, j+1} - u_{i+1/2, j}}{\Delta y} + \frac{v_{i+1, j+1/2} - v_{i, j+1/2}}{\Delta x} \right)}{\Delta y} \\
& - \frac{\mu_{xy, i+1/2, j-1/2} \left( \frac{u_{i+1/2, j} - u_{i+1/2, j-1}}{\Delta y} + \frac{v_{i+1, j-1/2} - v_{i, j-1/2}}{\Delta x} \right)}{\Delta y}
\end{aligned}
$$

Afin de préserver la continuité du flux tangentiel à travers les interfaces diphasiques présentant un fort saut de viscosité ($\mu_1 / \mu_2 \gg 1$), la viscosité au sommet $\mu_{xy}$ est calculée par **moyenne harmonique** des quatre mailles adjacentes :
$$
\mu_{xy, i+1/2, j+1/2} = \frac{4}{\frac{1}{\mu_{i, j}} + \frac{1}{\mu_{i+1, j}} + \frac{1}{\mu_{i, j+1}} + \frac{1}{\mu_{i+1, j+1}}}
$$

### 4.1.2. Localisation dans le Code Source

- **C++17 :**
  - Fichier : `COVOF_2D-C++/src/NavierStokesSolver.cpp`
  - Méthode : `NavierStokesSolver::computeViscousStress(const FlowFields& f, FlowFields& out)`
  - Lignes de code : Lignes 142 à 215.
- **Fortran 2018 :**
  - Fichier : `COVOF_2D-FORTRAN/src/module_navier_stokes.f90`
  - Routine : `subroutine ns_compute_viscous_stress(fields, mesh)`
  - Lignes de code : Lignes 198 à 280.

### 4.1.3. Analyse de l'Implémentation

```fortran
! Calcul de la moyenne harmonique aux sommets
do j = 1, ny + 1
    do i = 1, nx + 1
        mu_sum = 1.0_dp / max(fields%mu(i-1, j-1), 1.0e-12_dp) + &
                 1.0_dp / max(fields%mu(i,   j-1), 1.0e-12_dp) + &
                 1.0_dp / max(fields%mu(i-1, j  ), 1.0e-12_dp) + &
                 1.0_dp / max(fields%mu(i,   j  ), 1.0e-12_dp)
        mu_xy(i, j) = 4.0_dp / mu_sum
    end do
end do
```

*Détails d'optimisation et robustesse :*
1. **Régularisation par seuil minimal :** L'usage de `max(mu, 1.0e-12)` prémunit l'algorithme contre toute division par zéro dans le cas où un fluide fictif non visqueux est simulé.
2. **Localité spatiale :** Le tableau `mu_xy` est précalculé une fois par pas de temps avant l'assemblage de l'opérateur et réutilisé pour les deux composantes de vitesse $u$ et $v$.

---

## 4.2. Advection Convective : Upwind, Centré, Newton et Adams-Bashforth 2

### 4.2.1. Équation Discrète Exacte

Le terme convectif $\nabla \cdot (\rho \mathbf{u} \mathbf{u})$ pour la quantité de mouvement en $x$ s'écrit sous forme conservative :
$$
\mathcal{C}_u(\mathbf{u})_{i+1/2, j} = \frac{(\rho u u)_{i+1, j} - (\rho u u)_{i, j}}{\Delta x} + \frac{(\rho v u)_{i+1/2, j+1/2} - (\rho v u)_{i+1/2, j-1/2}}{\Delta y}
$$

COVOF 2D implémente quatre discrétisations du flux convectif :
1. **Centré d'ordre 2 (*Centered*) :**
   $$u_{i, j} = \frac{u_{i+1/2, j} + u_{i-1/2, j}}{2}, \quad u_{i+1/2, j+1/2} = \frac{u_{i+1/2, j} + u_{i+1/2, j+1}}{2}$$
2. **Amont du 1er ordre (*Upwind*) :**
   $$u_{i, j} = \begin{cases} u_{i-1/2, j} & \text{si } F_{i, j} \ge 0 \\ u_{i+1/2, j} & \text{si } F_{i, j} < 0 \end{cases}$$
3. **Linéarisation de Newton implicite (*Newton Linearization*) :**
   Lors de l'itération non linéaire $k+1$, le produit convectif est linéarisé par Taylor :
   $$\mathbf{u}^{k+1} \cdot \nabla \mathbf{u}^{k+1} \approx \mathbf{u}^k \cdot \nabla \mathbf{u}^{k+1} + \mathbf{u}^{k+1} \cdot \nabla \mathbf{u}^k - \mathbf{u}^k \cdot \nabla \mathbf{u}^k$$
4. **Schéma explicite d'Adams-Bashforth d'ordre 2 (AB2) :**
   $$\mathcal{C}(\mathbf{u})^{n+1/2} = \frac{3}{2} \mathcal{C}(\mathbf{u}^n) - \frac{1}{2} \mathcal{C}(\mathbf{u}^{n-1})$$

### 4.2.2. Localisation dans le Code Source

- **C++17 :**
  - Fichier : `COVOF_2D-C++/src/NavierStokesSolver.cpp`
  - Méthode : `NavierStokesSolver::computeConvectionTerm(const FlowFields& f, std::vector<double>& cu, std::vector<double>& cv)`
  - Lignes de code : Lignes 230 à 345.
- **Fortran 2018 :**
  - Fichier : `COVOF_2D-FORTRAN/src/module_navier_stokes.f90`
  - Routine : `subroutine ns_compute_convection(fields, mesh, cu, cv, scheme_type)`
  - Lignes de code : Lignes 310 à 425.

### 4.2.3. Analyse de l'Implémentation

Dans le code C++, l'aiguillage entre les schémas s'effectue sans surcoût de virtualité par branchement statique :

```cpp
if (config_.convective_scheme == ConvectionScheme::UPWIND) {
    double u_face_east = (u_conv > 0.0) ? u[IDX(i, j)] : u[IDX(i+1, j)];
    double u_face_west = (u_conv_w > 0.0) ? u[IDX(i-1, j)] : u[IDX(i, j)];
    cu[IDX(i, j)] = (u_conv * u_face_east - u_conv_w * u_face_west) * idx;
} else if (config_.convective_scheme == ConvectionScheme::CENTERED) {
    double u_face_east = 0.5 * (u[IDX(i, j)] + u[IDX(i+1, j)]);
    double u_face_west = 0.5 * (u[IDX(i-1, j)] + u[IDX(i, j)]);
    cu[IDX(i, j)] = (u_conv * u_face_east - u_conv_w * u_face_west) * idx;
}
```

---

## 4.3. Pénalisation par Lagrangien Augmenté : AL, AAL et C-AAL

### 4.3.1. Équation Discrète Exacte

La contrainte d'incompressibilité ou de compressibilité acoustique est imposée via l'algorithme du Lagrangien Augmenté :
$$
\mathcal{A}_{AL}(\mathbf{u}) = \mathcal{A}_{\text{NS}}(\mathbf{u}) - r_{\text{penalty}} \nabla (\nabla \cdot \mathbf{u})
$$

Trois variantes distinctes sont implémentées :
1. **AL Standard :**
   Pénalisation de la divergence pure :
   $$- r \nabla (\nabla \cdot \mathbf{u}^{n+1})$$
   avec mise à jour d'Uzawa pour la pression : $p^{k+1} = p^k - r \nabla \cdot \mathbf{u}^{k+1}$.
2. **AAL Incompressible (*Adaptive Augmented Lagrangian*) :**
   Introduction d'un terme mémoire $r_{\text{aal}}$ stocké au centre des mailles :
   $$r_{\text{aal}}^{n+1} = r_{\text{aal}}^n - \gamma (\nabla \cdot \mathbf{u}^{n+1})$$
   Ce terme absorbe l'erreur de divergence accumulée et autorise l'utilisation de pas de temps plus grands ($CFL > 1$) sans perte de conservation de la masse.
3. **C-AAL Compressible (*Compressible All-Mach Augmented Lagrangian*) :**
   Couplage direct avec l'équation d'état acoustique :
   $$\nabla \cdot \mathbf{u} + \frac{1}{\rho c_s^2} \frac{\partial p}{\partial t} = 0$$
   Le terme de pénalisation discrétisé s'enrichit de la célérité du son $c_s$ :
   $$\mathcal{A}_{C-AAL}(\mathbf{u}) = \mathcal{A}(\mathbf{u}) - r \nabla \left( \nabla \cdot \mathbf{u} + \frac{p^{n+1} - p^n}{\rho c_s^2 \Delta t} \right)$$

### 4.3.2. Localisation dans le Code Source

- **C++17 :**
  - Fichier : `COVOF_2D-C++/src/NavierStokesSolver.cpp`
  - Méthodes : `applyAugmentedLagrangianPenalty()`, `updateUzawaPressure()`
  - Lignes de code : Lignes 410 à 495.
- **Fortran 2018 :**
  - Fichier : `COVOF_2D-FORTRAN/src/module_navier_stokes.f90`
  - Routines : `subroutine ns_apply_augmented_lagrangian(...)`, `subroutine ns_update_pressure_uzawa(...)`
  - Lignes de code : Lignes 512 à 620.

### 4.3.3. Analyse de l'Implémentation

```fortran
! Application de l'opérateur -r * grad(div(u)) sur la composante u
do j = 2, ny - 1
    do i = 2, nx
        ! Calcul de la divergence au centre des mailles Est (i) et Ouest (i-1)
        div_east = (u(i+1, j) - u(i, j)) * idx + (v(i, j+1) - v(i, j)) * idy
        div_west = (u(i, j) - u(i-1, j)) * idx + (v(i-1, j+1) - v(i-1, j)) * idy
        
        if (config%use_caal) then
            div_east = div_east + (p(i, j) - p_old(i, j)) / (rho(i, j) * cs2 * dt)
            div_west = div_west + (p(i-1, j) - p_old(i-1, j)) / (rho(i-1, j) * cs2 * dt)
        end if
        
        grad_div = (div_east - div_west) * idx
        Ax_u(i, j) = Ax_u(i, j) - config%r_penalty * grad_div
    end do
end do
```

---

## 4.4. Traitement des Obstacles et Milieux Poreux par Pénalisation de Darcy-Brinkman

### 4.4.1. Équation Discrète Exacte

Pour modéliser des parois rigides immergées ou des matrices poreuses sans recourir à des maillages non structurés déformables, l'équation de Darcy-Brinkman introduit une force de traînée volumique :
$$
\mathbf{F}_{\text{Brinkman}} = - \frac{\mu}{K} \chi_T (\mathbf{u} - \mathbf{u}_{\text{solid}})
$$
où :
- $K$ est la perméabilité intrinsèque du milieu ($K \to 0$ pour un solide parfait, $K \to \infty$ pour un fluide libre).
- $\chi_T \in [0, 1]$ est la fonction caractéristique d'occupation de l'obstacle.

Dans l'assemblage implicite du système linéaire, ce terme ajoute une contribution positive directe sur la diagonale principale :
$$
\left( \frac{\rho}{\Delta t} + \frac{\mu_{i+1/2, j}}{K} \chi_{T, i+1/2, j} \right) u_{i+1/2, j}^{n+1} = \text{RHS}_u + \frac{\mu_{i+1/2, j}}{K} \chi_{T, i+1/2, j} u_{\text{solid}}
$$

### 4.4.2. Localisation dans le Code Source

- **C++17 :**
  - Fichier : `COVOF_2D-C++/src/NavierStokesSolver.cpp`
  - Méthode : `NavierStokesSolver::applyBrinkmanPenalization(...)`
  - Lignes de code : Lignes 510 à 560.
- **Fortran 2018 :**
  - Fichier : `COVOF_2D-FORTRAN/src/module_navier_stokes.f90`
  - Routine : `subroutine ns_apply_brinkman_drag(...)`
  - Lignes de code : Lignes 640 à 690.

---

## 4.5. Suivi d'Interface VOF-PLIC et Advection de Weymouth-Yue

### 4.5.1. Équation Discrète Exacte

L'advection de la fraction volumique $C \in [0, 1]$ est régie par l'équation de transport conservative :
$$
\frac{\partial C}{\partial t} + \nabla \cdot (C \mathbf{u}) = C (\nabla \cdot \mathbf{u})
$$

Pour éliminer toute diffusion numérique d'interface, l'algorithme géométrique en deux étapes se décompose comme suit :
1. **Reconstruction PLIC (*Piecewise Linear Interface Calculation*) :**
   Dans chaque maille mixte ($0 < C_{i,j} < 1$), l'interface est approchée par un segment de droite :
   $$\mathbf{n} \cdot \mathbf{x} + \alpha = 0$$
   où le vecteur normal unitaire $\mathbf{n} = -\frac{\nabla C}{\|\nabla C\|}$ est évalué par l'opérateur de Youngs ou ELVIRA (pochoir $3 \times 3$), et $\alpha$ est calculé par inversion analytique de la fraction volumique sous le plan.
2. **Advection Conservatrice de Weymouth-Yue (2010) :**
   L'advection dimensionnelle alternée (balayage $X$ puis $Y$) calcule le volume géométrique de fluide transitant à travers chaque face :
   $$C_{i,j}^{n+1} = C_{i,j}^n - \frac{\Delta t}{\Delta x} \left( F_{i+1/2, j}^{\text{WY}} - F_{i-1/2, j}^{\text{WY}} \right) - \frac{\Delta t}{\Delta y} \left( G_{i, j+1/2}^{\text{WY}} - G_{i, j-1/2}^{\text{WY}} \right) + C_{i,j}^n \Delta t (\nabla \cdot \mathbf{u})_{i,j}$$
   Le dernier terme compense exactement la divergence numérique locale et assure que $C \in [0, 1]$ rigoureusement sans écrêtage artificiel.

### 4.5.2. Localisation dans le Code Source

- **C++17 :**
  - Fichiers : `COVOF_2D-C++/src/InterfaceReconstruction.cpp`, `COVOF_2D-C++/src/InterfaceTracker.cpp`
  - Méthodes : `InterfaceReconstruction::reconstructPLIC()`, `InterfaceTracker::advectWeymouthYue()`
  - Lignes de code : Lignes 45 à 220 (`InterfaceReconstruction.cpp`) et Lignes 80 à 310 (`InterfaceTracker.cpp`).
- **Fortran 2018 :**
  - Fichiers : `COVOF_2D-FORTRAN/src/module_interface_reconstruction.f90`, `COVOF_2D-FORTRAN/src/module_interface_tracking.f90`
  - Routines : `subroutine plic_reconstruct_2d(...)`, `subroutine vof_advect_weymouth_yue(...)`
  - Lignes de code : Lignes 60 à 250 (`module_interface_reconstruction.f90`) et Lignes 95 à 380 (`module_interface_tracking.f90`).

---

## 4.6. Tension de Surface CSF (*Continuum Surface Force*) et Courbure $\kappa$

### 4.6.1. Équation Discrète Exacte

La force de tension de surface volumique introduite par Brackbill et al. s'écrit :
$$
\mathbf{F}_{\text{CSF}} = \sigma \kappa (\mathbf{x}) \nabla C
$$
où :
- $\sigma$ est le coefficient de tension superficielle ($\text{N}\cdot\text{m}^{-1}$).
- $\kappa$ est la courbure moyenne locale :
  $$\kappa = - \nabla \cdot \mathbf{n} = - \left( \frac{\partial n_x}{\partial x} + \frac{\partial n_y}{\partial y} \right)$$
- $\mathbf{n} = \frac{\nabla \tilde{C}}{\|\nabla \tilde{C}\|}$ est la normale calculée sur un champ de fraction volumique lissé par filtre gaussien ou moyenne pondérée :
  $$\tilde{C}_{i,j} = \frac{1}{16} \begin{bmatrix} 1 & 2 & 1 \\ 2 & 4 & 2 \\ 1 & 2 & 1 \end{bmatrix} * C_{i,j}$$

Sur la grille MAC, la projection de la force sur la face horizontale $u_{i+1/2, j}$ utilise la moyenne de courbure et le gradient discret de $C$ :
$$
F_{\text{CSF}, x, i+1/2, j} = \sigma \left( \frac{\kappa_{i,j} + \kappa_{i+1, j}}{2} \right) \left( \frac{C_{i+1, j} - C_{i, j}}{\Delta x} \right)
$$

### 4.6.2. Localisation dans le Code Source

- **C++17 :**
  - Fichier : `COVOF_2D-C++/src/InterfaceReconstruction.cpp`
  - Méthode : `InterfaceReconstruction::computeCurvature(const std::vector<double>& c, std::vector<double>& kappa)`
  - Lignes de code : Lignes 240 à 315.
- **Fortran 2018 :**
  - Fichier : `COVOF_2D-FORTRAN/src/module_interface_reconstruction.f90`
  - Routine : `subroutine compute_curvature_csf(c, kappa, mesh)`
  - Lignes de code : Lignes 270 à 365.

---

## 4.7. Préconditionneur Triangulaire par Blocs avec Approximation SBP (BT-SBP)

### 4.7.1. Équation Discrète Exacte

Le système couplé complet de Navier-Stokes discrétisé se met sous la forme bloc :
$$
\begin{bmatrix} A_u & 0 & G_x \\ 0 & A_v & G_y \\ D_x & D_y & -C_p \end{bmatrix} \begin{bmatrix} \mathbf{u} \\ \mathbf{v} \\ \mathbf{p} \end{bmatrix} = \begin{bmatrix} \mathbf{f}_u \\ \mathbf{f}_v \\ \mathbf{f}_p \end{bmatrix}
$$
où $A_u, A_v$ sont les blocs cinétiques d'advection-diffusion, $G = [G_x, G_y]^T$ est le gradient discret, $D = [D_x, D_y]$ est la divergence discrète, et $C_p$ représente la matrice de compressibilité éventuelle (nulle en incompressible pur).

Le préconditionneur bloc triangulaire supérieur BT-SBP applique l'opérateur :
$$
\mathcal{P}_{BT}^{-1} = \begin{bmatrix} \hat{A}_u^{-1} & 0 & -\hat{A}_u^{-1} G_x \hat{S}^{-1} \\ 0 & \hat{A}_v^{-1} & -\hat{A}_v^{-1} G_y \hat{S}^{-1} \\ 0 & 0 & -\hat{S}^{-1} \end{bmatrix}
$$
où $\hat{A}_u = \text{diag}(A_u)$ et le complément de Schur $\hat{S} \approx D \hat{A}^{-1} G + C_p$ est approché par un opérateur de Poisson enrichi pour la pression :
$$
\hat{S} \approx \nabla \cdot \left( \frac{\Delta t}{\rho} \nabla \right)
$$
Ce complément de Schur est résolu de manière approchée par quelques cycles de multigrille ou de relaxations SSOR.

### 4.7.2. Localisation dans le Code Source

- **C++17 :**
  - Fichier : `COVOF_2D-C++/src/Preconditioner.cpp`
  - Méthode : `BlockTriangularPreconditioner::apply(...)`
  - Lignes de code : Lignes 185 à 290.
- **Fortran 2018 :**
  - Fichier : `COVOF_2D-FORTRAN/src/module_preconditioners.f90`
  - Routine : `subroutine precond_apply_bt_sbp(r_u, r_v, r_p, z_u, z_v, z_p, mesh, config)`
  - Lignes de code : Lignes 210 à 335.
