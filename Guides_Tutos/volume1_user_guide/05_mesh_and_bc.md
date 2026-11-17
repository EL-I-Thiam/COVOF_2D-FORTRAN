# Géométrie, Maillages et Conditions aux Limites

## Génération et Topologie du Maillage Cartésien MAC

Le solveur **COVOF 2D** opère sur des domaines rectangulaires $[0, L_x] \times [0, L_y]$ discrétisés par un maillage cartésien structuré uniforme selon la topologie décalée MAC (*Marker-and-Cell*).

### Définition Géométrique du Domaine (`mesh`)

Le domaine de calcul est entièrement spécifié dans le bloc `mesh` du fichier JSON par quatre paramètres exacts :
- `mesh.nx` (`integer`, défaut `64`) : nombre global de cellules de calcul selon l'axe horizontal $X$.
- `mesh.ny` (`integer`, défaut `64`) : nombre global de cellules de calcul selon l'axe vertical $Y$.
- `mesh.lx` (`double`, défaut `1.0` m) : longueur physique totale du domaine selon $X$.
- `mesh.ly` (`double`, défaut `1.0` m) : hauteur physique totale du domaine selon $Y$.

Les pas d'espace constants sont calculés à l'initialisation selon :

$$\Delta x = \frac{L_x}{N_x}, \qquad \Delta y = \frac{L_y}{N_y}$$

### Coordonnées des Centres, Faces et Sommets

Pour une maille repérée localement par $(i, j)$ avec décalage global `(i_offset, j_offset)` en C++ (base 0) ou `(i_start_global - 1, j_start_global - 1)` en Fortran (base 1), les fonctions d'accès physiques retournent exactement :

- **Centres de cellules (Pression $p$, Fraction VOF $C$, $C_s$, $\widetilde{C}$, $\rho$, $\mu$, $T$, $\kappa$, $\chi_T$) :**
  $$x_c(i) = \left(i_{\text{offset}} + i - \frac{1}{2}\right)\Delta x, \qquad y_c(j) = \left(j_{\text{offset}} + j - \frac{1}{2}\right)\Delta y$$
- **Faces verticales (Vitesse $u$, $\rho_u$, $F_{\sigma, x}$) :**
  $$x_u(i) = \left(i_{\text{offset}} + i - 1\right)\Delta x, \qquad y_c(j) = \left(j_{\text{offset}} + j - \frac{1}{2}\right)\Delta y$$
  La face `u(1, j)` coïncide avec le bord gauche $x = 0$ du sous-domaine (ou du domaine global), et la face `u(nx_local + 1, j)` avec le bord droit.
- **Faces horizontales (Vitesse $v$, $\rho_v$, $F_{\sigma, y}$) :**
  $$x_c(i) = \left(i_{\text{offset}} + i - \frac{1}{2}\right)\Delta x, \qquad y_v(j) = \left(j_{\text{offset}} + j - 1\right)\Delta y$$
  La face `v(i, 1)` coïncide avec le bord inférieur $y = 0$, et `v(i, ny_local + 1)` avec le bord supérieur.
- **Sommets / Nœuds (Viscosité harmonique nodale $\mu_{\text{node}}$) :**
  $$\big(x_u(i), \, y_v(j)\big)$$

---

## Décomposition de Domaine et Parallélisme Cartésien 2D (MPI pur)

Pour l'exécution distribuée, COVOF 2D partitionne la grille globale $(N_x, N_y)$ en $P = P_x \times P_y$ sous-domaines rectangulaires gérés par des processus MPI indépendants (**MPI pur**, sans OpenMP).

```text
       +-----------------------+-----------------------+
       |  Rang MPI (0, 1)      |  Rang MPI (1, 1)      |
       |  is_boundary_left     |  is_boundary_right    |
       |  is_boundary_top      |  is_boundary_top      |
Ny / 2 +-----------------------+-----------------------+
       |  Rang MPI (0, 0)      |  Rang MPI (1, 0)      |
       |  is_boundary_left     |  is_boundary_right    |
       |  is_boundary_bottom   |  is_boundary_bottom   |
       +-----------------------+-----------------------+
       0                    Nx / 2                     Nx
```

### Topologie Virtuelle Cartésienne (`MPI_Cart_create`)

1. **Factorisation automatique et gestion des restes :**
   `MPI_Dims_create(num_procs, 2, dims)` détermine la grille $(P_x, P_y)$ la plus équilibrée. Si $N_x$ n'est pas divisible par $P_x$, les `nx_global % px` premiers rangs reçoivent une maille supplémentaire (`Mesh2D::init` / `init_mesh`).
2. **Périodicité et voisins cardinaux :**
   Les drapeaux `periodic_x` et `periodic_y` sont activés automatiquement si `left`/`right` ou `bottom`/`top` ont pour type `"PERIODIC"`. `MPI_Cart_shift` identifie les 4 voisins directs (`rank_left`, `rank_right`, `rank_bottom`, `rank_top` en C++ ; `neighbor_west`, `neighbor_east`, `neighbor_south`, `neighbor_north` en Fortran), valant `MPI_PROC_NULL` sur les parois physiques non périodiques.

### Couche Unique de Cellules Fantômes (*Halos*)

Chaque sous-domaine alloue exactement **1 couche de cellules fantômes** de chaque côté :
- Dimensions allouées pour les scalaires centrés : `nx_ghost = nx_local + 2`, `ny_ghost = ny_local + 2`, soit des indices locaux valides $i \in [0, N_{x,\text{loc}}+1]$ et $j \in [0, N_{y,\text{loc}}+1]$.
- Zone intérieure utile : $i \in [1, N_{x,\text{loc}}]$, $j \in [1, N_{y,\text{loc}}]$.
- L'échange non bloquant (`MpiComm::sync_cell_scalar`, `sync_u`, `sync_v` en C++ ; `exchange_halos_scalar`, `exchange_halos_velocity` en Fortran) procède en deux passes séquentielles (`MPI_Irecv`/`MPI_Isend` + `MPI_Waitall`) : d'abord selon $X$, puis selon $Y$ sur toute la largeur `[0 .. nx_local + 1]` afin de **propager exactement les 4 coins diagonaux** sans nécessiter de communications diagonales.

---

## Catalogue Exhaustif des Conditions aux Limites (CL)

Les quatre frontières du domaine (`left`, `right`, `bottom`, `top` dans `boundary_conditions`) sont configurées indépendamment. L'analyseur JSON reconnaît **quatre types hydrodynamiques fondamentaux (`BCType`)**, complétés par le mot-clé `profile` pour les entrées spécialisées et `T_type` pour la thermique.

:::note
Contrairement à certaines conventions génériques, les mots-clés `NO_SLIP`, `PRESSURE`, `PRESSURE_OUTLET`, `OUTFLOW`, `SYMMETRY` ou `PISTON` **n'existent pas** dans les analyseurs `Config.cpp` et `module_config.f90` (toute chaîne inconnue rebascule silencieusement sur `DIRICHLET`). Les comportements correspondants s'expriment exactement à l'aide des quatre types ci-dessous et de `profile`.
:::

### 1. Condition de Dirichlet (`"DIRICHLET"`, défaut)

Impose le vecteur vitesse $\mathbf{u} = (u_{\text{bc}}, v_{\text{bc}})^T$ sur la frontière (`u`, `v`).

- **Paroi solide immobile (*No-Slip*) :** `"type": "DIRICHLET", "u": 0.0, "v": 0.0`.
- **Paroi mobile (ex. cavité entraînée) :** `"type": "DIRICHLET", "u": 1.0, "v": 0.0` sur `"top"`.
- **Piston mobile compressible :** `"type": "DIRICHLET", "u": 0.0, "v": v_piston` sur `"bottom"`.

#### Implémentation Discrète MAC (`MpiComm.cpp` / `module_mpi_comm.f90`)
- **Composante normale sur face frontière :** imposée directement sur la face et sa fantôme externe (ex. sur `left` : $u_{1,j} = u_{0,j} = u_{\text{inlet}}$ ; sur `top` : $v_{i, N_y+1} = v_{i, N_y+2} = v_{\text{top}}$).
- **Composante tangentielle :**
  - Sur les parois horizontales (`bottom` et `top`), **extrapolation quadratique d'ordre 2 à 3 points** (exacte pour les profils paraboliques de Poiseuille lorsque $N_{y,\text{loc}} \ge 2$) :
    $$u_{i, 0} = \frac{8 u_{\text{bot}} - 6 u_{i, 1} + u_{i, 2}}{3}, \qquad u_{i, N_y+1} = \frac{8 u_{\text{top}} - 6 u_{i, N_y} + u_{i, N_y-1}}{3}$$
    *(avec repli linéaire $2 u_{\text{wall}} - u_{\text{int}}$ si $N_{y,\text{loc}} < 2$)*.
  - Sur les parois verticales (`left` et `right`), réflexion linéaire d'ordre 2 :
    $$v_{0, j} = 2 v_{\text{left}} - v_{1, j}, \qquad v_{N_x+1, j} = 2 v_{\text{right}} - v_{N_x, j}$$
- **Pression et VOF :** Neumann homogène sur les parois ($p_{\text{ghost}} = p_{\text{int}}$, $C_{\text{ghost}} = C_{\text{int}}$).

#### Profils Spécialisés (`boundary_conditions.<face>.profile`)

| Valeur de `profile` | Face Supportée | Paramètres Associés | Formulation Mathématique Exacte |
| :--- | :--- | :--- | :--- |
| `"UNIFORM"` (défaut) | Toutes | `u`, `v` | $u = u_{\text{bc}}, \; v = v_{\text{bc}}$ uniforme sur toute la face |
| `"PARABOLIC"` | `left` | `u_mean` (C++), `u` (Fortran) | $u(y) = 6 \, u_{\text{mean}} \frac{y (L_y - y)}{L_y^2}$ (profil de Poiseuille établi, $\bar{u} = u_{\text{mean}}$) |
| `"POISEUILLE_STEP"` | `left` | `u_mean`, `y_step` (C++) / `obs_y1` (Fortran) | $u(y) = 0$ si $y < y_{\text{step}}$ ; $u(y) = 6 \, u_{\text{mean}} \frac{(y - y_{\text{step}})(L_y - y)}{(L_y - y_{\text{step}})^2}$ si $y \ge y_{\text{step}}$ |
| `"JET_INLET"` | `top` | `v`, `jet_center`, `jet_width` | Injecteur centré en $x_c$ sur $|x - x_c| \le W/2$ avec $v_{\text{in,eff}} = v_{\text{top}} \frac{W}{D_{\text{disc}}}$, $C_{\text{ghost}} = 0$ (Fluide 2), $u=0$ ; et **débit compensatoire d'évacuation du gaz** hors buse $v_{\text{out}} = -v_{\text{top}} \frac{W}{L_x - D_{\text{disc}}}$ assurant $\int_{\partial\Omega} \mathbf{u}\cdot\mathbf{n}\,ds = 0$ |

### 2. Condition de Neumann (`"NEUMANN"`)

Prescrit une dérivée normale nulle pour la vitesse, la pression et la fraction volumique :

$$\frac{\partial u}{\partial n} = 0, \qquad \frac{\partial v}{\partial n} = 0, \qquad \frac{\partial p}{\partial n} = 0, \qquad \frac{\partial C}{\partial n} = 0$$

- **Implémentation :** Copie directe de l'intérieur vers la couche fantôme (ex. sur `right` : $u_{N_x+2, j} = u_{N_x+1, j}$, $v_{N_x+1, j} = v_{N_x, j}$, $p_{N_x+1, j} = p_{N_x, j}$).
- **Usage :** Sortie libre (*outflow*) en aval de canal, marche descendante ou sillage de cylindre.

### 3. Glissement Parfait / Symétrie (`"FREE_SLIP"`)

Impose l'imperméabilité normale et l'absence de cisaillement tangentiel :

$$\mathbf{u} \cdot \mathbf{n} = 0, \qquad \frac{\partial (\mathbf{u} \cdot \boldsymbol{\tau})}{\partial n} = 0$$

- **Implémentation :** Sur `bottom` par exemple, $v_{i, 1} = v_{i, 0} = 0$ (imperméabilité stricte) et $u_{i, 0} = u_{i, 1}$ (cisaillement nul).
- **Usage :** Plans de symétrie axiale, parois sans couche limite pour la rupture de barrage ou l'ascension de bulle, surface supérieure libre idéalisée. Notez la graphie exacte `"FREE_SLIP"` (avec tiret bas).

### 4. Périodicité Spatiale (`"PERIODIC"`)

Relie la frontière Ouest à la frontière Est (`left` et `right` en `"PERIODIC"`) ou Sud à Nord (`bottom` et `top` en `"PERIODIC"`). Les communications MPI connectent directement les sous-domaines extrêmes via `periods` dans `MPI_Cart_create` (ou copient directement sur 1 rang), et la migration des marqueurs CVI applique le repliement géométrique $\pm L_x$ ou $\pm L_y$.
- **Usage :** Nappe liquide cisaillée de Kelvin-Helmholtz (`LIQUID_SHEET`), instabilités de Rayleigh-Taylor périodiques latéralement, canaux infinis.

---

## Conditions aux Limites Thermiques et Comportement des Marqueurs

### Conditions Thermiques (`T_type`, `T`)

Indépendamment de la condition hydrodynamique, chaque bord supporte `T_type` :

| `T_type` | Équation Continue | Implémentation Fantôme (`apply_thermal_boundary_conditions`) |
| :--- | :--- | :--- |
| `"DIRICHLET"` | $T = T_{\text{bc}}$ sur $\partial\Omega$ | $T_{\text{ghost}} = 2 \, T_{\text{bc}} - T_{\text{interior}}$ (précision ordre 2) |
| `"NEUMANN"` (défaut) | $\frac{\partial T}{\partial n} = 0$ (adiabatique) | $T_{\text{ghost}} = T_{\text{interior}}$ |
| `"PERIODIC"` (C++) | $T(0,y) = T(L_x, y)$ | Synchronisation directe MPI ou copie si `dims == 1` |

### Traitement des Frontières pour les Marqueurs CVI (`InterfaceMarkers`)

Lors de l'advection RK2-Heun (`advect_rk2`), les positions prédites $\mathbf{x}_p^*$ sont traitées de manière cohérente avec les frontières physiques :
- Sur bord `"PERIODIC"` : repliement modulo $L_x$ ou $L_y$ dans $[0, L_x[ \times [0, L_y[$.
- Sur bord physique non périodique : maintien à $\epsilon_b \Delta x$ à l'intérieur des parois lors de la sous-étape d'évaluation de vitesse (réflexion sur parois `FREE_SLIP`, suppression ou sortie en zone pure hors domaine, et ré-ensemencement continu à la buse pour `JET_INLET`).

---

## Modélisation Géométrique des Obstacles Immergés (`obstacles`)

COVOF 2D représente les corps solides ou milieux poreux fixes par le tableau JSON `obstacles` (C++ supporte un nombre quelconque d'obstacles dans `std::vector<Obstacle>` ; Fortran lit le premier élément du tableau `"obstacles": [...]` via `json_get_obstacle`). Deux formes géométriques exactes sont implémentées :

### 1. Obstacle Rectangulaire (`"BOX"`)

Défini par les bornes $[x_0, x_1] \times [y_0, y_1]$ (`x0`, `y0`, `x1`, `y1`) et sa perméabilité `permeability` :

$$\Omega_s = \left\{ (x, y) \in \Omega \;\middle|\; x_0 - 10^{-9} \le x \le x_1 + 10^{-9} \;\text{ et }\; y_0 - 10^{-9} \le y \le y_1 + 10^{-9} \right\}$$

Usage typique : marche descendante (*backward-facing step*), chicane, plaque mince verticale, bloc poreux rectangulaire.

### 2. Obstacle Cylindrique Circulaire (`"CYLINDER"`)

Défini par son centre $(x_c, y_c)$ (`xc`, `yc`), son rayon $R$ (`radius`) et sa perméabilité `permeability` :

$$\Omega_s = \left\{ (x, y) \in \Omega \;\middle|\; (x - x_c)^2 + (y - y_c)^2 \le R^2 \right\}$$

Usage typique : écoulement autour d'un cylindre rigide (allée de von Kármán à $Re = 100$) ou matrice cylindrique poreuse.

Dans les deux cas, l'appartenance est testée indépendamment au centre de chaque face $u$ (`is_u_face_in_obstacle`), au centre de chaque face $v$ (`is_v_face_in_obstacle`) et au centre de cellule pour l'échelle AAL (`is_cell_in_obstacle`). Si `permeability <= 1.0e-12`, la face est traitée comme un **solide rigide adhérent** ($u = 0$ ou $v = 0$) ; si `permeability > 1.0e-12`, le terme de traînée de Darcy-Brinkman $\frac{\mu}{\mathcal{K}}\Delta x \Delta y$ est ajouté à la matrice et à la diagonale du préconditionneur.
