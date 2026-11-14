# Géométrie, Maillages et Conditions aux Limites

## Génération et Topologie du Maillage Cartésien MAC

Le solveur **COVOF 2D** opère sur des domaines géométriques rectangulaires discrétisés par un maillage cartésien structuré uniforme selon la topologie décalée MAC (*Marker-and-Cell*).

### Définition Géométrique du Domaine

Le domaine physique de calcul $\Omega = [x_{\text{min}}, x_{\text{max}}] \times [y_{\text{min}}, y_{\text{max}}]$ est caractérisé dans le bloc `mesh` du fichier JSON par :
- Les dimensions physiques : $L_x = x_{\text{max}} - x_{\text{min}}$ (`mesh.lx`) et $L_y = y_{\text{max}} - y_{\text{min}}$ (`mesh.ly`) $[\text{m}]$.
- L'origine des coordonnées cartésiennes : $(x_{\text{min}}, y_{\text{min}})$ (`mesh.x_min`, `mesh.y_min`) $[\text{m}]$ (défaut : $0.0, 0.0$).
- Le nombre total de subdivisions spatiales globales : $N_x$ mailles selon l'axe horizontal $X$ (`mesh.nx`) et $N_y$ mailles selon l'axe vertical $Y$ (`mesh.ny`).

Les pas d'espace constants du maillage sont calculés immédiatement selon :

$$\Delta x = \frac{L_x}{N_x}, \quad \Delta y = \frac{L_y}{N_y}$$

### Coordonnées Spatiales des Nœuds, Faces et Centres

Pour une maille repérée par les indices discrets globaux $(i, j)$ où $i \in [1, N_x]$ et $j \in [1, N_y]$ :
- **Coordonnées du centre de maille (Pression, VOF, Température) :**
  $$x_{c, i} = x_{\text{min}} + \left( i - \frac{1}{2} \right) \Delta x, \quad y_{c, j} = y_{\text{min}} + \left( j - \frac{1}{2} \right) \Delta y$$
- **Coordonnées des faces verticales (Vitesse $u$) :**
  $$x_{f, i} = x_{\text{min}} + (i - 1) \Delta x, \quad y_{c, j} = y_{\text{min}} + \left( j - \frac{1}{2} \right) \Delta y$$
  Il y a au total $(N_x + 1) \times N_y$ faces de vitesse horizontale dans le domaine complet.
- **Coordonnées des faces horizontales (Vitesse $v$) :**
  $$x_{c, i} = x_{\text{min}} + \left( i - \frac{1}{2} \right) \Delta x, \quad y_{f, j} = y_{\text{min}} + (j - 1) \Delta y$$
  Il y a au total $N_x \times (N_y + 1)$ faces de vitesse verticale dans le domaine complet.

---

## Décomposition de Domaine et Parallélisme Cartésien 2D (MPI)

Pour permettre l'exécution sur architectures distribuées (grappes de calcul et supercalculateurs), COVOF 2D partitionne le maillage global selon une grille 2D de processus MPI communicants.

```
       +-----------------------+-----------------------+
       |                       |                       |
       |     Rang MPI (0, 1)   |     Rang MPI (1, 1)   |
       |     Sous-domaine      |     Sous-domaine      |
       |                       |                       |
Ny / 2 +-----------------------+-----------------------+
       |                       |                       |
       |     Rang MPI (0, 0)   |     Rang MPI (1, 0)   |
       |     Sous-domaine      |     Sous-domaine      |
       |                       |                       |
       +-----------------------+-----------------------+
                   0                    Nx / 2                  Nx
```

### Topologie Virtuelle Cartésienne (`MPI_Cart_create`)

1. **Calcul automatique des dimensions de partitionnement :**
   À partir du nombre total de processus alloués $N_{\text{procs}}$, COVOF 2D détermine la décomposition optimale en $P_x$ tranches horizontales et $P_y$ tranches verticales via `MPI_Dims_create(n_procs, 2, dims)` :
   $$N_{\text{procs}} = P_x \times P_y \quad \text{avec} \quad \frac{P_x}{P_y} \approx \frac{N_x}{N_y}$$
2. **Identification des rangs voisins :**
   Grâce à `MPI_Cart_shift`, chaque processus détermine les rangs de ses quatre voisins directs dans la grille :
   - Voisin Ouest (`rank_west`) à la coordonnée $(p_x - 1, p_y)$
   - Voisin Est (`rank_east`) à la coordonnée $(p_x + 1, p_y)$
   - Voisin Sud (`rank_south`) à la coordonnée $(p_x, p_y - 1)$
   - Voisin Nord (`rank_north`) à la coordonnée $(p_x, p_y + 1)$
   Si un processus est situé sur une frontière physique sans condition périodique, le rang voisin correspondant est positionné à `MPI_PROC_NULL`, désactivant automatiquement les transferts réseau sans conditionnement lourd.

### Gestion des Cellules Fantômes (*Halos / Ghost Cells*)

Chaque sous-domaine local alloue une couche périphérique de cellules de garde :
- `nx_local` $\times$ `ny_local` : mailles intérieures possédées en propre par le processus.
- `nx_ghost = nx_local + 2` et `ny_ghost = ny_local + 2` : taille mémoire totale incluant les cellules fantômes d'indice $0$ et `nx_local + 1`.

À chaque étape de calcul nécessitant des dérivées spatiales ou des flux aux interfaces de sous-domaines :
1. Les données des mailles intérieures de bordure sont emballées dans des tampons contigus (*packing*).
2. Les communications asynchrones non bloquantes `MPI_Isend` et `MPI_Irecv` sont initiées vers les voisins E, W, N, S.
3. Après synchronisation par `MPI_Waitall`, les valeurs reçues sont déballées dans les cellules fantômes (*unpacking*).

---

## Catalogue Exhaustif des Conditions aux Limites (CL)

COVOF 2D dispose d'un catalogue complet de conditions aux limites physiques et numériques applicables indépendamment sur les quatre frontières cardinales du domaine : `left` (Ouest, $x = x_{\text{min}}$), `right` (Est, $x = x_{\text{max}}$), `bottom` (Sud, $y = y_{\text{min}}$) et `top` (Nord, $y = y_{\text{max}}$).

### 1. Vitesse Imposée (`DIRICHLET`)

La condition de Dirichlet prescrit explicitement la valeur vectorielle de la vitesse à la frontière :

$$\mathbf{u} = \mathbf{u}_{\text{wall}} = (u_{\text{bc}}, v_{\text{bc}})^T \quad \text{sur } \partial \Omega$$

#### Implémentation Numérique Décalée
- **Pour la composante normale :**
  La vitesse normale coïncide géométriquement avec la face frontière. Elle est imposée directement sans interpolation :
  $$u_{1/2, j} = u_{\text{bc}} \quad \text{(sur frontière gauche)}$$
- **Pour la composante tangentielle :**
  La vitesse tangentielle $v$ est localisée à une distance $\frac{\Delta x}{2}$ de la paroi. Pour imposer la valeur $v_{\text{bc}}$ sur la frontière physique avec une précision d'ordre 2, la cellule fantôme d'indice $i=0$ est extrapolée linéairement :
  $$\frac{v_{0, j+1/2} + v_{1, j+1/2}}{2} = v_{\text{bc}} \implies v_{0, j+1/2} = 2 v_{\text{bc}} - v_{1, j+1/2}$$

#### Cas d'Usage Typiques
- **Paroi solide immobile sans glissement (*No-Slip Wall*) :** $u_{\text{bc}} = 0, v_{\text{bc}} = 0$.
- **Paroi mobile de cavité entraînée (*Lid-Driven Cavity*) :** sur la frontière supérieure, $u_{\text{bc}} = U_{\text{lid}} = 1.0\,\text{m}\cdot\text{s}^{-1}, v_{\text{bc}} = 0$.
- **Profil d'entrée parabolique de Poiseuille :** $u(y) = 4 U_{\max} \frac{y (L_y - y)}{L_y^2}$.

### 2. Gradient Nul / Sortie Libre (`NEUMANN`)

La condition de Neumann prescrit un flux nul (dérivée normale nulle) des grandeurs de vitesse à travers la frontière :

$$\frac{\partial \mathbf{u}}{\partial \mathbf{n}} = \mathbf{0} \quad \iff \quad \frac{\partial u}{\partial n} = 0, \quad \frac{\partial v}{\partial n} = 0$$

#### Implémentation Numérique Décalée
Sur la frontière droite ($x = x_{\text{max}}$) :
$$u_{N_x+1/2, j} = u_{N_x-1/2, j}, \quad v_{N_x+1, j+1/2} = v_{N_x, j+1/2}$$
Pour la pression, la condition de Neumann homogène $\frac{\partial p}{\partial n} = 0$ est naturellement imposée pour assurer la compatibilité avec l'équation de quantité de mouvement.

#### Cas d'Usage Typiques
- **Sortie libre convective (*Outflow / Open Boundary*) :** Permet à l'écoulement et aux tourbillons de s'évacuer du domaine de calcul sans réflexion d'ondes numériques.

### 3. Périodicité Spatiale (`PERIODIC`)

La condition de périodicité relie topologiquement les frontières opposées du domaine pour simuler un milieu infini ou un écoulement en boucle fermée :

$$\mathbf{u}(x_{\text{min}}, y) = \mathbf{u}(x_{\text{max}}, y), \quad p(x_{\text{min}}, y) = p(x_{\text{max}}, y)$$

#### Implémentation Numérique Décalée
Les cellules fantômes de la frontière Ouest reçoivent les données des dernières cellules intérieures de la frontière Est, et réciproquement. Dans le module MPI, le communicateur cartésien active le drapeau de périodicité `periods = [true, false]`, de sorte que les processus de rang $0$ communiquent directement avec les processus de rang $P_x - 1$.

#### Cas d'Usage Typiques
- **Vortex de Taylor-Green et instabilités de cisaillement.**
- **Canaux infinis sous gradient de pression moteur.**

### 4. Glissement Parfait sans Frottement (`FREESLIP`)

La condition de glissement parfait impose une vitesse normale nulle (imperméabilité) et un cisaillement tangentiel nul (aucun frottement visqueux) :

$$\mathbf{u} \cdot \mathbf{n} = 0, \quad \frac{\partial (\mathbf{u} \times \mathbf{n})}{\partial \mathbf{n}} = \mathbf{0}$$

#### Implémentation Numérique Décalée
Sur la frontière inférieure ($y = y_{\text{min}}$, normale selon $y$) :
- Vitesse normale nulle sur la face frontière : $v_{i, 1/2} = 0$.
- Extrapolation symétrique de la vitesse tangentielle dans la maille fantôme :
  $$u_{i+1/2, 0} = u_{i+1/2, 1} \implies \left( \frac{\partial u}{\partial y} \right)_{y=y_{\text{min}}} = 0$$

#### Cas d'Usage Typiques
- **Plans de symétrie géométrique et hydrodynamique.**
- **Surfaces libres idéalisées et fonds marins sans couche limite visqueuse.**

---

## Conditions aux Limites Thermiques et VOF

En complément des conditions hydrodynamiques, les frontières physiques accueillent les conditions sur les scalaires passifs et actifs :

| Condition aux Limites | Mot-Clé JSON | Équation Continue | Implémentation Fantôme |
| :--- | :--- | :--- | :--- |
| **Température Imposée** | `THERMAL_DIRICHLET` | $T = T_{\text{wall}}$ | $T_{\text{ghost}} = 2 T_{\text{wall}} - T_{\text{interior}}$ |
| **Paroi Adiabatique** | `THERMAL_NEUMANN` | $\frac{\partial T}{\partial n} = 0$ | $T_{\text{ghost}} = T_{\text{interior}}$ |
| **Interface Imperméable** | `VOF_WALL` | $\mathbf{u} \cdot \mathbf{n}_{\text{wall}} = 0$ | Flux VOF nul à travers la face |
| **Sortie de Phase Libre** | `VOF_NEUMANN` | $\frac{\partial C}{\partial n} = 0$ | $C_{\text{ghost}} = C_{\text{interior}}$ |

---

## Modélisation Géométrique des Obstacles Immergés

COVOF 2D intègre un générateur géométrique interne d'obstacles solides immergés pénalisés par la méthode de Darcy-Brinkman (`obstacles.enable_obstacles = true`). Deux formes canoniques sont disponibles :

### 1. Obstacle Cylindrique (`CYLINDER`)

L'obstacle est défini par le centre de son disque $(x_c, y_c)$ (`obs_xc`, `obs_yc`) et son rayon $R$ (`obs_radius`) :

$$\Omega_s = \left\{ (x, y) \in \Omega \;\middle|\; (x - x_c)^2 + (y - y_c)^2 \le R^2 \right\}$$

- Pour chaque face verticale $u_{i+1/2, j}$, le test d'appartenance `is_u_face_in_obstacle` vérifie si le centre de face $(x_{f, i+1/2}, y_{c, j})$ est à l'intérieur du cylindre.
- Si le point est immergé et que la perméabilité `obs_perm` $\le 10^{-12}$, la composante de vitesse est directement annulée ($u = 0$, Dirichlet solide) ou pénalisée par le coefficient de Brinkman $\frac{\mu}{K}$.

### 2. Obstacle Rectangulaire / Plaque (`RECTANGLE`)

L'obstacle est défini par sa boîte englobante $[x_0, x_1] \times [y_0, y_1]$ (`obs_x0`, `obs_x1`, `obs_y0`, `obs_y1`) :

$$\Omega_s = \left\{ (x, y) \in \Omega \;\middle|\; x_0 \le x \le x_1 \quad \text{et} \quad y_0 \le y \le y_1 \right\}$$

Cette géométrie permet de représenter facilement des marches descendantes (*backward-facing step*), des chicanes, des parois internes minces ou des piliers carrés.
