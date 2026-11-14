# Structures de Données, Mémoire et Parallélisme HPC

La conception des structures de données internes d'un solveur CFD destiné aux supercalculateurs modernes requiert un compromis rigoureux entre l'expressivité algorithmique et l'efficacité matérielle. À l'échelle d'un nœud de calcul doté de plusieurs dizaines de cœurs CPU et de mémoires hiérarchiques (caches L1, L2, L3 et contrôleurs DDR/HBM NUMA), la performance brute est gouvernée par le débit mémoire (*memory bandwidth bound*), la localité spatiale et temporelle des caches, et l'alignement des données pour les registres vectoriels SIMD (AVX-512, ARM SVE).

Ce chapitre détaille l'organisation mémoire interne, la gestion des cellules fantômes, la décomposition de domaine distribuée MPI et les optimisations sans allocation (*zero-allocation*) mises en œuvre dans COVOF 2D.

---

## 3.1. Représentation Interne du Maillage et des Champs

### 3.1.1. Topologie de Grille MAC Décalée (*Marker-And-Cell*)

La grille MAC décalée répartit les variables hydrodynamiques sur différentes entités géométriques de la maille afin d'éviter tout découplage en damier (*checkerboard effect*) de la pression :
- **Centres de maille $(i, j)$ :** Stockent les variables scalaires : pression $p_{i,j}$, fraction volumique $C_{i,j}$, température $T_{i,j}$, masse volumique $\rho_{i,j}$, viscosité dynamique $\mu_{i,j}$, porosité $\chi_{T,i,j}$, et résidu de compressibilité $r_{\text{aal},i,j}$.
- **Faces verticales $(i \pm \frac{1}{2}, j)$ :** Stockent la composante horizontale de vitesse $u_{i-\frac{1}{2}, j}$ et les flux de masse horizontaux.
- **Faces horizontales $(i, j \pm \frac{1}{2})$ :** Stockent la composante verticale de vitesse $v_{i, j-\frac{1}{2}}$ et les flux de masse verticaux.
- **Nœuds / Sommets $(i \pm \frac{1}{2}, j \pm \frac{1}{2})$ :** Stockent la viscosité dynamique moyennée harmoniquement $\mu_{xy, i+\frac{1}{2}, j+\frac{1}{2}}$, le gradient croisé de vitesse et la vorticité $\omega_{z}$.

```
                 (i, j+1/2) : v(i, j+1/2)
            +-------------------------------+
            |                               |
            |                               |
 (i-1/2, j) |             (i, j)            | (i+1/2, j)
   u(i-1/2) |        p, C, T, rho, mu       |   u(i+1/2)
            |                               |
            |                               |
            +-------------------------------+
                 (i, j-1/2) : v(i, j-1/2)
   (i-1/2, j-1/2) : mu_xy, omega
```

### 3.1.2. Choix de Conception : Structure de Tableaux (SoA) vs Tableau de Structures (AoS)

Pour le stockage des champs discrets, deux paradigmes s'opposent :
1. **AoS (*Array of Structures*) :** Une structure unique `Cell` contenant tous les scalaires et composantes vectorielles, regroupée dans un tableau global `std::vector<Cell> cells`. Ce modèle favorise l'encapsulation orientée objet mais détruit la performance vectorielle lors de l'application des pochoirs de diffusion (les octets consécutifs en mémoire appartiennent à des variables différentes, saturant inutilement les lignes de cache L1).
2. **SoA (*Structure of Arrays*) :** Chaque grandeur physique est stockée dans son propre tableau contigu 1D indépendant en mémoire vive.

**Dans COVOF 2D, l'architecture SoA est strictement adoptée** pour les performances HPC :
- Côté C++17 : la structure `FlowFields` rassemble des vecteurs contigus `std::vector<double> u, v, p, rho, mu, c, T, chi_T, r_aal`.
- Côté Fortran 2018 : le type dérivé `flow_fields_t` encapsule des tableaux allouables 2D contigus `real(c_double), allocatable :: u(:, :), v(:, :), p(:, :)`.

### 3.1.3. Alignement Mémoire et Vectorisation SIMD (AVX2 / AVX-512)

Les processeurs x86_64 actuels exploitent des registres vectoriels SIMD de 256 bits (AVX2, traitant 4 réels 64 bits par cycle) ou 512 bits (AVX-512, traitant 8 réels 64 bits par cycle). Lorsque l'adresse de début d'un tableau est un multiple de 64 octets, les instructions de chargement vectoriel aligné (`_mm512_load_pd` / `vmovapd`) s'exécutent avec une latence minimale, évitant le franchissement de frontières de lignes de cache (*cache-line splits*).

- En C++17 :
  ```cpp
  // Allocateur aligné sur 64 octets
  template <typename T>
  struct AlignedAllocator {
      using value_type = T;
      T* allocate(std::size_t n) {
          void* ptr = nullptr;
          if (posix_memalign(&ptr, 64, n * sizeof(T)) != 0) throw std::bad_alloc();
          return static_cast<T*>(ptr);
      }
      void deallocate(T* p, std::size_t) noexcept { free(p); }
  };
  using AlignedVectorDouble = std::vector<double, AlignedAllocator<double>>;
  ```

- En Fortran 2018 :
  L'utilisation des directives du compilateur garantit l'alignement et force la vectorisation :
  ```fortran
  !DIR$ ATTRIBUTES ALIGN : 64 :: u, v, p
  !DIR$ ASSUME_ALIGNED u: 64, v: 64, p: 64
  !$omp simd aligned(u, v, p : 64)
  do i = 1, nx
      p(i, j) = p(i, j) - dt * rho(i, j) * div(i, j)
  end do
  ```

---

## 3.2. Gestion des Cellules Fantômes (*Halo / Ghost Cells*)

Afin de permettre l'évaluation locale des opérateurs différentiels spatiaux (gradients, divergences, laplaciens) sans test conditionnel de bordure dans les boucles de calcul critiques, chaque sous-domaine local est bordé d'une couche périphérique de cellules fictives appelées **cellules fantômes** (*ghost cells* ou *halos*).

### 3.2.1. Épaisseur des Couches Fantômes ($N_{\text{ghost}}$)

L'épaisseur requise $N_{\text{ghost}}$ dépend directement de la largeur du pochoir numérique (*stencil*) le plus étendu utilisé dans la simulation :

| Opérateur Numérique | Pochoir Spatial | Épaisseur Requise $N_{\text{ghost}}$ | Remarques |
| :--- | :--- | :---: | :--- |
| **Diffusion Visqueuse d'ordre 2** | Pochoir 5 points compact | 1 couche | Évaluation standard de $\Delta u$ et $\Delta v$. |
| **Convection Centrée / Upwind ordre 1** | Pochoir 3 points | 1 couche | Calcul standard des flux aux faces. |
| **Convection QUICK / Adams-Bashforth 2** | Pochoir décentré amont 3 points | 2 couches | Nécessite la valeur de deux mailles en amont. |
| **VOF PLIC (Youngs / ELVIRA)** | Pochoir $3 \times 3$ mailles | 2 couches | Calcul des gradients de fraction volumique $\nabla C$. |
| **Tension de Surface CSF ($\kappa$)** | Pochoir $5 \times 5$ mailles | 2 couches | Dérivation seconde pour la courbure de surface $\kappa = -\nabla \cdot \mathbf{n}$. |

Dans COVOF 2D, l'épaisseur standard par défaut est configurée à **$N_{\text{ghost}} = 2$ couches** de mailles tout autour du domaine local, ce qui autorise sans réallocation l'exécution de l'ensemble des schémas d'ordre élevé et des traceurs d'interface.

### 3.2.2. Topologie d'Indexation Locale

Pour un sous-domaine de dimension utile $N_x^{\text{loc}} \times N_y^{\text{loc}}$ mailles, la dimension totale mémoire du tableau alloué est :
$$N_{x,\text{total}} = N_x^{\text{loc}} + 2 N_{\text{ghost}}, \quad N_{y,\text{total}} = N_y^{\text{loc}} + 2 N_{\text{ghost}}$$

L'espace indiciel se décompose ainsi :
- **Zone utile (cœur de calcul) :**
  - C++ (0-based) : $i \in [N_{\text{ghost}}, N_{\text{ghost}} + N_x^{\text{loc}} - 1]$, $j \in [N_{\text{ghost}}, N_{\text{ghost}} + N_y^{\text{loc}} - 1]$.
  - Fortran (1-based) : $i \in [N_{\text{ghost}} + 1, N_{\text{ghost}} + N_x^{\text{loc}}]$, $j \in [N_{\text{ghost}} + 1, N_{\text{ghost}} + N_y^{\text{loc}}]$.
- **Cellules fantômes Ouest / Est :** $i \in [0, N_{\text{ghost}} - 1]$ et $i \in [N_{\text{ghost}} + N_x^{\text{loc}}, N_{x,\text{total}} - 1]$ (en C++).
- **Cellules fantômes Sud / Nord :** $j \in [0, N_{\text{ghost}} - 1]$ et $j \in [N_{\text{ghost}} + N_y^{\text{loc}}, N_{y,\text{total}} - 1]$ (en C++).

---

## 3.3. Parallélisme Distribué MPI : Topologie Cartésienne et Halos

Le passage à l'échelle massive sur supercalculateur repose sur le découpage du domaine spatial global $\Omega$ en une grille bidimensionnelle de $P_x \times P_y = P$ sous-domaines réguliers.

### 3.3.1. Création de la Topologie Virtuelle (`MPI_Cart_create`)

COVOF 2D utilise les topologies cartésiennes communicantes standard MPI. La fonction `MPI_Cart_create` permet au gestionnaire MPI matériel (sur réseaux InfiniBand, Omni-Path ou Slingshot) d'optimiser le placement des processus sur les processeurs physiques afin de minimiser les sauts réseau (*network hops*) :

```cpp
// C++17 : Initialisation cartésienne dans src/MpiComm.cpp
int dims[2] = {px, py};
int periods[2] = {periodic_x ? 1 : 0, periodic_y ? 1 : 0};
int reorder = 1;

MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cart_comm_);
MPI_Comm_rank(cart_comm_, &cart_rank_);
MPI_Cart_coords(cart_comm_, cart_rank_, 2, coords_);

// Recherche des 4 voisins immédiats
MPI_Cart_shift(cart_comm_, 0, 1, &rank_west_, &rank_east_);
MPI_Cart_shift(cart_comm_, 1, 1, &rank_south_, &rank_north_);
```

### 3.3.2. Communications Point-à-Point Non Bloquantes et Recouvrement

Pour transférer les mailles fantômes entre sous-domaines adjacents sans bloquer les cœurs de calcul, COVOF 2D utilise des requêtes non bloquantes `MPI_Isend` et `MPI_Irecv`.

Afin d'éviter la copie intermédiaire manuelle des colonnes de mailles (non contiguës en mémoire C++), un type dérivé de vecteur MPI (`MPI_Type_vector`) est créé pour les frontières Est et Ouest :

```cpp
// Création du type colonne pour le halo vertical
MPI_Datatype col_type;
MPI_Type_vector(ny_loc, nghost, nxt, MPI_DOUBLE, &col_type);
MPI_Type_commit(&col_type);
```

Le protocole d'échange synchrone sécurisé se déroule en deux étapes :
1. **Échange horizontal (Ouest $\leftrightarrow$ Est) :** Envoi des bordures utiles aux voisins et réception dans les halos correspondants.
2. **Échange vertical (Sud $\leftrightarrow$ Nord) :** Envoi des bordures utiles incluant les halos Ouest/Est déjà mis à jour. Cette séquence sérielle garantit le remplissage automatique et exact des quatre coins (*corner cells*) sans requêtes diagonales supplémentaires.

```
 Étape 1 : Halos Verticaux              Étape 2 : Halos Horizontaux
      (Ouest <--> Est)                      (Sud <--> Nord)
    +---+-----------+---+                 +---+-----------+---+
    | H |           | H |                 | H |     H     | H | <-- Coins reçus
    | a |  Domaine  | a |                 +---+-----------+---+
    | l |   Utile   | l |                 |   |  Domaine  |   |
    | o |   Local   | o |    ------>      |   |   Utile   |   |
    |   |           |   |                 +---+-----------+---+
    | W |           | E |                 | H |     H     | H | <-- Coins reçus
    +---+-----------+---+                 +---+-----------+---+
```

---

## 3.4. Élimination des Allocations Dynamiques dans les Boucles de Krylov

L'une des causes majeures de chute de performance dans les codes C++ orientés objet réside dans l'appel intempestif à l'allocateur mémoire du système (`malloc` / `free`, `new` / `delete`) à l'intérieur de la boucle temporelle ou au fil des itérations d'un solveur de Krylov. Chaque réallocation déclenche des verrous système (*kernel locks*) et fragmente le tas (*heap*).

### 3.4.1. Protocole de Préallocation Statique (*Zero-Allocation Guarantee*)

Dans COVOF 2D :
1. **Toutes les structures de travail sont allouées une seule fois lors de la phase d'initialisation (`Simulation::init()`).**
2. Les solveurs itératifs (BiCGSTAB2, BiCGSTAB, GMRES, MINRES, CG) disposent d'un ensemble de vecteurs de travail temporaires pré-dimensionnés membres de classe :
   - `std::vector<double> r0_, r_, p_, v_, s_, t_, z_, y_;` (pour BiCGSTAB)
   - `std::vector<double> v_p1_, v_p2_, w_, p1_, p2_;` (pour MINRES)
3. La routine d'application de l'opérateur couplé `ns_apply_system_matrix` reçoit impérativement les buffers de sortie préalloués en paramètres de référence :
   ```cpp
   void NavierStokesSolver::applySystemMatrix(const double* __restrict__ x_in,
                                              double* __restrict__ Ax_out);
   ```
4. L'utilisation du mot-clé `__restrict__` informe le compilateur C++ ou Fortran que les pointeurs ne se chevauchent pas en mémoire (*no aliasing*), autorisant la réorganisation optimale des instructions de chargement et de stockage vectoriel.

---

## 3.5. Structure de Stockage Creux CSR (*Compressed Sparse Row*)

Bien que les opérateurs de Navier-Stokes soient principalement évalués sans matrice (*matrix-free*), certaines étapes (telles que le solveur direct LU LAPACK ou les préconditionneurs algébriques ILU0) requièrent l'assemblage explicite de la matrice d'advection-diffusion-pression sous format creux CSR.

### 3.5.1. Représentation CSR à Trois Tableaux

Une matrice creuse $A \in \mathbb{R}^{M \times N}$ comportant $N_{nz}$ éléments non nuls est représentée par trois tableaux contigus :
1. `values` (taille $N_{nz}$, type `double`) : Stocke la valeur numérique de chaque entrée non nulle parcourue ligne par ligne.
2. `col_indices` (taille $N_{nz}$, type `int`) : Stocke l'indice de colonne $j \in [0, N-1]$ correspondant à chaque valeur.
3. `row_offsets` (taille $M+1$, type `int`) : Stocke l'indice de départ dans les tableaux `values` et `col_indices` de la première entrée non nulle de chaque ligne $i$. Le nombre d'éléments non nuls sur la ligne $i$ est donné par $\text{row\_offsets}[i+1] - \text{row\_offsets}[i]$.

### 3.5.2. Produit Matrice-Vecteur Creux SpMV Parallèle

Le produit matrice creuse - vecteur plein $\mathbf{y} = A \cdot \mathbf{x}$ s'écrit de manière compacte et hautement parallélisable en OpenMP :

```cpp
void CSRMatrix::spmv(const double* __restrict__ x, double* __restrict__ y) const {
    #pragma omp parallel for schedule(guided)
    for (int i = 0; i < num_rows_; ++i) {
        double sum = 0.0;
        const int row_start = row_offsets_[i];
        const int row_end   = row_offsets_[i + 1];
        #pragma omp simd reduction(+:sum)
        for (int k = row_start; k < row_end; ++k) {
            sum += values_[k] * x[col_indices_[k]];
        }
        y[i] = sum;
    }
}
```

La gémellité stricte de cette routine est assurée en Fortran 2018 au sein du module `module_csr_matrix.f90` via la sous-routine `csr_spmv_apply`, garantissant l'équivalence parfaite des calculs de résidus sur les maillages non structurés ou les matrices couplées assemblées.
