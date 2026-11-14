# Interopérabilité C++ $\leftrightarrow$ Fortran

L'un des défis majeurs dans le développement de codes de calcul scientifique modernes réside dans la cohabitation harmonieuse entre des langages à haut niveau d'abstraction (tels que C++17) et des langages hautement optimisés pour le calcul matriciel intensif (tels que Fortran 2018). Dans l'écosystème COVOF 2D, l'interopérabilité ne constitue pas un simple mécanisme d'appel ponctuel, mais le pivot assurant la stricte conformité sémantique, la parité des résultats à la précision machine près ($\epsilon_{\text{mach}} \approx 10^{-16}$) et la portabilité des bibliothèques de calcul haute performance (HPC).

Ce chapitre explicite l'ensemble des mécanismes formels, des normes de liaison binaire (ABI), des conventions d'ordonnancement mémoire et des protocoles de communication inter-langages mis en œuvre au sein du solveur.

---

## 2.1. Standards et Conventions d'Interface : `iso_c_binding` et `extern "C"`

Le standard Fortran 2003 (consolidé en Fortran 2008 et 2018) a introduit le module intrinsèque `iso_c_binding`. Ce module normalise formellement l'interopérabilité binaire directe entre Fortran et le langage C standard (ISO C). Dès lors que le langage C++17 utilise la directive `extern "C"`, l'incompatibilité historique due à la décoration des noms de fonctions C++ (*name mangling*) est complètement éliminée, permettant un appel direct sans coût d'indirection à la vitesse native du processeur.

### 2.1.1. Directive `extern "C"` côté C++17

En C++, le polymorphisme, les espaces de nommage (*namespaces*) et la surcharge des fonctions modifient le symbole binaire exporté dans la table des symboles du fichier objet (`.o`). La directive `extern "C"` ordonne au compilateur d'adopter la convention d'appel et la table des symboles C standard :

```cpp
// include/covof/InteropBindings.hpp
#pragma once
#include <cstddef>
#include <cstdint>

#ifdef __cplusplus
extern "C" {
#endif

// Prototypes des passerelles C exportées vers Fortran
void covof_cpp_init_mpi(int* f_comm, int* rank, int* size, int* ierr);
void covof_cpp_apply_pressure_laplacian(const double* p_in, double* p_out,
                                        const int* nx, const int* ny,
                                        const double* dx, const double* dy,
                                        int* ierr);

#ifdef __cplusplus
}
#endif
```

### 2.1.2. Module Intrinsèque `iso_c_binding` côté Fortran 2018

Côté Fortran, chaque sous-programme ou fonction destiné à franchir la frontière inter-langage est enrichi de l'attribut `bind(C, name="...")`. Cet attribut garantit que le compilateur Fortran n'ajoute aucun soulignement supplémentaire (*underscoring*) et applique la convention d'appel standard C (registres d'arguments et nettoyage de la pile conformes à l'ABI System V AMD64 / ARM64 AAPCS) :

```fortran
! src/module_interop_bindings.f90
module module_interop_bindings
    use, intrinsic :: iso_c_binding
    implicit none
    private

    public :: covof_fortran_apply_diffusion

contains

    subroutine covof_fortran_apply_diffusion(u, lap_u, nx, ny, dx, dy, ierr) &
               bind(C, name="covof_fortran_apply_diffusion")
        integer(c_int), intent(in), value               :: nx, ny
        real(c_double), intent(in), value               :: dx, dy
        real(c_double), dimension(nx, ny), intent(in)   :: u
        real(c_double), dimension(nx, ny), intent(out)  :: lap_u
        integer(c_int), intent(out)                     :: ierr

        integer(c_int) :: i, j
        real(c_double) :: idx2, idy2

        ierr = 0
        idx2 = 1.0_c_double / (dx * dx)
        idy2 = 1.0_c_double / (dy * dy)

        !$omp parallel do collapse(2) private(i, j)
        do j = 2, ny - 1
            do i = 2, nx - 1
                lap_u(i, j) = (u(i+1, j) - 2.0_c_double * u(i, j) + u(i-1, j)) * idx2 + &
                              (u(i, j+1) - 2.0_c_double * u(i, j) + u(i, j-1)) * idy2
            end do
        end do
        !$omp end parallel do
    end subroutine covof_fortran_apply_diffusion

end module module_interop_bindings
```

---

## 2.2. Table Complète de Correspondance des Types et Structures

Pour éviter tout écrasement mémoire ou décalage d'octets (*memory misalignment / padding*), chaque type élémentaire manipulé dans COVOF 2D dispose d'une correspondance biunivoque stricte.

### 2.2.1. Types Scalaires Primitifs

Le tableau ci-dessous recense la bijection exacte entre les types de base C++17 et Fortran 2018 :

| Grandeur Physique / Concept | Type C++17 (`<cstdint>`, `<cstddef>`) | Taille Binaire | Type Fortran 2018 (`iso_c_binding`) | Précision / Intervalle |
| :--- | :--- | :--- | :--- | :--- |
| Champ réel (vitesse, pression) | `double` | 64 bits (8 octets) | `real(c_double)` | IEEE 754 double précision ($\approx 16$ chiffres significatifs) |
| Champ simple (option basse mémoire) | `float` | 32 bits (4 octets) | `real(c_float)` | IEEE 754 simple précision ($\approx 7$ chiffres significatifs) |
| Indices de mailles et tailles locales | `int` ou `int32_t` | 32 bits (4 octets) | `integer(c_int)` | Entier signé standard ($-2 \times 10^9$ à $+2 \times 10^9$) |
| Indices de solveurs creux / NDOF global | `int64_t` | 64 bits (8 octets) | `integer(c_int64_t)` | Entier signé étendu ($-9 \times 10^{18}$ à $+9 \times 10^{18}$) |
| Drapeaux booléens et options | `bool` | 8 bits (1 octet) | `logical(c_bool)` | Booléen C (`0` ou `1`) |
| Caractère ASCII / Code de contrôle | `char` | 8 bits (1 octet) | `character(kind=c_char)` | Caractère ISO-8859 / UTF-8 |
| Pointeur générique non typé | `void*` ou `uintptr_t` | 64 bits (8 octets) | `type(c_ptr)` | Adresse mémoire virtuelle absolue |
| Pointeur de fonction (callback) | `void (*func)(...)` | 64 bits (8 octets) | `type(c_funptr)` | Pointeur d'exécution indirecte |

### 2.2.2. Structures Dérivées et Alignement Mémoire (`bind(C)`)

Lorsqu'un type composé ou une structure de données doit être transmise entre C++ et Fortran, la structure Fortran **doit impérativement** comporter l'attribut `bind(C)`, tandis que la structure C++ doit utiliser un agencement de types POD (*Plain Old Data*) exempt de fonctions virtuelles (pas de table `vtable`) et respecter l'alignement naturel des membres :

#### Définition C++17

```cpp
// include/covof/InteropTypes.hpp
#pragma once
#include <cstdint>

struct C_MeshMetrics {
    int32_t nx;
    int32_t ny;
    int32_t nghost;
    double  dx;
    double  dy;
    double  x_min;
    double  y_min;
};
```

#### Définition Fortran 2018 Correspondante

```fortran
! src/module_interop_types.f90
module module_interop_types
    use, intrinsic :: iso_c_binding
    implicit none

    type, bind(C) :: c_mesh_metrics_t
        integer(c_int32_t) :: nx
        integer(c_int32_t) :: ny
        integer(c_int32_t) :: nghost
        real(c_double)     :: dx
        real(c_double)     :: dy
        real(c_double)     :: x_min
        real(c_double)     :: y_min
    end type c_mesh_metrics_t

end module module_interop_types
```

### 2.2.3. Conversion de Pointeurs Opaques : `c_f_pointer` et `c_loc`

Dans les architectures hybrides avancées, Fortran peut allouer un bloc contigu de mémoire et le transmettre à C++ sous forme de `type(c_ptr)`, ou inversement, C++ peut transmettre un pointeur `double*` brut à Fortran. Pour manipuler ce pointeur en Fortran sous forme d'un tableau indicé multidimensionnel sans recopie mémoire, on utilise la routine intrinsèque `c_f_pointer` :

```fortran
subroutine process_cpp_buffer(raw_ptr, nx, ny) bind(C, name="process_cpp_buffer")
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), value, intent(in) :: raw_ptr
    integer(c_int), value, intent(in) :: nx, ny

    ! Pointeur Fortran local pointant sur le buffer C++
    real(c_double), pointer, dimension(:, :) :: field

    ! Établissement de la vue de tableau Fortran sans copie
    call c_f_pointer(raw_ptr, field, [nx, ny])

    ! Modification directe en mémoire partagée
    field(1, 1) = 0.0_c_double
end subroutine process_cpp_buffer
```

---

## 2.3. Gestion de la Disposition Mémoire (Row-Major vs Column-Major)

La divergence la plus critique entre C++ et Fortran réside dans la convention d'ordonnancement des éléments matriciels multidimensionnels dans l'espace d'adressage linéaire de la mémoire vive (RAM).

### 2.3.1. Théorie de l'Adressage Linéaire

Soit une matrice bidimensionnelle ou un champ scalaire discrétisé comportant $N_x$ mailles en direction horizontale et $N_y$ mailles en direction verticale.

1. **Convention C++ (Ligne majeure / *Row-Major Order*) :**
   Le deuxième indice (direction $x$, colonnes) varie le plus rapidement dans l'espace mémoire :
   $$\text{Offset}_{\text{C++}}(i, j) = j \times N_x + i$$
   Deux mailles consécutives sur la même ligne physique $(i, j)$ et $(i+1, j)$ sont situées à des adresses contiguës en mémoire ($\Delta \text{Adresse} = 8 \text{ octets}$).

2. **Convention Fortran (Colonne majeure / *Column-Major Order*) :**
   Le premier indice (direction $x$, lignes géométriques du tableau) varie le plus rapidement :
   $$\text{Offset}_{\text{Fortran}}(i, j) = (j - 1) \times N_x + (i - 1)$$
   (avec convention d'indexation 1-based par défaut). Deux mailles consécutives $(i, j)$ et $(i+1, j)$ sont contiguës en mémoire.

### 2.3.2. Règle d'Or de Parité Mémoire dans COVOF 2D

Afin de garantir que les noyaux de calcul Fortran et C++ s'exécutent avec des boucles vectorisables sans permutation de dimensions (*stride-1 access*), **la topologie mémoire dans COVOF 2D est unifiée sur la direction $x$ rapide**.

- En C++, les champs sont alloués sous forme d'un tableau unidimensionnel aplati `std::vector<double>` de taille totale $(N_x + 2 N_{\text{ghost}}) \times (N_y + 2 N_{\text{ghost}})$, accédé par la macro d'indexation en ligne majeure :
  $$\text{IDX}(i, j) = j \times N_{x,\text{total}} + i$$
  Dans une boucle de calcul, l'indice $i$ est donc placé dans la boucle la plus interne :
  ```cpp
  for (int j = 1; j <= ny; ++j) {
      #pragma omp simd
      for (int i = 1; i <= nx; ++i) {
          lap[j * nxt + i] = ...;
      }
  }
  ```

- En Fortran, le tableau 2D `field(nx_total, ny_total)` possède son premier indice $i$ aligné en colonne majeure, ce qui correspond **exactement au même agencement physique en mémoire** :
  ```fortran
  do j = 1, ny
      !$omp simd
      do i = 1, nx
          lap(i, j) = ...
      end do
  end do
  ```

| Représentation Logique | Indexation C++ | Indexation Fortran | Adresse Linéaire Relative |
| :---: | :---: | :---: | :---: |
| Maille $(i, j)$ | `field[j * nxt + i]` | `field(i + 1, j + 1)` | Base $+ (j \times N_{xt} + i) \times 8$ |
| Voisin Est $(i+1, j)$ | `field[j * nxt + (i + 1)]` | `field(i + 2, j + 1)` | Base $+ (j \times N_{xt} + i + 1) \times 8$ |
| Voisin Nord $(i, j+1)$ | `field[(j + 1) * nxt + i]` | `field(i + 1, j + 2)` | Base $+ ((j + 1) \times N_{xt} + i) \times 8$ |

Cette parfaite concordance d'agencement physique permet de passer directement l'adresse `&field[0]` du C++ au tableau Fortran `dimension(nxt, nyt)` sans transposer la moindre donnée.

---

## 2.4. Gestion des Chaînes de Caractères, Entrées/Sorties et Erreurs

### 2.4.1. Passage des Chaînes de Caractères

En langage C/C++, une chaîne de caractères est un pointeur vers un tableau d'octets terminé par le caractère nul sentinelle `\0` (`ASCII 0`). En Fortran standard, une chaîne `character(len=N)` est un bloc de longueur fixe complété par des espaces blancs (*blank padding*) sans caractère nul final.

Pour garantir l'interopérabilité sans écrasement de mémoire :
1. Toute chaîne transmise par C++ vers Fortran doit être terminée par `\0`.
2. Le sous-programme Fortran convertit le pointeur C en chaîne native grâce à une routine de normalisation :

```fortran
subroutine c_string_to_fortran(c_str_ptr, f_str)
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), intent(in)            :: c_str_ptr
    character(len=*), intent(out)       :: f_str

    character(kind=c_char), pointer    :: char_array(:)
    integer                            :: idx, max_len

    f_str = ""
    if (.not. c_associated(c_str_ptr)) return

    call c_f_pointer(c_str_ptr, char_array, [1024])
    idx = 1
    max_len = len(f_str)
    do while (idx <= max_len .and. char_array(idx) /= c_null_char)
        f_str(idx:idx) = char_array(idx)
        idx = idx + 1
    end do
end subroutine c_string_to_fortran
```

### 2.4.2. Synchronisation des Unités Logiques d'E/S

L'utilisation simultanée des flux C++ (`std::cout`, `std::cerr`) et des unités logiques Fortran (`write(*, *)`, `write(6, *)`) sur les clusters HPC peut provoquer des désynchronisations et des entrelacements de lignes dans les fichiers de logs Slurm / PBS.

Les règles de synchronisation suivantes sont obligatoires :
1. Avant tout appel à un noyau de calcul Fortran, vider systématiquement le tampon C++ :
   ```cpp
   std::cout << std::flush;
   std::cerr << std::flush;
   ```
2. Côté Fortran, forcer la purge du flux de sortie standard via l'instruction intrinsèque `flush` :
   ```fortran
   write(6, '(A, I6, A, ES12.4)') " [Iter] ", iter, " Residue: ", res_norm
   flush(6)
   ```

### 2.4.3. Propagation Sécurisée des Codes d'Erreur

Le lancement d'une exception C++ (`throw std::runtime_error(...)`) à travers la frontière d'un appel de fonction Fortran provoque un arrêt immédiat du programme (*undefined behavior*, `std::terminate`), car la pile d'exécution Fortran ne supporte pas le déroulement de pile C++ (*stack unwinding*). Réciproquement, un `stop` ou `error stop` Fortran interrompt brutalement l'exécutable MPI sans donner à C++ l'opportunité de fermer les fichiers VTK ou d'invoquer `MPI_Finalize()`.

Par conséquent, **toutes les interfaces inter-langages utilisent une variable de retour de statut scalaire `ierr`** :

```fortran
subroutine covof_execute_step(ierr) bind(C, name="covof_execute_step")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), intent(out) :: ierr

    ierr = 0
    ! Détection d'une anomalie physique ou numérique (ex: divergence CFL)
    if (cfl_max > 10.0_c_double .or. isnan(cfl_max)) then
        ierr = 101 ! Code d'erreur divergence CFL
        return
    end if
end subroutine covof_execute_step
```

Côté C++, le code de retour est testé immédiatement après l'appel, et une exception contrôlée est levée dans le contexte C++ :

```cpp
int ierr = 0;
covof_execute_step(&ierr);
if (ierr != 0) {
    throw std::runtime_error("Erreur critique detectee dans le sous-programme Fortran. Code: " + std::to_string(ierr));
}
```

---

## 2.5. Catalogue des Points d'Entrée et Passerelles HPC

Les interfaces partagées de bas niveau couvrent les communicateurs distribués MPI et les opérations d'algèbre linéaire BLAS / LAPACK.

### 2.5.1. Traduction des Communicateurs MPI (`MPI_Comm_c2f` / `MPI_Comm_f2c`)

Les communicateurs MPI sont représentés par des types différents :
- En C/C++ : `MPI_Comm` (généralement un pointeur ou un typedef opaque `struct ompi_communicator_t*` sous OpenMPI, ou un entier `int` sous MPICH).
- En Fortran : `integer` (type entier 32 bits standard).

Le standard MPI impose l'usage des fonctions de conversion pour garantir la parfaite portabilité entre implémentations :

```cpp
// Conversion d'un communicateur C++ en entier Fortran
MPI_Comm cpp_cart_comm = ...;
MPI_Fint f_comm = MPI_Comm_c2f(cpp_cart_comm);

// Passage de f_comm à Fortran
fortran_setup_cartesian_comm(&f_comm);
```

Et réciproquement en Fortran :

```fortran
subroutine fortran_setup_cartesian_comm(f_comm) bind(C, name="fortran_setup_cartesian_comm")
    use, intrinsic :: iso_c_binding
    use mpi
    implicit none
    integer(c_int), intent(in), value :: f_comm
    integer :: local_comm, ierr, rank

    local_comm = f_comm
    call MPI_Comm_rank(local_comm, rank, ierr)
end subroutine fortran_setup_cartesian_comm
```

### 2.5.2. Tableau Récapitulatif des Fonctions Jumelles C++ / Fortran

Le tableau suivant dresse l'inventaire exhaustif de la gémellité des opérateurs critiques du solveur :

| Opérateur Numérique | Routine Fortran 2018 | Méthode C++17 | Signature d'Appel / Invariants |
| :--- | :--- | :--- | :--- |
| **Opérateur $A \cdot \mathbf{x}$ (Navier-Stokes)** | `module_navier_stokes:ns_apply_system_matrix` | `NavierStokesSolver::applySystemMatrix` | `(x_in, Ax_out, fields, mesh, config)` |
| **Calcul Divergence $\nabla \cdot \mathbf{u}$** | `module_navier_stokes:compute_divergence` | `NavierStokesSolver::computeDivergence` | `(u, v, div_out, mesh, dx, dy)` |
| **Terme Source Gravité / Flottabilité** | `module_navier_stokes:add_boussinesq_buoyancy` | `NavierStokesSolver::addBuoyancyForce` | `(fy, T, rho, T_ref, beta, gy)` |
| **Pénalisation de Brinkman** | `module_navier_stokes:apply_brinkman_drag` | `NavierStokesSolver::applyBrinkmanPenalization` | `(u, v, chi_mask, K_perm, dt)` |
| **Reconstruction Normale PLIC** | `module_interface_reconstruction:compute_plic_normal` | `InterfaceReconstruction::computeNormals` | `(C, nx, ny, alpha, curvature)` |
| **Advection VOF Weymouth-Yue** | `module_interface_tracking:advect_vof_weymouth_yue` | `InterfaceTracker::advectConservative` | `(C, u, v, C_new, dt, mesh)` |
| **Préconditionneur SSOR** | `module_preconditioners:precond_apply_ssor` | `SSORPreconditioner::apply` | `(r_in, z_out, omega, diag, mesh)` |
| **Préconditionneur ILU(0)** | `module_preconditioners:precond_apply_ilu0` | `ILU0Preconditioner::apply` | `(r_in, z_out, l_mat, u_mat, mesh)` |
| **Préconditionneur BT-SBP** | `module_preconditioners:precond_apply_bt_sbp` | `BlockTriangularPreconditioner::apply` | `(r_u, r_v, r_p, z_u, z_v, z_p, mesh)` |
| **Solveur BiCGSTAB2** | `module_linear_solvers:solve_bicgstab2` | `BiCGSTAB2Solver::solve` | `(A_op, b, x, tol, max_iter, stats)` |
| **Solveur MINRES** | `module_linear_solvers:solve_minres` | `MINRESSolver::solve` | `(A_op, b, x, tol, max_iter, stats)` |
| **Échanges Mailles Fantômes** | `module_mpi:mpi_exchange_halos_2d` | `MpiComm::exchangeHalos` | `(field, nxt, nyt, nghost, tag)` |
