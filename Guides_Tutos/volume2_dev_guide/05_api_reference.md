# Référence Exhaustive de l'API Interne

Ce chapitre constitue le manuel de référence technique détaillé des interfaces logicielles (API) des deux branches du solveur COVOF 2D. Il recense l'intégralité des classes C++17 et des modules Fortran 2018 avec leurs membres, leurs constructeurs, leurs signatures de méthodes et leurs spécificateurs d'intention.

---

## 5.1. Référence Complète de l'API C++17

L'architecture C++17 s'articule autour de l'espace de nommage global `covof`.

### 5.1.1. Classe `Mesh2D` (`include/covof/Mesh.hpp`)

Modélise la grille cartésienne 2D décalée MAC et encapsule la métrique géométrique locale et globale.

- **Constructeur & Destructeur :**
  - `Mesh2D(const Config& config, const MpiComm& mpi_comm);`
    *Rôle :* Alloue et calcule les métriques spatiales, les coordonnées 1D des faces et centres, et initialise les bornes locales selon le rang cartésien.
  - `~Mesh2D() noexcept = default;`

- **Membres de données publics / accesseurs :**
  - `int nx_local() const noexcept;` : Nombre de mailles utiles locales en direction $x$.
  - `int ny_local() const noexcept;` : Nombre de mailles utiles locales en direction $y$.
  - `int nx_global() const noexcept;` : Nombre de mailles globales en direction $x$.
  - `int ny_global() const noexcept;` : Nombre de mailles globales en direction $y$.
  - `int nghost() const noexcept;` : Épaisseur de la couche de mailles fantômes ($N_{\text{ghost}} = 2$).
  - `double dx() const noexcept;` : Pas d'espace horizontal uniforme $\Delta x$.
  - `double dy() const noexcept;` : Pas d'espace vertical uniforme $\Delta y$.
  - `double x_min_local() const noexcept;` : Abscisse minimale du sous-domaine local.
  - `double y_min_local() const noexcept;` : Ordonnée minimale du sous-domaine local.
  - `const std::vector<double>& x_centers() const noexcept;` : Coordonnées $x$ des centres de mailles.
  - `const std::vector<double>& y_centers() const noexcept;` : Coordonnées $y$ des centres de mailles.
  - `const std::vector<double>& x_faces() const noexcept;` : Coordonnées $x$ des faces verticales.
  - `const std::vector<double>& y_faces() const noexcept;` : Coordonnées $y$ des faces horizontales.

### 5.1.2. Structure `FlowFields` (`include/covof/Field.hpp`)

Stocke l'ensemble des champs scalaires et vectoriels discrets au format SoA (*Structure of Arrays*).

- **Champs de Données Publics (`std::vector<double>`) :**
  - `u` : Composante horizontale de vitesse sur les faces verticales $(i+1/2, j)$, dimension $(N_{xt} + 1) \times N_{yt}$.
  - `v` : Composante verticale de vitesse sur les faces horizontales $(i, j+1/2)$, dimension $N_{xt} \times (N_{yt} + 1)$.
  - `p` : Champ scalaire de pression aux centres de mailles $(i, j)$, dimension $N_{xt} \times N_{yt}$.
  - `rho` : Masse volumique aux centres de mailles, dimension $N_{xt} \times N_{yt}$.
  - `mu` : Viscosité dynamique aux centres de mailles, dimension $N_{xt} \times N_{yt}$.
  - `c` : Fraction volumique VOF aux centres de mailles ($C \in [0, 1]$), dimension $N_{xt} \times N_{yt}$.
  - `T` : Température passive ou active aux centres de mailles, dimension $N_{xt} \times N_{yt}$.
  - `chi_T` : Fonction caractéristique de masque poreux/solide de Brinkman ($\chi_T \in [0, 1]$).
  - `r_aal` : Variable mémoire de pénalisation adaptative (AAL / C-AAL).

- **Méthodes Utilitaires :**
  - `void allocate(int nxt, int nyt);` : Alloue l'ensemble des tampons mémoires contigus avec alignement 64 octets.
  - `void reset();` : Réinitialise tous les champs à zéro.

### 5.1.3. Classe `MpiComm` (`include/covof/MpiComm.hpp`)

Gère la topologie cartésienne distribuée et les échanges de mailles fantômes.

- **Méthodes Publiques :**
  - `void init(int argc, char** argv);` : Initialise `MPI_Init_thread` avec support multithread `MPI_THREAD_FUNNELED`.
  - `void createCartesianTopology(int px, int py, bool periodic_x, bool periodic_y);` : Invoque `MPI_Cart_create`.
  - `void exchangeHalos(std::vector<double>& field, int nxt, int nyt, int nghost);` : Réalise les échanges non bloquants avec types dérivés de vecteurs.
  - `double globalMin(double local_val) const;` : Réduction globale `MPI_Allreduce` avec opérateur `MPI_MIN`.
  - `double globalMax(double local_val) const;` : Réduction globale `MPI_Allreduce` avec opérateur `MPI_MAX`.
  - `double globalSum(double local_val) const;` : Réduction globale `MPI_Allreduce` avec opérateur `MPI_SUM`.
  - `int rank() const noexcept;` : Rang cartésien du processus courant.
  - `int size() const noexcept;` : Nombre total de processus MPI alloués.

### 5.1.4. Classe `NavierStokesSolver` (`include/covof/NavierStokesSolver.hpp`)

Noyau numérique principal gérant l'intégration temporelle et spatiale des équations de Navier-Stokes.

- **Méthodes Publiques :**
  - `NavierStokesSolver(const Config& cfg, const Mesh2D& mesh, MpiComm& mpi);`
  - `void assembleRhs(const FlowFields& f, double dt, std::vector<double>& rhs_u, std::vector<double>& rhs_v);`
  - `void applySystemMatrix(const double* x_in, double* Ax_out);` : Opérateur couplé global sans matrice.
  - `void computeDivergence(const std::vector<double>& u, const std::vector<double>& v, std::vector<double>& div);`
  - `void updateUzawaPressure(FlowFields& f, double dt, double r_penalty);`
  - `void applyBrinkmanPenalization(FlowFields& f, double dt);`
  - `double computeMaxCFL(const FlowFields& f, double dt) const;`

### 5.1.5. Hiérarchie des Solveurs Linéaires (`include/covof/LinearSolver.hpp`)

- **Classe de Base Abstraite `LinearSolverBase` :**
  - `virtual void solve(const LinearOperator& A, const std::vector<double>& b, std::vector<double>& x, LinearSolverStats& stats) = 0;`
  - `virtual void setPreconditioner(std::shared_ptr<PreconditionerBase> precond);`

- **Classes Dérivées :**
  - `class BiCGSTAB2Solver : public LinearSolverBase;` : Biconjugué Gradient Stabilisé d'ordre 2 avec minimisation polynomiale.
  - `class BiCGSTABSolver : public LinearSolverBase;` : Biconjugué Gradient Stabilisé standard.
  - `class CGSolver : public LinearSolverBase;` : Gradient Conjugué symétrique défini positif.
  - `class GMRESSolver : public LinearSolverBase;` : GMRES redémarré avec orthogonalisation de Gram-Schmidt modifiée.
  - `class MINRESSolver : public LinearSolverBase;` : Résidu Minimal pour systèmes symétriques indéfinis.
  - `class DirectSolver : public LinearSolverBase;` : Solveur direct creux basé sur la factorisation LU LAPACK (`dgetrf`).

### 5.1.6. Hiérarchie des Préconditionneurs (`include/covof/Preconditioner.hpp`)

- **Classe de Base Abstraite `PreconditionerBase` :**
  - `virtual void setup(const LinearOperator& A) = 0;`
  - `virtual void apply(const double* r_in, double* z_out) const = 0;`

- **Classes Dérivées :**
  - `class JacobiPreconditioner : public PreconditionerBase;` : Inverse de la diagonale principale.
  - `class SSORPreconditioner : public PreconditionerBase;` : Sur-relaxation successive symétrique avec paramètre $\omega \in ]0, 2[$.
  - `class ILU0Preconditioner : public PreconditionerBase;` : Factorisation incomplète sans remplissage $ILU(0)$ sur format CSR.
  - `class BlockTriangularPreconditioner : public PreconditionerBase;` : Bloc triangulaire BT-SBP avec solveur de Schur approché.
  - `class CustomPreconditioner : public PreconditionerBase;` : Préconditionneur polynomial personnalisé.

---

## 5.2. Référence Complète de l'API Fortran 2018

L'architecture Fortran 2018 s'articule autour de 14 modules autonomes.

### 5.2.1. Module `module_precision.f90`

Définit les paramètres de types normalisés via `iso_fortran_env` et `iso_c_binding` :
- `integer, parameter :: dp = real64` : Précision double IEEE 754 (8 octets).
- `integer, parameter :: sp = real32` : Précision simple IEEE 754 (4 octets).
- `integer, parameter :: i4 = int32`  : Entier standard 32 bits.
- `integer, parameter :: i8 = int64`  : Entier grand volume 64 bits.

### 5.2.2. Module `module_mesh.f90`

Encapsule le type dérivé `mesh_2d_t` et ses sous-programmes de gestion :

```fortran
type, public :: mesh_2d_t
    integer(i4) :: nx_loc, ny_loc
    integer(i4) :: nx_glob, ny_glob
    integer(i4) :: nghost
    real(dp)    :: dx, dy
    real(dp)    :: x_min, y_min
    real(dp), allocatable :: xc(:), yc(:)
    real(dp), allocatable :: xf(:), yf(:)
contains
    procedure, public :: init => mesh_init_2d
    procedure, public :: destroy => mesh_destroy_2d
end type mesh_2d_t
```

- **Sous-routines associées :**
  - `subroutine mesh_init_2d(this, cfg, mpi_comm, ierr)`
    - `class(mesh_2d_t), intent(inout) :: this`
    - `type(config_t), intent(in) :: cfg`
    - `type(mpi_comm_t), intent(in) :: mpi_comm`
    - `integer(i4), intent(out) :: ierr`
  - `subroutine mesh_destroy_2d(this)`
    - `class(mesh_2d_t), intent(inout) :: this`

### 5.2.3. Module `module_flow_fields.f90`

Définit le type dérivé `flow_fields_t` pour le stockage des champs :

```fortran
type, public :: flow_fields_t
    real(dp), allocatable :: u(:, :)
    real(dp), allocatable :: v(:, :)
    real(dp), allocatable :: p(:, :)
    real(dp), allocatable :: rho(:, :)
    real(dp), allocatable :: mu(:, :)
    real(dp), allocatable :: c(:, :)
    real(dp), allocatable :: T(:, :)
    real(dp), allocatable :: chi_T(:, :)
    real(dp), allocatable :: r_aal(:, :)
contains
    procedure, public :: allocate_fields => fields_allocate_2d
    procedure, public :: deallocate_fields => fields_deallocate_2d
end type flow_fields_t
```

### 5.2.4. Module `module_navier_stokes.f90`

Noyau algorithmique des équations hydrodynamiques :

- `subroutine ns_apply_system_matrix(x_in, Ax_out, fields, mesh, cfg)`
  - `real(dp), dimension(:), intent(in) :: x_in`
  - `real(dp), dimension(:), intent(out) :: Ax_out`
  - `type(flow_fields_t), intent(in) :: fields`
  - `type(mesh_2d_t), intent(in) :: mesh`
  - `type(config_t), intent(in) :: cfg`
- `subroutine ns_compute_viscous_stress(fields, mesh, tau_u, tau_v)`
  - `type(flow_fields_t), intent(in) :: fields`
  - `type(mesh_2d_t), intent(in) :: mesh`
  - `real(dp), dimension(:, :), intent(out) :: tau_u, tau_v`
- `subroutine ns_compute_divergence(u, v, div, mesh)`
  - `real(dp), dimension(:, :), intent(in) :: u, v`
  - `real(dp), dimension(:, :), intent(out) :: div`
  - `type(mesh_2d_t), intent(in) :: mesh`
- `subroutine ns_apply_brinkman_drag(fields, dt, mesh, cfg)`
  - `type(flow_fields_t), intent(inout) :: fields`
  - `real(dp), intent(in) :: dt`
  - `type(mesh_2d_t), intent(in) :: mesh`
  - `type(config_t), intent(in) :: cfg`

### 5.2.5. Module `module_linear_solvers.f90`

Fournit l'implémentation vectorisée des méthodes de sous-espace de Krylov :

- `subroutine solve_bicgstab2(apply_A, apply_P, b, x, tol, max_iter, stats, comm)`
  - `procedure(linear_operator_interface) :: apply_A`
  - `procedure(precond_operator_interface) :: apply_P`
  - `real(dp), dimension(:), intent(in) :: b`
  - `real(dp), dimension(:), intent(inout) :: x`
  - `real(dp), intent(in) :: tol`
  - `integer(i4), intent(in) :: max_iter`
  - `type(solver_stats_t), intent(out) :: stats`
  - `integer(i4), intent(in) :: comm`
- `subroutine solve_minres(apply_A, apply_P, b, x, tol, max_iter, stats, comm)`
  - Mêmes spécificateurs d'intention, optimisé pour les matrices symétriques régies par l'algorithme de Lanczos.
- `subroutine solve_gmres(apply_A, apply_P, b, x, restart, tol, max_iter, stats, comm)`
  - `integer(i4), intent(in) :: restart` : Dimension de la base de Krylov orthogonale avant redémarrage.

### 5.2.6. Module `module_preconditioners.f90`

Contient les routines de préconditionnement :
- `subroutine precond_apply_ssor(r, z, omega, diag, mesh)`
  - `real(dp), dimension(:), intent(in) :: r`
  - `real(dp), dimension(:), intent(out) :: z`
  - `real(dp), intent(in) :: omega`
  - `real(dp), dimension(:), intent(in) :: diag`
  - `type(mesh_2d_t), intent(in) :: mesh`
- `subroutine precond_apply_ilu0(r, z, l_val, u_val, diag, row_ptr, col_ind)`
  - `real(dp), dimension(:), intent(in) :: r`
  - `real(dp), dimension(:), intent(out) :: z`
  - `real(dp), dimension(:), intent(in) :: l_val, u_val, diag`
  - `integer(i4), dimension(:), intent(in) :: row_ptr, col_ind`
- `subroutine precond_apply_bt_sbp(r_u, r_v, r_p, z_u, z_v, z_p, mesh, cfg)`
  - Application du préconditionneur bloc triangulaire avec résolution approchée du complément de Schur.
