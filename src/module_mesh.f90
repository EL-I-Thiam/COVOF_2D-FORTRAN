!===============================================================================
! Module        : module_mesh
! Description   : Définition de la topologie cartésienne 2D, gestion du
!                 partitionnement spatial multi-processus OpenMPI, calcul des
!                 coordonnées métriques (faces décalées et centres de cellules)
!                 et détection des inclusions géométriques d'obstacles.
! Auteur        : Elhadji Ibrahima Thiam
! Standard      : Fortran 2018 / MPI-3
!===============================================================================
module module_mesh
  ! Importation du descripteur de double précision
  use module_precision, only: dp
  ! Interdiction stricte du typage implicite pour garantir l'intégrité mémoire
  implicit none
  ! Par défaut, tous les membres et routines sont privés
  private
  ! Inclusion des définitions MPI requises pour les communicateurs cartésiens
  include 'mpif.h'

  ! Export explicite des types dérivés et fonctions publiques
  public :: mesh_2d_t, init_mesh, destroy_mesh
  public :: x_face, y_face, x_center, y_center, is_in_obstacle

  !-----------------------------------------------------------------------------
  ! Structure de données représentant le maillage 2D décomposé
  !-----------------------------------------------------------------------------
  type :: mesh_2d_t
     ! Dimensions globales du domaine de calcul (nombre total de mailles)
     integer  :: nx_global = 64, ny_global = 64
     ! Longueurs physiques totales du domaine en mètres [m]
     real(dp) :: lx = 1.0_dp, ly = 1.0_dp
     ! Pas d'espace uniformes selon les directions X et Y : dx = lx/nx, dy = ly/ny
     real(dp) :: dx = 0.0_dp, dy = 0.0_dp

     ! --- Topologie Cartésienne MPI et Rangs de Voisinage ---
     integer  :: rank = 0                 ! Rang du sous-domaine local courant
     integer  :: size = 1                 ! Nombre total de sous-domaines (processus)
     integer  :: comm = MPI_COMM_WORLD    ! Communicateur initial parent
     integer  :: cart_comm = MPI_COMM_NULL! Communicateur cartésien bidimensionnel dédié
     integer  :: px = 1, py = 1           ! Nombre de processus selon X et selon Y
     integer  :: coord_x = 0, coord_y = 0 ! Coordonnées cartésiennes du sous-domaine dans la grille MPI
     integer  :: neighbor_west = MPI_PROC_NULL  ! Rang du voisin Ouest (gauche, x-)
     integer  :: neighbor_east = MPI_PROC_NULL  ! Rang du voisin Est (droite, x+)
     integer  :: neighbor_south = MPI_PROC_NULL ! Rang du voisin Sud (bas, y-)
     integer  :: neighbor_north = MPI_PROC_NULL ! Rang du voisin Nord (haut, y+)

     ! --- Indicateurs Booléens de Frontières Physiques du Domaine Global ---
     logical  :: is_boundary_left   = .false. ! Vrai si le sous-domaine touche la paroi gauche x = 0
     logical  :: is_boundary_right  = .false. ! Vrai si le sous-domaine touche la paroi droite x = lx
     logical  :: is_boundary_bottom = .false. ! Vrai si le sous-domaine touche la paroi inférieure y = 0
     logical  :: is_boundary_top    = .false. ! Vrai si le sous-domaine touche la paroi supérieure y = ly

     ! --- Partitionnement Spatial Local ---
     integer  :: nx_local = 64, ny_local = 64 ! Nombre de cellules locales attribuées à ce processus
     integer  :: i_start_global = 1           ! Indice global de la première cellule locale en X (base 1)
     integer  :: j_start_global = 1           ! Indice global de la première cellule locale en Y (base 1)
  end type mesh_2d_t

contains

  !=============================================================================
  ! Routine : init_mesh
  ! Rôle    : Calcule la factorisation optimale du nombre de cœurs en une grille
  !           2D (px, py), crée la topologie cartésienne MPI avec périodicité,
  !           identifie les 4 voisins immédiats et répartit équitablement les mailles.
  !=============================================================================
  subroutine init_mesh(mesh, comm_world, nx, ny, lx, ly, periodic_x, periodic_y)
    type(mesh_2d_t), intent(out) :: mesh       ! Structure de maillage initialisée en sortie
    integer, intent(in)          :: comm_world ! Communicateur parent initial
    integer, intent(in)          :: nx, ny     ! Résolution globale souhaitée en cellules
    real(dp), intent(in)         :: lx, ly     ! Longueurs physiques du domaine rectangulaire
    logical, intent(in)          :: periodic_x ! Vrai si conditions périodiques en X
    logical, intent(in)          :: periodic_y ! Vrai si conditions périodiques en Y

    ! Variables de travail locales pour les primitives cartésiennes MPI
    integer :: dims(2), periods(2), coords(2)
    integer :: ierr, rank_in, size_in
    integer :: rem_x, rem_y

    ! 1. Enregistrement des dimensions physiques et de la taille du communicateur
    mesh%comm = comm_world
    call MPI_Comm_rank(comm_world, rank_in, ierr)
    call MPI_Comm_size(comm_world, size_in, ierr)
    mesh%rank = rank_in
    mesh%size = size_in

    mesh%nx_global = nx
    mesh%ny_global = ny
    mesh%lx = lx
    mesh%ly = ly
    ! Calcul des pas de discrétisation spatiale
    mesh%dx = lx / real(nx, dp)
    mesh%dy = ly / real(ny, dp)

    ! 2. Détermination automatique de la décomposition 2D optimale (px, py)
    dims = [0, 0]
    call MPI_Dims_create(mesh%size, 2, dims, ierr)
    mesh%px = dims(1)
    mesh%py = dims(2)

    ! Définition des drapeaux de périodicité pour la création de la topologie MPI
    periods(1) = merge(1, 0, periodic_x)
    periods(2) = merge(1, 0, periodic_y)

    ! 3. Création du communicateur virtuel cartésien MPI 2D
    call MPI_Cart_create(comm_world, 2, dims, periods, 1, mesh%cart_comm, ierr)
    call MPI_Comm_rank(mesh%cart_comm, mesh%rank, ierr)
    call MPI_Cart_coords(mesh%cart_comm, mesh%rank, 2, coords, ierr)
    mesh%coord_x = coords(1)
    mesh%coord_y = coords(2)

    ! 4. Recherche des rangs des 4 voisins adjacents par décalage unitaire cartésien
    ! Décalage direction 0 (axe X) : voisin Ouest (-1) et voisin Est (+1)
    call MPI_Cart_shift(mesh%cart_comm, 0, 1, mesh%neighbor_west, mesh%neighbor_east, ierr)
    ! Décalage direction 1 (axe Y) : voisin Sud (-1) et voisin Nord (+1)
    call MPI_Cart_shift(mesh%cart_comm, 1, 1, mesh%neighbor_south, mesh%neighbor_north, ierr)

    ! 5. Répartition équitable des mailles le long de l'axe X (gestion du reste mod)
    mesh%nx_local = nx / mesh%px
    rem_x = mod(nx, mesh%px)
    if (mesh%coord_x < rem_x) then
       ! Les premiers rem_x processus reçoivent une cellule supplémentaire
       mesh%nx_local = mesh%nx_local + 1
       mesh%i_start_global = mesh%coord_x * mesh%nx_local + 1
    else
       mesh%i_start_global = rem_x * (mesh%nx_local + 1) + (mesh%coord_x - rem_x) * mesh%nx_local + 1
    end if

    ! 6. Répartition équitable des mailles le long de l'axe Y
    mesh%ny_local = ny / mesh%py
    rem_y = mod(ny, mesh%py)
    if (mesh%coord_y < rem_y) then
       mesh%ny_local = mesh%ny_local + 1
       mesh%j_start_global = mesh%coord_y * mesh%ny_local + 1
    else
       mesh%j_start_global = rem_y * (mesh%ny_local + 1) + (mesh%coord_y - rem_y) * mesh%ny_local + 1
    end if

    ! 7. Détection des frontières physiques externes non périodiques
    mesh%is_boundary_left   = (mesh%neighbor_west == MPI_PROC_NULL)
    mesh%is_boundary_right  = (mesh%neighbor_east == MPI_PROC_NULL)
    mesh%is_boundary_bottom = (mesh%neighbor_south == MPI_PROC_NULL)
    mesh%is_boundary_top    = (mesh%neighbor_north == MPI_PROC_NULL)
  end subroutine init_mesh

  !=============================================================================
  ! Routine : destroy_mesh
  ! Rôle    : Libère le communicateur cartésien alloué auprès du runtime MPI.
  !=============================================================================
  subroutine destroy_mesh(mesh)
    type(mesh_2d_t), intent(inout) :: mesh
    integer :: ierr
    if (mesh%cart_comm /= MPI_COMM_NULL .and. mesh%cart_comm /= mesh%comm) then
       call MPI_Comm_free(mesh%cart_comm, ierr)
    end if
  end subroutine destroy_mesh

  !=============================================================================
  ! Fonction pure : x_face
  ! Rôle          : Calcule l'abscisse physique absolue de la face verticale i-1/2.
  !=============================================================================
  pure real(dp) function x_face(mesh, i)
    type(mesh_2d_t), intent(in) :: mesh
    integer, intent(in)          :: i
    x_face = real(mesh%i_start_global - 1 + i - 1, dp) * mesh%dx
  end function x_face

  !=============================================================================
  ! Fonction pure : y_face
  ! Rôle          : Calcule l'ordonnée physique absolue de la face horizontale j-1/2.
  !=============================================================================
  pure real(dp) function y_face(mesh, j)
    type(mesh_2d_t), intent(in) :: mesh
    integer, intent(in)          :: j
    y_face = real(mesh%j_start_global - 1 + j - 1, dp) * mesh%dy
  end function y_face

  !=============================================================================
  ! Fonction pure : x_center
  ! Rôle          : Calcule l'abscisse du barycentre géométrique de la cellule i.
  !=============================================================================
  pure real(dp) function x_center(mesh, i)
    type(mesh_2d_t), intent(in) :: mesh
    integer, intent(in)          :: i
    x_center = (real(mesh%i_start_global - 1 + i - 1, dp) + 0.5_dp) * mesh%dx
  end function x_center

  !=============================================================================
  ! Fonction pure : y_center
  ! Rôle          : Calcule l'ordonnée du barycentre géométrique de la cellule j.
  !=============================================================================
  pure real(dp) function y_center(mesh, j)
    type(mesh_2d_t), intent(in) :: mesh
    integer, intent(in)          :: j
    y_center = (real(mesh%j_start_global - 1 + j - 1, dp) + 0.5_dp) * mesh%dy
  end function y_center

  !=============================================================================
  ! Fonction pure : is_in_obstacle
  ! Rôle          : Teste l'appartenance d'un point (x, y) à un obstacle immergé
  !                 (rectangle ou cylindre) pour la pénalisation de Brinkman.
  !=============================================================================
  pure logical function is_in_obstacle(x, y, obs_type, x0, y0, x1, y1, xc, yc, r)
    real(dp), intent(in) :: x, y        ! Coordonnées du point à tester
    integer, intent(in)  :: obs_type    ! Type géométrique : 1=boîte, 2=cylindre
    real(dp), intent(in) :: x0, y0      ! Coin inférieur gauche du rectangle
    real(dp), intent(in) :: x1, y1      ! Coin supérieur droit du rectangle
    real(dp), intent(in) :: xc, yc      ! Centre du disque cylindrique
    real(dp), intent(in) :: r           ! Rayon du disque

    is_in_obstacle = .false.
    if (obs_type == 1) then ! Obstacle rectangulaire (ex: marche descendante)
       if (x >= x0 - 1.0e-9_dp .and. x <= x1 + 1.0e-9_dp .and. &
           y >= y0 - 1.0e-9_dp .and. y <= y1 + 1.0e-9_dp) then
          is_in_obstacle = .true.
       end if
    else if (obs_type == 2) then ! Obstacle cylindrique circulaire
       if ((x - xc)**2 + (y - yc)**2 <= r**2) then
          is_in_obstacle = .true.
       end if
    end if
  end function is_in_obstacle

end module module_mesh
