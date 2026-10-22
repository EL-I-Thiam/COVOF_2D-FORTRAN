!===============================================================================
! Module        : module_fields
! Description   : Définition et gestion de la mémoire des champs hydrodynamiques,
!                 thermiques et diphasiques sur la grille décalée MAC (Harlow-Welch).
!                 Assure l'allocation contiguë avec cellules fantômes pour les
!                 échanges de halos MPI et garantit la politique Zero-Allocation
!                 au sein de la boucle en temps.
! Auteur        : Elhadji Ibrahima Thiam
! Standard      : Fortran 2018
!===============================================================================
module module_fields
  ! Importation du descripteur de double précision
  use module_precision, only: dp
  ! Importation du type dérivé de maillage pour connaître les dimensions locales
  use module_mesh, only: mesh_2d_t
  ! Typage strict
  implicit none
  ! Encapsulation privée
  private

  ! Exportation publique du type conteneur et des routines associées
  public :: flow_fields_t, allocate_fields, deallocate_fields, reset_fields

  !-----------------------------------------------------------------------------
  ! Structure de données contenant l'ensemble des champs physiques du sous-domaine
  !-----------------------------------------------------------------------------
  type :: flow_fields_t
     ! Nombre de cellules physiques locales (hors cellules fantômes)
     integer :: nx = 0, ny = 0

     ! --- Champs Hydrodynamiques Discrétisés sur Grille Décalée (MAC) ---
     ! Indices : (0:nx+1, 0:ny+1) pour inclure 1 rangée de cellules fantômes
     real(dp), allocatable :: u(:, :)        ! Vitesse axiale aux faces verticales (i-1/2, j) [m/s]
     real(dp), allocatable :: v(:, :)        ! Vitesse transversale aux faces horizontales (i, j-1/2) [m/s]
     real(dp), allocatable :: p(:, :)        ! Pression statique scalaire au centre des mailles (i, j) [Pa]
     real(dp), allocatable :: div(:, :)      ! Divergence discrète locale de la vitesse (nabla . u) [1/s]

     ! --- Propriétés Thermodynamiques et Équation d'Énergie ---
     real(dp), allocatable :: C(:, :)        ! Fraction volumique VOF du fluide 1 (C in [0, 1])
     real(dp), allocatable :: T(:, :)        ! Température thermodynamique au centre des mailles [K]
     real(dp), allocatable :: rho(:, :)      ! Masse volumique effective au centre des cellules [kg/m^3]
     real(dp), allocatable :: rho_u(:, :)    ! Masse volumique interpolée aux faces u [kg/m^3]
     real(dp), allocatable :: rho_v(:, :)    ! Masse volumique interpolée aux faces v [kg/m^3]
     real(dp), allocatable :: mu(:, :)       ! Viscosité dynamique au centre des mailles [Pa.s]
     real(dp), allocatable :: mu_node(:, :)  ! Viscosité harmonique aux sommets/nœuds (i+1/2, j+1/2) [Pa.s]

     ! --- Forces Capillaires et Tension de Surface (Modèle CSF Brackbill) ---
     real(dp), allocatable :: fx_sigma(:, :) ! Force capillaire volumique selon X aux faces u [N/m^3]
     real(dp), allocatable :: fy_sigma(:, :) ! Force capillaire volumique selon Y aux faces v [N/m^3]
     real(dp), allocatable :: curvature(:, :)! Courbure locale d'interface kappa = -nabla . n [1/m]

     ! --- Masque d'Obstacle Immergé (Pénalisation de Brinkman) ---
     integer,  allocatable :: solid_mask(:, :)! 1 si la cellule est dans le solide, 0 en fluide pur

     ! --- Indicateurs et Diagnostics de Validation ---
     real(dp) :: mean_nusselt = 0.0_dp       ! Nombre de Nusselt moyen mesuré à la paroi chaude

  contains
     ! Méthodes liées au type dérivé pour syntaxe orientée objet
     procedure :: allocate => fields_allocate_mesh
     procedure :: free     => fields_free
  end type flow_fields_t

contains

  !=============================================================================
  ! Routine liée : fields_allocate_mesh
  ! Rôle         : Alloue les tableaux à partir d'une instance de mesh_2d_t.
  !=============================================================================
  subroutine fields_allocate_mesh(this, mesh)
    class(flow_fields_t), intent(inout) :: this
    type(mesh_2d_t), intent(in)         :: mesh
    call allocate_fields(this, mesh%nx_local, mesh%ny_local)
  end subroutine fields_allocate_mesh

  !=============================================================================
  ! Routine liée : fields_free
  ! Rôle         : Désalloue proprement l'ensemble des tableaux du conteneur.
  !=============================================================================
  subroutine fields_free(this)
    class(flow_fields_t), intent(inout) :: this
    call deallocate_fields(this)
  end subroutine fields_free

  !=============================================================================
  ! Routine : allocate_fields
  ! Rôle    : Alloue de façon contiguë en mémoire chaque champ 2D avec ses halos
  !           fantômes (indices de 0 à nx+1 et de 0 à ny+1) et initialise à zéro.
  !=============================================================================
  subroutine allocate_fields(fields, nx, ny)
    type(flow_fields_t), intent(inout) :: fields
    integer, intent(in)                :: nx, ny

    ! Enregistrement des dimensions locales
    fields%nx = nx
    fields%ny = ny

    ! Allocation dynamique unique (zéro allocation subséquente dans la boucle temporelle)
    allocate(fields%u(0:nx+1, 0:ny+1))
    allocate(fields%v(0:nx+1, 0:ny+1))
    allocate(fields%p(0:nx+1, 0:ny+1))
    allocate(fields%div(0:nx+1, 0:ny+1))
    allocate(fields%C(0:nx+1, 0:ny+1))
    allocate(fields%T(0:nx+1, 0:ny+1))
    allocate(fields%rho(0:nx+1, 0:ny+1))
    allocate(fields%rho_u(0:nx+1, 0:ny+1))
    allocate(fields%rho_v(0:nx+1, 0:ny+1))
    allocate(fields%mu(0:nx+1, 0:ny+1))
    allocate(fields%mu_node(0:nx+1, 0:ny+1))
    allocate(fields%fx_sigma(0:nx+1, 0:ny+1))
    allocate(fields%fy_sigma(0:nx+1, 0:ny+1))
    allocate(fields%curvature(0:nx+1, 0:ny+1))
    allocate(fields%solid_mask(0:nx+1, 0:ny+1))

    ! Remise à zéro des valeurs initiales
    call reset_fields(fields)
  end subroutine allocate_fields

  !=============================================================================
  ! Routine : deallocate_fields
  ! Rôle    : Désalloue chaque tableau et remet les compteurs de dimension à zéro.
  !=============================================================================
  subroutine deallocate_fields(fields)
    type(flow_fields_t), intent(inout) :: fields

    if (allocated(fields%u))          deallocate(fields%u)
    if (allocated(fields%v))          deallocate(fields%v)
    if (allocated(fields%p))          deallocate(fields%p)
    if (allocated(fields%div))        deallocate(fields%div)
    if (allocated(fields%C))          deallocate(fields%C)
    if (allocated(fields%T))          deallocate(fields%T)
    if (allocated(fields%rho))        deallocate(fields%rho)
    if (allocated(fields%rho_u))      deallocate(fields%rho_u)
    if (allocated(fields%rho_v))      deallocate(fields%rho_v)
    if (allocated(fields%mu))         deallocate(fields%mu)
    if (allocated(fields%mu_node))    deallocate(fields%mu_node)
    if (allocated(fields%fx_sigma))   deallocate(fields%fx_sigma)
    if (allocated(fields%fy_sigma))   deallocate(fields%fy_sigma)
    if (allocated(fields%curvature))  deallocate(fields%curvature)
    if (allocated(fields%solid_mask)) deallocate(fields%solid_mask)

    fields%nx = 0
    fields%ny = 0
  end subroutine deallocate_fields

  !=============================================================================
  ! Routine : reset_fields
  ! Rôle    : Initialise les champs à des valeurs physiques par défaut neutres.
  !=============================================================================
  subroutine reset_fields(fields)
    type(flow_fields_t), intent(inout) :: fields

    fields%u            = 0.0_dp     ! Vitesse axiale initiale nulle
    fields%v            = 0.0_dp     ! Vitesse transversale initiale nulle
    fields%p            = 0.0_dp     ! Pression de référence nulle
    fields%div          = 0.0_dp     ! Divergence initiale nulle
    fields%C            = 1.0_dp     ! Fluide 1 pur par défaut
    fields%T            = 0.0_dp     ! Température de référence nulle
    fields%rho          = 1.0_dp     ! Masse volumique unitaire par défaut
    fields%rho_u        = 1.0_dp
    fields%rho_v        = 1.0_dp
    fields%mu           = 1.0e-3_dp  ! Viscosité dynamique standard de l'eau (1 mPa.s)
    fields%mu_node      = 1.0e-3_dp
    fields%fx_sigma     = 0.0_dp     ! Tension superficielle nulle
    fields%fy_sigma     = 0.0_dp
    fields%curvature    = 0.0_dp
    fields%solid_mask   = 0          ! 100% domaine fluide sans obstacle
    fields%mean_nusselt = 0.0_dp
  end subroutine reset_fields

end module module_fields
