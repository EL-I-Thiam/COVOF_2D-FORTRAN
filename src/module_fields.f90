module module_fields
  use module_precision, only: dp
  use module_mesh, only: mesh_2d_t
  implicit none
  private

  public :: flow_fields_t, allocate_fields, deallocate_fields, reset_fields

  type :: flow_fields_t
     ! Dimensions locales
     integer :: nx = 0, ny = 0

     ! Champs hydrodynamiques MAC (avec une couche de cellules fantômes 0..nx+1, 0..ny+1)
     real(dp), allocatable :: u(:, :)        ! Vitesse aux faces x (i-1/2, j)
     real(dp), allocatable :: v(:, :)        ! Vitesse aux faces y (i, j-1/2)
     real(dp), allocatable :: p(:, :)        ! Pression aux centres (i, j)
     real(dp), allocatable :: div(:, :)      ! Divergence discrète locale

     ! Scalaires et grandeurs thermodynamiques
     real(dp), allocatable :: C(:, :)        ! Fraction volumique VOF [0, 1]
     real(dp), allocatable :: T(:, :)        ! Température thermodynamique [K]
     real(dp), allocatable :: rho(:, :)      ! Masse volumique effective
     real(dp), allocatable :: rho_u(:, :)    ! Masse volumique aux faces u
     real(dp), allocatable :: rho_v(:, :)    ! Masse volumique aux faces v
     real(dp), allocatable :: mu(:, :)       ! Viscosité dynamique cellule
     real(dp), allocatable :: mu_node(:, :)  ! Viscosité dynamique nodale harmonique

     ! Forces capillaires CSF et courbure
     real(dp), allocatable :: fx_sigma(:, :)
     real(dp), allocatable :: fy_sigma(:, :)
     real(dp), allocatable :: curvature(:, :)

     ! Masque de pénalisation solide de Brinkman (1 = solide, 0 = fluide)
     integer,  allocatable :: solid_mask(:, :)

     ! Diagnostic thermique
     real(dp) :: mean_nusselt = 0.0_dp
  contains
     procedure :: allocate => fields_allocate_mesh
     procedure :: free => fields_free
  end type flow_fields_t

contains

  subroutine fields_allocate_mesh(this, mesh)
    class(flow_fields_t), intent(inout) :: this
    type(mesh_2d_t), intent(in) :: mesh
    call allocate_fields(this, mesh%nx_local, mesh%ny_local)
  end subroutine fields_allocate_mesh

  subroutine fields_free(this)
    class(flow_fields_t), intent(inout) :: this
    call deallocate_fields(this)
  end subroutine fields_free

  subroutine allocate_fields(fields, nx, ny)
    type(flow_fields_t), intent(inout) :: fields
    integer, intent(in) :: nx, ny

    fields%nx = nx
    fields%ny = ny

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

    call reset_fields(fields)
  end subroutine allocate_fields

  subroutine deallocate_fields(fields)
    type(flow_fields_t), intent(inout) :: fields

    if (allocated(fields%u))         deallocate(fields%u)
    if (allocated(fields%v))         deallocate(fields%v)
    if (allocated(fields%p))         deallocate(fields%p)
    if (allocated(fields%div))       deallocate(fields%div)
    if (allocated(fields%C))         deallocate(fields%C)
    if (allocated(fields%T))         deallocate(fields%T)
    if (allocated(fields%rho))       deallocate(fields%rho)
    if (allocated(fields%rho_u))     deallocate(fields%rho_u)
    if (allocated(fields%rho_v))     deallocate(fields%rho_v)
    if (allocated(fields%mu))        deallocate(fields%mu)
    if (allocated(fields%mu_node))   deallocate(fields%mu_node)
    if (allocated(fields%fx_sigma))  deallocate(fields%fx_sigma)
    if (allocated(fields%fy_sigma))  deallocate(fields%fy_sigma)
    if (allocated(fields%curvature)) deallocate(fields%curvature)
    if (allocated(fields%solid_mask))deallocate(fields%solid_mask)
    fields%nx = 0
    fields%ny = 0
  end subroutine deallocate_fields

  subroutine reset_fields(fields)
    type(flow_fields_t), intent(inout) :: fields

    fields%u = 0.0_dp
    fields%v = 0.0_dp
    fields%p = 0.0_dp
    fields%div = 0.0_dp
    fields%C = 1.0_dp
    fields%T = 0.0_dp
    fields%rho = 1.0_dp
    fields%rho_u = 1.0_dp
    fields%rho_v = 1.0_dp
    fields%mu = 1.0e-3_dp
    fields%mu_node = 1.0e-3_dp
    fields%fx_sigma = 0.0_dp
    fields%fy_sigma = 0.0_dp
    fields%curvature = 0.0_dp
    fields%solid_mask = 0
    fields%mean_nusselt = 0.0_dp
  end subroutine reset_fields

end module module_fields
