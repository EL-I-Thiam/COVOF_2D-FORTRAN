module module_constants
  use module_precision, only: dp
  implicit none
  private

  ! Mathematical & physical constants
  real(dp), parameter, public :: PI = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter, public :: EPS_SMALL = 1.0e-30_dp
  real(dp), parameter, public :: PERMEABILITY_SOLID_DEFAULT = 1.0e-40_dp

  ! Boundary condition types
  integer, parameter, public :: BC_DIRICHLET = 1
  integer, parameter, public :: BC_NEUMANN   = 2
  integer, parameter, public :: BC_PERIODIC  = 3
  integer, parameter, public :: BC_FREESLIP  = 4

  ! Thermal boundary condition types
  integer, parameter, public :: THERMAL_DIRICHLET = 1
  integer, parameter, public :: THERMAL_NEUMANN   = 2

  ! Navier-Stokes solver formulation types
  integer, parameter, public :: SOLVER_AUGMENTED_LAGRANGIAN = 1
  integer, parameter, public :: SOLVER_MONOLITHIC           = 2

  ! Linear solver algorithms
  integer, parameter, public :: LINEAR_BICGSTAB2 = 1
  integer, parameter, public :: LINEAR_BICGSTAB  = 2
  integer, parameter, public :: LINEAR_CG        = 3
  integer, parameter, public :: LINEAR_GMRES     = 4
  integer, parameter, public :: LINEAR_DIRECT    = 5

  ! Preconditioner types
  integer, parameter, public :: PRECOND_NONE   = 0
  integer, parameter, public :: PRECOND_JACOBI = 1
  integer, parameter, public :: PRECOND_SSOR   = 2
  integer, parameter, public :: PRECOND_ILU0   = 3

  ! Convective scheme types
  integer, parameter, public :: CONV_UPWIND       = 1
  integer, parameter, public :: CONV_CENTERED     = 2
  integer, parameter, public :: CONV_TVD_SUPERBEE = 3

  ! Viscosity models
  integer, parameter, public :: VISC_CONSTANT    = 1
  integer, parameter, public :: VISC_EXPONENTIAL = 2

  ! Initial condition types
  integer, parameter, public :: IC_UNIFORM          = 1
  integer, parameter, public :: IC_POISEUILLE       = 2
  integer, parameter, public :: IC_VORTEX_SHEAR     = 3
  integer, parameter, public :: IC_RAYLEIGH_TAYLOR  = 4
  integer, parameter, public :: IC_DROPLET          = 5
  integer, parameter, public :: IC_BUBBLE           = 6
  integer, parameter, public :: IC_DAM_BREAK        = 7
  integer, parameter, public :: IC_LAPLACIAN_TEST   = 8

  ! Immersed boundary obstacle shapes
  integer, parameter, public :: OBS_NONE     = 0
  integer, parameter, public :: OBS_BOX      = 1
  integer, parameter, public :: OBS_CYLINDER = 2

end module module_constants
