module module_config
  use module_precision, only: dp
  use module_constants
  use module_json
  implicit none
  private

  public :: bc_t, config_t, load_config_from_json

  type :: bc_t
     integer  :: btype = BC_DIRICHLET
     real(dp) :: u = 0.0_dp
     real(dp) :: v = 0.0_dp
     integer  :: thermal_type = THERMAL_NEUMANN
     real(dp) :: T = 0.0_dp
  end type bc_t

  type :: config_t
     character(len=64) :: case_name = "fugu2d_case"
     character(len=32) :: mode = "NAVIER_STOKES"

     ! Mesh parameters
     integer  :: nx = 64, ny = 64
     real(dp) :: lx = 1.0_dp, ly = 1.0_dp

     ! Physical models & flags
     logical  :: enable_two_phase = .false.
     logical  :: enable_surface_tension = .false.
     logical  :: enable_convection = .true.
     integer  :: convective_scheme = CONV_UPWIND
     integer  :: viscosity_model = VISC_CONSTANT
     real(dp) :: visc_mu0 = 1.0e-3_dp, visc_a = 2.0_dp, visc_e = 0.01_dp
     real(dp) :: sigma = 0.0_dp
     real(dp) :: gx = 0.0_dp, gy = 0.0_dp
     real(dp) :: rho1 = 1.0_dp, mu1 = 1.0e-3_dp
     real(dp) :: rho2 = 1.0_dp, mu2 = 1.0e-3_dp

     ! Thermal parameters
     logical  :: enable_energy = .false.
     real(dp) :: Ra = 1.0e4_dp, Pr = 0.71_dp, T_ref = 0.0_dp, beta = 3.41e-3_dp
     real(dp) :: Cp = 1.0_dp, k_thermal = 0.014_dp

     ! Time integration
     real(dp) :: dt = 0.01_dp, t_final = 1.0_dp
     logical  :: adaptive_dt = .false.
     real(dp) :: cfl_max = 0.5_dp

     ! Boundary conditions
     type(bc_t) :: bc_left, bc_right, bc_bottom, bc_top

     ! Immersed boundary obstacle
     logical  :: obs_enabled = .false.
     integer  :: obs_type = OBS_NONE
     real(dp) :: obs_x0 = 0.0_dp, obs_y0 = 0.0_dp, obs_x1 = 0.0_dp, obs_y1 = 0.0_dp
     real(dp) :: obs_xc = 0.0_dp, obs_yc = 0.0_dp, obs_radius = 0.0_dp
     real(dp) :: obs_perm = 1.0e-40_dp

     ! Initial conditions
     integer  :: ic_type = IC_UNIFORM
     real(dp) :: ic_u = 0.0_dp, ic_v = 0.0_dp
     real(dp) :: ic_center_x = 0.5_dp, ic_center_y = 0.5_dp, ic_radius = 0.25_dp
     real(dp) :: ic_box_x0 = 0.0_dp, ic_box_y0 = 0.0_dp, ic_box_x1 = 0.5_dp, ic_box_y1 = 0.5_dp
     real(dp) :: ic_pert_amp = 0.05_dp, ic_pert_k = 1.0_dp

     ! Solvers
     integer  :: solver_type = SOLVER_AUGMENTED_LAGRANGIAN
     integer  :: linear_solver = LINEAR_BICGSTAB2
     integer  :: precond_type = PRECOND_JACOBI
     real(dp) :: tol = 1.0e-5_dp
     integer  :: max_iter = 300
     real(dp) :: r_al = 100.0_dp
     integer  :: uzawa_iter = 1
     logical  :: use_precond = .true.
     real(dp) :: steady_state_tol = 1.0e-5_dp

     ! Outputs
     character(len=128) :: output_dir = "output"
     integer  :: vtk_interval = 50
     integer  :: log_interval = 10
     logical  :: write_csv = .true.
  end type config_t

contains

  subroutine parse_bc(buf, prefix, bc)
    character(len=*), intent(in) :: buf, prefix
    type(bc_t), intent(out) :: bc
    character(len=32) :: btype_str, ttype_str

    call json_get_string(buf, trim(prefix) // ".type", btype_str, "DIRICHLET")
    if (trim(btype_str) == "DIRICHLET") then
       bc%btype = BC_DIRICHLET
    else if (trim(btype_str) == "NEUMANN") then
       bc%btype = BC_NEUMANN
    else if (trim(btype_str) == "PERIODIC") then
       bc%btype = BC_PERIODIC
    else if (trim(btype_str) == "FREE_SLIP") then
       bc%btype = BC_FREESLIP
    else
       bc%btype = BC_DIRICHLET
    end if

    call json_get_real(buf, trim(prefix) // ".u", bc%u, 0.0_dp)
    call json_get_real(buf, trim(prefix) // ".u_mean", bc%u, bc%u)
    call json_get_real(buf, trim(prefix) // ".u_val", bc%u, bc%u)
    call json_get_real(buf, trim(prefix) // ".v", bc%v, 0.0_dp)
    call json_get_real(buf, trim(prefix) // ".v_val", bc%v, bc%v)

    call json_get_string(buf, trim(prefix) // ".T_type", ttype_str, "NEUMANN")
    if (trim(ttype_str) == "DIRICHLET") then
       bc%thermal_type = THERMAL_DIRICHLET
    else
       bc%thermal_type = THERMAL_NEUMANN
    end if
    call json_get_real(buf, trim(prefix) // ".T", bc%T, 0.0_dp)
    call json_get_real(buf, trim(prefix) // ".T_val", bc%T, bc%T)
  end subroutine parse_bc

  subroutine load_config_from_json(filename, cfg, ierr)
    character(len=*), intent(in) :: filename
    type(config_t), intent(out) :: cfg
    integer, intent(out), optional :: ierr

    character(len=:), allocatable :: buf
    character(len=64) :: str_val
    integer :: stat

    stat = 0
    call json_read_file(filename, buf, stat)
    if (present(ierr)) ierr = stat
    if (stat /= 0) return

    ! Case name & mode
    call json_get_string(buf, "case_name", cfg%case_name, "fugu_case")
    call json_get_string(buf, "output.case_name", cfg%case_name, cfg%case_name)
    call json_get_string(buf, "mode", cfg%mode, "NAVIER_STOKES")

    ! Mesh
    call json_get_int(buf, "mesh.nx", cfg%nx, 64)
    call json_get_int(buf, "mesh.ny", cfg%ny, 64)
    call json_get_real(buf, "mesh.lx", cfg%lx, 1.0_dp)
    call json_get_real(buf, "mesh.ly", cfg%ly, 1.0_dp)

    ! Physics
    call json_get_bool(buf, "physics.enable_two_phase", cfg%enable_two_phase, .false.)
    call json_get_bool(buf, "physics.enable_surface_tension", cfg%enable_surface_tension, .false.)
    call json_get_bool(buf, "physics.enable_convection", cfg%enable_convection, .true.)
    call json_get_real(buf, "physics.sigma", cfg%sigma, 0.0_dp)
    call json_get_real(buf, "physics.gx", cfg%gx, 0.0_dp)
    call json_get_real(buf, "physics.gy", cfg%gy, 0.0_dp)

    call json_get_string(buf, "physics.convective_scheme", str_val, "UPWIND")
    if (trim(str_val) == "CENTERED") then
       cfg%convective_scheme = CONV_CENTERED
    else if (trim(str_val) == "TVD_SUPERBEE") then
       cfg%convective_scheme = CONV_TVD_SUPERBEE
    else
       cfg%convective_scheme = CONV_UPWIND
    end if

    call json_get_string(buf, "physics.viscosity_model", str_val, "CONSTANT")
    if (trim(str_val) == "EXPONENTIAL") then
       cfg%viscosity_model = VISC_EXPONENTIAL
       call json_get_real(buf, "physics.visc_mu0", cfg%visc_mu0, 1.0e-3_dp)
       call json_get_real(buf, "physics.visc_a", cfg%visc_a, 2.0_dp)
       call json_get_real(buf, "physics.visc_e", cfg%visc_e, 0.01_dp)
    else
       cfg%viscosity_model = VISC_CONSTANT
    end if

    call json_get_real(buf, "physics.fluid1.rho", cfg%rho1, 1.0_dp)
    call json_get_real(buf, "physics.fluid1.mu", cfg%mu1, 1.0e-3_dp)
    call json_get_real(buf, "physics.fluid2.rho", cfg%rho2, cfg%rho1)
    call json_get_real(buf, "physics.fluid2.mu", cfg%mu2, cfg%mu1)

    ! Thermal
    call json_get_bool(buf, "physics.enable_energy", cfg%enable_energy, .false.)
    call json_get_bool(buf, "thermal.enable_energy", cfg%enable_energy, cfg%enable_energy)
    call json_get_real(buf, "physics.Ra", cfg%Ra, 1.0e4_dp)
    call json_get_real(buf, "physics.Pr", cfg%Pr, 0.71_dp)
    call json_get_real(buf, "physics.T_ref", cfg%T_ref, 0.0_dp)
    call json_get_real(buf, "physics.beta", cfg%beta, 3.41e-3_dp)
    call json_get_real(buf, "physics.fluid1.Cp", cfg%Cp, 1.0_dp)
    call json_get_real(buf, "physics.fluid1.k_thermal", cfg%k_thermal, 0.014_dp)

    ! Time
    call json_get_real(buf, "time.dt", cfg%dt, 0.01_dp)
    call json_get_real(buf, "time.t_final", cfg%t_final, 1.0_dp)
    call json_get_bool(buf, "time.adaptive_dt", cfg%adaptive_dt, .false.)
    call json_get_real(buf, "time.cfl_max", cfg%cfl_max, 0.5_dp)
    call json_get_real(buf, "steady_state_tol", cfg%steady_state_tol, 1.0e-5_dp)

    ! Boundary conditions
    call parse_bc(buf, "boundary_conditions.left", cfg%bc_left)
    call parse_bc(buf, "boundary_conditions.right", cfg%bc_right)
    call parse_bc(buf, "boundary_conditions.bottom", cfg%bc_bottom)
    call parse_bc(buf, "boundary_conditions.top", cfg%bc_top)

    ! Obstacles
    call json_get_obstacle(buf, cfg%obs_type, cfg%obs_x0, cfg%obs_y0, cfg%obs_x1, cfg%obs_y1, &
                           cfg%obs_xc, cfg%obs_yc, cfg%obs_radius, cfg%obs_perm, cfg%obs_enabled)

    ! Initial condition
    call json_get_string(buf, "initial_condition.type", str_val, "UNIFORM")
    if (trim(str_val) == "UNIFORM") then
       cfg%ic_type = IC_UNIFORM
    else if (trim(str_val) == "POISEUILLE") then
       cfg%ic_type = IC_POISEUILLE
    else if (trim(str_val) == "VORTEX_SHEAR") then
       cfg%ic_type = IC_VORTEX_SHEAR
    else if (trim(str_val) == "RAYLEIGH_TAYLOR") then
       cfg%ic_type = IC_RAYLEIGH_TAYLOR
    else if (trim(str_val) == "DROPLET") then
       cfg%ic_type = IC_DROPLET
    else if (trim(str_val) == "BUBBLE") then
       cfg%ic_type = IC_BUBBLE
    else if (trim(str_val) == "DAM_BREAK") then
       cfg%ic_type = IC_DAM_BREAK
    else if (trim(str_val) == "LAPLACIAN_TEST") then
       cfg%ic_type = IC_LAPLACIAN_TEST
    else if (trim(str_val) == "SINGLE_PHASE") then
       cfg%ic_type = IC_UNIFORM
    end if

    call json_get_real(buf, "initial_condition.u", cfg%ic_u, 0.0_dp)
    call json_get_real(buf, "initial_condition.v", cfg%ic_v, 0.0_dp)
    call json_get_real(buf, "initial_condition.center_x", cfg%ic_center_x, 0.5_dp)
    call json_get_real(buf, "initial_condition.center_y", cfg%ic_center_y, 0.5_dp)
    call json_get_real(buf, "initial_condition.radius", cfg%ic_radius, 0.25_dp)
    call json_get_real(buf, "initial_condition.box_x0", cfg%ic_box_x0, 0.0_dp)
    call json_get_real(buf, "initial_condition.box_y0", cfg%ic_box_y0, 0.0_dp)
    call json_get_real(buf, "initial_condition.box_x1", cfg%ic_box_x1, 0.5_dp)
    call json_get_real(buf, "initial_condition.box_y1", cfg%ic_box_y1, 0.5_dp)
    call json_get_real(buf, "initial_condition.pert_amp", cfg%ic_pert_amp, 0.05_dp)
    call json_get_real(buf, "initial_condition.pert_k", cfg%ic_pert_k, 1.0_dp)

    ! Solver
    call json_get_string(buf, "solver.type", str_val, "AUGMENTED_LAGRANGIAN")
    if (trim(str_val) == "CONJUGATE_GRADIENT") then
       cfg%solver_type = SOLVER_AUGMENTED_LAGRANGIAN
       cfg%linear_solver = LINEAR_CG
    else
       cfg%solver_type = SOLVER_AUGMENTED_LAGRANGIAN
       cfg%linear_solver = LINEAR_BICGSTAB2
    end if

    call json_get_string(buf, "solver.linear_solver", str_val, "BICGSTAB2")
    if (trim(str_val) == "BICGSTAB") then
       cfg%linear_solver = LINEAR_BICGSTAB
    else if (trim(str_val) == "CG" .or. trim(str_val) == "CONJUGATE_GRADIENT") then
       cfg%linear_solver = LINEAR_CG
    else if (trim(str_val) == "GMRES") then
       cfg%linear_solver = LINEAR_GMRES
    end if

    call json_get_real(buf, "solver.tol", cfg%tol, 1.0e-5_dp)
    call json_get_int(buf, "solver.max_iter", cfg%max_iter, 300)
    call json_get_real(buf, "solver.r_al", cfg%r_al, 100.0_dp)
    call json_get_int(buf, "solver.uzawa_iter", cfg%uzawa_iter, 1)
    call json_get_bool(buf, "solver.use_precond", cfg%use_precond, .true.)
    call json_get_real(buf, "steady_state_tol", cfg%steady_state_tol, cfg%steady_state_tol)
    call json_get_real(buf, "solver.steady_state_tol", cfg%steady_state_tol, cfg%steady_state_tol)

    ! Output
    call json_get_string(buf, "output.output_dir", cfg%output_dir, "output")
    call json_get_string(buf, "output.directory", cfg%output_dir, cfg%output_dir)
    call json_get_int(buf, "output.vtk_interval", cfg%vtk_interval, 50)
    call json_get_int(buf, "output.log_interval", cfg%log_interval, 10)
    call json_get_bool(buf, "output.write_csv_history", cfg%write_csv, .true.)
  end subroutine load_config_from_json

end module module_config
