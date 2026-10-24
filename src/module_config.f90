! [Procédure/Module] module module_config
module module_config
  use module_precision, only: dp
  use module_constants
  use module_json
  implicit none
  private

  public :: bc_t, config_t, load_config_from_json

  type :: bc_t
     ! Déclaration de variable locale / paramètre de calcul
     integer  :: btype = BC_DIRICHLET
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: u = 0.0_dp
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: v = 0.0_dp
     ! Profil spécifique pour la frontière ("UNIFORM", "PARABOLIC", "JET_INLET")
     character(len=32) :: profile = "UNIFORM"
     ! Largeur de l'injecteur pour profil JET_INLET [m] (Cas 3.3 Fluids 2021)
     real(dp) :: jet_width = 0.08_dp
     ! Position centrale de l'injecteur selon X [m]
     real(dp) :: jet_center = 0.5_dp
     ! Déclaration de variable locale / paramètre de calcul
     integer  :: thermal_type = THERMAL_NEUMANN
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: T = 0.0_dp
  end type bc_t

  type :: config_t
     ! Déclaration de variable locale / paramètre de calcul
     character(len=64) :: case_name = "covof2d_case"
     ! Déclaration de variable locale / paramètre de calcul
     character(len=32) :: mode = "NAVIER_STOKES"

     ! Mesh parameters
     integer  :: nx = 64, ny = 64
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: lx = 1.0_dp, ly = 1.0_dp

     ! Physical models & flags
     logical  :: enable_two_phase = .false.
     ! Déclaration de variable locale / paramètre de calcul
     logical  :: enable_surface_tension = .false.
     ! Déclaration de variable locale / paramètre de calcul
     logical  :: enable_convection = .true.
     ! Déclaration de variable locale / paramètre de calcul
     integer  :: convective_scheme = CONV_UPWIND
     ! Déclaration de variable locale / paramètre de calcul
     integer  :: viscosity_model = VISC_CONSTANT
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: visc_mu0 = 1.0e-3_dp, visc_a = 2.0_dp, visc_e = 0.01_dp
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: sigma = 0.0_dp
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: gx = 0.0_dp, gy = 0.0_dp
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: rho1 = 1.0_dp, mu1 = 1.0e-3_dp
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: rho2 = 1.0_dp, mu2 = 1.0e-3_dp

     ! VOF and interface tracking
     character(len=32) :: reconstruction_method = "PLIC"
     ! Déclaration de variable locale / paramètre de calcul
     character(len=32) :: advection_method = "VOF_TVD"
     ! Déclaration de variable locale / paramètre de calcul
     character(len=32) :: tvd_limiter = "SUPERBEE"
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: cfl_vof = 0.4_dp

     ! Thermal parameters
     logical  :: enable_energy = .false.
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: Ra = 1.0e4_dp, Pr = 0.71_dp, T_ref = 0.0_dp, beta = 3.41e-3_dp
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: Cp = 1.0_dp, k_thermal = 0.014_dp

     ! Time integration
     real(dp) :: dt = 0.01_dp, t_final = 1.0_dp
     ! Déclaration de variable locale / paramètre de calcul
     logical  :: adaptive_dt = .false.
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: cfl_max = 0.5_dp

     ! Boundary conditions
     type(bc_t) :: bc_left, bc_right, bc_bottom, bc_top

     ! Immersed boundary obstacle
     logical  :: obs_enabled = .false.
     ! Déclaration de variable locale / paramètre de calcul
     integer  :: obs_type = OBS_NONE
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: obs_x0 = 0.0_dp, obs_y0 = 0.0_dp, obs_x1 = 0.0_dp, obs_y1 = 0.0_dp
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: obs_xc = 0.0_dp, obs_yc = 0.0_dp, obs_radius = 0.0_dp
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: obs_perm = 1.0e-40_dp

     ! Initial conditions
     integer  :: ic_type = IC_UNIFORM
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: ic_u = 0.0_dp, ic_v = 0.0_dp
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: ic_center_x = 0.5_dp, ic_center_y = 0.5_dp, ic_radius = 0.25_dp
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: ic_box_x0 = 0.0_dp, ic_box_y0 = 0.0_dp, ic_box_x1 = 0.5_dp, ic_box_y1 = 0.5_dp
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: ic_pert_amp = 0.05_dp, ic_pert_k = 1.0_dp
     ! Déclaration de l'épaisseur initiale de nappe liquide (Cas 3.2 Fluids 2021)
     real(dp) :: ic_thickness = 0.0003_dp
     ! Déclaration de la largeur initiale de jet liquide (Cas 3.3 Fluids 2021)
     real(dp) :: ic_jet_width = 0.08_dp

     ! Solvers
     integer  :: solver_type = SOLVER_AUGMENTED_LAGRANGIAN
     ! Déclaration de variable locale / paramètre de calcul
     integer  :: linear_solver = LINEAR_BICGSTAB2
     ! Déclaration de variable locale / paramètre de calcul
     integer  :: precond_type = PRECOND_JACOBI
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: tol = 1.0e-5_dp
     ! Déclaration de variable locale / paramètre de calcul
     integer  :: max_iter = 300
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: r_al = 100.0_dp
     ! Déclaration de variable locale / paramètre de calcul
     integer  :: uzawa_iter = 1
     ! Déclaration de variable locale / paramètre de calcul
     logical  :: use_precond = .true.
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: steady_state_tol = 1.0e-5_dp

     ! Outputs
     character(len=128) :: output_dir = "output"
     ! Déclaration de variable locale / paramètre de calcul
     integer  :: vtk_interval = 50
     ! Déclaration de variable locale / paramètre de calcul
     integer  :: log_interval = 10
     ! Déclaration de variable locale / paramètre de calcul
     logical  :: write_csv = .true.
  end type config_t

contains

  ! [Procédure/Module] subroutine parse_bc
  subroutine parse_bc(buf, prefix, bc)
    ! Déclaration de variable locale / paramètre de calcul
    character(len=*), intent(in) :: buf, prefix
    ! Déclaration de variable locale / paramètre de calcul
    type(bc_t), intent(out) :: bc
    ! Déclaration de variable locale / paramètre de calcul
    character(len=32) :: btype_str, ttype_str

    ! Appel de sous-programme modulaire
    call json_get_string(buf, trim(prefix) // ".type", btype_str, "DIRICHLET")
    ! Évaluation de la condition logique
    if (trim(btype_str) == "DIRICHLET") then
       ! Affectation et calcul numérique
       bc%btype = BC_DIRICHLET
    ! Branche conditionnelle alternative
    else if (trim(btype_str) == "NEUMANN") then
       ! Affectation et calcul numérique
       bc%btype = BC_NEUMANN
    ! Branche conditionnelle alternative
    else if (trim(btype_str) == "PERIODIC") then
       ! Affectation et calcul numérique
       bc%btype = BC_PERIODIC
    ! Branche conditionnelle alternative
    else if (trim(btype_str) == "FREE_SLIP") then
       ! Affectation et calcul numérique
       bc%btype = BC_FREESLIP
    ! Branche par défaut si condition non satisfaite
    else
       ! Affectation et calcul numérique
       bc%btype = BC_DIRICHLET
    end if

    ! Appel de sous-programme modulaire
    call json_get_real(buf, trim(prefix) // ".u", bc%u, 0.0_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, trim(prefix) // ".u_mean", bc%u, bc%u)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, trim(prefix) // ".u_val", bc%u, bc%u)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, trim(prefix) // ".v", bc%v, 0.0_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, trim(prefix) // ".v_val", bc%v, bc%v)
    ! Lecture du profil de frontière (Cas 3.3 Fluids 2021)
    call json_get_string(buf, trim(prefix) // ".profile", bc%profile, "UNIFORM")
    ! Lecture de la largeur de l'injecteur
    call json_get_real(buf, trim(prefix) // ".jet_width", bc%jet_width, 0.08_dp)
    ! Lecture du centre de l'injecteur
    call json_get_real(buf, trim(prefix) // ".jet_center", bc%jet_center, 0.5_dp)

    ! Appel de sous-programme modulaire
    call json_get_string(buf, trim(prefix) // ".T_type", ttype_str, "NEUMANN")
    ! Évaluation de la condition logique
    if (trim(ttype_str) == "DIRICHLET") then
       ! Affectation et calcul numérique
       bc%thermal_type = THERMAL_DIRICHLET
    ! Branche par défaut si condition non satisfaite
    else
       ! Affectation et calcul numérique
       bc%thermal_type = THERMAL_NEUMANN
    end if
    ! Appel de sous-programme modulaire
    call json_get_real(buf, trim(prefix) // ".T", bc%T, 0.0_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, trim(prefix) // ".T_val", bc%T, bc%T)
  end subroutine parse_bc

  ! [Procédure/Module] subroutine load_config_from_json
  subroutine load_config_from_json(filename, cfg, ierr)
    ! Déclaration de variable locale / paramètre de calcul
    character(len=*), intent(in) :: filename
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(out) :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(out), optional :: ierr

    ! Déclaration de variable locale / paramètre de calcul
    character(len=:), allocatable :: buf
    ! Déclaration de variable locale / paramètre de calcul
    character(len=64) :: str_val
    ! Déclaration de variable locale / paramètre de calcul
    integer :: stat

    ! Affectation et calcul numérique
    stat = 0
    ! Appel de sous-programme modulaire
    call json_read_file(filename, buf, stat)
    ! Évaluation de la norme du résidu ou de l'erreur numérique
    if (present(ierr)) ierr = stat
    ! Affectation et calcul numérique
    if (stat /= 0) return

    ! Case name & mode
    call json_get_string(buf, "case_name", cfg%case_name, "covof_case")
    ! Appel de sous-programme modulaire
    call json_get_string(buf, "output.case_name", cfg%case_name, cfg%case_name)
    ! Appel de sous-programme modulaire
    call json_get_string(buf, "mode", cfg%mode, "NAVIER_STOKES")

    ! Mesh
    call json_get_int(buf, "mesh.nx", cfg%nx, 64)
    ! Appel de sous-programme modulaire
    call json_get_int(buf, "mesh.ny", cfg%ny, 64)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "mesh.lx", cfg%lx, 1.0_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "mesh.ly", cfg%ly, 1.0_dp)

    ! Physics
    call json_get_bool(buf, "physics.enable_two_phase", cfg%enable_two_phase, .false.)
    ! Appel de sous-programme modulaire
    call json_get_bool(buf, "physics.enable_surface_tension", cfg%enable_surface_tension, .false.)
    ! Appel de sous-programme modulaire
    call json_get_bool(buf, "physics.enable_convection", cfg%enable_convection, .true.)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "physics.sigma", cfg%sigma, 0.0_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "physics.gx", cfg%gx, 0.0_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "physics.gy", cfg%gy, 0.0_dp)

    ! Appel de sous-programme modulaire
    call json_get_string(buf, "physics.convective_scheme", str_val, "UPWIND")
    ! Évaluation de la condition logique
    if (trim(str_val) == "CENTERED") then
       ! Affectation et calcul numérique
       cfg%convective_scheme = CONV_CENTERED
    ! Branche conditionnelle alternative
    else if (trim(str_val) == "TVD_SUPERBEE") then
       ! Affectation et calcul numérique
       cfg%convective_scheme = CONV_TVD_SUPERBEE
    ! Branche par défaut si condition non satisfaite
    else
       ! Affectation et calcul numérique
       cfg%convective_scheme = CONV_UPWIND
    end if

    ! Appel de sous-programme modulaire
    call json_get_string(buf, "physics.viscosity_model", str_val, "CONSTANT")
    ! Évaluation de la condition logique
    if (trim(str_val) == "EXPONENTIAL") then
       ! Affectation et calcul numérique
       cfg%viscosity_model = VISC_EXPONENTIAL
       ! Appel de sous-programme modulaire
       call json_get_real(buf, "physics.visc_mu0", cfg%visc_mu0, 1.0e-3_dp)
       ! Appel de sous-programme modulaire
       call json_get_real(buf, "physics.visc_a", cfg%visc_a, 2.0_dp)
       ! Appel de sous-programme modulaire
       call json_get_real(buf, "physics.visc_e", cfg%visc_e, 0.01_dp)
    ! Branche par défaut si condition non satisfaite
    else
       ! Affectation et calcul numérique
       cfg%viscosity_model = VISC_CONSTANT
    end if

    ! Appel de sous-programme modulaire
    call json_get_real(buf, "physics.fluid1.rho", cfg%rho1, 1.0_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "physics.fluid1.mu", cfg%mu1, 1.0e-3_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "physics.fluid2.rho", cfg%rho2, cfg%rho1)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "physics.fluid2.mu", cfg%mu2, cfg%mu1)

    ! Thermal
    call json_get_bool(buf, "physics.enable_energy", cfg%enable_energy, .false.)
    ! Appel de sous-programme modulaire
    call json_get_bool(buf, "thermal.enable_energy", cfg%enable_energy, cfg%enable_energy)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "physics.Ra", cfg%Ra, 1.0e4_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "physics.Pr", cfg%Pr, 0.71_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "physics.T_ref", cfg%T_ref, 0.0_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "physics.beta", cfg%beta, 3.41e-3_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "physics.fluid1.Cp", cfg%Cp, 1.0_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "physics.fluid1.k_thermal", cfg%k_thermal, 0.014_dp)

    ! Time
    call json_get_real(buf, "time.dt", cfg%dt, 0.01_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "time.t_final", cfg%t_final, 1.0_dp)
    ! Appel de sous-programme modulaire
    call json_get_bool(buf, "time.adaptive_dt", cfg%adaptive_dt, .false.)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "time.cfl_max", cfg%cfl_max, 0.5_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "steady_state_tol", cfg%steady_state_tol, 1.0e-5_dp)

    ! Boundary conditions
    call parse_bc(buf, "boundary_conditions.left", cfg%bc_left)
    ! Appel de sous-programme modulaire
    call parse_bc(buf, "boundary_conditions.right", cfg%bc_right)
    ! Appel de sous-programme modulaire
    call parse_bc(buf, "boundary_conditions.bottom", cfg%bc_bottom)
    ! Appel de sous-programme modulaire
    call parse_bc(buf, "boundary_conditions.top", cfg%bc_top)

    ! Obstacles
    call json_get_obstacle(buf, cfg%obs_type, cfg%obs_x0, cfg%obs_y0, cfg%obs_x1, cfg%obs_y1, &
                           cfg%obs_xc, cfg%obs_yc, cfg%obs_radius, cfg%obs_perm, cfg%obs_enabled)

    ! Initial condition
    call json_get_string(buf, "initial_condition.type", str_val, "UNIFORM")
    ! Évaluation de la condition logique
    if (trim(str_val) == "UNIFORM") then
       ! Affectation et calcul numérique
       cfg%ic_type = IC_UNIFORM
    ! Branche conditionnelle alternative
    else if (trim(str_val) == "POISEUILLE") then
       ! Affectation et calcul numérique
       cfg%ic_type = IC_POISEUILLE
    ! Branche conditionnelle alternative
    else if (trim(str_val) == "VORTEX_SHEAR") then
       ! Affectation et calcul numérique
       cfg%ic_type = IC_VORTEX_SHEAR
    ! Branche conditionnelle alternative
    else if (trim(str_val) == "RAYLEIGH_TAYLOR") then
       ! Affectation et calcul numérique
       cfg%ic_type = IC_RAYLEIGH_TAYLOR
    ! Branche conditionnelle alternative
    else if (trim(str_val) == "DROPLET") then
       ! Affectation et calcul numérique
       cfg%ic_type = IC_DROPLET
    ! Branche conditionnelle alternative pour bulle ou cylindre immergé
    else if (trim(str_val) == "BUBBLE" .or. trim(str_val) == "CYLINDER") then
       ! Affectation et calcul numérique
       cfg%ic_type = IC_BUBBLE
    ! Branche conditionnelle alternative
    else if (trim(str_val) == "DAM_BREAK") then
       ! Affectation et calcul numérique
       cfg%ic_type = IC_DAM_BREAK
    ! Branche conditionnelle alternative pour nappe liquide (Cas 3.2 Fluids 2021)
    else if (trim(str_val) == "LIQUID_SHEET") then
       ! Affectation et calcul numérique
       cfg%ic_type = IC_LIQUID_SHEET
    ! Branche conditionnelle alternative pour jet liquide visqueux (Cas 3.3 Fluids 2021)
    else if (trim(str_val) == "VISCOUS_JET" .or. trim(str_val) == "JET_BUCKLING") then
       ! Affectation et calcul numérique
       cfg%ic_type = IC_VISCOUS_JET
    ! Branche conditionnelle alternative
    else if (trim(str_val) == "LAPLACIAN_TEST") then
       ! Affectation et calcul numérique
       cfg%ic_type = IC_LAPLACIAN_TEST
    ! Branche conditionnelle alternative
    else if (trim(str_val) == "SINGLE_PHASE") then
       ! Affectation et calcul numérique
       cfg%ic_type = IC_UNIFORM
    end if

    ! Appel de sous-programme modulaire
    call json_get_real(buf, "initial_condition.u", cfg%ic_u, 0.0_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "initial_condition.v", cfg%ic_v, 0.0_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "initial_condition.center_x", cfg%ic_center_x, 0.5_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "initial_condition.center_y", cfg%ic_center_y, 0.5_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "initial_condition.radius", cfg%ic_radius, 0.25_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "initial_condition.box_x0", cfg%ic_box_x0, 0.0_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "initial_condition.box_y0", cfg%ic_box_y0, 0.0_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "initial_condition.box_x1", cfg%ic_box_x1, 0.5_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "initial_condition.box_y1", cfg%ic_box_y1, 0.5_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "initial_condition.pert_amp", cfg%ic_pert_amp, 0.05_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "initial_condition.pert_k", cfg%ic_pert_k, 1.0_dp)
    ! Lecture de l'épaisseur de nappe liquide (Cas 3.2 Fluids 2021)
    call json_get_real(buf, "initial_condition.thickness", cfg%ic_thickness, 0.0003_dp)
    ! Lecture de la largeur de jet liquide (Cas 3.3 Fluids 2021)
    call json_get_real(buf, "initial_condition.jet_width", cfg%ic_jet_width, 0.08_dp)

    ! Solver
    call json_get_string(buf, "solver.type", str_val, "AUGMENTED_LAGRANGIAN")
    ! Évaluation de la condition logique
    if (trim(str_val) == "CONJUGATE_GRADIENT") then
       ! Affectation et calcul numérique
       cfg%solver_type = SOLVER_AUGMENTED_LAGRANGIAN
       ! Affectation et calcul numérique
       cfg%linear_solver = LINEAR_CG
    ! Branche par défaut si condition non satisfaite
    else
       ! Affectation et calcul numérique
       cfg%solver_type = SOLVER_AUGMENTED_LAGRANGIAN
       ! Affectation et calcul numérique
       cfg%linear_solver = LINEAR_BICGSTAB2
    end if

    ! Appel de sous-programme modulaire
    call json_get_string(buf, "solver.linear_solver", str_val, "BICGSTAB2")
    ! Évaluation de la condition logique
    if (trim(str_val) == "BICGSTAB") then
       ! Affectation et calcul numérique
       cfg%linear_solver = LINEAR_BICGSTAB
    ! Branche conditionnelle alternative
    else if (trim(str_val) == "CG" .or. trim(str_val) == "CONJUGATE_GRADIENT") then
       ! Affectation et calcul numérique
       cfg%linear_solver = LINEAR_CG
    ! Branche conditionnelle alternative
    else if (trim(str_val) == "GMRES") then
       ! Affectation et calcul numérique
       cfg%linear_solver = LINEAR_GMRES
    end if

    ! Appel de sous-programme modulaire
    call json_get_real(buf, "solver.tol", cfg%tol, 1.0e-5_dp)
    ! Appel de sous-programme modulaire
    call json_get_int(buf, "solver.max_iter", cfg%max_iter, 300)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "solver.r_al", cfg%r_al, 100.0_dp)
    ! Appel de sous-programme modulaire
    call json_get_int(buf, "solver.uzawa_iter", cfg%uzawa_iter, 1)
    ! Appel de sous-programme modulaire
    call json_get_bool(buf, "solver.use_precond", cfg%use_precond, .true.)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "steady_state_tol", cfg%steady_state_tol, cfg%steady_state_tol)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "solver.steady_state_tol", cfg%steady_state_tol, cfg%steady_state_tol)

    ! VOF parameters
    call json_get_string(buf, "vof.advection_method", cfg%advection_method, "WEYMOUTH_YUE")
    ! Appel de sous-programme modulaire
    call json_get_string(buf, "vof.reconstruction_method", cfg%reconstruction_method, "PLIC")
    ! Appel de sous-programme modulaire
    call json_get_string(buf, "vof.tvd_limiter", cfg%tvd_limiter, "NONE")
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "vof.cfl_vof", cfg%cfl_vof, 0.5_dp)

    ! Cohérence stricte : le schéma algébrique TVD ne s'accompagne d'aucune reconstruction géométrique,
    ! tandis que les schémas géométriques (WEYMOUTH_YUE, CIAM) n'utilisent pas de limiteur TVD.
    if (trim(cfg%advection_method) == "VOF_TVD" .or. trim(cfg%advection_method) == "TVD") then
       cfg%reconstruction_method = "NONE"
       if (trim(cfg%tvd_limiter) == "NONE" .or. len_trim(cfg%tvd_limiter) == 0) then
          cfg%tvd_limiter = "SUPERBEE"
       end if
    else
       if (trim(cfg%reconstruction_method) == "NONE" .or. len_trim(cfg%reconstruction_method) == 0) then
          cfg%reconstruction_method = "PLIC"
       end if
       cfg%tvd_limiter = "NONE"
    end if

    ! Output
    call json_get_string(buf, "output.output_dir", cfg%output_dir, "output")
    ! Appel de sous-programme modulaire
    call json_get_string(buf, "output.directory", cfg%output_dir, cfg%output_dir)
    ! Appel de sous-programme modulaire
    call json_get_int(buf, "output.vtk_interval", cfg%vtk_interval, 50)
    ! Appel de sous-programme modulaire
    call json_get_int(buf, "output.log_interval", cfg%log_interval, 10)
    ! Appel de sous-programme modulaire
    call json_get_bool(buf, "output.write_csv_history", cfg%write_csv, .true.)
  end subroutine load_config_from_json

end module module_config
