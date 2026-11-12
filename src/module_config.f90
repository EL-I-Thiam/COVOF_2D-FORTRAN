! [Procédure/Module] module module_config
module module_config
  use module_precision, only: dp
  use module_constants
  use module_json
  use module_eos, only: eos_parse_type
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
     character(len=64) :: test_case = "CAS_3_1"
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
     integer  :: convective_linearization = CONV_LIN_PICARD
     logical  :: enable_svof_smoothing = .false.
     real(dp) :: svof_filter_width = 1.0_dp
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

     ! Formulation All-Mach / Compressible (Fluids 2024, El Ouafa et al.)
     logical  :: compressible = .false.
     character(len=32) :: flow_regime = "INCOMPRESSIBLE"
     logical  :: enable_ssp_rk3_momentum = .false.
     logical  :: enable_diffusion_smoothing = .false.
     integer  :: diffusion_smoothing_steps = 3

     ! Fluid 1 EOS
     integer  :: eos1_type = EOS_INCOMPRESSIBLE
     real(dp) :: chi_T1 = 0.0_dp, gamma1 = 1.4_dp, r_gas1 = 287.05_dp, cv1 = 717.5_dp
     real(dp) :: p_inf1 = 0.0_dp, b_nasg1 = 0.0_dp, q_nasg1 = 0.0_dp

     ! Fluid 2 EOS
     integer  :: eos2_type = EOS_INCOMPRESSIBLE
     real(dp) :: chi_T2 = 0.0_dp, gamma2 = 1.4_dp, r_gas2 = 287.05_dp, cv2 = 717.5_dp
     real(dp) :: p_inf2 = 0.0_dp, b_nasg2 = 0.0_dp, q_nasg2 = 0.0_dp

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
     integer  :: time_scheme = TIME_SCHEME_EULER
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
     real(dp) :: ic_p0 = 0.0_dp

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
     logical  :: adaptative_al = .false.     ! Activer le Lagrangien augmenté adaptatif r(t, M) (Éq. 14)
     real(dp) :: k_aal         = 100.0_dp    ! Facteur d'amplification K in [10, 1000] (Éq. 14)
     real(dp) :: uzawa_tol     = 1.0e-12_dp  ! Seuil d'arrêt d'Uzawa sur ||div(u^{*,m})||_inf (Éq. 10)
     real(dp) :: l0_ref        = 0.0_dp      ! Échelle de longueur L0 [m] (0 = auto)
     real(dp) :: t0_ref        = 0.0_dp      ! Échelle de temps t0 [s] (0 = auto = dt)
     real(dp) :: u0_ref        = 0.0_dp      ! Échelle de vitesse u0 [m/s] (0 = auto)
     real(dp) :: p0_ref        = 0.0_dp      ! Échelle de pression p0 [Pa] (0 = auto)
     ! Déclaration de variable locale / paramètre de calcul
     integer  :: uzawa_iter = 1
     ! Déclaration de variable locale / paramètre de calcul
     logical  :: use_precond = .true.
     ! Déclaration de variable locale / paramètre de calcul
     real(dp) :: steady_state_tol = 0.0_dp
     ! Terme source pour le mode LAPLACIAN
     real(dp) :: source_term = 1.0_dp

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
    call json_get_string(buf, "output.case_name", cfg%case_name, cfg%case_name)
    call json_get_string(buf, "test_case", cfg%test_case, "CAS_3_1")
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

    call json_get_string(buf, "physics.convective_linearization", str_val, "PICARD")
    if (trim(str_val) == "PICARD") then
       call json_get_string(buf, "convective_linearization", str_val, "PICARD")
    end if
    if (trim(str_val) == "ADAMS_BASHFORTH" .or. trim(str_val) == "AB2") then
       cfg%convective_linearization = CONV_LIN_ADAMS_BASHFORTH
    else if (trim(str_val) == "NEWTON") then
       cfg%convective_linearization = CONV_LIN_NEWTON
    else
       cfg%convective_linearization = CONV_LIN_PICARD
    end if

    call json_get_bool(buf, "physics.enable_svof_smoothing", cfg%enable_svof_smoothing, .false.)
    if (.not. cfg%enable_svof_smoothing) then
       call json_get_bool(buf, "enable_svof_smoothing", cfg%enable_svof_smoothing, .false.)
    end if
    call json_get_real(buf, "physics.svof_filter_width", cfg%svof_filter_width, 1.0_dp)
    if (cfg%svof_filter_width == 1.0_dp) then
       call json_get_real(buf, "svof_filter_width", cfg%svof_filter_width, 1.0_dp)
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

    ! Formulation Compressible / All-Mach (Fluids 2024)
    call json_get_bool(buf, "physics.compressible", cfg%compressible, .false.)
    call json_get_string(buf, "physics.flow_regime", str_val, "INCOMPRESSIBLE")
    cfg%flow_regime = trim(str_val)
    if (trim(str_val) == "COMPRESSIBLE" .or. trim(str_val) == "ALL_MACH") then
       cfg%compressible = .true.
    end if
    call json_get_bool(buf, "physics.enable_ssp_rk3_momentum", cfg%enable_ssp_rk3_momentum, cfg%compressible)
    call json_get_bool(buf, "physics.enable_diffusion_smoothing", cfg%enable_diffusion_smoothing, .false.)
    call json_get_int(buf, "physics.diffusion_smoothing_steps", cfg%diffusion_smoothing_steps, 3)

    ! EOS Fluide 1
    call json_get_string(buf, "physics.fluid1.eos", str_val, "INCOMPRESSIBLE")
    cfg%eos1_type = eos_parse_type(str_val)
    call json_get_real(buf, "physics.fluid1.chi_T", cfg%chi_T1, 0.0_dp)
    call json_get_real(buf, "physics.fluid1.gamma", cfg%gamma1, 1.4_dp)
    call json_get_real(buf, "physics.fluid1.r_gas", cfg%r_gas1, 287.05_dp)
    call json_get_real(buf, "physics.fluid1.cv", cfg%cv1, 717.5_dp)
    call json_get_real(buf, "physics.fluid1.p_inf", cfg%p_inf1, 0.0_dp)
    call json_get_real(buf, "physics.fluid1.b_nasg", cfg%b_nasg1, 0.0_dp)
    call json_get_real(buf, "physics.fluid1.q_nasg", cfg%q_nasg1, 0.0_dp)

    ! EOS Fluide 2
    call json_get_string(buf, "physics.fluid2.eos", str_val, "INCOMPRESSIBLE")
    cfg%eos2_type = eos_parse_type(str_val)
    call json_get_real(buf, "physics.fluid2.chi_T", cfg%chi_T2, 0.0_dp)
    call json_get_real(buf, "physics.fluid2.gamma", cfg%gamma2, 1.4_dp)
    call json_get_real(buf, "physics.fluid2.r_gas", cfg%r_gas2, 287.05_dp)
    call json_get_real(buf, "physics.fluid2.cv", cfg%cv2, 717.5_dp)
    call json_get_real(buf, "physics.fluid2.p_inf", cfg%p_inf2, 0.0_dp)
    call json_get_real(buf, "physics.fluid2.b_nasg", cfg%b_nasg2, 0.0_dp)
    call json_get_real(buf, "physics.fluid2.q_nasg", cfg%q_nasg2, 0.0_dp)

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
    call json_get_string(buf, "time.time_scheme", str_val, "EULER")
    if (trim(str_val) == "EULER") then
       call json_get_string(buf, "time_scheme", str_val, "EULER")
    end if
    if (trim(str_val) == "GEAR" .or. trim(str_val) == "BDF2") then
       cfg%time_scheme = TIME_SCHEME_GEAR
    else
       cfg%time_scheme = TIME_SCHEME_EULER
    end if
    call json_get_real(buf, "time.dt", cfg%dt, 0.01_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "time.t_final", cfg%t_final, 1.0_dp)
    ! Appel de sous-programme modulaire
    call json_get_bool(buf, "time.adaptive_dt", cfg%adaptive_dt, .false.)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "time.cfl_max", cfg%cfl_max, 0.5_dp)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "steady_state_tol", cfg%steady_state_tol, 0.0_dp)

    ! Terme source pour le mode LAPLACIAN
    call json_get_real(buf, "physics.source_term", cfg%source_term, 1.0_dp)
    call json_get_real(buf, "physics.f", cfg%source_term, cfg%source_term)
    call json_get_real(buf, "source_term", cfg%source_term, cfg%source_term)
    call json_get_real(buf, "f", cfg%source_term, cfg%source_term)
    call json_get_real(buf, "laplacian.source_term", cfg%source_term, cfg%source_term)
    call json_get_real(buf, "laplacian.f", cfg%source_term, cfg%source_term)

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
    call json_get_string(buf, "initial_condition.type", str_val, "")
    if (len_trim(str_val) == 0) then
       call json_get_string(buf, "initial_conditions.type", str_val, "UNIFORM")
    end if
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
    ! Branche conditionnelle alternative pour cylindre immergé (Cas 3.1 Fluids 2021)
    else if (trim(str_val) == "CYLINDER") then
       ! Affectation et calcul numérique
       cfg%ic_type = IC_CYLINDER
    ! Branche conditionnelle alternative pour bulle
    else if (trim(str_val) == "BUBBLE") then
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
    call json_get_real(buf, "initial_conditions.u", cfg%ic_u, cfg%ic_u)
    call json_get_real(buf, "initial_condition.v", cfg%ic_v, 0.0_dp)
    call json_get_real(buf, "initial_conditions.v", cfg%ic_v, cfg%ic_v)
    call json_get_real(buf, "initial_condition.center_x", cfg%ic_center_x, 0.5_dp)
    call json_get_real(buf, "initial_conditions.center_x", cfg%ic_center_x, cfg%ic_center_x)
    call json_get_real(buf, "initial_condition.center_y", cfg%ic_center_y, 0.5_dp)
    call json_get_real(buf, "initial_conditions.center_y", cfg%ic_center_y, cfg%ic_center_y)
    call json_get_real(buf, "initial_condition.radius", cfg%ic_radius, 0.25_dp)
    call json_get_real(buf, "initial_conditions.radius", cfg%ic_radius, cfg%ic_radius)
    call json_get_real(buf, "initial_condition.box_x0", cfg%ic_box_x0, 0.0_dp)
    call json_get_real(buf, "initial_conditions.box_x0", cfg%ic_box_x0, cfg%ic_box_x0)
    call json_get_real(buf, "initial_condition.box_y0", cfg%ic_box_y0, 0.0_dp)
    call json_get_real(buf, "initial_conditions.box_y0", cfg%ic_box_y0, cfg%ic_box_y0)
    call json_get_real(buf, "initial_condition.box_x1", cfg%ic_box_x1, 0.5_dp)
    call json_get_real(buf, "initial_conditions.box_x1", cfg%ic_box_x1, cfg%ic_box_x1)
    call json_get_real(buf, "initial_condition.box_y1", cfg%ic_box_y1, 0.5_dp)
    call json_get_real(buf, "initial_conditions.box_y1", cfg%ic_box_y1, cfg%ic_box_y1)
    call json_get_real(buf, "initial_condition.pert_amp", cfg%ic_pert_amp, 0.05_dp)
    call json_get_real(buf, "initial_conditions.pert_amp", cfg%ic_pert_amp, cfg%ic_pert_amp)
    call json_get_real(buf, "initial_condition.pert_k", cfg%ic_pert_k, 1.0_dp)
    call json_get_real(buf, "initial_conditions.pert_k", cfg%ic_pert_k, cfg%ic_pert_k)
    call json_get_real(buf, "initial_condition.thickness", cfg%ic_thickness, 0.0003_dp)
    call json_get_real(buf, "initial_conditions.thickness", cfg%ic_thickness, cfg%ic_thickness)
    call json_get_real(buf, "initial_condition.jet_width", cfg%ic_jet_width, 0.08_dp)
    call json_get_real(buf, "initial_conditions.jet_width", cfg%ic_jet_width, cfg%ic_jet_width)
    call json_get_real(buf, "initial_condition.p_init", cfg%ic_p0, 0.0_dp)
    call json_get_real(buf, "initial_condition.p0", cfg%ic_p0, cfg%ic_p0)
    call json_get_real(buf, "initial_conditions.p_init", cfg%ic_p0, cfg%ic_p0)
    call json_get_real(buf, "initial_conditions.p0", cfg%ic_p0, cfg%ic_p0)

    ! Solver
    call json_get_string(buf, "solver.type", str_val, "AUGMENTED_LAGRANGIAN")
    ! Évaluation de la condition logique
    if (trim(str_val) == "CONJUGATE_GRADIENT") then
       ! Affectation et calcul numérique
       cfg%solver_type = SOLVER_AUGMENTED_LAGRANGIAN
       ! Affectation et calcul numérique
       cfg%linear_solver = LINEAR_CG
    else if (trim(str_val) == "MONOLITHIC" .or. trim(str_val) == "FULLY_COUPLED" .or. trim(str_val) == "FC") then
       ! Affectation et calcul numérique
       cfg%solver_type = SOLVER_MONOLITHIC
       ! Affectation et calcul numérique
       cfg%linear_solver = LINEAR_BICGSTAB2
    else if (trim(str_val) == "STANDARD_PROJECTION" .or. trim(str_val) == "PROJECTION") then
       ! Affectation et calcul numérique
       cfg%solver_type = SOLVER_STANDARD_PROJECTION
       ! Affectation et calcul numérique
       cfg%linear_solver = LINEAR_BICGSTAB2
    else if (trim(str_val) == "C-AAL" .or. trim(str_val) == "CAAL" .or. &
             trim(str_val) == "COMPRESSIBLE_AAL") then
       ! Formulation C-AAL Compressible Bas-Mach (Fluids 2024, El Ouafa et al.)
       cfg%solver_type = SOLVER_CAAL
       cfg%adaptative_al = .true.
       cfg%compressible = .true.
       cfg%linear_solver = LINEAR_BICGSTAB2
    else if (trim(str_val) == "ADAPTATIVE_AUGMENTED_LAGRANGIAN" .or. &
             trim(str_val) == "ADAPTIVE_AUGMENTED_LAGRANGIAN" .or. &
             trim(str_val) == "AAL") then
       ! Formulation AAL Incompressible pure (Vincent et al. 2004)
       cfg%solver_type = SOLVER_AAL
       cfg%adaptative_al = .true.
       cfg%linear_solver = LINEAR_BICGSTAB2
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
    else if (trim(str_val) == "MINRES") then
       ! Solveur MINRES (Paige & Saunders 1975)
       cfg%linear_solver = LINEAR_MINRES
    ! Branche conditionnelle alternative
    else if (trim(str_val) == "DIRECT" .or. trim(str_val) == "MUMPS") then
       ! Affectation et calcul numérique
       cfg%linear_solver = LINEAR_DIRECT
    end if

    ! Appel de sous-programme modulaire
    call json_get_real(buf, "solver.tol", cfg%tol, 1.0e-5_dp)
    ! Appel de sous-programme modulaire
    call json_get_int(buf, "solver.max_iter", cfg%max_iter, 300)
    ! Appel de sous-programme modulaire
    call json_get_real(buf, "solver.r_al", cfg%r_al, 100.0_dp)
    call json_get_bool(buf, "solver.adaptative_al", cfg%adaptative_al, cfg%adaptative_al)
    call json_get_bool(buf, "solver.adaptative_r", cfg%adaptative_al, cfg%adaptative_al)
    call json_get_bool(buf, "solver.adaptative_r_al", cfg%adaptative_al, cfg%adaptative_al)
    call json_get_real(buf, "solver.k_aal", cfg%k_aal, cfg%k_aal)
    call json_get_real(buf, "solver.r_factor", cfg%k_aal, cfg%k_aal)
    call json_get_real(buf, "solver.uzawa_tol", cfg%uzawa_tol, cfg%uzawa_tol)
    call json_get_real(buf, "solver.l0_ref", cfg%l0_ref, cfg%l0_ref)
    call json_get_real(buf, "solver.t0_ref", cfg%t0_ref, cfg%t0_ref)
    call json_get_real(buf, "solver.u0_ref", cfg%u0_ref, cfg%u0_ref)
    call json_get_real(buf, "solver.p0_ref", cfg%p0_ref, cfg%p0_ref)
    ! Appel de sous-programme modulaire
    call json_get_int(buf, "solver.uzawa_iter", cfg%uzawa_iter, 1)
    ! Appel de sous-programme modulaire
    call json_get_bool(buf, "solver.use_precond", cfg%use_precond, .true.)
    ! Appel de sous-programme modulaire
    call json_get_string(buf, "solver.precond_type", str_val, "")
    if (len_trim(str_val) == 0) then
       call json_get_string(buf, "solver.preconditioner", str_val, "JACOBI")
    end if
    if (trim(str_val) == "JACOBI" .or. trim(str_val) == "DIAGONAL") then
       cfg%precond_type = PRECOND_JACOBI
    else if (trim(str_val) == "SSOR" .or. trim(str_val) == "GAUSS_SEIDEL") then
       cfg%precond_type = PRECOND_SSOR
    else if (trim(str_val) == "ILU0" .or. trim(str_val) == "ILU") then
       cfg%precond_type = PRECOND_ILU0
    else if (trim(str_val) == "BLOCK_TRIANGULAR" .or. trim(str_val) == "SCHUR" .or. &
             trim(str_val) == "BT_SBP" .or. trim(str_val) == "BT-SBP" .or. &
             trim(str_val) == "BLOCK" .or. trim(str_val) == "BLOCK_PRECONDITIONER") then
       cfg%precond_type = PRECOND_BLOCK_TRIANGULAR
    else if (trim(str_val) == "CUSTOM") then
       cfg%precond_type = PRECOND_CUSTOM
    else if (trim(str_val) == "NONE") then
       cfg%precond_type = PRECOND_NONE
    end if
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
