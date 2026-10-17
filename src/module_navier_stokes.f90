!===============================================================================
! Module: module_navier_stokes
! Description: Solveur Navier-Stokes 2D incompressible diphasique et thermique
!              selon l'architecture Fortran 2018 de FUGU 2D.
!              - Discrétisation spatiale MAC décalée (Harlow & Welch 1965)
!              - Couplage Vitesse-Pression par Lagrangien Augmenté (Fortin & Glowinski)
!              - Algorithme d'Uzawa pour la pression
!              - Viscosités harmonique aux nœuds et arithmétique aux mailles
!              - Équation de l'énergie thermique couplée (Boussinesq)
!              - Pénalisation volumique d'obstacles (Brinkman)
!===============================================================================
module module_navier_stokes
  use module_precision, only: dp, ip
  use module_constants
  use module_config, only: config_t
  use module_mesh
  use module_fields
  use module_mpi_comm
  implicit none
  private

  public :: navier_stokes_t
  public :: is_u_face_in_obstacle, is_v_face_in_obstacle

  type :: navier_stokes_t
    integer :: n_owned_u = 0
    integer :: n_owned_v = 0
    integer :: n_owned_p = 0
    integer :: n_owned_total = 0

    ! Diagonales précalculées pour le préconditionneur Jacobi
    real(dp), allocatable :: diag_u(:)
    real(dp), allocatable :: diag_v(:)
    real(dp), allocatable :: diag_p(:)
  contains
    procedure :: init => ns_init
    procedure :: pack_vector => ns_pack_vector
    procedure :: unpack_vector => ns_unpack_vector
    procedure :: apply_system_matrix => ns_apply_system_matrix
    procedure :: assemble_rhs => ns_assemble_rhs
    procedure :: update_pressure_al => ns_update_pressure_al
    procedure :: assemble_preconditioner_blocks => ns_assemble_preconditioner_blocks
    procedure :: solve_energy => ns_solve_energy
    procedure :: compute_mean_nusselt_hot_wall => ns_compute_mean_nusselt_hot_wall
  end type navier_stokes_t

contains

  !-----------------------------------------------------------------------------
  ! Détection si une face de vitesse normale u(i, j) est dans un obstacle solide
  !-----------------------------------------------------------------------------
  pure function is_u_face_in_obstacle(mesh, cfg, i, j) result(in_obs)
    type(mesh_2d_t), intent(in)   :: mesh
    type(config_t), intent(in)    :: cfg
    integer, intent(in)           :: i, j
    logical                       :: in_obs
    real(dp) :: x_f, y_c

    in_obs = .false.
    if (.not. cfg%obs_enabled) return
    x_f = x_face(mesh, i)
    y_c = y_center(mesh, j)
    in_obs = is_in_obstacle(x_f, y_c, cfg%obs_type, cfg%obs_x0, cfg%obs_y0, &
                            cfg%obs_x1, cfg%obs_y1, cfg%obs_xc, cfg%obs_yc, cfg%obs_radius)
  end function is_u_face_in_obstacle

  !-----------------------------------------------------------------------------
  ! Détection si une face de vitesse normale v(i, j) est dans un obstacle solide
  !-----------------------------------------------------------------------------
  pure function is_v_face_in_obstacle(mesh, cfg, i, j) result(in_obs)
    type(mesh_2d_t), intent(in)   :: mesh
    type(config_t), intent(in)    :: cfg
    integer, intent(in)           :: i, j
    logical                       :: in_obs
    real(dp) :: x_c, y_f

    in_obs = .false.
    if (.not. cfg%obs_enabled) return
    x_c = x_center(mesh, i)
    y_f = y_face(mesh, j)
    in_obs = is_in_obstacle(x_c, y_f, cfg%obs_type, cfg%obs_x0, cfg%obs_y0, &
                            cfg%obs_x1, cfg%obs_y1, cfg%obs_xc, cfg%obs_yc, cfg%obs_radius)
  end function is_v_face_in_obstacle

  !-----------------------------------------------------------------------------
  ! Initialisation du solveur Navier-Stokes et comptage des DDLs possédés
  !-----------------------------------------------------------------------------
  subroutine ns_init(this, mesh, cfg)
    class(navier_stokes_t), intent(inout) :: this
    type(mesh_2d_t), intent(in)           :: mesh
    type(config_t), intent(in)            :: cfg
    integer :: u_nx, v_ny

    u_nx = mesh%nx_local
    if (mesh%is_boundary_right) u_nx = u_nx + 1
    this%n_owned_u = u_nx * mesh%ny_local

    v_ny = mesh%ny_local
    if (mesh%is_boundary_top) v_ny = v_ny + 1
    this%n_owned_v = mesh%nx_local * v_ny

    this%n_owned_p = mesh%nx_local * mesh%ny_local

    if (cfg%solver_type == SOLVER_AUGMENTED_LAGRANGIAN) then
      this%n_owned_total = this%n_owned_u + this%n_owned_v
    else
      this%n_owned_total = this%n_owned_u + this%n_owned_v + this%n_owned_p
    end if

    if (allocated(this%diag_u)) deallocate(this%diag_u)
    if (allocated(this%diag_v)) deallocate(this%diag_v)
    if (allocated(this%diag_p)) deallocate(this%diag_p)

    allocate(this%diag_u(this%n_owned_u))
    allocate(this%diag_v(this%n_owned_v))
    allocate(this%diag_p(this%n_owned_p))

    this%diag_u = 1.0_dp
    this%diag_v = 1.0_dp
    this%diag_p = 1.0_dp
  end subroutine ns_init

  !-----------------------------------------------------------------------------
  ! Empaquette les champs physiques 2D dans un vecteur 1D contigu
  !-----------------------------------------------------------------------------
  subroutine ns_pack_vector(this, mesh, cfg, fields, vec)
    class(navier_stokes_t), intent(in) :: this
    type(mesh_2d_t), intent(in)        :: mesh
    type(config_t), intent(in)         :: cfg
    type(flow_fields_t), intent(in)    :: fields
    real(dp), intent(out)              :: vec(:)

    integer :: i, j, dof, u_nx, v_ny, v_offset, p_offset

    u_nx = mesh%nx_local
    if (mesh%is_boundary_right) u_nx = u_nx + 1

    ! 1. Pack u
    do j = 1, mesh%ny_local
      do i = 1, u_nx
        dof = (j - 1) * u_nx + i
        vec(dof) = fields%u(i, j)
      end do
    end do

    ! 2. Pack v
    v_offset = this%n_owned_u
    v_ny = mesh%ny_local
    if (mesh%is_boundary_top) v_ny = v_ny + 1

    do j = 1, v_ny
      do i = 1, mesh%nx_local
        dof = v_offset + (j - 1) * mesh%nx_local + i
        vec(dof) = fields%v(i, j)
      end do
    end do

    ! 3. Pack p (si monolithique)
    if (cfg%solver_type == SOLVER_MONOLITHIC) then
      p_offset = this%n_owned_u + this%n_owned_v
      do j = 1, mesh%ny_local
        do i = 1, mesh%nx_local
          dof = p_offset + (j - 1) * mesh%nx_local + i
          vec(dof) = fields%p(i, j)
        end do
      end do
    end if
  end subroutine ns_pack_vector

  !-----------------------------------------------------------------------------
  ! Dépaquette le vecteur 1D dans les champs 2D et synchronise les halos MPI
  !-----------------------------------------------------------------------------
  subroutine ns_unpack_vector(this, mesh, cfg, vec, fields, homogeneous)
    class(navier_stokes_t), intent(in) :: this
    type(mesh_2d_t), intent(in)        :: mesh
    type(config_t), intent(in)         :: cfg
    real(dp), intent(in)               :: vec(:)
    type(flow_fields_t), intent(inout) :: fields
    logical, intent(in)                :: homogeneous

    integer :: i, j, dof, u_nx, v_ny, v_offset, p_offset

    u_nx = mesh%nx_local
    if (mesh%is_boundary_right) u_nx = u_nx + 1

    ! 1. Unpack u
    do j = 1, mesh%ny_local
      do i = 1, u_nx
        dof = (j - 1) * u_nx + i
        fields%u(i, j) = vec(dof)
      end do
    end do

    ! 2. Unpack v
    v_offset = this%n_owned_u
    v_ny = mesh%ny_local
    if (mesh%is_boundary_top) v_ny = v_ny + 1

    do j = 1, v_ny
      do i = 1, mesh%nx_local
        dof = v_offset + (j - 1) * mesh%nx_local + i
        fields%v(i, j) = vec(dof)
      end do
    end do

    ! 3. Unpack p (si monolithique)
    if (cfg%solver_type == SOLVER_MONOLITHIC) then
      p_offset = this%n_owned_u + this%n_owned_v
      do j = 1, mesh%ny_local
        do i = 1, mesh%nx_local
          dof = p_offset + (j - 1) * mesh%nx_local + i
          fields%p(i, j) = vec(dof)
        end do
      end do
      call sync_cell_scalar(mesh, fields%p)
    end if

    ! Synchronisation des mailles fantômes et conditions aux limites
    call sync_u(mesh, fields%u)
    call sync_v(mesh, fields%v)
    call apply_boundary_conditions(mesh, cfg, fields, homogeneous)
  end subroutine ns_unpack_vector

  !-----------------------------------------------------------------------------
  ! Application de l'opérateur matriciel complet y = A * x
  !-----------------------------------------------------------------------------
  subroutine ns_apply_system_matrix(this, mesh, cfg, fields_old, dt, x, y)
    class(navier_stokes_t), intent(in) :: this
    type(mesh_2d_t), intent(in)        :: mesh
    type(config_t), intent(in)         :: cfg
    type(flow_fields_t), intent(in)    :: fields_old
    real(dp), intent(in)               :: dt
    real(dp), intent(in)               :: x(:)
    real(dp), intent(out)              :: y(:)

    type(flow_fields_t) :: tf
    real(dp) :: dx, dy, dvol
    integer :: i, j, dof, u_nx, v_ny, v_offset, p_offset
    real(dp) :: div_cell(0:mesh%nx_local+1, 0:mesh%ny_local+1)
    real(dp) :: rho_u, rho_v, diag_bc, val
    real(dp) :: u_n, u_n_E, u_n_W, v_n_N, v_n_S
    real(dp) :: conv_E, conv_W, conv_N, conv_S
    real(dp) :: mu_cell_E, mu_cell_W, diff_elong_E, diff_elong_W
    real(dp) :: mu_node_N, mu_node_S, shear_N, shear_S
    real(dp) :: grad_div_x, grad_div_y, grad_p
    real(dp) :: v_n, v_n_S_val
    real(dp) :: mu_cell_N, mu_cell_S, diff_elong_N, diff_elong_S
    real(dp) :: mu_node_E, mu_node_W, shear_E, shear_W
    real(dp) :: div_val, diag_p_term
    logical  :: is_wall_bot, is_wall_top, in_obs

    dx = mesh%dx
    dy = mesh%dy
    dvol = dx * dy

    u_nx = mesh%nx_local
    if (mesh%is_boundary_right) u_nx = u_nx + 1

    v_ny = mesh%ny_local
    if (mesh%is_boundary_top) v_ny = v_ny + 1

    call tf%allocate(mesh)
    call this%unpack_vector(mesh, cfg, x, tf, .true.)

    ! Précalcul de la divergence sur chaque cellule locale + halo
    do j = 0, mesh%ny_local + 1
      do i = 0, mesh%nx_local + 1
        div_cell(i, j) = (tf%u(i + 1, j) - tf%u(i, j)) / dx + &
                         (tf%v(i, j + 1) - tf%v(i, j)) / dy
      end do
    end do

    ! =========================================================================
    ! 1. Équation de Quantité de Mouvement en u (Faces verticales)
    ! =========================================================================
    do j = 1, mesh%ny_local
      do i = 1, u_nx
        dof = (j - 1) * u_nx + i
        rho_u = fields_old%rho_u(i, j)
        diag_bc = (rho_u / dt) * dvol

        ! Condition aux limites physique sur les parois gauche / droite
        if (i == 1 .and. mesh%is_boundary_left) then
          if (cfg%bc_left%btype == BC_NEUMANN) then
            y(dof) = diag_bc * (tf%u(1, j) - tf%u(2, j))
          else
            y(dof) = diag_bc * tf%u(1, j)
          end if
          cycle
        end if

        if (i == u_nx .and. mesh%is_boundary_right) then
          if (cfg%bc_right%btype == BC_NEUMANN) then
            y(dof) = diag_bc * (tf%u(u_nx, j) - tf%u(u_nx - 1, j))
          else
            y(dof) = diag_bc * tf%u(u_nx, j)
          end if
          cycle
        end if

        ! Face dans un obstacle solide (adhérence nulle)
        if (is_u_face_in_obstacle(mesh, cfg, i, j)) then
          y(dof) = diag_bc * tf%u(i, j)
          cycle
        end if

        ! Terme instationnaire
        val = (rho_u / dt) * tf%u(i, j) * dvol

        ! Terme convectif semi-implicite F_uu
        if (cfg%enable_convection) then
          u_n = fields_old%u(i, j)
          u_n_E = fields_old%u(i + 1, j)
          u_n_W = fields_old%u(i - 1, j)

          v_n_N = 0.5_dp * (fields_old%v(i - 1, j + 1) + fields_old%v(i, j + 1))
          v_n_S = 0.5_dp * (fields_old%v(i - 1, j) + fields_old%v(i, j))

          conv_E = 0.5_dp * (u_n + u_n_E) * 0.5_dp * (tf%u(i, j) + tf%u(i + 1, j)) * dy
          conv_W = 0.5_dp * (u_n_W + u_n) * 0.5_dp * (tf%u(i - 1, j) + tf%u(i, j)) * dy
          conv_N = v_n_N * 0.5_dp * (tf%u(i, j) + tf%u(i, j + 1)) * dx
          conv_S = v_n_S * 0.5_dp * (tf%u(i, j - 1) + tf%u(i, j)) * dx

          val = val + rho_u * (conv_E - conv_W + conv_N - conv_S)
        end if

        ! Termes visqueux implicites d'ordre 2
        mu_cell_E = fields_old%mu(i, j)
        mu_cell_W = fields_old%mu(i - 1, j)
        diff_elong_E = 2.0_dp * mu_cell_E * (tf%u(i + 1, j) - tf%u(i, j)) / dx * dy
        diff_elong_W = 2.0_dp * mu_cell_W * (tf%u(i, j) - tf%u(i - 1, j)) / dx * dy
        val = val - (diff_elong_E - diff_elong_W)

        mu_node_N = fields_old%mu_node(i, j + 1)
        mu_node_S = fields_old%mu_node(i, j)

        shear_N = mu_node_N * ((tf%u(i, j + 1) - tf%u(i, j)) / dy + &
                               (tf%v(i, j + 1) - tf%v(i - 1, j + 1)) / dx) * dx
        shear_S = mu_node_S * ((tf%u(i, j) - tf%u(i, j - 1)) / dy + &
                               (tf%v(i, j) - tf%v(i - 1, j)) / dx) * dx
        val = val - (shear_N - shear_S)

        ! Lagrangien augmenté : - r * d(div)/dx * dvol
        grad_div_x = (div_cell(i, j) - div_cell(i - 1, j)) * dy
        val = val - cfg%r_al * grad_div_x

        ! Gradient de pression discret B_u^T (si monolithique)
        if (cfg%solver_type == SOLVER_MONOLITHIC) then
          grad_p = (tf%p(i, j) - tf%p(i - 1, j)) * dy
          val = val + grad_p
        end if

        y(dof) = val
      end do
    end do

    ! =========================================================================
    ! 2. Équation de Quantité de Mouvement en v (Faces horizontales)
    ! =========================================================================
    v_offset = this%n_owned_u
    do j = 1, v_ny
      do i = 1, mesh%nx_local
        dof = v_offset + (j - 1) * mesh%nx_local + i
        rho_v = fields_old%rho_v(i, j)
        diag_bc = (rho_v / dt) * dvol

        is_wall_bot = (j == 1 .and. mesh%is_boundary_bottom .and. &
                       (cfg%bc_bottom%btype == BC_DIRICHLET .or. cfg%bc_bottom%btype == BC_FREESLIP))
        is_wall_top = (j == v_ny .and. mesh%is_boundary_top .and. &
                       (cfg%bc_top%btype == BC_DIRICHLET .or. cfg%bc_top%btype == BC_FREESLIP))
        in_obs = is_v_face_in_obstacle(mesh, cfg, i, j)

        if (is_wall_bot .or. is_wall_top .or. in_obs) then
          y(dof) = diag_bc * tf%v(i, j)
          cycle
        end if

        val = (rho_v / dt) * tf%v(i, j) * dvol

        if (cfg%enable_convection) then
          v_n = fields_old%v(i, j)
          v_n_N = fields_old%v(i, j + 1)
          v_n_S_val = fields_old%v(i, j - 1)

          u_n_E = 0.5_dp * (fields_old%u(i + 1, j - 1) + fields_old%u(i + 1, j))
          u_n_W = 0.5_dp * (fields_old%u(i, j - 1) + fields_old%u(i, j))

          conv_E = u_n_E * 0.5_dp * (tf%v(i, j) + tf%v(i + 1, j)) * dy
          conv_W = u_n_W * 0.5_dp * (tf%v(i - 1, j) + tf%v(i, j)) * dy
          conv_N = 0.5_dp * (v_n + v_n_N) * 0.5_dp * (tf%v(i, j) + tf%v(i, j + 1)) * dx
          conv_S = 0.5_dp * (v_n_S_val + v_n) * 0.5_dp * (tf%v(i, j - 1) + tf%v(i, j)) * dx

          val = val + rho_v * (conv_E - conv_W + conv_N - conv_S)
        end if

        mu_cell_N = fields_old%mu(i, j)
        mu_cell_S = fields_old%mu(i, j - 1)
        diff_elong_N = 2.0_dp * mu_cell_N * (tf%v(i, j + 1) - tf%v(i, j)) / dy * dx
        diff_elong_S = 2.0_dp * mu_cell_S * (tf%v(i, j) - tf%v(i, j - 1)) / dy * dx
        val = val - (diff_elong_N - diff_elong_S)

        mu_node_E = fields_old%mu_node(i + 1, j)
        mu_node_W = fields_old%mu_node(i, j)

        shear_E = mu_node_E * ((tf%v(i + 1, j) - tf%v(i, j)) / dx + &
                               (tf%u(i + 1, j) - tf%u(i + 1, j - 1)) / dy) * dy
        shear_W = mu_node_W * ((tf%v(i, j) - tf%v(i - 1, j)) / dx + &
                               (tf%u(i, j) - tf%u(i, j - 1)) / dy) * dy
        val = val - (shear_E - shear_W)

        grad_div_y = (div_cell(i, j) - div_cell(i, j - 1)) * dx
        val = val - cfg%r_al * grad_div_y

        if (cfg%solver_type == SOLVER_MONOLITHIC) then
          grad_p = (tf%p(i, j) - tf%p(i, j - 1)) * dx
          val = val + grad_p
        end if

        y(dof) = val
      end do
    end do

    ! =========================================================================
    ! 3. Équation de Continuité (si monolithique régularisé)
    ! =========================================================================
    if (cfg%solver_type == SOLVER_MONOLITHIC) then
      p_offset = this%n_owned_u + this%n_owned_v
      do j = 1, mesh%ny_local
        do i = 1, mesh%nx_local
          dof = p_offset + (j - 1) * mesh%nx_local + i
          div_val = div_cell(i, j) * dvol
          diag_p_term = (dvol / cfg%r_al) * tf%p(i, j)
          y(dof) = div_val - diag_p_term
        end do
      end do
    end if

    call tf%free()
  end subroutine ns_apply_system_matrix

  !-----------------------------------------------------------------------------
  ! Assemblage du second membre rhs pour le système linéaire
  !-----------------------------------------------------------------------------
  subroutine ns_assemble_rhs(this, mesh, cfg, fields, dt, rhs)
    class(navier_stokes_t), intent(in) :: this
    type(mesh_2d_t), intent(in)        :: mesh
    type(config_t), intent(in)         :: cfg
    type(flow_fields_t), intent(in)    :: fields
    real(dp), intent(in)               :: dt
    real(dp), intent(out)              :: rhs(:)

    real(dp) :: dx, dy, dvol
    integer :: i, j, dof, u_nx, v_ny, v_offset, p_offset
    real(dp) :: rho_u, rho_v, diag_bc, val, u_bc, v_bc, y_pos, h_step, grad_p
    real(dp) :: mu_N, mu_S, mu_W, mu_E, T_v

    dx = mesh%dx
    dy = mesh%dy
    dvol = dx * dy

    u_nx = mesh%nx_local
    if (mesh%is_boundary_right) u_nx = u_nx + 1

    v_ny = mesh%ny_local
    if (mesh%is_boundary_top) v_ny = v_ny + 1

    rhs = 0.0_dp

    ! =========================================================================
    ! 1. Second membre pour u
    ! =========================================================================
    do j = 1, mesh%ny_local
      do i = 1, u_nx
        dof = (j - 1) * u_nx + i
        rho_u = fields%rho_u(i, j)
        diag_bc = (rho_u / dt) * dvol

        if (is_u_face_in_obstacle(mesh, cfg, i, j)) then
          rhs(dof) = 0.0_dp
          cycle
        end if

        if (i == 1 .and. mesh%is_boundary_left) then
          u_bc = 0.0_dp
          if (cfg%bc_left%btype == BC_DIRICHLET) then
            if (cfg%obs_enabled .and. cfg%obs_type == OBS_BOX) then
              y_pos = y_center(mesh, j)
              if (y_pos >= cfg%obs_y1) then
                h_step = mesh%ly - cfg%obs_y1
                u_bc = 6.0_dp * cfg%bc_left%u * (y_pos - cfg%obs_y1) * (mesh%ly - y_pos) / (h_step * h_step)
              end if
            else
              u_bc = cfg%bc_left%u
            end if
          end if
          rhs(dof) = diag_bc * u_bc
          cycle
        end if

        if (i == u_nx .and. mesh%is_boundary_right) then
          u_bc = 0.0_dp
          if (cfg%bc_right%btype == BC_DIRICHLET) then
            u_bc = cfg%bc_right%u
          end if
          rhs(dof) = diag_bc * u_bc
          cycle
        end if

        val = (rho_u / dt) * fields%u(i, j) * dvol
        val = val + rho_u * cfg%gx * dvol
        val = val + fields%fx_sigma(i, j) * dvol

        if (cfg%solver_type == SOLVER_AUGMENTED_LAGRANGIAN) then
          grad_p = (fields%p(i, j) - fields%p(i - 1, j)) * dy
          val = val - grad_p
        end if

        if (j == mesh%ny_local .and. mesh%is_boundary_top .and. cfg%bc_top%btype == BC_DIRICHLET) then
          mu_N = fields%mu_node(i, j + 1)
          val = val + mu_N * (2.0_dp * cfg%bc_top%u / dy) * dx
        end if
        if (j == 1 .and. mesh%is_boundary_bottom .and. cfg%bc_bottom%btype == BC_DIRICHLET) then
          mu_S = fields%mu_node(i, j)
          val = val + mu_S * (2.0_dp * cfg%bc_bottom%u / dy) * dx
        end if

        rhs(dof) = val
      end do
    end do

    ! =========================================================================
    ! 2. Second membre pour v
    ! =========================================================================
    v_offset = this%n_owned_u
    do j = 1, v_ny
      do i = 1, mesh%nx_local
        dof = v_offset + (j - 1) * mesh%nx_local + i
        rho_v = fields%rho_v(i, j)
        diag_bc = (rho_v / dt) * dvol

        if (is_v_face_in_obstacle(mesh, cfg, i, j)) then
          rhs(dof) = 0.0_dp
          cycle
        end if

        if (j == 1 .and. mesh%is_boundary_bottom .and. &
            (cfg%bc_bottom%btype == BC_DIRICHLET .or. cfg%bc_bottom%btype == BC_FREESLIP)) then
          v_bc = 0.0_dp
          if (cfg%bc_bottom%btype == BC_DIRICHLET) v_bc = cfg%bc_bottom%v
          rhs(dof) = diag_bc * v_bc
          cycle
        end if

        if (j == v_ny .and. mesh%is_boundary_top .and. &
            (cfg%bc_top%btype == BC_DIRICHLET .or. cfg%bc_top%btype == BC_FREESLIP)) then
          v_bc = 0.0_dp
          if (cfg%bc_top%btype == BC_DIRICHLET) v_bc = cfg%bc_top%v
          rhs(dof) = diag_bc * v_bc
          cycle
        end if

        val = (rho_v / dt) * fields%v(i, j) * dvol
        val = val + rho_v * cfg%gy * dvol

        if (cfg%enable_energy) then
          T_v = 0.5_dp * (fields%T(i, j) + fields%T(i, j - 1))
          val = val + rho_v * cfg%beta * (T_v - cfg%T_ref) * cfg%gy * dvol
        end if

        val = val + fields%fy_sigma(i, j) * dvol

        if (cfg%solver_type == SOLVER_AUGMENTED_LAGRANGIAN) then
          grad_p = (fields%p(i, j) - fields%p(i, j - 1)) * dx
          val = val - grad_p
        end if

        if (i == 1 .and. mesh%is_boundary_left .and. cfg%bc_left%btype == BC_DIRICHLET) then
          mu_W = fields%mu_node(i, j)
          val = val + mu_W * (2.0_dp * cfg%bc_left%v / dx) * dy
        end if
        if (i == mesh%nx_local .and. mesh%is_boundary_right .and. cfg%bc_right%btype == BC_DIRICHLET) then
          mu_E = fields%mu_node(i + 1, j)
          val = val + mu_E * (2.0_dp * cfg%bc_right%v / dx) * dy
        end if

        rhs(dof) = val
      end do
    end do

    ! =========================================================================
    ! 3. Second membre pour p (si monolithique)
    ! =========================================================================
    if (cfg%solver_type == SOLVER_MONOLITHIC) then
      p_offset = this%n_owned_u + this%n_owned_v
      do j = 1, mesh%ny_local
        do i = 1, mesh%nx_local
          dof = p_offset + (j - 1) * mesh%nx_local + i
          rhs(dof) = - (dvol / cfg%r_al) * fields%p(i, j)
        end do
      end do
    end if
  end subroutine ns_assemble_rhs

  !-----------------------------------------------------------------------------
  ! Mise à jour de la pression par l'étape d'Uzawa du Lagrangien augmenté :
  !   p^{m+1} = p^m - r_al * div(u^{m+1})
  !-----------------------------------------------------------------------------
  subroutine ns_update_pressure_al(this, mesh, cfg, fields)
    class(navier_stokes_t), intent(in) :: this
    type(mesh_2d_t), intent(in)        :: mesh
    type(config_t), intent(in)         :: cfg
    type(flow_fields_t), intent(inout) :: fields

    integer :: i, j
    real(dp) :: dx, dy, div_val

    dx = mesh%dx
    dy = mesh%dy

    do j = 1, mesh%ny_local
      do i = 1, mesh%nx_local
        div_val = (fields%u(i + 1, j) - fields%u(i, j)) / dx + &
                  (fields%v(i, j + 1) - fields%v(i, j)) / dy
        fields%p(i, j) = fields%p(i, j) - cfg%r_al * div_val
      end do
    end do

    call sync_cell_scalar(mesh, fields%p)
  end subroutine ns_update_pressure_al

  !-----------------------------------------------------------------------------
  ! Assemblage des diagonales pour le préconditionnement Jacobi
  !-----------------------------------------------------------------------------
  subroutine ns_assemble_preconditioner_blocks(this, mesh, cfg, fields, dt)
    class(navier_stokes_t), intent(inout) :: this
    type(mesh_2d_t), intent(in)           :: mesh
    type(config_t), intent(in)            :: cfg
    type(flow_fields_t), intent(in)       :: fields
    real(dp), intent(in)                  :: dt

    real(dp) :: dx, dy, dvol, rho_u, rho_v, mu_c, d_val
    integer :: i, j, dof, u_nx, v_ny
    logical :: is_wall

    dx = mesh%dx
    dy = mesh%dy
    dvol = dx * dy

    u_nx = mesh%nx_local
    if (mesh%is_boundary_right) u_nx = u_nx + 1

    v_ny = mesh%ny_local
    if (mesh%is_boundary_top) v_ny = v_ny + 1

    ! 1. Diagonale Bloc u
    do j = 1, mesh%ny_local
      do i = 1, u_nx
        dof = (j - 1) * u_nx + i
        rho_u = fields%rho_u(i, j)
        is_wall = (i == 1 .and. mesh%is_boundary_left) .or. &
                  (i == u_nx .and. mesh%is_boundary_right) .or. &
                  is_u_face_in_obstacle(mesh, cfg, i, j)

        if (is_wall) then
          this%diag_u(dof) = (rho_u / dt) * dvol
          cycle
        end if

        mu_c = fields%mu(i, j)
        d_val = (rho_u / dt) * dvol + 4.0_dp * mu_c * (dy / dx) + 2.0_dp * mu_c * (dx / dy) + &
                2.0_dp * cfg%r_al * (dy / dx)
        if (d_val > 1.0e-15_dp) then
          this%diag_u(dof) = d_val
        else
          this%diag_u(dof) = 1.0_dp
        end if
      end do
    end do

    ! 2. Diagonale Bloc v
    do j = 1, v_ny
      do i = 1, mesh%nx_local
        dof = (j - 1) * mesh%nx_local + i
        rho_v = fields%rho_v(i, j)
        is_wall = (j == 1 .and. mesh%is_boundary_bottom .and. &
                   (cfg%bc_bottom%btype == BC_DIRICHLET .or. cfg%bc_bottom%btype == BC_FREESLIP)) .or. &
                  (j == v_ny .and. mesh%is_boundary_top .and. &
                   (cfg%bc_top%btype == BC_DIRICHLET .or. cfg%bc_top%btype == BC_FREESLIP)) .or. &
                  is_v_face_in_obstacle(mesh, cfg, i, j)

        if (is_wall) then
          this%diag_v(dof) = (rho_v / dt) * dvol
          cycle
        end if

        mu_c = fields%mu(i, j)
        d_val = (rho_v / dt) * dvol + 4.0_dp * mu_c * (dx / dy) + 2.0_dp * mu_c * (dy / dx) + &
                2.0_dp * cfg%r_al * (dx / dy)
        if (d_val > 1.0e-15_dp) then
          this%diag_v(dof) = d_val
        else
          this%diag_v(dof) = 1.0_dp
        end if
      end do
    end do

    ! 3. Diagonale Bloc p (si monolithique)
    if (cfg%solver_type == SOLVER_MONOLITHIC) then
      this%diag_p = dvol / cfg%r_al
    end if
  end subroutine ns_assemble_preconditioner_blocks

  !-----------------------------------------------------------------------------
  ! Résolution de l'énergie thermique (Boussinesq)
  !-----------------------------------------------------------------------------
  subroutine ns_solve_energy(this, mesh, cfg, fields, dt)
    class(navier_stokes_t), intent(in) :: this
    type(mesh_2d_t), intent(in)        :: mesh
    type(config_t), intent(in)         :: cfg
    type(flow_fields_t), intent(inout) :: fields
    real(dp), intent(in)               :: dt

    integer :: nx, ny, n_cells, i, j, c_dof, iter, max_it
    real(dp) :: dx, dy, dvol, rho, Cp, k_th
    real(dp) :: u_e, u_w, v_n, v_s
    real(dp) :: d_e, d_w, d_n, d_s, conv_e, conv_w, conv_n, conv_s
    real(dp) :: a_c, val
    real(dp), allocatable :: a_diag(:), a_east(:), a_west(:), a_north(:), a_south(:)
    real(dp), allocatable :: rhs(:), T_sol(:), r(:), p(:), v_vec(:), s_vec(:), t_vec(:)
    real(dp), allocatable :: r0_hat(:), y_prec(:), z_prec(:)
    real(dp) :: rho_old, alpha, omega, rho_new, r0v, s_norm, tt, ts, r_norm, res_init

    if (.not. cfg%enable_energy) return

    nx = mesh%nx_local
    ny = mesh%ny_local
    n_cells = nx * ny

    dx = mesh%dx
    dy = mesh%dy
    dvol = dx * dy
    rho = fields%rho(1, 1)
    Cp = cfg%Cp
    k_th = cfg%k_thermal

    allocate(a_diag(n_cells), a_east(n_cells), a_west(n_cells), a_north(n_cells), a_south(n_cells))
    allocate(rhs(n_cells), T_sol(n_cells), r(n_cells), p(n_cells), v_vec(n_cells), s_vec(n_cells), t_vec(n_cells))
    allocate(r0_hat(n_cells), y_prec(n_cells), z_prec(n_cells))

    a_diag = 0.0_dp
    a_east = 0.0_dp
    a_west = 0.0_dp
    a_north = 0.0_dp
    a_south = 0.0_dp
    rhs = 0.0_dp

    do j = 1, ny
      do i = 1, nx
        c_dof = (j - 1) * nx + i
        T_sol(c_dof) = fields%T(i, j)

        a_c = (rho * Cp / dt) * dvol
        rhs(c_dof) = (rho * Cp / dt) * dvol * fields%T(i, j)

        u_e = fields%u(i + 1, j)
        d_e = k_th * dy / dx
        conv_e = rho * Cp * u_e * dy * 0.5_dp

        u_w = fields%u(i, j)
        d_w = k_th * dy / dx
        conv_w = rho * Cp * u_w * dy * 0.5_dp

        v_n = fields%v(i, j + 1)
        d_n = k_th * dx / dy
        conv_n = rho * Cp * v_n * dx * 0.5_dp

        v_s = fields%v(i, j)
        d_s = k_th * dx / dy
        conv_s = rho * Cp * v_s * dx * 0.5_dp

        if (i < nx) then
          a_c = a_c + d_e + conv_e
          a_east(c_dof) = -d_e + conv_e
        else if (mesh%is_boundary_right) then
          if (cfg%bc_right%thermal_type == THERMAL_DIRICHLET) then
            a_c = a_c + 2.0_dp * d_e
            rhs(c_dof) = rhs(c_dof) + 2.0_dp * d_e * cfg%bc_right%T
          end if
        else
          a_c = a_c + d_e + conv_e
          rhs(c_dof) = rhs(c_dof) - (-d_e + conv_e) * fields%T(nx + 1, j)
        end if

        if (i > 1) then
          a_c = a_c + d_w - conv_w
          a_west(c_dof) = -d_w - conv_w
        else if (mesh%is_boundary_left) then
          if (cfg%bc_left%thermal_type == THERMAL_DIRICHLET) then
            a_c = a_c + 2.0_dp * d_w
            rhs(c_dof) = rhs(c_dof) + 2.0_dp * d_w * cfg%bc_left%T
          end if
        else
          a_c = a_c + d_w - conv_w
          rhs(c_dof) = rhs(c_dof) - (-d_w - conv_w) * fields%T(0, j)
        end if

        if (j < ny) then
          a_c = a_c + d_n + conv_n
          a_north(c_dof) = -d_n + conv_n
        else if (mesh%is_boundary_top) then
          if (cfg%bc_top%thermal_type == THERMAL_DIRICHLET) then
            a_c = a_c + 2.0_dp * d_n
            rhs(c_dof) = rhs(c_dof) + 2.0_dp * d_n * cfg%bc_top%T
          end if
        else
          a_c = a_c + d_n + conv_n
          rhs(c_dof) = rhs(c_dof) - (-d_n + conv_n) * fields%T(i, ny + 1)
        end if

        if (j > 1) then
          a_c = a_c + d_s - conv_s
          a_south(c_dof) = -d_s - conv_s
        else if (mesh%is_boundary_bottom) then
          if (cfg%bc_bottom%thermal_type == THERMAL_DIRICHLET) then
            a_c = a_c + 2.0_dp * d_s
            rhs(c_dof) = rhs(c_dof) + 2.0_dp * d_s * cfg%bc_bottom%T
          end if
        else
          a_c = a_c + d_s - conv_s
          rhs(c_dof) = rhs(c_dof) - (-d_s - conv_s) * fields%T(i, 0)
        end if

        a_diag(c_dof) = a_c
      end do
    end do

    ! Résolution itérative BiCGSTAB avec préconditionnement Jacobi
    do j = 1, ny
      do i = 1, nx
        c_dof = (j - 1) * nx + i
        val = a_diag(c_dof) * T_sol(c_dof)
        if (i < nx) val = val + a_east(c_dof) * T_sol(c_dof + 1)
        if (i > 1)  val = val + a_west(c_dof) * T_sol(c_dof - 1)
        if (j < ny) val = val + a_north(c_dof) * T_sol(c_dof + nx)
        if (j > 1)  val = val + a_south(c_dof) * T_sol(c_dof - nx)
        r(c_dof) = rhs(c_dof) - val
      end do
    end do

    r0_hat = r
    p = r
    rho_old = 1.0_dp
    alpha = 1.0_dp
    omega = 1.0_dp

    res_init = sqrt(sum(r*r))
    if (res_init > 1.0e-13_dp) then
      max_it = 500
      do iter = 1, max_it
        rho_new = sum(r0_hat * r)
        if (abs(rho_new) < 1.0e-30_dp) exit

        if (iter > 1) then
          p = r + (rho_new / rho_old) * (alpha / omega) * (p - omega * v_vec)
        else
          p = r
        end if

        y_prec = p / a_diag

        do j = 1, ny
          do i = 1, nx
            c_dof = (j - 1) * nx + i
            val = a_diag(c_dof) * y_prec(c_dof)
            if (i < nx) val = val + a_east(c_dof) * y_prec(c_dof + 1)
            if (i > 1)  val = val + a_west(c_dof) * y_prec(c_dof - 1)
            if (j < ny) val = val + a_north(c_dof) * y_prec(c_dof + nx)
            if (j > 1)  val = val + a_south(c_dof) * y_prec(c_dof - nx)
            v_vec(c_dof) = val
          end do
        end do

        r0v = sum(r0_hat * v_vec)
        if (abs(r0v) < 1.0e-30_dp) exit
        alpha = rho_new / r0v

        s_vec = r - alpha * v_vec
        s_norm = sqrt(sum(s_vec * s_vec))
        if (s_norm / res_init < 1.0e-8_dp) then
          T_sol = T_sol + alpha * y_prec
          exit
        end if

        z_prec = s_vec / a_diag

        do j = 1, ny
          do i = 1, nx
            c_dof = (j - 1) * nx + i
            val = a_diag(c_dof) * z_prec(c_dof)
            if (i < nx) val = val + a_east(c_dof) * z_prec(c_dof + 1)
            if (i > 1)  val = val + a_west(c_dof) * z_prec(c_dof - 1)
            if (j < ny) val = val + a_north(c_dof) * z_prec(c_dof + nx)
            if (j > 1)  val = val + a_south(c_dof) * z_prec(c_dof - nx)
            t_vec(c_dof) = val
          end do
        end do

        tt = sum(t_vec * t_vec)
        ts = sum(t_vec * s_vec)
        if (tt > 1.0e-30_dp) then
          omega = ts / tt
        else
          omega = 0.0_dp
        end if

        T_sol = T_sol + alpha * y_prec + omega * z_prec
        r = s_vec - omega * t_vec
        rho_old = rho_new

        r_norm = sqrt(sum(r * r))
        if (r_norm / res_init < 1.0e-8_dp) exit
      end do
    end if

    do j = 1, ny
      do i = 1, nx
        fields%T(i, j) = T_sol((j - 1) * nx + i)
      end do
    end do

    call sync_cell_scalar(mesh, fields%T)
    call apply_thermal_boundary_conditions(mesh, cfg, fields)

    deallocate(a_diag, a_east, a_west, a_north, a_south)
    deallocate(rhs, T_sol, r, p, v_vec, s_vec, t_vec, r0_hat, y_prec, z_prec)
  end subroutine ns_solve_energy

  !-----------------------------------------------------------------------------
  ! Calcul du nombre de Nusselt moyen sur la paroi chaude (gauche)
  !-----------------------------------------------------------------------------
  function ns_compute_mean_nusselt_hot_wall(this, mesh, cfg, fields) result(nu_mean)
    class(navier_stokes_t), intent(in) :: this
    type(mesh_2d_t), intent(in)        :: mesh
    type(config_t), intent(in)         :: cfg
    type(flow_fields_t), intent(in)    :: fields
    real(dp)                           :: nu_mean

    integer :: j, count_loc
    real(dp) :: local_nu_sum, global_nu_sum, global_count
    real(dp) :: T_w, T_1, delta_T, nu_local

    local_nu_sum = 0.0_dp
    count_loc = 0

    if (mesh%is_boundary_left) then
      delta_T = abs(cfg%bc_left%T - cfg%bc_right%T)
      if (delta_T < 1.0e-12_dp) delta_T = 1.0_dp

      do j = 1, mesh%ny_local
        T_w = cfg%bc_left%T
        T_1 = fields%T(1, j)
        nu_local = 2.0_dp * (T_w - T_1) / mesh%dx * (mesh%lx / delta_T)
        local_nu_sum = local_nu_sum + nu_local
        count_loc = count_loc + 1
      end do
    end if

    global_nu_sum = global_sum(mesh, local_nu_sum)
    global_count = global_sum(mesh, real(count_loc, dp))

    if (global_count > 0.0_dp) then
      nu_mean = global_nu_sum / global_count
    else
      nu_mean = 1.0_dp
    end if
  end function ns_compute_mean_nusselt_hot_wall

end module module_navier_stokes
