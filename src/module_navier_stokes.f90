!===============================================================================
! Module: module_navier_stokes
! Description: Solveur Navier-Stokes 2D incompressible diphasique et thermique
!              selon l'architecture Fortran 2018 de COVOF 2D.
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
    ! Déclaration de variable locale / paramètre de calcul
    integer :: n_owned_u = 0
    ! Déclaration de variable locale / paramètre de calcul
    integer :: n_owned_v = 0
    ! Déclaration de variable locale / paramètre de calcul
    integer :: n_owned_p = 0
    ! Déclaration de variable locale / paramètre de calcul
    integer :: n_owned_total = 0

    ! Diagonales précalculées pour le préconditionneur Jacobi
    real(dp), allocatable :: diag_u(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: diag_v(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: diag_p(:)

    ! Complément de Schur 5-points pour le préconditionneur par blocs triangulaires (M. El Ouafa 2022/2023)
    real(dp), allocatable :: schur_aE(:, :)
    real(dp), allocatable :: schur_aW(:, :)
    real(dp), allocatable :: schur_aN(:, :)
    real(dp), allocatable :: schur_aS(:, :)
    real(dp), allocatable :: schur_aP(:, :)

    ! Stencils 5-points pour blocs vitesse Fu et Fv (M. El Ouafa 2022/2023)
    real(dp), allocatable :: u_aE(:), u_aW(:), u_aN(:), u_aS(:), u_aP(:)
    real(dp), allocatable :: v_aE(:), v_aW(:), v_aN(:), v_aS(:), v_aP(:)
  contains
    ! Affectation et calcul numérique
    procedure :: init => ns_init
    ! Affectation et calcul numérique
    procedure :: pack_vector => ns_pack_vector
    ! Affectation et calcul numérique
    procedure :: unpack_vector => ns_unpack_vector
    ! Affectation et calcul numérique
    procedure :: apply_system_matrix => ns_apply_system_matrix
    ! Affectation et calcul numérique
    procedure :: assemble_rhs => ns_assemble_rhs
    ! Évaluation de la norme du résidu ou de l'erreur numérique
    procedure :: update_pressure_al => ns_update_pressure_al
    ! Affectation et calcul numérique
    procedure :: assemble_preconditioner_blocks => ns_assemble_preconditioner_blocks
    ! Affectation et calcul numérique
    procedure :: solve_energy => ns_solve_energy
    procedure :: advect_ssp_rk3_momentum => ns_advect_ssp_rk3_momentum
    procedure :: solve_total_energy => ns_solve_total_energy
    ! Affectation et calcul numérique
    procedure :: compute_mean_nusselt_hot_wall => ns_compute_mean_nusselt_hot_wall
  end type navier_stokes_t

contains

  !-----------------------------------------------------------------------------
  ! Détection si une face de vitesse normale u(i, j) est dans un obstacle solide
  !-----------------------------------------------------------------------------
  pure function is_u_face_in_obstacle(mesh, cfg, i, j) result(in_obs)
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)   :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)    :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(in)           :: i, j
    ! Déclaration de variable locale / paramètre de calcul
    logical                       :: in_obs
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: x_f, y_c

    ! Affectation et calcul numérique
    in_obs = .false.
    if (.not. cfg%obs_enabled) return
    ! Affectation et calcul numérique
    x_f = x_face(mesh, i)
    ! Affectation et calcul numérique
    y_c = y_center(mesh, j)
    ! Affectation et calcul numérique
    in_obs = is_in_obstacle(x_f, y_c, cfg%obs_type, cfg%obs_x0, cfg%obs_y0, &
                            cfg%obs_x1, cfg%obs_y1, cfg%obs_xc, cfg%obs_yc, cfg%obs_radius)
  end function is_u_face_in_obstacle

  !-----------------------------------------------------------------------------
  ! Détection si une face de vitesse normale v(i, j) est dans un obstacle solide
  !-----------------------------------------------------------------------------
  pure function is_v_face_in_obstacle(mesh, cfg, i, j) result(in_obs)
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)   :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)    :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(in)           :: i, j
    ! Déclaration de variable locale / paramètre de calcul
    logical                       :: in_obs
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: x_c, y_f

    ! Affectation et calcul numérique
    in_obs = .false.
    if (.not. cfg%obs_enabled) return
    ! Affectation et calcul numérique
    x_c = x_center(mesh, i)
    ! Affectation et calcul numérique
    y_f = y_face(mesh, j)
    ! Affectation et calcul numérique
    in_obs = is_in_obstacle(x_c, y_f, cfg%obs_type, cfg%obs_x0, cfg%obs_y0, &
                            cfg%obs_x1, cfg%obs_y1, cfg%obs_xc, cfg%obs_yc, cfg%obs_radius)
  end function is_v_face_in_obstacle

  !-----------------------------------------------------------------------------
  ! Initialisation du solveur Navier-Stokes et comptage des DDLs possédés
  !-----------------------------------------------------------------------------
  subroutine ns_init(this, mesh, cfg)
    ! Déclaration de variable locale / paramètre de calcul
    class(navier_stokes_t), intent(inout) :: this
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)           :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)            :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    integer :: u_nx, v_ny

    ! Affectation et calcul numérique
    u_nx = mesh%nx_local
    ! Affectation et calcul numérique
    if (mesh%is_boundary_right) u_nx = u_nx + 1
    ! Affectation et calcul numérique
    this%n_owned_u = u_nx * mesh%ny_local

    ! Affectation et calcul numérique
    v_ny = mesh%ny_local
    ! Affectation et calcul numérique
    if (mesh%is_boundary_top) v_ny = v_ny + 1
    ! Affectation et calcul numérique
    this%n_owned_v = mesh%nx_local * v_ny

    ! Affectation et calcul numérique
    this%n_owned_p = mesh%nx_local * mesh%ny_local

    ! Évaluation de la condition logique
    if (cfg%solver_type == SOLVER_AUGMENTED_LAGRANGIAN .or. &
        cfg%solver_type == SOLVER_STANDARD_PROJECTION) then
      ! Affectation et calcul numérique
      this%n_owned_total = this%n_owned_u + this%n_owned_v
    ! Branche par défaut si condition non satisfaite
    else
      ! Affectation et calcul numérique
      this%n_owned_total = this%n_owned_u + this%n_owned_v + this%n_owned_p
    end if

    if (allocated(this%diag_u)) deallocate(this%diag_u)
    if (allocated(this%diag_v)) deallocate(this%diag_v)
    if (allocated(this%diag_p)) deallocate(this%diag_p)
    if (allocated(this%schur_aE)) deallocate(this%schur_aE)
    if (allocated(this%schur_aW)) deallocate(this%schur_aW)
    if (allocated(this%schur_aN)) deallocate(this%schur_aN)
    if (allocated(this%schur_aS)) deallocate(this%schur_aS)
    if (allocated(this%schur_aP)) deallocate(this%schur_aP)
    if (allocated(this%u_aE)) deallocate(this%u_aE)
    if (allocated(this%u_aW)) deallocate(this%u_aW)
    if (allocated(this%u_aN)) deallocate(this%u_aN)
    if (allocated(this%u_aS)) deallocate(this%u_aS)
    if (allocated(this%u_aP)) deallocate(this%u_aP)
    if (allocated(this%v_aE)) deallocate(this%v_aE)
    if (allocated(this%v_aW)) deallocate(this%v_aW)
    if (allocated(this%v_aN)) deallocate(this%v_aN)
    if (allocated(this%v_aS)) deallocate(this%v_aS)
    if (allocated(this%v_aP)) deallocate(this%v_aP)

    ! Allocation dynamique de tableau contigu en mémoire
    allocate(this%diag_u(this%n_owned_u))
    allocate(this%diag_v(this%n_owned_v))
    allocate(this%diag_p(this%n_owned_p))

    allocate(this%schur_aE(mesh%nx_local, mesh%ny_local))
    allocate(this%schur_aW(mesh%nx_local, mesh%ny_local))
    allocate(this%schur_aN(mesh%nx_local, mesh%ny_local))
    allocate(this%schur_aS(mesh%nx_local, mesh%ny_local))
    allocate(this%schur_aP(mesh%nx_local, mesh%ny_local))

    allocate(this%u_aE(this%n_owned_u), this%u_aW(this%n_owned_u), &
             this%u_aN(this%n_owned_u), this%u_aS(this%n_owned_u), &
             this%u_aP(this%n_owned_u))
    allocate(this%v_aE(this%n_owned_v), this%v_aW(this%n_owned_v), &
             this%v_aN(this%n_owned_v), this%v_aS(this%n_owned_v), &
             this%v_aP(this%n_owned_v))

    ! Affectation et calcul numérique
    this%diag_u = 1.0_dp
    this%diag_v = 1.0_dp
    this%diag_p = 1.0_dp

    this%schur_aE = 0.0_dp
    this%schur_aW = 0.0_dp
    this%schur_aN = 0.0_dp
    this%schur_aS = 0.0_dp
    this%schur_aP = 1.0_dp

    this%u_aE = 0.0_dp; this%u_aW = 0.0_dp; this%u_aN = 0.0_dp; this%u_aS = 0.0_dp; this%u_aP = 1.0_dp
    this%v_aE = 0.0_dp; this%v_aW = 0.0_dp; this%v_aN = 0.0_dp; this%v_aS = 0.0_dp; this%v_aP = 1.0_dp
  end subroutine ns_init

  !-----------------------------------------------------------------------------
  ! Empaquette les champs physiques 2D dans un vecteur 1D contigu
  !-----------------------------------------------------------------------------
  subroutine ns_pack_vector(this, mesh, cfg, fields, vec)
    ! Déclaration de variable locale / paramètre de calcul
    class(navier_stokes_t), intent(in) :: this
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)        :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)         :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(in)    :: fields
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(out)              :: vec(:)

    ! Déclaration de variable locale / paramètre de calcul
    integer :: i, j, dof, u_nx, v_ny, v_offset, p_offset

    ! Affectation et calcul numérique
    u_nx = mesh%nx_local
    ! Affectation et calcul numérique
    if (mesh%is_boundary_right) u_nx = u_nx + 1

    ! 1. Pack u
    do j = 1, mesh%ny_local
      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 1, u_nx
        ! Affectation et calcul numérique
        dof = (j - 1) * u_nx + i
        ! Affectation et calcul numérique
        vec(dof) = fields%u(i, j)
      end do
    end do

    ! 2. Pack v
    v_offset = this%n_owned_u
    ! Affectation et calcul numérique
    v_ny = mesh%ny_local
    ! Affectation et calcul numérique
    if (mesh%is_boundary_top) v_ny = v_ny + 1

    ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
    do j = 1, v_ny
      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 1, mesh%nx_local
        ! Affectation et calcul numérique
        dof = v_offset + (j - 1) * mesh%nx_local + i
        ! Affectation et calcul numérique
        vec(dof) = fields%v(i, j)
      end do
    end do

    ! 3. Pack p (si monolithique)
    if (cfg%solver_type == SOLVER_MONOLITHIC) then
      ! Affectation et calcul numérique
      p_offset = this%n_owned_u + this%n_owned_v
      ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
      do j = 1, mesh%ny_local
        ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
        do i = 1, mesh%nx_local
          ! Affectation et calcul numérique
          dof = p_offset + (j - 1) * mesh%nx_local + i
          ! Affectation et calcul numérique
          vec(dof) = fields%p(i, j)
        end do
      end do
    end if
  end subroutine ns_pack_vector

  !-----------------------------------------------------------------------------
  ! Dépaquette le vecteur 1D dans les champs 2D et synchronise les halos MPI
  !-----------------------------------------------------------------------------
  subroutine ns_unpack_vector(this, mesh, cfg, vec, fields, homogeneous)
    ! Déclaration de variable locale / paramètre de calcul
    class(navier_stokes_t), intent(in) :: this
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)        :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)         :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)               :: vec(:)
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(inout) :: fields
    ! Déclaration de variable locale / paramètre de calcul
    logical, intent(in)                :: homogeneous

    ! Déclaration de variable locale / paramètre de calcul
    integer :: i, j, dof, u_nx, v_ny, v_offset, p_offset

    ! Affectation et calcul numérique
    u_nx = mesh%nx_local
    ! Affectation et calcul numérique
    if (mesh%is_boundary_right) u_nx = u_nx + 1

    ! 1. Unpack u
    do j = 1, mesh%ny_local
      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 1, u_nx
        ! Affectation et calcul numérique
        dof = (j - 1) * u_nx + i
        ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
        fields%u(i, j) = vec(dof)
      end do
    end do

    ! 2. Unpack v
    v_offset = this%n_owned_u
    ! Affectation et calcul numérique
    v_ny = mesh%ny_local
    ! Affectation et calcul numérique
    if (mesh%is_boundary_top) v_ny = v_ny + 1

    ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
    do j = 1, v_ny
      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 1, mesh%nx_local
        ! Affectation et calcul numérique
        dof = v_offset + (j - 1) * mesh%nx_local + i
        ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
        fields%v(i, j) = vec(dof)
      end do
    end do

    ! 3. Unpack p (si monolithique)
    if (cfg%solver_type == SOLVER_MONOLITHIC) then
      ! Affectation et calcul numérique
      p_offset = this%n_owned_u + this%n_owned_v
      ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
      do j = 1, mesh%ny_local
        ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
        do i = 1, mesh%nx_local
          ! Affectation et calcul numérique
          dof = p_offset + (j - 1) * mesh%nx_local + i
          ! Actualisation de la pression statique au centre de la maille
          fields%p(i, j) = vec(dof)
        end do
      end do
      ! Appel de sous-programme modulaire
      call sync_cell_scalar(mesh, fields%p)
    end if

    ! Synchronisation des mailles fantômes et conditions aux limites
    call sync_u(mesh, fields%u)
    ! Appel de sous-programme modulaire
    call sync_v(mesh, fields%v)
    ! Appel de sous-programme modulaire
    call apply_boundary_conditions(mesh, cfg, fields, homogeneous)
  end subroutine ns_unpack_vector

  !-----------------------------------------------------------------------------
  ! Application de l'opérateur matriciel complet y = A * x
  !-----------------------------------------------------------------------------
  subroutine ns_apply_system_matrix(this, mesh, cfg, fields_old, dt, x, y)
    ! Déclaration de variable locale / paramètre de calcul
    class(navier_stokes_t), intent(in) :: this
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)        :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)         :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(in)    :: fields_old
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)               :: dt
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)               :: x(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(out)              :: y(:)

    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t) :: tf
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: dx, dy, dvol
    ! Déclaration de variable locale / paramètre de calcul
    integer :: i, j, dof, u_nx, v_ny, v_offset, p_offset
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: div_cell(0:mesh%nx_local+1, 0:mesh%ny_local+1)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: rho_u, rho_v, diag_bc, val
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: u_n, u_n_E, u_n_W, v_n_N, v_n_S
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: conv_E, conv_W, conv_N, conv_S
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: mu_cell_E, mu_cell_W, diff_elong_E, diff_elong_W
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: mu_node_N, mu_node_S, shear_N, shear_S
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: grad_div_x, grad_div_y, grad_p
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: v_n, v_n_S_val
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: mu_cell_N, mu_cell_S, diff_elong_N, diff_elong_S
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: mu_node_E, mu_node_W, shear_E, shear_W
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: div_val, diag_p_term
    real(dp) :: chi, u_c, v_c, dp_dx, dp_dy
    ! Déclaration de variable locale / paramètre de calcul
    logical  :: is_wall_bot, is_wall_top, in_obs

    ! Calcul de la métrique spatiale élémentaire
    dx = mesh%dx
    ! Calcul de la métrique spatiale élémentaire
    dy = mesh%dy
    ! Affectation et calcul numérique
    dvol = dx * dy

    ! Affectation et calcul numérique
    u_nx = mesh%nx_local
    ! Affectation et calcul numérique
    if (mesh%is_boundary_right) u_nx = u_nx + 1

    ! Affectation et calcul numérique
    v_ny = mesh%ny_local
    ! Affectation et calcul numérique
    if (mesh%is_boundary_top) v_ny = v_ny + 1

    ! Appel de sous-programme modulaire
    call tf%allocate(mesh)
    ! Appel de sous-programme modulaire
    call this%unpack_vector(mesh, cfg, x, tf, .true.)

    ! Précalcul de la divergence sur chaque cellule locale + halo (dans les bornes allouées 0:nx+1, 0:ny+1)
    div_cell = 0.0_dp
    do j = 0, mesh%ny_local
      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 0, mesh%nx_local
        ! Affectation et calcul numérique
        div_cell(i, j) = (tf%u(i + 1, j) - tf%u(i, j)) / dx + &
                         (tf%v(i, j + 1) - tf%v(i, j)) / dy
      end do
    end do

    ! =========================================================================
    ! 1. Équation de Quantité de Mouvement en u (Faces verticales)
    ! =========================================================================
    do j = 1, mesh%ny_local
      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 1, u_nx
        ! Affectation et calcul numérique
        dof = (j - 1) * u_nx + i
        ! Calcul de la masse volumique locale ou interpolée aux faces
        rho_u = fields_old%rho_u(i, j)
        ! Affectation et calcul numérique
        diag_bc = (rho_u / dt) * dvol

        ! Condition aux limites physique sur les parois gauche / droite
        if (i == 1 .and. mesh%is_boundary_left) then
          ! Évaluation de la condition logique
          if (cfg%bc_left%btype == BC_NEUMANN) then
            ! Affectation et calcul numérique
            y(dof) = diag_bc * (tf%u(1, j) - tf%u(2, j))
          ! Branche par défaut si condition non satisfaite
          else
            ! Affectation et calcul numérique
            y(dof) = diag_bc * tf%u(1, j)
          end if
          cycle
        end if

        ! Évaluation de la condition logique
        if (i == u_nx .and. mesh%is_boundary_right) then
          ! Évaluation de la condition logique
          if (cfg%bc_right%btype == BC_NEUMANN) then
            ! Affectation et calcul numérique
            y(dof) = diag_bc * (tf%u(u_nx, j) - tf%u(u_nx - 1, j))
          ! Branche par défaut si condition non satisfaite
          else
            ! Affectation et calcul numérique
            y(dof) = diag_bc * tf%u(u_nx, j)
          end if
          cycle
        end if

        ! Face dans un obstacle solide (adhérence nulle)
        if (is_u_face_in_obstacle(mesh, cfg, i, j)) then
          ! Affectation et calcul numérique
          y(dof) = diag_bc * tf%u(i, j)
          cycle
        end if

        ! Terme instationnaire
        val = (rho_u / dt) * tf%u(i, j) * dvol

        ! Terme convectif semi-implicite F_uu (forme advective non-conservative a densite harmonique)
        if (cfg%enable_convection) then
          u_n = fields_old%u(i, j)
          u_n_E = fields_old%u(i + 1, j)
          u_n_W = fields_old%u(i - 1, j)

          v_n_N = 0.5_dp * (fields_old%v(i - 1, j + 1) + fields_old%v(i, j + 1))
          v_n_S = 0.5_dp * (fields_old%v(i - 1, j) + fields_old%v(i, j))

          mu_cell_E = (2.0_dp * rho_u * fields_old%rho_u(i + 1, j)) / max(1.0e-12_dp, rho_u + fields_old%rho_u(i + 1, j))
          mu_cell_W = (2.0_dp * rho_u * fields_old%rho_u(i - 1, j)) / max(1.0e-12_dp, rho_u + fields_old%rho_u(i - 1, j))
          if (j < mesh%ny_local) then
            mu_node_N = (2.0_dp * rho_u * fields_old%rho_u(i, j + 1)) / max(1.0e-12_dp, rho_u + fields_old%rho_u(i, j + 1))
          else
            mu_node_N = rho_u
          end if
          if (j > 1) then
            mu_node_S = (2.0_dp * rho_u * fields_old%rho_u(i, j - 1)) / max(1.0e-12_dp, rho_u + fields_old%rho_u(i, j - 1))
          else
            mu_node_S = rho_u
          end if

          conv_E = 0.5_dp * (u_n + u_n_E) * 0.5_dp * (tf%u(i + 1, j) - tf%u(i, j)) * dy
          conv_W = 0.5_dp * (u_n_W + u_n) * 0.5_dp * (tf%u(i, j) - tf%u(i - 1, j)) * dy
          conv_N = v_n_N * 0.5_dp * (tf%u(i, j + 1) - tf%u(i, j)) * dx
          conv_S = v_n_S * 0.5_dp * (tf%u(i, j) - tf%u(i, j - 1)) * dx

          val = val + (mu_cell_E * conv_E + mu_cell_W * conv_W + mu_node_N * conv_N + mu_node_S * conv_S)
        end if

        ! Termes visqueux implicites d'ordre 2
        mu_cell_E = fields_old%mu(i, j)
        ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
        mu_cell_W = fields_old%mu(i - 1, j)
        ! Affectation et calcul numérique
        diff_elong_E = 2.0_dp * mu_cell_E * (tf%u(i + 1, j) - tf%u(i, j)) / dx * dy
        ! Affectation et calcul numérique
        diff_elong_W = 2.0_dp * mu_cell_W * (tf%u(i, j) - tf%u(i - 1, j)) / dx * dy
        ! Affectation et calcul numérique
        val = val - (diff_elong_E - diff_elong_W)

        ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
        mu_node_N = fields_old%mu_node(i, j + 1)
        ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
        mu_node_S = fields_old%mu_node(i, j)

        ! Affectation et calcul numérique
        shear_N = mu_node_N * ((tf%u(i, j + 1) - tf%u(i, j)) / dy + &
                               (tf%v(i, j + 1) - tf%v(i - 1, j + 1)) / dx) * dx
        ! Affectation et calcul numérique
        shear_S = mu_node_S * ((tf%u(i, j) - tf%u(i, j - 1)) / dy + &
                               (tf%v(i, j) - tf%v(i - 1, j)) / dx) * dx
        ! Affectation et calcul numérique
        val = val - (shear_N - shear_S)

        ! Lagrangien augmenté : - r * d(div)/dx * dvol (uniquement en Augmented Lagrangian)
        if (cfg%solver_type == SOLVER_AUGMENTED_LAGRANGIAN) then
          grad_div_x = (div_cell(i, j) - div_cell(i - 1, j)) * dy
          val = val - cfg%r_al * grad_div_x
        end if

        ! Gradient de pression discret B_u^T (si monolithique)
        if (cfg%solver_type == SOLVER_MONOLITHIC) then
          ! Affectation et calcul numérique
          grad_p = (tf%p(i, j) - tf%p(i - 1, j)) * dy
          ! Affectation et calcul numérique
          val = val + grad_p
        end if

        ! Affectation et calcul numérique
        y(dof) = val
      end do
    end do

    ! =========================================================================
    ! 2. Équation de Quantité de Mouvement en v (Faces horizontales)
    ! =========================================================================
    v_offset = this%n_owned_u
    ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
    do j = 1, v_ny
      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 1, mesh%nx_local
        ! Affectation et calcul numérique
        dof = v_offset + (j - 1) * mesh%nx_local + i
        ! Calcul de la masse volumique locale ou interpolée aux faces
        rho_v = fields_old%rho_v(i, j)
        ! Affectation et calcul numérique
        diag_bc = (rho_v / dt) * dvol

        ! Affectation et calcul numérique
        is_wall_bot = (j == 1 .and. mesh%is_boundary_bottom .and. &
                       ! Affectation et calcul numérique
                       (cfg%bc_bottom%btype == BC_DIRICHLET .or. cfg%bc_bottom%btype == BC_FREESLIP))
        ! Affectation et calcul numérique
        is_wall_top = (j == v_ny .and. mesh%is_boundary_top .and. &
                       ! Affectation et calcul numérique
                       (cfg%bc_top%btype == BC_DIRICHLET .or. cfg%bc_top%btype == BC_FREESLIP))
        ! Affectation et calcul numérique
        in_obs = is_v_face_in_obstacle(mesh, cfg, i, j)

        ! Évaluation de la condition logique
        if (is_wall_bot .or. is_wall_top .or. in_obs) then
          ! Affectation et calcul numérique
          y(dof) = diag_bc * tf%v(i, j)
          cycle
        end if

        ! Affectation et calcul numérique
        val = (rho_v / dt) * tf%v(i, j) * dvol

        ! Évaluation de la condition logique
        if (cfg%enable_convection) then
          v_n = fields_old%v(i, j)
          v_n_N = fields_old%v(i, j + 1)
          v_n_S_val = fields_old%v(i, j - 1)

          u_n_E = 0.5_dp * (fields_old%u(i + 1, j - 1) + fields_old%u(i + 1, j))
          u_n_W = 0.5_dp * (fields_old%u(i, j - 1) + fields_old%u(i, j))

          if (i < mesh%nx_local) then
            mu_node_E = (2.0_dp * rho_v * fields_old%rho_v(i + 1, j)) / max(1.0e-12_dp, rho_v + fields_old%rho_v(i + 1, j))
          else
            mu_node_E = rho_v
          end if
          if (i > 1) then
            mu_node_W = (2.0_dp * rho_v * fields_old%rho_v(i - 1, j)) / max(1.0e-12_dp, rho_v + fields_old%rho_v(i - 1, j))
          else
            mu_node_W = rho_v
          end if
          mu_cell_N = (2.0_dp * rho_v * fields_old%rho_v(i, j + 1)) / max(1.0e-12_dp, rho_v + fields_old%rho_v(i, j + 1))
          mu_cell_S = (2.0_dp * rho_v * fields_old%rho_v(i, j - 1)) / max(1.0e-12_dp, rho_v + fields_old%rho_v(i, j - 1))

          conv_E = u_n_E * 0.5_dp * (tf%v(i + 1, j) - tf%v(i, j)) * dy
          conv_W = u_n_W * 0.5_dp * (tf%v(i, j) - tf%v(i - 1, j)) * dy
          conv_N = 0.5_dp * (v_n + v_n_N) * 0.5_dp * (tf%v(i, j + 1) - tf%v(i, j)) * dx
          conv_S = 0.5_dp * (v_n_S_val + v_n) * 0.5_dp * (tf%v(i, j) - tf%v(i, j - 1)) * dx

          val = val + (mu_node_E * conv_E + mu_node_W * conv_W + mu_cell_N * conv_N + mu_cell_S * conv_S)
        end if

        ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
        mu_cell_N = fields_old%mu(i, j)
        ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
        mu_cell_S = fields_old%mu(i, j - 1)
        ! Affectation et calcul numérique
        diff_elong_N = 2.0_dp * mu_cell_N * (tf%v(i, j + 1) - tf%v(i, j)) / dy * dx
        ! Affectation et calcul numérique
        diff_elong_S = 2.0_dp * mu_cell_S * (tf%v(i, j) - tf%v(i, j - 1)) / dy * dx
        ! Affectation et calcul numérique
        val = val - (diff_elong_N - diff_elong_S)

        ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
        mu_node_E = fields_old%mu_node(i + 1, j)
        ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
        mu_node_W = fields_old%mu_node(i, j)

        ! Affectation et calcul numérique
        shear_E = mu_node_E * ((tf%v(i + 1, j) - tf%v(i, j)) / dx + &
                               (tf%u(i + 1, j) - tf%u(i + 1, j - 1)) / dy) * dy
        ! Affectation et calcul numérique
        shear_W = mu_node_W * ((tf%v(i, j) - tf%v(i - 1, j)) / dx + &
                               (tf%u(i, j) - tf%u(i, j - 1)) / dy) * dy
        ! Affectation et calcul numérique
        val = val - (shear_E - shear_W)

        ! Lagrangien augmenté : - r * d(div)/dy * dvol (uniquement en Augmented Lagrangian)
        if (cfg%solver_type == SOLVER_AUGMENTED_LAGRANGIAN) then
          grad_div_y = (div_cell(i, j) - div_cell(i, j - 1)) * dx
          val = val - cfg%r_al * grad_div_y
        end if

        ! Évaluation de la condition logique
        if (cfg%solver_type == SOLVER_MONOLITHIC) then
          ! Affectation et calcul numérique
          grad_p = (tf%p(i, j) - tf%p(i, j - 1)) * dx
          ! Affectation et calcul numérique
          val = val + grad_p
        end if

        ! Affectation et calcul numérique
        y(dof) = val
      end do
    end do

    ! =========================================================================
    ! 3. Équation de Continuité (si monolithique régularisé ou compressible)
    ! =========================================================================
    if (cfg%solver_type == SOLVER_MONOLITHIC) then
      ! Affectation et calcul numérique
      p_offset = this%n_owned_u + this%n_owned_v
      ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
      do j = 1, mesh%ny_local
        ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
        do i = 1, mesh%nx_local
          ! Affectation et calcul numérique
          dof = p_offset + (j - 1) * mesh%nx_local + i
          ! Affectation et calcul numérique
          div_val = div_cell(i, j) * dvol

          if (cfg%compressible) then
             chi = fields_old%chi_T(i, j)
             div_val = div_val + (chi / dt) * tf%p(i, j) * dvol
             u_c = 0.5_dp * (tf%u(i, j) + tf%u(i + 1, j))
             v_c = 0.5_dp * (tf%v(i, j) + tf%v(i, j + 1))
             dp_dx = (fields_old%p(i + 1, j) - fields_old%p(i - 1, j)) / (2.0_dp * dx)
             dp_dy = (fields_old%p(i, j + 1) - fields_old%p(i, j - 1)) / (2.0_dp * dy)
             div_val = div_val + chi * (u_c * dp_dx + v_c * dp_dy) * dvol
          end if

          y(dof) = div_val
        end do
      end do
    end if

    ! Appel de sous-programme modulaire
    call tf%free()
  end subroutine ns_apply_system_matrix

  !-----------------------------------------------------------------------------
  ! Assemblage du second membre rhs pour le système linéaire
  !-----------------------------------------------------------------------------
  subroutine ns_assemble_rhs(this, mesh, cfg, fields, dt, rhs)
    ! Déclaration de variable locale / paramètre de calcul
    class(navier_stokes_t), intent(in) :: this
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)        :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)         :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(in)    :: fields
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)               :: dt
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(out)              :: rhs(:)

    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: dx, dy, dvol
    ! Déclaration de variable locale / paramètre de calcul
    integer :: i, j, dof, u_nx, v_ny, v_offset, p_offset
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: rho_u, rho_v, diag_bc, val, u_bc, v_bc, y_pos, h_step, grad_p
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: mu_N, mu_S, mu_W, mu_E, T_v
    ! Déclaration des variables pour l'injecteur JET_INLET (Cas 3.3 Fluids 2021)
    real(dp) :: xc_jet, half_w_jet

    ! Calcul de la métrique spatiale élémentaire
    dx = mesh%dx
    ! Calcul de la métrique spatiale élémentaire
    dy = mesh%dy
    ! Affectation et calcul numérique
    dvol = dx * dy

    ! Affectation et calcul numérique
    u_nx = mesh%nx_local
    ! Affectation et calcul numérique
    if (mesh%is_boundary_right) u_nx = u_nx + 1

    ! Affectation et calcul numérique
    v_ny = mesh%ny_local
    ! Affectation et calcul numérique
    if (mesh%is_boundary_top) v_ny = v_ny + 1

    ! Affectation et calcul numérique
    rhs = 0.0_dp

    ! =========================================================================
    ! 1. Second membre pour u
    ! =========================================================================
    do j = 1, mesh%ny_local
      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 1, u_nx
        ! Affectation et calcul numérique
        dof = (j - 1) * u_nx + i
        ! Calcul de la masse volumique locale ou interpolée aux faces
        rho_u = fields%rho_u(i, j)
        ! Affectation et calcul numérique
        diag_bc = (rho_u / dt) * dvol

        ! Évaluation de la condition logique
        if (is_u_face_in_obstacle(mesh, cfg, i, j)) then
          ! Affectation et calcul numérique
          rhs(dof) = 0.0_dp
          cycle
        end if

        ! Évaluation de la condition logique
        if (i == 1 .and. mesh%is_boundary_left) then
          ! Affectation et calcul numérique
          u_bc = 0.0_dp
          ! Évaluation de la condition logique
          if (cfg%bc_left%btype == BC_DIRICHLET) then
            ! Évaluation de la condition logique
            if (cfg%obs_enabled .and. cfg%obs_type == OBS_BOX) then
              ! Affectation et calcul numérique
              y_pos = y_center(mesh, j)
              ! Évaluation de la condition logique
              if (y_pos >= cfg%obs_y1) then
                ! Affectation et calcul numérique
                h_step = mesh%ly - cfg%obs_y1
                ! Affectation et calcul numérique
                u_bc = 6.0_dp * cfg%bc_left%u * (y_pos - cfg%obs_y1) * (mesh%ly - y_pos) / (h_step * h_step)
              end if
            ! Branche par défaut si condition non satisfaite
            else
              ! Affectation et calcul numérique
              u_bc = cfg%bc_left%u
            end if
          end if
          ! Affectation et calcul numérique
          rhs(dof) = diag_bc * u_bc
          cycle
        end if

        ! Évaluation de la condition logique
        if (i == u_nx .and. mesh%is_boundary_right) then
          ! Affectation et calcul numérique
          u_bc = 0.0_dp
          ! Évaluation de la condition logique
          if (cfg%bc_right%btype == BC_DIRICHLET) then
            ! Affectation et calcul numérique
            u_bc = cfg%bc_right%u
          end if
          ! Affectation et calcul numérique
          rhs(dof) = diag_bc * u_bc
          cycle
        end if

        ! Affectation et calcul numérique
        if (cfg%enable_ssp_rk3_momentum .and. allocated(fields%rhou_star)) then
          val = (fields%rhou_star(i, j) / dt) * dvol
        else
          val = (rho_u / dt) * fields%u(i, j) * dvol
        end if
        ! Affectation et calcul numérique
        val = val + rho_u * cfg%gx * dvol
        ! Affectation et calcul numérique
        val = val + fields%fx_sigma(i, j) * dvol

        ! Évaluation de la condition logique
        if (cfg%solver_type == SOLVER_AUGMENTED_LAGRANGIAN) then
          ! Affectation et calcul numérique
          grad_p = (fields%p(i, j) - fields%p(i - 1, j)) * dy
          ! Affectation et calcul numérique
          val = val - grad_p
        end if

        ! Évaluation de la condition logique
        if (j == mesh%ny_local .and. mesh%is_boundary_top .and. cfg%bc_top%btype == BC_DIRICHLET) then
          ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
          mu_N = fields%mu_node(i, j + 1)
          ! Affectation et calcul numérique
          val = val + mu_N * (2.0_dp * cfg%bc_top%u / dy) * dx
        end if
        ! Évaluation de la condition logique
        if (j == 1 .and. mesh%is_boundary_bottom .and. cfg%bc_bottom%btype == BC_DIRICHLET) then
          ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
          mu_S = fields%mu_node(i, j)
          ! Affectation et calcul numérique
          val = val + mu_S * (2.0_dp * cfg%bc_bottom%u / dy) * dx
        end if

        ! Affectation et calcul numérique
        rhs(dof) = val
      end do
    end do

    ! =========================================================================
    ! 2. Second membre pour v
    ! =========================================================================
    v_offset = this%n_owned_u
    ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
    do j = 1, v_ny
      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 1, mesh%nx_local
        ! Affectation et calcul numérique
        dof = v_offset + (j - 1) * mesh%nx_local + i
        ! Calcul de la masse volumique locale ou interpolée aux faces
        rho_v = fields%rho_v(i, j)
        ! Affectation et calcul numérique
        diag_bc = (rho_v / dt) * dvol

        ! Évaluation de la condition logique
        if (is_v_face_in_obstacle(mesh, cfg, i, j)) then
          ! Affectation et calcul numérique
          rhs(dof) = 0.0_dp
          cycle
        end if

        ! Affectation et calcul numérique
        if (j == 1 .and. mesh%is_boundary_bottom .and. &
            ! Affectation et calcul numérique
            (cfg%bc_bottom%btype == BC_DIRICHLET .or. cfg%bc_bottom%btype == BC_FREESLIP)) then
          ! Affectation et calcul numérique
          v_bc = 0.0_dp
          ! Affectation et calcul numérique
          if (cfg%bc_bottom%btype == BC_DIRICHLET) v_bc = cfg%bc_bottom%v
          ! Affectation et calcul numérique
          rhs(dof) = diag_bc * v_bc
          cycle
        end if

        ! Affectation et calcul numérique
        if (j == v_ny .and. mesh%is_boundary_top .and. &
            ! Affectation et calcul numérique
            (cfg%bc_top%btype == BC_DIRICHLET .or. cfg%bc_top%btype == BC_FREESLIP)) then
          ! Affectation et calcul numérique
          v_bc = 0.0_dp
          ! Évaluation du profil d'injecteur JET_INLET (Cas 3.3 Fluids 2021)
          if (trim(cfg%bc_top%profile) == "JET_INLET") then
             xc_jet = merge(cfg%bc_top%jet_center, 0.5_dp * mesh%lx, cfg%bc_top%jet_center > 0.0_dp)
             half_w_jet = merge(0.5_dp * cfg%bc_top%jet_width, 0.04_dp, cfg%bc_top%jet_width > 0.0_dp)
             T_v = 0.0_dp
             do p_offset = 1, mesh%nx_global
                if (abs((real(p_offset, dp) - 0.5_dp) * mesh%dx - xc_jet) <= half_w_jet) T_v = T_v + mesh%dx
             end do
             if (T_v < 1.0e-12_dp) T_v = 2.0_dp * half_w_jet
             if (abs(x_center(mesh, i) - xc_jet) <= half_w_jet) then
                v_bc = cfg%bc_top%v * (2.0_dp * half_w_jet) / T_v
             else
                v_bc = -cfg%bc_top%v * (2.0_dp * half_w_jet) / (mesh%lx - T_v)
             end if
          else if (cfg%bc_top%btype == BC_DIRICHLET) then
             v_bc = cfg%bc_top%v
          end if
          ! Affectation et calcul numérique
          rhs(dof) = diag_bc * v_bc
          cycle
        end if

        ! Affectation et calcul numérique
        if (cfg%enable_ssp_rk3_momentum .and. allocated(fields%rhov_star)) then
          val = (fields%rhov_star(i, j) / dt) * dvol
        else
          val = (rho_v / dt) * fields%v(i, j) * dvol
        end if
        ! Affectation et calcul numérique
        val = val + rho_v * cfg%gy * dvol

        ! Évaluation de la condition logique
        if (cfg%enable_energy) then
          ! Affectation et calcul numérique
          T_v = 0.5_dp * (fields%T(i, j) + fields%T(i, j - 1))
          ! Affectation et calcul numérique
          val = val + rho_v * cfg%beta * (T_v - cfg%T_ref) * cfg%gy * dvol
        end if

        ! Affectation et calcul numérique
        val = val + fields%fy_sigma(i, j) * dvol

        ! Évaluation de la condition logique
        if (cfg%solver_type == SOLVER_AUGMENTED_LAGRANGIAN) then
          ! Affectation et calcul numérique
          grad_p = (fields%p(i, j) - fields%p(i, j - 1)) * dx
          ! Affectation et calcul numérique
          val = val - grad_p
        end if

        ! Évaluation de la condition logique
        if (i == 1 .and. mesh%is_boundary_left .and. cfg%bc_left%btype == BC_DIRICHLET) then
          ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
          mu_W = fields%mu_node(i, j)
          ! Affectation et calcul numérique
          val = val + mu_W * (2.0_dp * cfg%bc_left%v / dx) * dy
        end if
        ! Évaluation de la condition logique
        if (i == mesh%nx_local .and. mesh%is_boundary_right .and. cfg%bc_right%btype == BC_DIRICHLET) then
          ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
          mu_E = fields%mu_node(i + 1, j)
          ! Affectation et calcul numérique
          val = val + mu_E * (2.0_dp * cfg%bc_right%v / dx) * dy
        end if

        ! Affectation et calcul numérique
        rhs(dof) = val
      end do
    end do

    ! =========================================================================
    ! 3. Second membre pour p (si monolithique)
    ! =========================================================================
    if (cfg%solver_type == SOLVER_MONOLITHIC) then
      ! Affectation et calcul numérique
      p_offset = this%n_owned_u + this%n_owned_v
      ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
      do j = 1, mesh%ny_local
        ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
        do i = 1, mesh%nx_local
          ! Affectation et calcul numérique
          dof = p_offset + (j - 1) * mesh%nx_local + i
          if (cfg%compressible) then
             rhs(dof) = (fields%chi_T(i, j) / dt) * fields%p(i, j) * dvol
          else
             rhs(dof) = 0.0_dp
          end if
        end do
      end do
    end if
  end subroutine ns_assemble_rhs

  !-----------------------------------------------------------------------------
  ! Mise à jour de la pression par l'étape d'Uzawa du Lagrangien augmenté :
  !   p^{m+1} = p^m - r_al * div(u^{m+1})
  !-----------------------------------------------------------------------------
  subroutine ns_update_pressure_al(this, mesh, cfg, fields)
    ! Déclaration de variable locale / paramètre de calcul
    class(navier_stokes_t), intent(in) :: this
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)        :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)         :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(inout) :: fields

    ! Déclaration de variable locale / paramètre de calcul
    integer :: i, j
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: dx, dy, div_val

    ! Calcul de la métrique spatiale élémentaire
    dx = mesh%dx
    ! Calcul de la métrique spatiale élémentaire
    dy = mesh%dy

    ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
    do j = 1, mesh%ny_local
      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 1, mesh%nx_local
        ! Affectation et calcul numérique
        div_val = (fields%u(i + 1, j) - fields%u(i, j)) / dx + &
                  (fields%v(i, j + 1) - fields%v(i, j)) / dy
        ! Actualisation de la pression statique au centre de la maille
        fields%p(i, j) = fields%p(i, j) - cfg%r_al * div_val
      end do
    end do

    ! Appel de sous-programme modulaire
    call sync_cell_scalar(mesh, fields%p)
  end subroutine ns_update_pressure_al

  !-----------------------------------------------------------------------------
  ! Assemblage des diagonales pour le préconditionnement Jacobi
  !-----------------------------------------------------------------------------
  subroutine ns_assemble_preconditioner_blocks(this, mesh, cfg, fields, dt)
    ! Déclaration de variable locale / paramètre de calcul
    class(navier_stokes_t), intent(inout) :: this
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)           :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)            :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(in)       :: fields
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)                  :: dt

    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: dx, dy, dvol, rho_u, rho_v, mu_c, d_val, al_u, al_v
    real(dp) :: Du, Dv, mu_e, mu_w, mu_n, mu_s, mu_avg
    ! Déclaration de variable locale / paramètre de calcul
    integer :: i, j, dof, u_nx, v_ny, dof_u, dof_v
    ! Déclaration de variable locale / paramètre de calcul
    logical :: is_wall

    ! Calcul de la métrique spatiale élémentaire
    dx = mesh%dx
    ! Calcul de la métrique spatiale élémentaire
    dy = mesh%dy
    ! Affectation et calcul numérique
    dvol = dx * dy

    ! Affectation et calcul numérique
    u_nx = mesh%nx_local
    ! Affectation et calcul numérique
    if (mesh%is_boundary_right) u_nx = u_nx + 1

    ! Affectation et calcul numérique
    v_ny = mesh%ny_local
    ! Affectation et calcul numérique
    if (mesh%is_boundary_top) v_ny = v_ny + 1

    ! 1. Diagonale et Stencil Bloc u
    do j = 1, mesh%ny_local
      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 1, u_nx
        ! Affectation et calcul numérique
        dof = (j - 1) * u_nx + i
        ! Calcul de la masse volumique locale ou interpolée aux faces
        rho_u = fields%rho_u(i, j)
        Du = (rho_u / dt) * dvol
        ! Affectation et calcul numérique
        is_wall = (i == 1 .and. mesh%is_boundary_left) .or. &
                  ! Affectation et calcul numérique
                  (i == u_nx .and. mesh%is_boundary_right) .or. &
                  is_u_face_in_obstacle(mesh, cfg, i, j)

        ! Évaluation de la condition logique
        if (is_wall) then
          this%diag_u(dof) = Du
          if (allocated(this%u_aP)) then
            this%u_aE(dof) = 0.0_dp
            this%u_aW(dof) = 0.0_dp
            this%u_aN(dof) = 0.0_dp
            this%u_aS(dof) = 0.0_dp
            this%u_aP(dof) = merge(Du, 1.0_dp, Du > 1.0e-30_dp)
          end if
          cycle
        end if

        ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
        mu_c = fields%mu(i, j)
        al_u = 0.0_dp
        if (cfg%solver_type == SOLVER_AUGMENTED_LAGRANGIAN) then
          al_u = 2.0_dp * cfg%r_al * (dy / dx)
        end if
        d_val = Du + 4.0_dp * mu_c * (dy / dx) + 2.0_dp * mu_c * (dx / dy) + al_u
        this%diag_u(dof) = merge(d_val, 1.0_dp, d_val > 1.0e-15_dp)

        if (allocated(this%u_aP)) then
          mu_e = fields%mu(i, j)
          mu_w = fields%mu(i - 1, j)
          mu_avg = 0.5_dp * (mu_e + mu_w)
          this%u_aE(dof) = 2.0_dp * mu_e * (dy / dx)
          this%u_aW(dof) = 2.0_dp * mu_w * (dy / dx)
          this%u_aN(dof) = mu_avg * (dx / dy)
          this%u_aS(dof) = mu_avg * (dx / dy)
          this%u_aP(dof) = Du + this%u_aE(dof) + this%u_aW(dof) + this%u_aN(dof) + this%u_aS(dof)
        end if
      end do
    end do

    ! 2. Diagonale et Stencil Bloc v
    do j = 1, v_ny
      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 1, mesh%nx_local
        ! Affectation et calcul numérique
        dof = (j - 1) * mesh%nx_local + i
        ! Calcul de la masse volumique locale ou interpolée aux faces
        rho_v = fields%rho_v(i, j)
        Dv = (rho_v / dt) * dvol
        ! Affectation et calcul numérique
        is_wall = (j == 1 .and. mesh%is_boundary_bottom .and. &
                   ! Affectation et calcul numérique
                   (cfg%bc_bottom%btype == BC_DIRICHLET .or. cfg%bc_bottom%btype == BC_FREESLIP)) .or. &
                  ! Affectation et calcul numérique
                  (j == v_ny .and. mesh%is_boundary_top .and. &
                   ! Affectation et calcul numérique
                   (cfg%bc_top%btype == BC_DIRICHLET .or. cfg%bc_top%btype == BC_FREESLIP)) .or. &
                  is_v_face_in_obstacle(mesh, cfg, i, j)

        ! Évaluation de la condition logique
        if (is_wall) then
          this%diag_v(dof) = Dv
          if (allocated(this%v_aP)) then
            this%v_aE(dof) = 0.0_dp
            this%v_aW(dof) = 0.0_dp
            this%v_aN(dof) = 0.0_dp
            this%v_aS(dof) = 0.0_dp
            this%v_aP(dof) = merge(Dv, 1.0_dp, Dv > 1.0e-30_dp)
          end if
          cycle
        end if

        ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
        mu_c = fields%mu(i, j)
        al_v = 0.0_dp
        if (cfg%solver_type == SOLVER_AUGMENTED_LAGRANGIAN) then
          al_v = 2.0_dp * cfg%r_al * (dx / dy)
        end if
        d_val = Dv + 4.0_dp * mu_c * (dx / dy) + 2.0_dp * mu_c * (dy / dx) + al_v
        this%diag_v(dof) = merge(d_val, 1.0_dp, d_val > 1.0e-15_dp)

        if (allocated(this%v_aP)) then
          mu_n = fields%mu(i, j)
          mu_s = fields%mu(i, j - 1)
          mu_avg = 0.5_dp * (mu_n + mu_s)
          this%v_aN(dof) = 2.0_dp * mu_n * (dx / dy)
          this%v_aS(dof) = 2.0_dp * mu_s * (dx / dy)
          this%v_aE(dof) = mu_avg * (dy / dx)
          this%v_aW(dof) = mu_avg * (dy / dx)
          this%v_aP(dof) = Dv + this%v_aN(dof) + this%v_aS(dof) + this%v_aE(dof) + this%v_aW(dof)
        end if
      end do
    end do

    ! 3. Diagonale Bloc p (si monolithique)
    if (cfg%solver_type == SOLVER_MONOLITHIC) then
      this%diag_p = 1.0_dp
    end if

    ! 4. Assemblage du complément de Schur pour préconditionneur par blocs (M. El Ouafa 2022/2023)
    if (allocated(this%schur_aP)) then
      this%schur_aE = 0.0_dp
      this%schur_aW = 0.0_dp
      this%schur_aN = 0.0_dp
      this%schur_aS = 0.0_dp

      do j = 1, mesh%ny_local
        do i = 1, mesh%nx_local
          ! Face Ouest (u_i, j)
          is_wall = (i == 1 .and. mesh%is_boundary_left .and. &
                     (cfg%bc_left%btype == BC_DIRICHLET .or. cfg%bc_left%btype == BC_FREESLIP)) .or. &
                    is_u_face_in_obstacle(mesh, cfg, i, j)
          if (.not. is_wall) then
            dof_u = (j - 1) * u_nx + i
            Du = this%u_aP(dof_u)
            if (abs(Du) < 1.0e-30_dp) Du = 1.0_dp
            this%schur_aW(i, j) = (dy * dy) / Du
          end if

          ! Face Est (u_i+1, j)
          is_wall = (i == mesh%nx_local .and. mesh%is_boundary_right .and. &
                     (cfg%bc_right%btype == BC_DIRICHLET .or. cfg%bc_right%btype == BC_FREESLIP)) .or. &
                    is_u_face_in_obstacle(mesh, cfg, i + 1, j)
          if (.not. is_wall) then
            dof_u = (j - 1) * u_nx + (i + 1)
            Du = this%u_aP(dof_u)
            if (abs(Du) < 1.0e-30_dp) Du = 1.0_dp
            this%schur_aE(i, j) = (dy * dy) / Du
          end if

          ! Face Sud (v_i, j)
          is_wall = (j == 1 .and. mesh%is_boundary_bottom .and. &
                     (cfg%bc_bottom%btype == BC_DIRICHLET .or. cfg%bc_bottom%btype == BC_FREESLIP)) .or. &
                    is_v_face_in_obstacle(mesh, cfg, i, j)
          if (.not. is_wall) then
            dof_v = (j - 1) * mesh%nx_local + i
            Dv = this%v_aP(dof_v)
            if (abs(Dv) < 1.0e-30_dp) Dv = 1.0_dp
            this%schur_aS(i, j) = (dx * dx) / Dv
          end if

          ! Face Nord (v_i, j+1)
          is_wall = (j == mesh%ny_local .and. mesh%is_boundary_top .and. &
                     (cfg%bc_top%btype == BC_DIRICHLET .or. cfg%bc_top%btype == BC_FREESLIP)) .or. &
                    is_v_face_in_obstacle(mesh, cfg, i, j + 1)
          if (.not. is_wall) then
            dof_v = j * mesh%nx_local + i
            Dv = this%v_aP(dof_v)
            if (abs(Dv) < 1.0e-30_dp) Dv = 1.0_dp
            this%schur_aN(i, j) = (dx * dx) / Dv
          end if

          this%schur_aP(i, j) = this%schur_aW(i, j) + this%schur_aE(i, j) + &
                                this%schur_aS(i, j) + this%schur_aN(i, j) + &
                                1.0e-12_dp * dvol
          if (cfg%compressible) then
             this%schur_aP(i, j) = this%schur_aP(i, j) + (fields%chi_T(i, j) / dt) * dvol
          end if
        end do
      end do
    end if
  end subroutine ns_assemble_preconditioner_blocks

  !-----------------------------------------------------------------------------
  ! Résolution de l'énergie thermique (Boussinesq)
  !-----------------------------------------------------------------------------
  subroutine ns_solve_energy(this, mesh, cfg, fields, dt)
    ! Déclaration de variable locale / paramètre de calcul
    class(navier_stokes_t), intent(in) :: this
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)        :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)         :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(inout) :: fields
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)               :: dt

    ! Déclaration de variable locale / paramètre de calcul
    integer :: nx, ny, n_cells, i, j, c_dof, iter, max_it
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: dx, dy, dvol, rho, Cp, k_th
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: u_e, u_w, v_n, v_s
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: d_e, d_w, d_n, d_s, conv_e, conv_w, conv_n, conv_s
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: a_c, val
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: a_diag(:), a_east(:), a_west(:), a_north(:), a_south(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: rhs(:), T_sol(:), r(:), p(:), v_vec(:), s_vec(:), t_vec(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: r0_hat(:), y_prec(:), z_prec(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: rho_old, alpha, omega, rho_new, r0v, s_norm, tt, ts, r_norm, res_init

    if (.not. cfg%enable_energy) return

    ! Affectation et calcul numérique
    nx = mesh%nx_local
    ! Affectation et calcul numérique
    ny = mesh%ny_local
    ! Affectation et calcul numérique
    n_cells = nx * ny

    ! Calcul de la métrique spatiale élémentaire
    dx = mesh%dx
    ! Calcul de la métrique spatiale élémentaire
    dy = mesh%dy
    ! Affectation et calcul numérique
    dvol = dx * dy
    ! Calcul de la masse volumique locale ou interpolée aux faces
    rho = fields%rho(1, 1)
    ! Affectation et calcul numérique
    Cp = cfg%Cp
    ! Affectation et calcul numérique
    k_th = cfg%k_thermal

    ! Allocation dynamique de tableau contigu en mémoire
    allocate(a_diag(n_cells), a_east(n_cells), a_west(n_cells), a_north(n_cells), a_south(n_cells))
    ! Allocation dynamique de tableau contigu en mémoire
    allocate(rhs(n_cells), T_sol(n_cells), r(n_cells), p(n_cells), v_vec(n_cells), s_vec(n_cells), t_vec(n_cells))
    ! Allocation dynamique de tableau contigu en mémoire
    allocate(r0_hat(n_cells), y_prec(n_cells), z_prec(n_cells))

    ! Affectation et calcul numérique
    a_diag = 0.0_dp
    ! Affectation et calcul numérique
    a_east = 0.0_dp
    ! Affectation et calcul numérique
    a_west = 0.0_dp
    ! Affectation et calcul numérique
    a_north = 0.0_dp
    ! Affectation et calcul numérique
    a_south = 0.0_dp
    ! Affectation et calcul numérique
    rhs = 0.0_dp

    ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
    do j = 1, ny
      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 1, nx
        ! Affectation et calcul numérique
        c_dof = (j - 1) * nx + i
        ! Affectation et calcul numérique
        T_sol(c_dof) = fields%T(i, j)

        ! Affectation et calcul numérique
        a_c = (rho * Cp / dt) * dvol
        ! Affectation et calcul numérique
        rhs(c_dof) = (rho * Cp / dt) * dvol * fields%T(i, j)

        ! Affectation et calcul numérique
        u_e = fields%u(i + 1, j)
        ! Affectation et calcul numérique
        d_e = k_th * dy / dx
        ! Affectation et calcul numérique
        conv_e = rho * Cp * u_e * dy * 0.5_dp

        ! Affectation et calcul numérique
        u_w = fields%u(i, j)
        ! Affectation et calcul numérique
        d_w = k_th * dy / dx
        ! Affectation et calcul numérique
        conv_w = rho * Cp * u_w * dy * 0.5_dp

        ! Affectation et calcul numérique
        v_n = fields%v(i, j + 1)
        ! Affectation et calcul numérique
        d_n = k_th * dx / dy
        ! Affectation et calcul numérique
        conv_n = rho * Cp * v_n * dx * 0.5_dp

        ! Affectation et calcul numérique
        v_s = fields%v(i, j)
        ! Affectation et calcul numérique
        d_s = k_th * dx / dy
        ! Affectation et calcul numérique
        conv_s = rho * Cp * v_s * dx * 0.5_dp

        ! Évaluation de la condition logique
        if (i < nx) then
          ! Affectation et calcul numérique
          a_c = a_c + d_e + conv_e
          ! Affectation et calcul numérique
          a_east(c_dof) = -d_e + conv_e
        ! Branche conditionnelle alternative
        else if (mesh%is_boundary_right) then
          ! Évaluation de la condition logique
          if (cfg%bc_right%thermal_type == THERMAL_DIRICHLET) then
            ! Affectation et calcul numérique
            a_c = a_c + 2.0_dp * d_e
            ! Affectation et calcul numérique
            rhs(c_dof) = rhs(c_dof) + 2.0_dp * d_e * cfg%bc_right%T
          end if
        ! Branche par défaut si condition non satisfaite
        else
          ! Affectation et calcul numérique
          a_c = a_c + d_e + conv_e
          ! Affectation et calcul numérique
          rhs(c_dof) = rhs(c_dof) - (-d_e + conv_e) * fields%T(nx + 1, j)
        end if

        ! Évaluation de la condition logique
        if (i > 1) then
          ! Affectation et calcul numérique
          a_c = a_c + d_w - conv_w
          ! Affectation et calcul numérique
          a_west(c_dof) = -d_w - conv_w
        ! Branche conditionnelle alternative
        else if (mesh%is_boundary_left) then
          ! Évaluation de la condition logique
          if (cfg%bc_left%thermal_type == THERMAL_DIRICHLET) then
            ! Affectation et calcul numérique
            a_c = a_c + 2.0_dp * d_w
            ! Affectation et calcul numérique
            rhs(c_dof) = rhs(c_dof) + 2.0_dp * d_w * cfg%bc_left%T
          end if
        ! Branche par défaut si condition non satisfaite
        else
          ! Affectation et calcul numérique
          a_c = a_c + d_w - conv_w
          ! Affectation et calcul numérique
          rhs(c_dof) = rhs(c_dof) - (-d_w - conv_w) * fields%T(0, j)
        end if

        ! Évaluation de la condition logique
        if (j < ny) then
          ! Affectation et calcul numérique
          a_c = a_c + d_n + conv_n
          ! Affectation et calcul numérique
          a_north(c_dof) = -d_n + conv_n
        ! Branche conditionnelle alternative
        else if (mesh%is_boundary_top) then
          ! Évaluation de la condition logique
          if (cfg%bc_top%thermal_type == THERMAL_DIRICHLET) then
            ! Affectation et calcul numérique
            a_c = a_c + 2.0_dp * d_n
            ! Affectation et calcul numérique
            rhs(c_dof) = rhs(c_dof) + 2.0_dp * d_n * cfg%bc_top%T
          end if
        ! Branche par défaut si condition non satisfaite
        else
          ! Affectation et calcul numérique
          a_c = a_c + d_n + conv_n
          ! Affectation et calcul numérique
          rhs(c_dof) = rhs(c_dof) - (-d_n + conv_n) * fields%T(i, ny + 1)
        end if

        ! Évaluation de la condition logique
        if (j > 1) then
          ! Affectation et calcul numérique
          a_c = a_c + d_s - conv_s
          ! Affectation et calcul numérique
          a_south(c_dof) = -d_s - conv_s
        ! Branche conditionnelle alternative
        else if (mesh%is_boundary_bottom) then
          ! Évaluation de la condition logique
          if (cfg%bc_bottom%thermal_type == THERMAL_DIRICHLET) then
            ! Affectation et calcul numérique
            a_c = a_c + 2.0_dp * d_s
            ! Affectation et calcul numérique
            rhs(c_dof) = rhs(c_dof) + 2.0_dp * d_s * cfg%bc_bottom%T
          end if
        ! Branche par défaut si condition non satisfaite
        else
          ! Affectation et calcul numérique
          a_c = a_c + d_s - conv_s
          ! Affectation et calcul numérique
          rhs(c_dof) = rhs(c_dof) - (-d_s - conv_s) * fields%T(i, 0)
        end if

        ! Affectation et calcul numérique
        a_diag(c_dof) = a_c
      end do
    end do

    ! Résolution itérative BiCGSTAB avec préconditionnement Jacobi
    do j = 1, ny
      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 1, nx
        ! Affectation et calcul numérique
        c_dof = (j - 1) * nx + i
        ! Affectation et calcul numérique
        val = a_diag(c_dof) * T_sol(c_dof)
        ! Affectation et calcul numérique
        if (i < nx) val = val + a_east(c_dof) * T_sol(c_dof + 1)
        ! Affectation et calcul numérique
        if (i > 1)  val = val + a_west(c_dof) * T_sol(c_dof - 1)
        ! Affectation et calcul numérique
        if (j < ny) val = val + a_north(c_dof) * T_sol(c_dof + nx)
        ! Affectation et calcul numérique
        if (j > 1)  val = val + a_south(c_dof) * T_sol(c_dof - nx)
        ! Affectation et calcul numérique
        r(c_dof) = rhs(c_dof) - val
      end do
    end do

    ! Affectation et calcul numérique
    r0_hat = r
    ! Affectation et calcul numérique
    p = r
    ! Calcul de la masse volumique locale ou interpolée aux faces
    rho_old = 1.0_dp
    ! Affectation et calcul numérique
    alpha = 1.0_dp
    ! Affectation et calcul numérique
    omega = 1.0_dp

    ! Évaluation de la norme du résidu ou de l'erreur numérique
    res_init = sqrt(sum(r*r))
    ! Évaluation de la condition logique
    if (res_init > 1.0e-13_dp) then
      ! Affectation et calcul numérique
      max_it = 500
      ! Boucle d'itération numérique
      do iter = 1, max_it
        ! Calcul de la masse volumique locale ou interpolée aux faces
        rho_new = sum(r0_hat * r)
        if (abs(rho_new) < 1.0e-30_dp) exit

        ! Évaluation de la condition logique
        if (iter > 1) then
          ! Affectation et calcul numérique
          p = r + (rho_new / rho_old) * (alpha / omega) * (p - omega * v_vec)
        ! Branche par défaut si condition non satisfaite
        else
          ! Affectation et calcul numérique
          p = r
        end if

        ! Affectation et calcul numérique
        y_prec = p / a_diag

        ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
        do j = 1, ny
          ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
          do i = 1, nx
            ! Affectation et calcul numérique
            c_dof = (j - 1) * nx + i
            ! Affectation et calcul numérique
            val = a_diag(c_dof) * y_prec(c_dof)
            ! Affectation et calcul numérique
            if (i < nx) val = val + a_east(c_dof) * y_prec(c_dof + 1)
            ! Affectation et calcul numérique
            if (i > 1)  val = val + a_west(c_dof) * y_prec(c_dof - 1)
            ! Affectation et calcul numérique
            if (j < ny) val = val + a_north(c_dof) * y_prec(c_dof + nx)
            ! Affectation et calcul numérique
            if (j > 1)  val = val + a_south(c_dof) * y_prec(c_dof - nx)
            ! Affectation et calcul numérique
            v_vec(c_dof) = val
          end do
        end do

        ! Affectation et calcul numérique
        r0v = sum(r0_hat * v_vec)
        if (abs(r0v) < 1.0e-30_dp) exit
        ! Affectation et calcul numérique
        alpha = rho_new / r0v

        ! Affectation et calcul numérique
        s_vec = r - alpha * v_vec
        ! Évaluation de la norme du résidu ou de l'erreur numérique
        s_norm = sqrt(sum(s_vec * s_vec))
        ! Évaluation de la condition logique
        if (s_norm / res_init < 1.0e-8_dp) then
          ! Affectation et calcul numérique
          T_sol = T_sol + alpha * y_prec
          exit
        end if

        ! Affectation et calcul numérique
        z_prec = s_vec / a_diag

        ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
        do j = 1, ny
          ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
          do i = 1, nx
            ! Affectation et calcul numérique
            c_dof = (j - 1) * nx + i
            ! Affectation et calcul numérique
            val = a_diag(c_dof) * z_prec(c_dof)
            ! Affectation et calcul numérique
            if (i < nx) val = val + a_east(c_dof) * z_prec(c_dof + 1)
            ! Affectation et calcul numérique
            if (i > 1)  val = val + a_west(c_dof) * z_prec(c_dof - 1)
            ! Affectation et calcul numérique
            if (j < ny) val = val + a_north(c_dof) * z_prec(c_dof + nx)
            ! Affectation et calcul numérique
            if (j > 1)  val = val + a_south(c_dof) * z_prec(c_dof - nx)
            ! Affectation et calcul numérique
            t_vec(c_dof) = val
          end do
        end do

        ! Affectation et calcul numérique
        tt = sum(t_vec * t_vec)
        ! Affectation et calcul numérique
        ts = sum(t_vec * s_vec)
        ! Évaluation de la condition logique
        if (tt > 1.0e-30_dp) then
          ! Affectation et calcul numérique
          omega = ts / tt
        ! Branche par défaut si condition non satisfaite
        else
          ! Affectation et calcul numérique
          omega = 0.0_dp
        end if

        ! Affectation et calcul numérique
        T_sol = T_sol + alpha * y_prec + omega * z_prec
        ! Affectation et calcul numérique
        r = s_vec - omega * t_vec
        ! Calcul de la masse volumique locale ou interpolée aux faces
        rho_old = rho_new

        ! Évaluation de la norme du résidu ou de l'erreur numérique
        r_norm = sqrt(sum(r * r))
        if (r_norm / res_init < 1.0e-8_dp) exit
      end do
    end if

    ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
    do j = 1, ny
      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 1, nx
        ! Évolution du champ scalaire de température
        fields%T(i, j) = T_sol((j - 1) * nx + i)
      end do
    end do

    ! Appel de sous-programme modulaire
    call sync_cell_scalar(mesh, fields%T)
    ! Appel de sous-programme modulaire
    call apply_thermal_boundary_conditions(mesh, cfg, fields)

    ! Libération de la mémoire allouée
    deallocate(a_diag, a_east, a_west, a_north, a_south)
    ! Libération de la mémoire allouée
    deallocate(rhs, T_sol, r, p, v_vec, s_vec, t_vec, r0_hat, y_prec, z_prec)
  end subroutine ns_solve_energy

  !-----------------------------------------------------------------------------
  ! Calcul du nombre de Nusselt moyen sur la paroi chaude (gauche)
  !-----------------------------------------------------------------------------
  function ns_compute_mean_nusselt_hot_wall(this, mesh, cfg, fields) result(nu_mean)
    ! Déclaration de variable locale / paramètre de calcul
    class(navier_stokes_t), intent(in) :: this
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)        :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)         :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(in)    :: fields
    ! Déclaration de variable locale / paramètre de calcul
    real(dp)                           :: nu_mean

    ! Déclaration de variable locale / paramètre de calcul
    integer :: j, count_loc
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: local_nu_sum, global_nu_sum, global_count
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: T_w, T_1, delta_T, nu_local

    ! Affectation et calcul numérique
    local_nu_sum = 0.0_dp
    ! Affectation et calcul numérique
    count_loc = 0

    ! Évaluation de la condition logique
    if (mesh%is_boundary_left) then
      ! Affectation et calcul numérique
      delta_T = abs(cfg%bc_left%T - cfg%bc_right%T)
      ! Affectation et calcul numérique
      if (delta_T < 1.0e-12_dp) delta_T = 1.0_dp

      ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
      do j = 1, mesh%ny_local
        ! Affectation et calcul numérique
        T_w = cfg%bc_left%T
        ! Affectation et calcul numérique
        T_1 = fields%T(1, j)
        ! Affectation et calcul numérique
        nu_local = 2.0_dp * (T_w - T_1) / mesh%dx * (mesh%lx / delta_T)
        ! Affectation et calcul numérique
        local_nu_sum = local_nu_sum + nu_local
        ! Affectation et calcul numérique
        count_loc = count_loc + 1
      end do
    end if

    ! Affectation et calcul numérique
    global_nu_sum = global_sum(mesh, local_nu_sum)
    ! Affectation et calcul numérique
    global_count = global_sum(mesh, real(count_loc, dp))

    ! Évaluation de la condition logique
    if (global_count > 0.0_dp) then
      ! Affectation et calcul numérique
      nu_mean = global_nu_sum / global_count
    ! Branche par défaut si condition non satisfaite
    else
      ! Affectation et calcul numérique
      nu_mean = 1.0_dp
    end if
  end function ns_compute_mean_nusselt_hot_wall

  !-----------------------------------------------------------------------------
  ! Transport inertiel SSP-RK3 préservant la quantité de mouvement (Fluids 2024, Étape 1)
  !-----------------------------------------------------------------------------
  subroutine ns_advect_ssp_rk3_momentum(this, mesh, cfg, fields, dt)
    class(navier_stokes_t), intent(in) :: this
    type(mesh_2d_t),        intent(in) :: mesh
    type(config_t),         intent(in) :: cfg
    type(flow_fields_t),    intent(inout) :: fields
    real(dp),               intent(in) :: dt

    integer  :: i, j, nx, ny, u_nx, v_ny
    real(dp) :: dx, dy
    real(dp), allocatable :: u_adv1(:, :), v_adv1(:, :), u_adv2(:, :), v_adv2(:, :)
    real(dp), allocatable :: rho_n(:, :), rho_1(:, :), rho_2(:, :), div_rho(:, :)
    real(dp), allocatable :: rhou_n(:, :), rhou_1(:, :), rhou_2(:, :), div_u_k(:, :)
    real(dp), allocatable :: rhov_n(:, :), rhov_1(:, :), rhov_2(:, :), div_v_k(:, :)
    real(dp), allocatable :: rhoe_n(:, :), rhoe_1(:, :), rhoe_2(:, :), div_e(:, :)
    real(dp) :: u_e, u_w, v_n, v_s, phi_e, phi_w, phi_n, phi_s
    real(dp) :: u_c, v_c, v_node_n, v_node_s, u_node_e, u_node_w, val_e, val_w, val_n, val_s

    nx = mesh%nx_local
    ny = mesh%ny_local
    dx = mesh%dx
    dy = mesh%dy
    u_nx = nx + merge(1, 0, mesh%is_boundary_right)
    v_ny = ny + merge(1, 0, mesh%is_boundary_top)

    allocate(u_adv1(0:nx+1, 0:ny+1), v_adv1(0:nx+1, 0:ny+1))
    allocate(u_adv2(0:nx+1, 0:ny+1), v_adv2(0:nx+1, 0:ny+1))

    ! Vitesses advectrices extrapolées
    u_adv1 = 2.0_dp * fields%u - fields%u_old
    v_adv1 = 2.0_dp * fields%v - fields%v_old
    u_adv2 = 1.5_dp * fields%u - 0.5_dp * fields%u_old
    v_adv2 = 1.5_dp * fields%v - 0.5_dp * fields%v_old

    ! 1. Advection SSP-RK3 pour rho
    allocate(rho_n(0:nx+1, 0:ny+1), rho_1(0:nx+1, 0:ny+1), rho_2(0:nx+1, 0:ny+1), div_rho(0:nx+1, 0:ny+1))
    rho_n = fields%rho

    ! Étape 1a
    div_rho = 0.0_dp
    do j = 1, ny
       do i = 1, nx
          u_e = fields%u(i + 1, j)
          u_w = fields%u(i, j)
          v_n = fields%v(i, j + 1)
          v_s = fields%v(i, j)
          phi_e = merge(rho_n(i, j), rho_n(i + 1, j), u_e >= 0.0_dp)
          phi_w = merge(rho_n(i - 1, j), rho_n(i, j), u_w >= 0.0_dp)
          phi_n = merge(rho_n(i, j), rho_n(i, j + 1), v_n >= 0.0_dp)
          phi_s = merge(rho_n(i, j - 1), rho_n(i, j), v_s >= 0.0_dp)
          div_rho(i, j) = (u_e * phi_e - u_w * phi_w) / dx + (v_n * phi_n - v_s * phi_s) / dy
       end do
    end do
    rho_1 = rho_n - dt * div_rho
    call exchange_halos_scalar(rho_1, mesh)

    ! Étape 1b
    div_rho = 0.0_dp
    do j = 1, ny
       do i = 1, nx
          u_e = u_adv1(i + 1, j)
          u_w = u_adv1(i, j)
          v_n = v_adv1(i, j + 1)
          v_s = v_adv1(i, j)
          phi_e = merge(rho_1(i, j), rho_1(i + 1, j), u_e >= 0.0_dp)
          phi_w = merge(rho_1(i - 1, j), rho_1(i, j), u_w >= 0.0_dp)
          phi_n = merge(rho_1(i, j), rho_1(i, j + 1), v_n >= 0.0_dp)
          phi_s = merge(rho_1(i, j - 1), rho_1(i, j), v_s >= 0.0_dp)
          div_rho(i, j) = (u_e * phi_e - u_w * phi_w) / dx + (v_n * phi_n - v_s * phi_s) / dy
       end do
    end do
    rho_2 = 0.75_dp * rho_n + 0.25_dp * rho_1 - 0.25_dp * dt * div_rho
    call exchange_halos_scalar(rho_2, mesh)

    ! Étape 1c
    div_rho = 0.0_dp
    do j = 1, ny
       do i = 1, nx
          u_e = u_adv2(i + 1, j)
          u_w = u_adv2(i, j)
          v_n = v_adv2(i, j + 1)
          v_s = v_adv2(i, j)
          phi_e = merge(rho_2(i, j), rho_2(i + 1, j), u_e >= 0.0_dp)
          phi_w = merge(rho_2(i - 1, j), rho_2(i, j), u_w >= 0.0_dp)
          phi_n = merge(rho_2(i, j), rho_2(i, j + 1), v_n >= 0.0_dp)
          phi_s = merge(rho_2(i, j - 1), rho_2(i, j), v_s >= 0.0_dp)
          div_rho(i, j) = (u_e * phi_e - u_w * phi_w) / dx + (v_n * phi_n - v_s * phi_s) / dy
       end do
    end do
    fields%rho = (1.0_dp / 3.0_dp) * rho_n + (2.0_dp / 3.0_dp) * rho_2 - (2.0_dp / 3.0_dp) * dt * div_rho
    call exchange_halos_scalar(fields%rho, mesh)

    ! Actualiser rho_u et rho_v
    do j = 1, ny
       do i = 1, nx + 1
          fields%rho_u(i, j) = 0.5_dp * (fields%rho(i - 1, j) + fields%rho(i, j))
       end do
    end do
    do j = 1, ny + 1
       do i = 1, nx
          fields%rho_v(i, j) = 0.5_dp * (fields%rho(i, j - 1) + fields%rho(i, j))
       end do
    end do

    ! 2. Advection SSP-RK3 pour rhou
    allocate(rhou_n(0:nx+1, 0:ny+1), rhou_1(0:nx+1, 0:ny+1), rhou_2(0:nx+1, 0:ny+1), div_u_k(0:nx+1, 0:ny+1))
    rhou_n = fields%rho_u * fields%u

    div_u_k = 0.0_dp
    do j = 1, ny
       do i = 1, u_nx
          u_c = 0.5_dp * (fields%u(i, j) + merge(fields%u(i + 1, j), fields%u(i, j), i < nx + 1))
          u_w = 0.5_dp * (fields%u(i, j) + merge(fields%u(i - 1, j), fields%u(i, j), i > 1))
          v_node_n = 0.5_dp * (fields%v(min(nx, i), j + 1) + fields%v(max(1, i - 1), j + 1))
          v_node_s = 0.5_dp * (fields%v(min(nx, i), j) + fields%v(max(1, i - 1), j))
          val_e = merge(rhou_n(i, j), rhou_n(min(nx + 1, i + 1), j), u_c >= 0.0_dp)
          val_w = merge(rhou_n(max(1, i - 1), j), rhou_n(i, j), u_w >= 0.0_dp)
          val_n = merge(rhou_n(i, j), rhou_n(i, min(ny, j + 1)), v_node_n >= 0.0_dp)
          val_s = merge(rhou_n(i, max(1, j - 1)), rhou_n(i, j), v_node_s >= 0.0_dp)
          div_u_k(i, j) = (u_c * val_e - u_w * val_w) / dx + (v_node_n * val_n - v_node_s * val_s) / dy
       end do
    end do
    rhou_1 = rhou_n - dt * div_u_k

    div_u_k = 0.0_dp
    do j = 1, ny
       do i = 1, u_nx
          u_c = 0.5_dp * (u_adv1(i, j) + merge(u_adv1(i + 1, j), u_adv1(i, j), i < nx + 1))
          u_w = 0.5_dp * (u_adv1(i, j) + merge(u_adv1(i - 1, j), u_adv1(i, j), i > 1))
          v_node_n = 0.5_dp * (v_adv1(min(nx, i), j + 1) + v_adv1(max(1, i - 1), j + 1))
          v_node_s = 0.5_dp * (v_adv1(min(nx, i), j) + v_adv1(max(1, i - 1), j))
          val_e = merge(rhou_1(i, j), rhou_1(min(nx + 1, i + 1), j), u_c >= 0.0_dp)
          val_w = merge(rhou_1(max(1, i - 1), j), rhou_1(i, j), u_w >= 0.0_dp)
          val_n = merge(rhou_1(i, j), rhou_1(i, min(ny, j + 1)), v_node_n >= 0.0_dp)
          val_s = merge(rhou_1(i, max(1, j - 1)), rhou_1(i, j), v_node_s >= 0.0_dp)
          div_u_k(i, j) = (u_c * val_e - u_w * val_w) / dx + (v_node_n * val_n - v_node_s * val_s) / dy
       end do
    end do
    rhou_2 = 0.75_dp * rhou_n + 0.25_dp * rhou_1 - 0.25_dp * dt * div_u_k

    div_u_k = 0.0_dp
    do j = 1, ny
       do i = 1, u_nx
          u_c = 0.5_dp * (u_adv2(i, j) + merge(u_adv2(i + 1, j), u_adv2(i, j), i < nx + 1))
          u_w = 0.5_dp * (u_adv2(i, j) + merge(u_adv2(i - 1, j), u_adv2(i, j), i > 1))
          v_node_n = 0.5_dp * (v_adv2(min(nx, i), j + 1) + v_adv2(max(1, i - 1), j + 1))
          v_node_s = 0.5_dp * (v_adv2(min(nx, i), j) + v_adv2(max(1, i - 1), j))
          val_e = merge(rhou_2(i, j), rhou_2(min(nx + 1, i + 1), j), u_c >= 0.0_dp)
          val_w = merge(rhou_2(max(1, i - 1), j), rhou_2(i, j), u_w >= 0.0_dp)
          val_n = merge(rhou_2(i, j), rhou_2(i, min(ny, j + 1)), v_node_n >= 0.0_dp)
          val_s = merge(rhou_2(i, max(1, j - 1)), rhou_2(i, j), v_node_s >= 0.0_dp)
          div_u_k(i, j) = (u_c * val_e - u_w * val_w) / dx + (v_node_n * val_n - v_node_s * val_s) / dy
       end do
    end do
    fields%rhou_star = (1.0_dp / 3.0_dp) * rhou_n + (2.0_dp / 3.0_dp) * rhou_2 - (2.0_dp / 3.0_dp) * dt * div_u_k

    ! 3. Advection SSP-RK3 pour rhov
    allocate(rhov_n(0:nx+1, 0:ny+1), rhov_1(0:nx+1, 0:ny+1), rhov_2(0:nx+1, 0:ny+1), div_v_k(0:nx+1, 0:ny+1))
    rhov_n = fields%rho_v * fields%v

    div_v_k = 0.0_dp
    do j = 1, v_ny
       do i = 1, nx
          v_c = 0.5_dp * (fields%v(i, j) + merge(fields%v(i, j + 1), fields%v(i, j), j < ny + 1))
          v_s = 0.5_dp * (fields%v(i, j) + merge(fields%v(i, j - 1), fields%v(i, j), j > 1))
          u_node_e = 0.5_dp * (fields%u(i + 1, min(ny, j)) + fields%u(i + 1, max(1, j - 1)))
          u_node_w = 0.5_dp * (fields%u(i, min(ny, j)) + fields%u(i, max(1, j - 1)))
          val_e = merge(rhov_n(max(1, i - 1), j), rhov_n(i, j), u_node_e >= 0.0_dp)
          val_w = merge(rhov_n(i, j), rhov_n(min(nx, i + 1), j), u_node_w >= 0.0_dp)
          val_n = merge(rhov_n(i, j), rhov_n(i, min(ny + 1, j + 1)), v_c >= 0.0_dp)
          val_s = merge(rhov_n(i, max(1, j - 1)), rhov_n(i, j), v_s >= 0.0_dp)
          div_v_k(i, j) = (u_node_e * val_w - u_node_w * val_e) / dx + (v_c * val_n - v_s * val_s) / dy
       end do
    end do
    rhov_1 = rhov_n - dt * div_v_k

    div_v_k = 0.0_dp
    do j = 1, v_ny
       do i = 1, nx
          v_c = 0.5_dp * (v_adv1(i, j) + merge(v_adv1(i, j + 1), v_adv1(i, j), j < ny + 1))
          v_s = 0.5_dp * (v_adv1(i, j) + merge(v_adv1(i, j - 1), v_adv1(i, j), j > 1))
          u_node_e = 0.5_dp * (u_adv1(i + 1, min(ny, j)) + u_adv1(i + 1, max(1, j - 1)))
          u_node_w = 0.5_dp * (u_adv1(i, min(ny, j)) + u_adv1(i, max(1, j - 1)))
          val_e = merge(rhov_1(max(1, i - 1), j), rhov_1(i, j), u_node_e >= 0.0_dp)
          val_w = merge(rhov_1(i, j), rhov_1(min(nx, i + 1), j), u_node_w >= 0.0_dp)
          val_n = merge(rhov_1(i, j), rhov_1(i, min(ny + 1, j + 1)), v_c >= 0.0_dp)
          val_s = merge(rhov_1(i, max(1, j - 1)), rhov_1(i, j), v_s >= 0.0_dp)
          div_v_k(i, j) = (u_node_e * val_w - u_node_w * val_e) / dx + (v_c * val_n - v_s * val_s) / dy
       end do
    end do
    rhov_2 = 0.75_dp * rhov_n + 0.25_dp * rhov_1 - 0.25_dp * dt * div_v_k

    div_v_k = 0.0_dp
    do j = 1, v_ny
       do i = 1, nx
          v_c = 0.5_dp * (v_adv2(i, j) + merge(v_adv2(i, j + 1), v_adv2(i, j), j < ny + 1))
          v_s = 0.5_dp * (v_adv2(i, j) + merge(v_adv2(i, j - 1), v_adv2(i, j), j > 1))
          u_node_e = 0.5_dp * (u_adv2(i + 1, min(ny, j)) + u_adv2(i + 1, max(1, j - 1)))
          u_node_w = 0.5_dp * (u_adv2(i, min(ny, j)) + u_adv2(i, max(1, j - 1)))
          val_e = merge(rhov_2(max(1, i - 1), j), rhov_2(i, j), u_node_e >= 0.0_dp)
          val_w = merge(rhov_2(i, j), rhov_2(min(nx, i + 1), j), u_node_w >= 0.0_dp)
          val_n = merge(rhov_2(i, j), rhov_2(i, min(ny + 1, j + 1)), v_c >= 0.0_dp)
          val_s = merge(rhov_2(i, max(1, j - 1)), rhov_2(i, j), v_s >= 0.0_dp)
          div_v_k(i, j) = (u_node_e * val_w - u_node_w * val_e) / dx + (v_c * val_n - v_s * val_s) / dy
       end do
    end do
    fields%rhov_star = (1.0_dp / 3.0_dp) * rhov_n + (2.0_dp / 3.0_dp) * rhov_2 - (2.0_dp / 3.0_dp) * dt * div_v_k

    ! 4. Advection SSP-RK3 pour rhoe si compressible
    if (cfg%compressible) then
       allocate(rhoe_n(0:nx+1, 0:ny+1), rhoe_1(0:nx+1, 0:ny+1), rhoe_2(0:nx+1, 0:ny+1), div_e(0:nx+1, 0:ny+1))
       rhoe_n = fields%rho * fields%energy

       div_e = 0.0_dp
       do j = 1, ny
          do i = 1, nx
             u_e = fields%u(i + 1, j)
             u_w = fields%u(i, j)
             v_n = fields%v(i, j + 1)
             v_s = fields%v(i, j)
             phi_e = merge(rhoe_n(i, j), rhoe_n(i + 1, j), u_e >= 0.0_dp)
             phi_w = merge(rhoe_n(i - 1, j), rhoe_n(i, j), u_w >= 0.0_dp)
             phi_n = merge(rhoe_n(i, j), rhoe_n(i, j + 1), v_n >= 0.0_dp)
             phi_s = merge(rhoe_n(i, j - 1), rhoe_n(i, j), v_s >= 0.0_dp)
             div_e(i, j) = (u_e * phi_e - u_w * phi_w) / dx + (v_n * phi_n - v_s * phi_s) / dy
          end do
       end do
       rhoe_1 = rhoe_n - dt * div_e
       call exchange_halos_scalar(rhoe_1, mesh)

       div_e = 0.0_dp
       do j = 1, ny
          do i = 1, nx
             u_e = u_adv1(i + 1, j)
             u_w = u_adv1(i, j)
             v_n = v_adv1(i, j + 1)
             v_s = v_adv1(i, j)
             phi_e = merge(rhoe_1(i, j), rhoe_1(i + 1, j), u_e >= 0.0_dp)
             phi_w = merge(rhoe_1(i - 1, j), rhoe_1(i, j), u_w >= 0.0_dp)
             phi_n = merge(rhoe_1(i, j), rhoe_1(i, j + 1), v_n >= 0.0_dp)
             phi_s = merge(rhoe_1(i, j - 1), rhoe_1(i, j), v_s >= 0.0_dp)
             div_e(i, j) = (u_e * phi_e - u_w * phi_w) / dx + (v_n * phi_n - v_s * phi_s) / dy
          end do
       end do
       rhoe_2 = 0.75_dp * rhoe_n + 0.25_dp * rhoe_1 - 0.25_dp * dt * div_e
       call exchange_halos_scalar(rhoe_2, mesh)

       div_e = 0.0_dp
       do j = 1, ny
          do i = 1, nx
             u_e = u_adv2(i + 1, j)
             u_w = u_adv2(i, j)
             v_n = v_adv2(i, j + 1)
             v_s = v_adv2(i, j)
             phi_e = merge(rhoe_2(i, j), rhoe_2(i + 1, j), u_e >= 0.0_dp)
             phi_w = merge(rhoe_2(i - 1, j), rhoe_2(i, j), u_w >= 0.0_dp)
             phi_n = merge(rhoe_2(i, j), rhoe_2(i, j + 1), v_n >= 0.0_dp)
             phi_s = merge(rhoe_2(i, j - 1), rhoe_2(i, j), v_s >= 0.0_dp)
             div_e(i, j) = (u_e * phi_e - u_w * phi_w) / dx + (v_n * phi_n - v_s * phi_s) / dy
          end do
       end do
       fields%rhoe_star = (1.0_dp / 3.0_dp) * rhoe_n + (2.0_dp / 3.0_dp) * rhoe_2 - (2.0_dp / 3.0_dp) * dt * div_e
       deallocate(rhoe_n, rhoe_1, rhoe_2, div_e)
    end if

    deallocate(u_adv1, v_adv1, u_adv2, v_adv2)
    deallocate(rho_n, rho_1, rho_2, div_rho)
    deallocate(rhou_n, rhou_1, rhou_2, div_u_k)
    deallocate(rhov_n, rhov_1, rhov_2, div_v_k)
  end subroutine ns_advect_ssp_rk3_momentum

  !-----------------------------------------------------------------------------
  ! Résolution de l'équation de l'énergie totale (Fluids 2024, Étape 7, Éq. 13)
  !-----------------------------------------------------------------------------
  subroutine ns_solve_total_energy(this, mesh, cfg, fields, dt)
    class(navier_stokes_t), intent(in) :: this
    type(mesh_2d_t),        intent(in) :: mesh
    type(config_t),         intent(in) :: cfg
    type(flow_fields_t),    intent(inout) :: fields
    real(dp),               intent(in) :: dt

    integer  :: i, j, nx, ny
    real(dp) :: dx, dy, u_e, u_w, v_n, v_s, p_e, p_w, p_n, p_s, div_vp
    real(dp) :: uc, vc, v_sq, grav_work, rhoe_star_val, rhoe_new, rho_new, e_new, u_int

    if (.not. cfg%compressible) return

    nx = mesh%nx_local
    ny = mesh%ny_local
    dx = mesh%dx
    dy = mesh%dy

    do j = 1, ny
       do i = 1, nx
          u_e = fields%u(i + 1, j)
          u_w = fields%u(i, j)
          v_n = fields%v(i, j + 1)
          v_s = fields%v(i, j)

          p_e = 0.5_dp * (fields%p(i, j) + fields%p(i + 1, j))
          p_w = 0.5_dp * (fields%p(i, j) + fields%p(i - 1, j))
          p_n = 0.5_dp * (fields%p(i, j) + fields%p(i, j + 1))
          p_s = 0.5_dp * (fields%p(i, j) + fields%p(i, j - 1))

          div_vp = (u_e * p_e - u_w * p_w) / dx + (v_n * p_n - v_s * p_s) / dy

          uc = 0.5_dp * (u_e + u_w)
          vc = 0.5_dp * (v_n + v_s)
          v_sq = uc * uc + vc * vc

          grav_work = fields%rho(i, j) * (cfg%gx * uc + cfg%gy * vc)

          rhoe_star_val = fields%rho(i, j) * fields%energy(i, j)
          if (allocated(fields%rhoe_star)) then
             rhoe_star_val = fields%rhoe_star(i, j)
          end if
          rhoe_new = rhoe_star_val + dt * (-div_vp + grav_work)

          rho_new = max(1.0e-6_dp, fields%rho(i, j))
          e_new = rhoe_new / rho_new

          fields%rho_e(i, j) = rhoe_new
          fields%energy(i, j) = e_new

          u_int = max(1.0_dp, e_new - 0.5_dp * v_sq)
          fields%T(i, j) = (1.4_dp - 1.0_dp) / 287.05_dp * u_int
       end do
    end do

    call exchange_halos_scalar(fields%energy, mesh)
    call exchange_halos_scalar(fields%T, mesh)
  end subroutine ns_solve_total_energy

end module module_navier_stokes
