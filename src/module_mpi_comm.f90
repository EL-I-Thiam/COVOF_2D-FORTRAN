! [Procédure/Module] module module_mpi_comm
module module_mpi_comm
  use module_precision, only: dp
  use module_constants
  use module_mesh, only: mesh_2d_t, x_center, y_center
  use module_fields, only: flow_fields_t
  use module_config, only: config_t
  implicit none
  private
  include 'mpif.h'

  public :: exchange_halos_scalar
  public :: exchange_halos_velocity
  public :: apply_velocity_bc
  public :: apply_thermal_bc
  public :: apply_vof_bc
  public :: global_max_val, global_min_val, global_sum_val
  public :: sync_u, sync_v, sync_cell_scalar
  public :: apply_boundary_conditions, apply_thermal_boundary_conditions
  public :: global_max, global_min, global_sum, global_norm2, global_dot_product

contains

  ! [Procédure/Module] subroutine exchange_halos_scalar
  subroutine exchange_halos_scalar(field, mesh)
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(inout) :: field(0:mesh%nx_local+1, 0:mesh%ny_local+1)

    ! Déclaration de variable locale / paramètre de calcul
    integer :: nx, ny, ierr
    ! Déclaration de variable locale / paramètre de calcul
    integer :: reqs(8), n_reqs
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: s_left(:), r_left(:), s_right(:), r_right(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: s_bot(:), r_bot(:), s_top(:), r_top(:)

    ! Affectation et calcul numérique
    if (mesh%size <= 1) return

    ! Affectation et calcul numérique
    nx = mesh%nx_local
    ! Affectation et calcul numérique
    ny = mesh%ny_local

    ! Allocation dynamique de tableau contigu en mémoire
    allocate(s_left(ny), r_left(ny), s_right(ny), r_right(ny))
    ! Allocation dynamique de tableau contigu en mémoire
    allocate(s_bot(nx), r_bot(nx), s_top(nx), r_top(nx))

    ! Affectation et calcul numérique
    s_left(:)  = field(1, 1:ny)
    ! Affectation et calcul numérique
    s_right(:) = field(nx, 1:ny)

    ! Affectation et calcul numérique
    n_reqs = 0
    ! Évaluation de la condition logique
    if (mesh%neighbor_west /= MPI_PROC_NULL) then
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Irecv(r_left, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 10, mesh%comm, reqs(n_reqs+1), ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Isend(s_left, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 11, mesh%comm, reqs(n_reqs+2), ierr)
       ! Affectation et calcul numérique
       n_reqs = n_reqs + 2
    end if
    ! Évaluation de la condition logique
    if (mesh%neighbor_east /= MPI_PROC_NULL) then
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Irecv(r_right, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 11, mesh%comm, reqs(n_reqs+1), ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Isend(s_right, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 10, mesh%comm, reqs(n_reqs+2), ierr)
       ! Affectation et calcul numérique
       n_reqs = n_reqs + 2
    end if

    ! Évaluation de la condition logique
    if (n_reqs > 0) then
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Waitall(n_reqs, reqs, MPI_STATUSES_IGNORE, ierr)
    end if

    ! Affectation et calcul numérique
    if (mesh%neighbor_west /= MPI_PROC_NULL) field(0, 1:ny) = r_left(:)
    ! Affectation et calcul numérique
    if (mesh%neighbor_east /= MPI_PROC_NULL) field(nx+1, 1:ny) = r_right(:)

    ! Affectation et calcul numérique
    s_bot(:) = field(1:nx, 1)
    ! Actualisation de la pression statique au centre de la maille
    s_top(:) = field(1:nx, ny)

    ! Affectation et calcul numérique
    n_reqs = 0
    ! Évaluation de la condition logique
    if (mesh%neighbor_south /= MPI_PROC_NULL) then
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Irecv(r_bot, nx, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 20, mesh%comm, reqs(n_reqs+1), ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Isend(s_bot, nx, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 21, mesh%comm, reqs(n_reqs+2), ierr)
       ! Affectation et calcul numérique
       n_reqs = n_reqs + 2
    end if
    ! Évaluation de la condition logique
    if (mesh%neighbor_north /= MPI_PROC_NULL) then
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Irecv(r_top, nx, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 21, mesh%comm, reqs(n_reqs+1), ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Isend(s_top, nx, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 20, mesh%comm, reqs(n_reqs+2), ierr)
       ! Affectation et calcul numérique
       n_reqs = n_reqs + 2
    end if

    ! Évaluation de la condition logique
    if (n_reqs > 0) then
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Waitall(n_reqs, reqs, MPI_STATUSES_IGNORE, ierr)
    end if

    ! Affectation et calcul numérique
    if (mesh%neighbor_south /= MPI_PROC_NULL) field(1:nx, 0) = r_bot(:)
    ! Affectation et calcul numérique
    if (mesh%neighbor_north /= MPI_PROC_NULL) field(1:nx, ny+1) = r_top(:)
  end subroutine exchange_halos_scalar

  ! [Procédure/Module] subroutine exchange_halos_velocity
  subroutine exchange_halos_velocity(u, v, mesh)
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(inout) :: u(0:mesh%nx_local+1, 0:mesh%ny_local+1)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(inout) :: v(0:mesh%nx_local+1, 0:mesh%ny_local+1)

    ! Déclaration de variable locale / paramètre de calcul
    integer :: nx, ny, ierr
    ! Déclaration de variable locale / paramètre de calcul
    integer :: reqs(8), n_reqs
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: su_l(:), ru_l(:), su_r(:), ru_r(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: sv_l(:), rv_l(:), sv_r(:), rv_r(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: su_b(:), ru_b(:), su_t(:), ru_t(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: sv_b(:), rv_b(:), sv_t(:), rv_t(:)

    ! Affectation et calcul numérique
    if (mesh%size <= 1) return

    ! Affectation et calcul numérique
    nx = mesh%nx_local
    ! Affectation et calcul numérique
    ny = mesh%ny_local

    ! Allocation dynamique de tableau contigu en mémoire
    allocate(su_l(ny), ru_l(ny), su_r(ny), ru_r(ny))
    ! Allocation dynamique de tableau contigu en mémoire
    allocate(sv_l(ny), rv_l(ny), sv_r(ny), rv_r(ny))

    ! Affectation et calcul numérique
    su_l(:) = u(1, 1:ny)
    ! Affectation et calcul numérique
    su_r(:) = u(nx, 1:ny)
    ! Affectation et calcul numérique
    sv_l(:) = v(1, 1:ny)
    ! Affectation et calcul numérique
    sv_r(:) = v(nx, 1:ny)

    ! Affectation et calcul numérique
    n_reqs = 0
    ! Évaluation de la condition logique
    if (mesh%neighbor_west /= MPI_PROC_NULL) then
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Irecv(ru_l, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 30, mesh%comm, reqs(n_reqs+1), ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Isend(su_l, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 31, mesh%comm, reqs(n_reqs+2), ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Irecv(rv_l, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 32, mesh%comm, reqs(n_reqs+3), ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Isend(sv_l, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 33, mesh%comm, reqs(n_reqs+4), ierr)
       ! Affectation et calcul numérique
       n_reqs = n_reqs + 4
    end if
    ! Évaluation de la condition logique
    if (mesh%neighbor_east /= MPI_PROC_NULL) then
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Irecv(ru_r, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 31, mesh%comm, reqs(n_reqs+1), ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Isend(su_r, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 30, mesh%comm, reqs(n_reqs+2), ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Irecv(rv_r, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 33, mesh%comm, reqs(n_reqs+3), ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Isend(sv_r, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 32, mesh%comm, reqs(n_reqs+4), ierr)
       ! Affectation et calcul numérique
       n_reqs = n_reqs + 4
    end if

    if (n_reqs > 0) call MPI_Waitall(n_reqs, reqs, MPI_STATUSES_IGNORE, ierr)

    ! Évaluation de la condition logique
    if (mesh%neighbor_west /= MPI_PROC_NULL) then
       ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
       u(0, 1:ny) = ru_l(:)
       ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
       v(0, 1:ny) = rv_l(:)
    end if
    ! Évaluation de la condition logique
    if (mesh%neighbor_east /= MPI_PROC_NULL) then
       ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
       u(nx+1, 1:ny) = ru_r(:)
       ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
       v(nx+1, 1:ny) = rv_r(:)
    end if

    ! Allocation dynamique de tableau contigu en mémoire
    allocate(su_b(nx+2), ru_b(nx+2), su_t(nx+2), ru_t(nx+2))
    ! Allocation dynamique de tableau contigu en mémoire
    allocate(sv_b(nx),   rv_b(nx),   sv_t(nx),   rv_t(nx))

    ! Affectation et calcul numérique
    su_b(:) = u(0:nx+1, 1)
    ! Affectation et calcul numérique
    su_t(:) = u(0:nx+1, ny)
    ! Affectation et calcul numérique
    sv_b(:) = v(1:nx, 1)
    ! Affectation et calcul numérique
    sv_t(:) = v(1:nx, ny)

    ! Affectation et calcul numérique
    n_reqs = 0
    ! Évaluation de la condition logique
    if (mesh%neighbor_south /= MPI_PROC_NULL) then
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Irecv(ru_b, nx+2, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 40, mesh%comm, reqs(n_reqs+1), ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Isend(su_b, nx+2, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 41, mesh%comm, reqs(n_reqs+2), ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Irecv(rv_b, nx,   MPI_DOUBLE_PRECISION, mesh%neighbor_south, 42, mesh%comm, reqs(n_reqs+3), ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Isend(sv_b, nx,   MPI_DOUBLE_PRECISION, mesh%neighbor_south, 43, mesh%comm, reqs(n_reqs+4), ierr)
       ! Affectation et calcul numérique
       n_reqs = n_reqs + 4
    end if
    ! Évaluation de la condition logique
    if (mesh%neighbor_north /= MPI_PROC_NULL) then
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Irecv(ru_t, nx+2, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 41, mesh%comm, reqs(n_reqs+1), ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Isend(su_t, nx+2, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 40, mesh%comm, reqs(n_reqs+2), ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Irecv(rv_t, nx,   MPI_DOUBLE_PRECISION, mesh%neighbor_north, 43, mesh%comm, reqs(n_reqs+3), ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Isend(sv_t, nx,   MPI_DOUBLE_PRECISION, mesh%neighbor_north, 42, mesh%comm, reqs(n_reqs+4), ierr)
       ! Affectation et calcul numérique
       n_reqs = n_reqs + 4
    end if

    if (n_reqs > 0) call MPI_Waitall(n_reqs, reqs, MPI_STATUSES_IGNORE, ierr)

    ! Évaluation de la condition logique
    if (mesh%neighbor_south /= MPI_PROC_NULL) then
       ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
       u(0:nx+1, 0) = ru_b(:)
       ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
       v(1:nx, 0)   = rv_b(:)
    end if
    ! Évaluation de la condition logique
    if (mesh%neighbor_north /= MPI_PROC_NULL) then
       ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
       u(0:nx+1, ny+1) = ru_t(:)
       ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
       v(1:nx, ny+1)   = rv_t(:)
    end if
  end subroutine exchange_halos_velocity

  ! [Procédure/Module] subroutine apply_velocity_bc
  subroutine apply_velocity_bc(fields, mesh, cfg, homogeneous)
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(inout) :: fields
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in) :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    logical, intent(in) :: homogeneous

    ! Déclaration de variable locale / paramètre de calcul
    integer :: nx, ny, i, j
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: u_l, v_l, u_r, v_r, u_b, v_b, u_t, v_t
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: y, h, u_inlet
    ! Déclaration des variables pour injecteur supérieur (Cas 3.3 Fluids 2021)
    real(dp) :: x_c, xc_top, half_w_top, v_out_top

    ! Affectation et calcul numérique
    nx = mesh%nx_local
    ! Affectation et calcul numérique
    ny = mesh%ny_local

    ! Affectation et calcul numérique
    u_l = merge(0.0_dp, cfg%bc_left%u, homogeneous)
    ! Affectation et calcul numérique
    v_l = merge(0.0_dp, cfg%bc_left%v, homogeneous)
    ! Affectation et calcul numérique
    u_r = merge(0.0_dp, cfg%bc_right%u, homogeneous)
    ! Affectation et calcul numérique
    v_r = merge(0.0_dp, cfg%bc_right%v, homogeneous)
    ! Affectation et calcul numérique
    u_b = merge(0.0_dp, cfg%bc_bottom%u, homogeneous)
    ! Affectation et calcul numérique
    v_b = merge(0.0_dp, cfg%bc_bottom%v, homogeneous)
    ! Affectation et calcul numérique
    u_t = merge(0.0_dp, cfg%bc_top%u, homogeneous)
    ! Affectation et calcul numérique
    v_t = merge(0.0_dp, cfg%bc_top%v, homogeneous)

    ! Traitement périodicité si 1 seul processus
    if (cfg%bc_left%btype == BC_PERIODIC .and. mesh%px == 1) then
       ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
       fields%u(0, 1:ny)    = fields%u(nx, 1:ny)
       ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
       fields%u(nx+1, 1:ny) = fields%u(1, 1:ny)
       ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
       fields%v(0, 1:ny)    = fields%v(nx, 1:ny)
       ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
       fields%v(nx+1, 1:ny) = fields%v(1, 1:ny)
    end if
    ! Évaluation de la condition logique
    if (cfg%bc_bottom%btype == BC_PERIODIC .and. mesh%py == 1) then
       ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
       fields%v(1:nx, 0)    = fields%v(1:nx, ny)
       ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
       fields%v(1:nx, ny+1) = fields%v(1:nx, 1)
       ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
       fields%u(0:nx+1, 0)    = fields%u(0:nx+1, ny)
       ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
       fields%u(0:nx+1, ny+1) = fields%u(0:nx+1, 1)
    end if

    ! Bord Ouest physique (x = 0)
    if (mesh%neighbor_west == MPI_PROC_NULL .and. cfg%bc_left%btype /= BC_PERIODIC) then
       ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
       do j = 1, ny
          ! Évaluation de la condition logique
          if (cfg%bc_left%btype == BC_DIRICHLET) then
             ! Affectation et calcul numérique
             u_inlet = u_l
             ! Évaluation de la condition logique
             if (.not. homogeneous .and. cfg%obs_enabled .and. cfg%obs_type == OBS_BOX) then
                ! Profil d'entrée marche descendante : parabole sur [y_step, ly]
                y = y_center(mesh, j)
                ! Évaluation de la condition logique
                if (y < cfg%obs_y1) then
                   ! Affectation et calcul numérique
                   u_inlet = 0.0_dp
                ! Branche par défaut si condition non satisfaite
                else
                   ! Affectation et calcul numérique
                   h = mesh%ly - cfg%obs_y1
                   ! Affectation et calcul numérique
                   u_inlet = 6.0_dp * cfg%bc_left%u * (y - cfg%obs_y1) * (mesh%ly - y) / (h * h)
                end if
             end if
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             if (.not. homogeneous) fields%u(1, j) = u_inlet
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%u(0, j) = u_inlet
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%v(0, j) = 2.0_dp * v_l - fields%v(1, j)
          ! Branche conditionnelle alternative
          else if (cfg%bc_left%btype == BC_FREESLIP) then
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             if (.not. homogeneous) fields%u(1, j) = 0.0_dp
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%u(0, j) = 0.0_dp
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%v(0, j) = fields%v(1, j)
          ! Branche conditionnelle alternative
          else if (cfg%bc_left%btype == BC_NEUMANN) then
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%u(0, j) = fields%u(1, j)
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%v(0, j) = fields%v(1, j)
          end if
       end do
    end if

    ! Bord Est physique (x = Lx)
    if (mesh%neighbor_east == MPI_PROC_NULL .and. cfg%bc_right%btype /= BC_PERIODIC) then
       ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
       do j = 1, ny
          ! Évaluation de la condition logique
          if (cfg%bc_right%btype == BC_DIRICHLET) then
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             if (.not. homogeneous) fields%u(nx+1, j) = u_r
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%u(nx+1, j) = u_r
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%v(nx+1, j) = 2.0_dp * v_r - fields%v(nx, j)
          ! Branche conditionnelle alternative
          else if (cfg%bc_right%btype == BC_FREESLIP) then
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             if (.not. homogeneous) fields%u(nx+1, j) = 0.0_dp
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%u(nx+1, j) = 0.0_dp
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%v(nx+1, j) = fields%v(nx, j)
          ! Branche conditionnelle alternative
          else if (cfg%bc_right%btype == BC_NEUMANN) then
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%u(nx+1, j) = fields%u(nx, j)
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%v(nx+1, j) = fields%v(nx, j)
          end if
       end do
    end if

    ! Bord Sud physique (y = 0)
    if (mesh%neighbor_south == MPI_PROC_NULL .and. cfg%bc_bottom%btype /= BC_PERIODIC) then
       ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
       do i = 1, nx
          ! Évaluation de la condition logique
          if (cfg%bc_bottom%btype == BC_DIRICHLET) then
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             if (.not. homogeneous) fields%v(i, 1) = v_b
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%v(i, 0) = v_b
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%u(i, 0) = 2.0_dp * u_b - fields%u(i, 1)
          ! Branche conditionnelle alternative
          else if (cfg%bc_bottom%btype == BC_FREESLIP) then
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             if (.not. homogeneous) fields%v(i, 1) = 0.0_dp
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%v(i, 0) = 0.0_dp
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%u(i, 0) = fields%u(i, 1)
          ! Branche conditionnelle alternative
          else if (cfg%bc_bottom%btype == BC_NEUMANN) then
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%u(i, 0) = fields%u(i, 1)
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%v(i, 0) = fields%v(i, 1)
          end if
       end do
    end if

    ! Bord Nord physique (y = Ly)
    if (mesh%neighbor_north == MPI_PROC_NULL .and. cfg%bc_top%btype /= BC_PERIODIC) then
       ! Traitement du profil d'injecteur JET_INLET (Cas 3.3 Fluids 2021)
       if (trim(cfg%bc_top%profile) == "JET_INLET") then
          xc_top = merge(cfg%bc_top%jet_center, 0.5_dp * mesh%lx, cfg%bc_top%jet_center > 0.0_dp)
          half_w_top = merge(0.5_dp * cfg%bc_top%jet_width, 0.04_dp, cfg%bc_top%jet_width > 0.0_dp)
          v_out_top = -v_t * (2.0_dp * half_w_top) / (mesh%lx - 2.0_dp * half_w_top)
          do i = 1, nx
             x_c = x_center(mesh, i)
             if (abs(x_c - xc_top) <= half_w_top) then
                if (.not. homogeneous) fields%v(i, ny+1) = v_t
                fields%u(i, ny+1) = 0.0_dp
             else
                if (.not. homogeneous) fields%v(i, ny+1) = v_out_top
                fields%u(i, ny+1) = fields%u(i, ny)
             end if
          end do
       else
          ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
          do i = 1, nx
             ! Évaluation de la condition logique
             if (cfg%bc_top%btype == BC_DIRICHLET) then
                ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
                if (.not. homogeneous) fields%v(i, ny+1) = v_t
                fields%u(i, ny+1) = 2.0_dp * u_t - fields%u(i, ny)
             ! Branche conditionnelle alternative
             else if (cfg%bc_top%btype == BC_FREESLIP) then
                ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
                if (.not. homogeneous) fields%v(i, ny+1) = 0.0_dp
                fields%u(i, ny+1) = fields%u(i, ny)
             ! Branche conditionnelle alternative
             else if (cfg%bc_top%btype == BC_NEUMANN) then
                ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
                fields%u(i, ny+1) = fields%u(i, ny)
                fields%v(i, ny+1) = fields%v(i, ny)
             end if
          end do
       end if
    end if
  end subroutine apply_velocity_bc

  ! [Procédure/Module] subroutine apply_thermal_bc
  subroutine apply_thermal_bc(fields, mesh, cfg)
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(inout) :: fields
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in) :: cfg

    ! Déclaration de variable locale / paramètre de calcul
    integer :: nx, ny, i, j
    ! Affectation et calcul numérique
    nx = mesh%nx_local
    ! Affectation et calcul numérique
    ny = mesh%ny_local

    ! Périodicité thermique
    if (cfg%bc_left%thermal_type == BC_PERIODIC .and. mesh%px == 1) then
       ! Évolution du champ scalaire de température
       fields%T(0, 1:ny)    = fields%T(nx, 1:ny)
       ! Évolution du champ scalaire de température
       fields%T(nx+1, 1:ny) = fields%T(1, 1:ny)
    end if
    ! Évaluation de la condition logique
    if (cfg%bc_bottom%thermal_type == BC_PERIODIC .and. mesh%py == 1) then
       ! Évolution du champ scalaire de température
       fields%T(1:nx, 0)    = fields%T(1:nx, ny)
       ! Évolution du champ scalaire de température
       fields%T(1:nx, ny+1) = fields%T(1:nx, 1)
    end if

    ! Bord Ouest physique (x = 0)
    if (mesh%neighbor_west == MPI_PROC_NULL) then
       ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
       do j = 1, ny
          ! Évaluation de la condition logique
          if (cfg%bc_left%thermal_type == THERMAL_DIRICHLET) then
             ! Évolution du champ scalaire de température
             fields%T(0, j) = 2.0_dp * cfg%bc_left%T - fields%T(1, j)
          ! Branche par défaut si condition non satisfaite
          else
             ! Évolution du champ scalaire de température
             fields%T(0, j) = fields%T(1, j)
          end if
       end do
    end if

    ! Bord Est physique (x = Lx)
    if (mesh%neighbor_east == MPI_PROC_NULL) then
       ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
       do j = 1, ny
          ! Évaluation de la condition logique
          if (cfg%bc_right%thermal_type == THERMAL_DIRICHLET) then
             ! Évolution du champ scalaire de température
             fields%T(nx+1, j) = 2.0_dp * cfg%bc_right%T - fields%T(nx, j)
          ! Branche par défaut si condition non satisfaite
          else
             ! Évolution du champ scalaire de température
             fields%T(nx+1, j) = fields%T(nx, j)
          end if
       end do
    end if

    ! Bord Sud physique (y = 0)
    if (mesh%neighbor_south == MPI_PROC_NULL) then
       ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
       do i = 1, nx
          ! Évaluation de la condition logique
          if (cfg%bc_bottom%thermal_type == THERMAL_DIRICHLET) then
             ! Évolution du champ scalaire de température
             fields%T(i, 0) = 2.0_dp * cfg%bc_bottom%T - fields%T(i, 1)
          ! Branche par défaut si condition non satisfaite
          else
             ! Évolution du champ scalaire de température
             fields%T(i, 0) = fields%T(i, 1)
          end if
       end do
    end if

    ! Bord Nord physique (y = Ly)
    if (mesh%neighbor_north == MPI_PROC_NULL) then
       ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
       do i = 1, nx
          ! Évaluation de la condition logique
          if (cfg%bc_top%thermal_type == THERMAL_DIRICHLET) then
             ! Évolution du champ scalaire de température
             fields%T(i, ny+1) = 2.0_dp * cfg%bc_top%T - fields%T(i, ny)
          ! Branche par défaut si condition non satisfaite
          else
             ! Évolution du champ scalaire de température
             fields%T(i, ny+1) = fields%T(i, ny)
          end if
       end do
    end if
  end subroutine apply_thermal_bc

  ! [Procédure/Module] subroutine apply_vof_bc
  subroutine apply_vof_bc(fields, mesh, cfg)
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(inout) :: fields
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in) :: cfg

    ! Déclaration de variable locale / paramètre de calcul
    integer :: nx, ny, i
    real(dp) :: x_c, xc_top, half_w_top
    ! Affectation et calcul numérique
    nx = mesh%nx_local
    ! Affectation et calcul numérique
    ny = mesh%ny_local

    ! Affectation et calcul numérique
    if (mesh%neighbor_west == MPI_PROC_NULL)  fields%C(0, 1:ny)    = fields%C(1, 1:ny)
    ! Affectation et calcul numérique
    if (mesh%neighbor_east == MPI_PROC_NULL)  fields%C(nx+1, 1:ny) = fields%C(nx, 1:ny)
    ! Affectation et calcul numérique
    if (mesh%neighbor_south == MPI_PROC_NULL) fields%C(1:nx, 0)    = fields%C(1:nx, 1)
    ! Affectation et calcul numérique pour bord Nord
    if (mesh%neighbor_north == MPI_PROC_NULL) then
       if (trim(cfg%bc_top%profile) == "JET_INLET") then
          xc_top = merge(cfg%bc_top%jet_center, 0.5_dp * mesh%lx, cfg%bc_top%jet_center > 0.0_dp)
          half_w_top = merge(0.5_dp * cfg%bc_top%jet_width, 0.04_dp, cfg%bc_top%jet_width > 0.0_dp)
          do i = 1, nx
             x_c = x_center(mesh, i)
             if (abs(x_c - xc_top) <= half_w_top) then
                fields%C(i, ny+1) = 0.0_dp ! Injection continue de liquide (Fluide 2)
             else
                fields%C(i, ny+1) = fields%C(i, ny) ! Évacuation libre du gaz
             end if
          end do
       else
          fields%C(1:nx, ny+1) = fields%C(1:nx, ny)
       end if
    end if
  end subroutine apply_vof_bc

  ! Déclaration de variable locale / paramètre de calcul
  real(dp) function global_max_val(local_val, comm)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in) :: local_val
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(in) :: comm
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: res
    ! Déclaration de variable locale / paramètre de calcul
    integer :: ierr
    ! Appel à la bibliothèque de communication parallèle MPI
    call MPI_Allreduce(local_val, res, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
    ! Affectation et calcul numérique
    global_max_val = res
  end function global_max_val

  ! Déclaration de variable locale / paramètre de calcul
  real(dp) function global_min_val(local_val, comm)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in) :: local_val
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(in) :: comm
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: res
    ! Déclaration de variable locale / paramètre de calcul
    integer :: ierr
    ! Appel à la bibliothèque de communication parallèle MPI
    call MPI_Allreduce(local_val, res, 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
    ! Affectation et calcul numérique
    global_min_val = res
  end function global_min_val

  ! Déclaration de variable locale / paramètre de calcul
  real(dp) function global_sum_val(local_val, comm)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in) :: local_val
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(in) :: comm
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: res
    ! Déclaration de variable locale / paramètre de calcul
    integer :: ierr
    ! Appel à la bibliothèque de communication parallèle MPI
    call MPI_Allreduce(local_val, res, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    ! Affectation et calcul numérique
    global_sum_val = res
  end function global_sum_val

  ! [Procédure/Module] subroutine sync_u
  subroutine sync_u(mesh, u)
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(inout) :: u(0:mesh%nx_local+1, 0:mesh%ny_local+1)
    ! Appel de sous-programme modulaire
    call exchange_halos_scalar(u, mesh)
  end subroutine sync_u

  ! [Procédure/Module] subroutine sync_v
  subroutine sync_v(mesh, v)
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(inout) :: v(0:mesh%nx_local+1, 0:mesh%ny_local+1)
    ! Appel de sous-programme modulaire
    call exchange_halos_scalar(v, mesh)
  end subroutine sync_v

  ! [Procédure/Module] subroutine sync_cell_scalar
  subroutine sync_cell_scalar(mesh, s)
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(inout) :: s(0:mesh%nx_local+1, 0:mesh%ny_local+1)
    ! Appel de sous-programme modulaire
    call exchange_halos_scalar(s, mesh)
  end subroutine sync_cell_scalar

  ! [Procédure/Module] subroutine apply_boundary_conditions
  subroutine apply_boundary_conditions(mesh, cfg, fields, homogeneous)
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in) :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(inout) :: fields
    ! Déclaration de variable locale / paramètre de calcul
    logical, intent(in) :: homogeneous
    ! Appel de sous-programme modulaire
    call apply_velocity_bc(fields, mesh, cfg, homogeneous)
  end subroutine apply_boundary_conditions

  ! [Procédure/Module] subroutine apply_thermal_boundary_conditions
  subroutine apply_thermal_boundary_conditions(mesh, cfg, fields)
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in) :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(inout) :: fields
    ! Appel de sous-programme modulaire
    call apply_thermal_bc(fields, mesh, cfg)
  end subroutine apply_thermal_boundary_conditions

  ! Déclaration de variable locale / paramètre de calcul
  real(dp) function global_max(mesh, local_val)
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in) :: local_val
    ! Affectation et calcul numérique
    global_max = global_max_val(local_val, mesh%comm)
  end function global_max

  ! Déclaration de variable locale / paramètre de calcul
  real(dp) function global_min(mesh, local_val)
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in) :: local_val
    ! Affectation et calcul numérique
    global_min = global_min_val(local_val, mesh%comm)
  end function global_min

  ! Déclaration de variable locale / paramètre de calcul
  real(dp) function global_sum(mesh, local_val)
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in) :: local_val
    ! Affectation et calcul numérique
    global_sum = global_sum_val(local_val, mesh%comm)
  end function global_sum

  ! Déclaration de variable locale / paramètre de calcul
  real(dp) function global_dot_product(mesh, x, y, n)
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in) :: x(:), y(:)
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(in) :: n
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: local_dot
    ! Affectation et calcul numérique
    local_dot = dot_product(x(1:n), y(1:n))
    ! Affectation et calcul numérique
    global_dot_product = global_sum_val(local_dot, mesh%comm)
  end function global_dot_product

  ! Déclaration de variable locale / paramètre de calcul
  real(dp) function global_norm2(mesh, x, n)
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in) :: x(:)
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(in) :: n
    ! Évaluation de la norme du résidu ou de l'erreur numérique
    global_norm2 = sqrt(global_dot_product(mesh, x, x, n))
  end function global_norm2

end module module_mpi_comm
