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

  subroutine exchange_halos_scalar(field, mesh)
    type(mesh_2d_t), intent(in) :: mesh
    real(dp), intent(inout) :: field(0:mesh%nx_local+1, 0:mesh%ny_local+1)

    integer :: nx, ny, ierr
    integer :: reqs(8), n_reqs
    real(dp), allocatable :: s_left(:), r_left(:), s_right(:), r_right(:)
    real(dp), allocatable :: s_bot(:), r_bot(:), s_top(:), r_top(:)

    if (mesh%size <= 1) return

    nx = mesh%nx_local
    ny = mesh%ny_local

    allocate(s_left(ny), r_left(ny), s_right(ny), r_right(ny))
    allocate(s_bot(nx), r_bot(nx), s_top(nx), r_top(nx))

    s_left(:)  = field(1, 1:ny)
    s_right(:) = field(nx, 1:ny)

    n_reqs = 0
    if (mesh%neighbor_west /= MPI_PROC_NULL) then
       call MPI_Irecv(r_left, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 10, mesh%comm, reqs(n_reqs+1), ierr)
       call MPI_Isend(s_left, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 11, mesh%comm, reqs(n_reqs+2), ierr)
       n_reqs = n_reqs + 2
    end if
    if (mesh%neighbor_east /= MPI_PROC_NULL) then
       call MPI_Irecv(r_right, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 11, mesh%comm, reqs(n_reqs+1), ierr)
       call MPI_Isend(s_right, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 10, mesh%comm, reqs(n_reqs+2), ierr)
       n_reqs = n_reqs + 2
    end if

    if (n_reqs > 0) then
       call MPI_Waitall(n_reqs, reqs, MPI_STATUSES_IGNORE, ierr)
    end if

    if (mesh%neighbor_west /= MPI_PROC_NULL) field(0, 1:ny) = r_left(:)
    if (mesh%neighbor_east /= MPI_PROC_NULL) field(nx+1, 1:ny) = r_right(:)

    s_bot(:) = field(1:nx, 1)
    s_top(:) = field(1:nx, ny)

    n_reqs = 0
    if (mesh%neighbor_south /= MPI_PROC_NULL) then
       call MPI_Irecv(r_bot, nx, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 20, mesh%comm, reqs(n_reqs+1), ierr)
       call MPI_Isend(s_bot, nx, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 21, mesh%comm, reqs(n_reqs+2), ierr)
       n_reqs = n_reqs + 2
    end if
    if (mesh%neighbor_north /= MPI_PROC_NULL) then
       call MPI_Irecv(r_top, nx, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 21, mesh%comm, reqs(n_reqs+1), ierr)
       call MPI_Isend(s_top, nx, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 20, mesh%comm, reqs(n_reqs+2), ierr)
       n_reqs = n_reqs + 2
    end if

    if (n_reqs > 0) then
       call MPI_Waitall(n_reqs, reqs, MPI_STATUSES_IGNORE, ierr)
    end if

    if (mesh%neighbor_south /= MPI_PROC_NULL) field(1:nx, 0) = r_bot(:)
    if (mesh%neighbor_north /= MPI_PROC_NULL) field(1:nx, ny+1) = r_top(:)
  end subroutine exchange_halos_scalar

  subroutine exchange_halos_velocity(u, v, mesh)
    type(mesh_2d_t), intent(in) :: mesh
    real(dp), intent(inout) :: u(0:mesh%nx_local+1, 0:mesh%ny_local+1)
    real(dp), intent(inout) :: v(0:mesh%nx_local+1, 0:mesh%ny_local+1)

    integer :: nx, ny, ierr
    integer :: reqs(8), n_reqs
    real(dp), allocatable :: su_l(:), ru_l(:), su_r(:), ru_r(:)
    real(dp), allocatable :: sv_l(:), rv_l(:), sv_r(:), rv_r(:)
    real(dp), allocatable :: su_b(:), ru_b(:), su_t(:), ru_t(:)
    real(dp), allocatable :: sv_b(:), rv_b(:), sv_t(:), rv_t(:)

    if (mesh%size <= 1) return

    nx = mesh%nx_local
    ny = mesh%ny_local

    allocate(su_l(ny), ru_l(ny), su_r(ny), ru_r(ny))
    allocate(sv_l(ny), rv_l(ny), sv_r(ny), rv_r(ny))

    su_l(:) = u(1, 1:ny)
    su_r(:) = u(nx, 1:ny)
    sv_l(:) = v(1, 1:ny)
    sv_r(:) = v(nx, 1:ny)

    n_reqs = 0
    if (mesh%neighbor_west /= MPI_PROC_NULL) then
       call MPI_Irecv(ru_l, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 30, mesh%comm, reqs(n_reqs+1), ierr)
       call MPI_Isend(su_l, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 31, mesh%comm, reqs(n_reqs+2), ierr)
       call MPI_Irecv(rv_l, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 32, mesh%comm, reqs(n_reqs+3), ierr)
       call MPI_Isend(sv_l, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 33, mesh%comm, reqs(n_reqs+4), ierr)
       n_reqs = n_reqs + 4
    end if
    if (mesh%neighbor_east /= MPI_PROC_NULL) then
       call MPI_Irecv(ru_r, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 31, mesh%comm, reqs(n_reqs+1), ierr)
       call MPI_Isend(su_r, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 30, mesh%comm, reqs(n_reqs+2), ierr)
       call MPI_Irecv(rv_r, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 33, mesh%comm, reqs(n_reqs+3), ierr)
       call MPI_Isend(sv_r, ny, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 32, mesh%comm, reqs(n_reqs+4), ierr)
       n_reqs = n_reqs + 4
    end if

    if (n_reqs > 0) call MPI_Waitall(n_reqs, reqs, MPI_STATUSES_IGNORE, ierr)

    if (mesh%neighbor_west /= MPI_PROC_NULL) then
       u(0, 1:ny) = ru_l(:)
       v(0, 1:ny) = rv_l(:)
    end if
    if (mesh%neighbor_east /= MPI_PROC_NULL) then
       u(nx+1, 1:ny) = ru_r(:)
       v(nx+1, 1:ny) = rv_r(:)
    end if

    allocate(su_b(nx+2), ru_b(nx+2), su_t(nx+2), ru_t(nx+2))
    allocate(sv_b(nx),   rv_b(nx),   sv_t(nx),   rv_t(nx))

    su_b(:) = u(0:nx+1, 1)
    su_t(:) = u(0:nx+1, ny)
    sv_b(:) = v(1:nx, 1)
    sv_t(:) = v(1:nx, ny)

    n_reqs = 0
    if (mesh%neighbor_south /= MPI_PROC_NULL) then
       call MPI_Irecv(ru_b, nx+2, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 40, mesh%comm, reqs(n_reqs+1), ierr)
       call MPI_Isend(su_b, nx+2, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 41, mesh%comm, reqs(n_reqs+2), ierr)
       call MPI_Irecv(rv_b, nx,   MPI_DOUBLE_PRECISION, mesh%neighbor_south, 42, mesh%comm, reqs(n_reqs+3), ierr)
       call MPI_Isend(sv_b, nx,   MPI_DOUBLE_PRECISION, mesh%neighbor_south, 43, mesh%comm, reqs(n_reqs+4), ierr)
       n_reqs = n_reqs + 4
    end if
    if (mesh%neighbor_north /= MPI_PROC_NULL) then
       call MPI_Irecv(ru_t, nx+2, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 41, mesh%comm, reqs(n_reqs+1), ierr)
       call MPI_Isend(su_t, nx+2, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 40, mesh%comm, reqs(n_reqs+2), ierr)
       call MPI_Irecv(rv_t, nx,   MPI_DOUBLE_PRECISION, mesh%neighbor_north, 43, mesh%comm, reqs(n_reqs+3), ierr)
       call MPI_Isend(sv_t, nx,   MPI_DOUBLE_PRECISION, mesh%neighbor_north, 42, mesh%comm, reqs(n_reqs+4), ierr)
       n_reqs = n_reqs + 4
    end if

    if (n_reqs > 0) call MPI_Waitall(n_reqs, reqs, MPI_STATUSES_IGNORE, ierr)

    if (mesh%neighbor_south /= MPI_PROC_NULL) then
       u(0:nx+1, 0) = ru_b(:)
       v(1:nx, 0)   = rv_b(:)
    end if
    if (mesh%neighbor_north /= MPI_PROC_NULL) then
       u(0:nx+1, ny+1) = ru_t(:)
       v(1:nx, ny+1)   = rv_t(:)
    end if
  end subroutine exchange_halos_velocity

  subroutine apply_velocity_bc(fields, mesh, cfg, homogeneous)
    type(flow_fields_t), intent(inout) :: fields
    type(mesh_2d_t), intent(in) :: mesh
    type(config_t), intent(in) :: cfg
    logical, intent(in) :: homogeneous

    integer :: nx, ny, i, j
    real(dp) :: u_l, v_l, u_r, v_r, u_b, v_b, u_t, v_t
    real(dp) :: y, h, u_inlet

    nx = mesh%nx_local
    ny = mesh%ny_local

    u_l = merge(0.0_dp, cfg%bc_left%u, homogeneous)
    v_l = merge(0.0_dp, cfg%bc_left%v, homogeneous)
    u_r = merge(0.0_dp, cfg%bc_right%u, homogeneous)
    v_r = merge(0.0_dp, cfg%bc_right%v, homogeneous)
    u_b = merge(0.0_dp, cfg%bc_bottom%u, homogeneous)
    v_b = merge(0.0_dp, cfg%bc_bottom%v, homogeneous)
    u_t = merge(0.0_dp, cfg%bc_top%u, homogeneous)
    v_t = merge(0.0_dp, cfg%bc_top%v, homogeneous)

    ! Traitement périodicité si 1 seul processus
    if (cfg%bc_left%btype == BC_PERIODIC .and. mesh%px == 1) then
       fields%u(0, 1:ny)    = fields%u(nx, 1:ny)
       fields%u(nx+1, 1:ny) = fields%u(1, 1:ny)
       fields%v(0, 1:ny)    = fields%v(nx, 1:ny)
       fields%v(nx+1, 1:ny) = fields%v(1, 1:ny)
    end if
    if (cfg%bc_bottom%btype == BC_PERIODIC .and. mesh%py == 1) then
       fields%v(1:nx, 0)    = fields%v(1:nx, ny)
       fields%v(1:nx, ny+1) = fields%v(1:nx, 1)
       fields%u(0:nx+1, 0)    = fields%u(0:nx+1, ny)
       fields%u(0:nx+1, ny+1) = fields%u(0:nx+1, 1)
    end if

    ! Bord Ouest physique (x = 0)
    if (mesh%neighbor_west == MPI_PROC_NULL .and. cfg%bc_left%btype /= BC_PERIODIC) then
       do j = 1, ny
          if (cfg%bc_left%btype == BC_DIRICHLET) then
             u_inlet = u_l
             if (.not. homogeneous .and. cfg%obs_enabled .and. cfg%obs_type == OBS_BOX) then
                ! Profil d'entrée marche descendante : parabole sur [y_step, ly]
                y = y_center(mesh, j)
                if (y < cfg%obs_y1) then
                   u_inlet = 0.0_dp
                else
                   h = mesh%ly - cfg%obs_y1
                   u_inlet = 6.0_dp * cfg%bc_left%u * (y - cfg%obs_y1) * (mesh%ly - y) / (h * h)
                end if
             end if
             if (.not. homogeneous) fields%u(1, j) = u_inlet
             fields%u(0, j) = u_inlet
             fields%v(0, j) = 2.0_dp * v_l - fields%v(1, j)
          else if (cfg%bc_left%btype == BC_FREESLIP) then
             if (.not. homogeneous) fields%u(1, j) = 0.0_dp
             fields%u(0, j) = 0.0_dp
             fields%v(0, j) = fields%v(1, j)
          else if (cfg%bc_left%btype == BC_NEUMANN) then
             fields%u(0, j) = fields%u(1, j)
             fields%v(0, j) = fields%v(1, j)
          end if
       end do
    end if

    ! Bord Est physique (x = Lx)
    if (mesh%neighbor_east == MPI_PROC_NULL .and. cfg%bc_right%btype /= BC_PERIODIC) then
       do j = 1, ny
          if (cfg%bc_right%btype == BC_DIRICHLET) then
             if (.not. homogeneous) fields%u(nx+1, j) = u_r
             fields%u(nx+1, j) = u_r
             fields%v(nx+1, j) = 2.0_dp * v_r - fields%v(nx, j)
          else if (cfg%bc_right%btype == BC_FREESLIP) then
             if (.not. homogeneous) fields%u(nx+1, j) = 0.0_dp
             fields%u(nx+1, j) = 0.0_dp
             fields%v(nx+1, j) = fields%v(nx, j)
          else if (cfg%bc_right%btype == BC_NEUMANN) then
             fields%u(nx+1, j) = fields%u(nx, j)
             fields%v(nx+1, j) = fields%v(nx, j)
          end if
       end do
    end if

    ! Bord Sud physique (y = 0)
    if (mesh%neighbor_south == MPI_PROC_NULL .and. cfg%bc_bottom%btype /= BC_PERIODIC) then
       do i = 1, nx
          if (cfg%bc_bottom%btype == BC_DIRICHLET) then
             if (.not. homogeneous) fields%v(i, 1) = v_b
             fields%v(i, 0) = v_b
             fields%u(i, 0) = 2.0_dp * u_b - fields%u(i, 1)
          else if (cfg%bc_bottom%btype == BC_FREESLIP) then
             if (.not. homogeneous) fields%v(i, 1) = 0.0_dp
             fields%v(i, 0) = 0.0_dp
             fields%u(i, 0) = fields%u(i, 1)
          else if (cfg%bc_bottom%btype == BC_NEUMANN) then
             fields%u(i, 0) = fields%u(i, 1)
             fields%v(i, 0) = fields%v(i, 1)
          end if
       end do
    end if

    ! Bord Nord physique (y = Ly)
    if (mesh%neighbor_north == MPI_PROC_NULL .and. cfg%bc_top%btype /= BC_PERIODIC) then
       do i = 1, nx
          if (cfg%bc_top%btype == BC_DIRICHLET) then
             if (.not. homogeneous) fields%v(i, ny+1) = v_t
             fields%v(i, ny+1) = v_t
             fields%u(i, ny+1) = 2.0_dp * u_t - fields%u(i, ny)
          else if (cfg%bc_top%btype == BC_FREESLIP) then
             if (.not. homogeneous) fields%v(i, ny+1) = 0.0_dp
             fields%v(i, ny+1) = 0.0_dp
             fields%u(i, ny+1) = fields%u(i, ny)
          else if (cfg%bc_top%btype == BC_NEUMANN) then
             fields%u(i, ny+1) = fields%u(i, ny)
             fields%v(i, ny+1) = fields%v(i, ny)
          end if
       end do
    end if
  end subroutine apply_velocity_bc

  subroutine apply_thermal_bc(fields, mesh, cfg)
    type(flow_fields_t), intent(inout) :: fields
    type(mesh_2d_t), intent(in) :: mesh
    type(config_t), intent(in) :: cfg

    integer :: nx, ny, i, j
    nx = mesh%nx_local
    ny = mesh%ny_local

    ! Périodicité thermique
    if (cfg%bc_left%thermal_type == BC_PERIODIC .and. mesh%px == 1) then
       fields%T(0, 1:ny)    = fields%T(nx, 1:ny)
       fields%T(nx+1, 1:ny) = fields%T(1, 1:ny)
    end if
    if (cfg%bc_bottom%thermal_type == BC_PERIODIC .and. mesh%py == 1) then
       fields%T(1:nx, 0)    = fields%T(1:nx, ny)
       fields%T(1:nx, ny+1) = fields%T(1:nx, 1)
    end if

    ! Bord Ouest physique (x = 0)
    if (mesh%neighbor_west == MPI_PROC_NULL) then
       do j = 1, ny
          if (cfg%bc_left%thermal_type == THERMAL_DIRICHLET) then
             fields%T(0, j) = 2.0_dp * cfg%bc_left%T - fields%T(1, j)
          else
             fields%T(0, j) = fields%T(1, j)
          end if
       end do
    end if

    ! Bord Est physique (x = Lx)
    if (mesh%neighbor_east == MPI_PROC_NULL) then
       do j = 1, ny
          if (cfg%bc_right%thermal_type == THERMAL_DIRICHLET) then
             fields%T(nx+1, j) = 2.0_dp * cfg%bc_right%T - fields%T(nx, j)
          else
             fields%T(nx+1, j) = fields%T(nx, j)
          end if
       end do
    end if

    ! Bord Sud physique (y = 0)
    if (mesh%neighbor_south == MPI_PROC_NULL) then
       do i = 1, nx
          if (cfg%bc_bottom%thermal_type == THERMAL_DIRICHLET) then
             fields%T(i, 0) = 2.0_dp * cfg%bc_bottom%T - fields%T(i, 1)
          else
             fields%T(i, 0) = fields%T(i, 1)
          end if
       end do
    end if

    ! Bord Nord physique (y = Ly)
    if (mesh%neighbor_north == MPI_PROC_NULL) then
       do i = 1, nx
          if (cfg%bc_top%thermal_type == THERMAL_DIRICHLET) then
             fields%T(i, ny+1) = 2.0_dp * cfg%bc_top%T - fields%T(i, ny)
          else
             fields%T(i, ny+1) = fields%T(i, ny)
          end if
       end do
    end if
  end subroutine apply_thermal_bc

  subroutine apply_vof_bc(fields, mesh, cfg)
    type(flow_fields_t), intent(inout) :: fields
    type(mesh_2d_t), intent(in) :: mesh
    type(config_t), intent(in) :: cfg

    integer :: nx, ny
    nx = mesh%nx_local
    ny = mesh%ny_local

    if (mesh%neighbor_west == MPI_PROC_NULL)  fields%C(0, 1:ny)    = fields%C(1, 1:ny)
    if (mesh%neighbor_east == MPI_PROC_NULL)  fields%C(nx+1, 1:ny) = fields%C(nx, 1:ny)
    if (mesh%neighbor_south == MPI_PROC_NULL) fields%C(1:nx, 0)    = fields%C(1:nx, 1)
    if (mesh%neighbor_north == MPI_PROC_NULL) fields%C(1:nx, ny+1) = fields%C(1:nx, ny)
  end subroutine apply_vof_bc

  real(dp) function global_max_val(local_val, comm)
    real(dp), intent(in) :: local_val
    integer, intent(in) :: comm
    real(dp) :: res
    integer :: ierr
    call MPI_Allreduce(local_val, res, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
    global_max_val = res
  end function global_max_val

  real(dp) function global_min_val(local_val, comm)
    real(dp), intent(in) :: local_val
    integer, intent(in) :: comm
    real(dp) :: res
    integer :: ierr
    call MPI_Allreduce(local_val, res, 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
    global_min_val = res
  end function global_min_val

  real(dp) function global_sum_val(local_val, comm)
    real(dp), intent(in) :: local_val
    integer, intent(in) :: comm
    real(dp) :: res
    integer :: ierr
    call MPI_Allreduce(local_val, res, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    global_sum_val = res
  end function global_sum_val

  subroutine sync_u(mesh, u)
    type(mesh_2d_t), intent(in) :: mesh
    real(dp), intent(inout) :: u(0:mesh%nx_local+1, 0:mesh%ny_local+1)
    call exchange_halos_scalar(u, mesh)
  end subroutine sync_u

  subroutine sync_v(mesh, v)
    type(mesh_2d_t), intent(in) :: mesh
    real(dp), intent(inout) :: v(0:mesh%nx_local+1, 0:mesh%ny_local+1)
    call exchange_halos_scalar(v, mesh)
  end subroutine sync_v

  subroutine sync_cell_scalar(mesh, s)
    type(mesh_2d_t), intent(in) :: mesh
    real(dp), intent(inout) :: s(0:mesh%nx_local+1, 0:mesh%ny_local+1)
    call exchange_halos_scalar(s, mesh)
  end subroutine sync_cell_scalar

  subroutine apply_boundary_conditions(mesh, cfg, fields, homogeneous)
    type(mesh_2d_t), intent(in) :: mesh
    type(config_t), intent(in) :: cfg
    type(flow_fields_t), intent(inout) :: fields
    logical, intent(in) :: homogeneous
    call apply_velocity_bc(fields, mesh, cfg, homogeneous)
  end subroutine apply_boundary_conditions

  subroutine apply_thermal_boundary_conditions(mesh, cfg, fields)
    type(mesh_2d_t), intent(in) :: mesh
    type(config_t), intent(in) :: cfg
    type(flow_fields_t), intent(inout) :: fields
    call apply_thermal_bc(fields, mesh, cfg)
  end subroutine apply_thermal_boundary_conditions

  real(dp) function global_max(mesh, local_val)
    type(mesh_2d_t), intent(in) :: mesh
    real(dp), intent(in) :: local_val
    global_max = global_max_val(local_val, mesh%comm)
  end function global_max

  real(dp) function global_min(mesh, local_val)
    type(mesh_2d_t), intent(in) :: mesh
    real(dp), intent(in) :: local_val
    global_min = global_min_val(local_val, mesh%comm)
  end function global_min

  real(dp) function global_sum(mesh, local_val)
    type(mesh_2d_t), intent(in) :: mesh
    real(dp), intent(in) :: local_val
    global_sum = global_sum_val(local_val, mesh%comm)
  end function global_sum

  real(dp) function global_dot_product(mesh, x, y, n)
    type(mesh_2d_t), intent(in) :: mesh
    real(dp), intent(in) :: x(:), y(:)
    integer, intent(in) :: n
    real(dp) :: local_dot
    local_dot = dot_product(x(1:n), y(1:n))
    global_dot_product = global_sum_val(local_dot, mesh%comm)
  end function global_dot_product

  real(dp) function global_norm2(mesh, x, n)
    type(mesh_2d_t), intent(in) :: mesh
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: n
    global_norm2 = sqrt(global_dot_product(mesh, x, x, n))
  end function global_norm2

end module module_mpi_comm
