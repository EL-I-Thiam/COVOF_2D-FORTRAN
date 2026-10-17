module module_mesh
  use module_precision, only: dp
  implicit none
  private
  include 'mpif.h'

  public :: mesh_2d_t, init_mesh, destroy_mesh
  public :: x_face, y_face, x_center, y_center, is_in_obstacle

  type :: mesh_2d_t
     integer  :: nx_global = 64, ny_global = 64
     real(dp) :: lx = 1.0_dp, ly = 1.0_dp
     real(dp) :: dx = 0.0_dp, dy = 0.0_dp

     ! MPI topology
     integer  :: rank = 0, size = 1
     integer  :: comm = MPI_COMM_WORLD
     integer  :: cart_comm = MPI_COMM_NULL
     integer  :: px = 1, py = 1
     integer  :: coord_x = 0, coord_y = 0
     integer  :: neighbor_west = MPI_PROC_NULL
     integer  :: neighbor_east = MPI_PROC_NULL
     integer  :: neighbor_south = MPI_PROC_NULL
     integer  :: neighbor_north = MPI_PROC_NULL

     ! Boundary flags
     logical  :: is_boundary_left = .false.
     logical  :: is_boundary_right = .false.
     logical  :: is_boundary_bottom = .false.
     logical  :: is_boundary_top = .false.

     ! Local partitioning
     integer  :: nx_local = 64, ny_local = 64
     integer  :: i_start_global = 1, j_start_global = 1
  end type mesh_2d_t

contains

  subroutine init_mesh(mesh, comm_world, nx, ny, lx, ly, periodic_x, periodic_y)
    type(mesh_2d_t), intent(out) :: mesh
    integer, intent(in) :: comm_world, nx, ny
    real(dp), intent(in) :: lx, ly
    logical, intent(in) :: periodic_x, periodic_y

    integer :: dims(2), periods(2), coords(2)
    integer :: ierr, rank_in, size_in
    integer :: rem_x, rem_y

    mesh%comm = comm_world
    call MPI_Comm_rank(comm_world, rank_in, ierr)
    call MPI_Comm_size(comm_world, size_in, ierr)
    mesh%rank = rank_in
    mesh%size = size_in

    mesh%nx_global = nx
    mesh%ny_global = ny
    mesh%lx = lx
    mesh%ly = ly
    mesh%dx = lx / real(nx, dp)
    mesh%dy = ly / real(ny, dp)

    ! Décomposition 2D cartésienne optimale
    dims = [0, 0]
    call MPI_Dims_create(mesh%size, 2, dims, ierr)
    mesh%px = dims(1)
    mesh%py = dims(2)

    periods(1) = merge(1, 0, periodic_x)
    periods(2) = merge(1, 0, periodic_y)

    call MPI_Cart_create(comm_world, 2, dims, periods, 1, mesh%cart_comm, ierr)
    call MPI_Comm_rank(mesh%cart_comm, mesh%rank, ierr)
    call MPI_Cart_coords(mesh%cart_comm, mesh%rank, 2, coords, ierr)
    mesh%coord_x = coords(1)
    mesh%coord_y = coords(2)

    call MPI_Cart_shift(mesh%cart_comm, 0, 1, mesh%neighbor_west, mesh%neighbor_east, ierr)
    call MPI_Cart_shift(mesh%cart_comm, 1, 1, mesh%neighbor_south, mesh%neighbor_north, ierr)

    ! Répartition équilibrée des cellules entre processus
    mesh%nx_local = nx / mesh%px
    rem_x = mod(nx, mesh%px)
    if (mesh%coord_x < rem_x) then
       mesh%nx_local = mesh%nx_local + 1
       mesh%i_start_global = mesh%coord_x * mesh%nx_local + 1
    else
       mesh%i_start_global = rem_x * (mesh%nx_local + 1) + (mesh%coord_x - rem_x) * mesh%nx_local + 1
    end if

    mesh%ny_local = ny / mesh%py
    rem_y = mod(ny, mesh%py)
    if (mesh%coord_y < rem_y) then
       mesh%ny_local = mesh%ny_local + 1
       mesh%j_start_global = mesh%coord_y * mesh%ny_local + 1
    else
       mesh%j_start_global = rem_y * (mesh%ny_local + 1) + (mesh%coord_y - rem_y) * mesh%ny_local + 1
    end if

    mesh%is_boundary_left   = (mesh%neighbor_west == MPI_PROC_NULL)
    mesh%is_boundary_right  = (mesh%neighbor_east == MPI_PROC_NULL)
    mesh%is_boundary_bottom = (mesh%neighbor_south == MPI_PROC_NULL)
    mesh%is_boundary_top    = (mesh%neighbor_north == MPI_PROC_NULL)
  end subroutine init_mesh

  subroutine destroy_mesh(mesh)
    type(mesh_2d_t), intent(inout) :: mesh
    integer :: ierr
    if (mesh%cart_comm /= MPI_COMM_NULL .and. mesh%cart_comm /= mesh%comm) then
       call MPI_Comm_free(mesh%cart_comm, ierr)
    end if
  end subroutine destroy_mesh

  pure real(dp) function x_face(mesh, i)
    type(mesh_2d_t), intent(in) :: mesh
    integer, intent(in) :: i
    x_face = real(mesh%i_start_global - 1 + i - 1, dp) * mesh%dx
  end function x_face

  pure real(dp) function y_face(mesh, j)
    type(mesh_2d_t), intent(in) :: mesh
    integer, intent(in) :: j
    y_face = real(mesh%j_start_global - 1 + j - 1, dp) * mesh%dy
  end function y_face

  pure real(dp) function x_center(mesh, i)
    type(mesh_2d_t), intent(in) :: mesh
    integer, intent(in) :: i
    x_center = (real(mesh%i_start_global - 1 + i - 1, dp) + 0.5_dp) * mesh%dx
  end function x_center

  pure real(dp) function y_center(mesh, j)
    type(mesh_2d_t), intent(in) :: mesh
    integer, intent(in) :: j
    y_center = (real(mesh%j_start_global - 1 + j - 1, dp) + 0.5_dp) * mesh%dy
  end function y_center

  pure logical function is_in_obstacle(x, y, obs_type, x0, y0, x1, y1, xc, yc, r)
    real(dp), intent(in) :: x, y
    integer, intent(in)  :: obs_type
    real(dp), intent(in) :: x0, y0, x1, y1, xc, yc, r

    is_in_obstacle = .false.
    if (obs_type == 1) then ! OBS_BOX
       if (x >= x0 - 1.0e-9_dp .and. x <= x1 + 1.0e-9_dp .and. &
           y >= y0 - 1.0e-9_dp .and. y <= y1 + 1.0e-9_dp) then
          is_in_obstacle = .true.
       end if
    else if (obs_type == 2) then ! OBS_CYLINDER
       if ((x - xc)**2 + (y - yc)**2 <= r**2) then
          is_in_obstacle = .true.
       end if
    end if
  end function is_in_obstacle

end module module_mesh
