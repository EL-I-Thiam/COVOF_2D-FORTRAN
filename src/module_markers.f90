!===============================================================================
! Module        : module_markers
! Description   : Gestion des marqueurs lagrangiens d'interface non connectés
!                 (Unconnected Marker-VOF) et interpolation conservative de
!                 vitesse (CVI - Conservative Velocity Interpolation) à
!                 divergence analytiquement nulle sur grille décalée MAC 2D.
!                 - Ensemencement initial et ré-ensemencement local depuis PLIC
!                 - Interpolation CVI (H(div)-conforme) et bilinéaire sur grille MAC
!                 - Intégration temporelle explicite RK2-Heun
!                 - Migration parallèle MPI en 2 passes (Est/Ouest puis Nord/Sud)
!                 - Estimation locale de la courbure et des normales par
!                   ajustement parabolique pondéré dans le repère tourné (tau, eta)
!                   avec filtrage angulaire pour ligaments et films minces
!                 - Export VTK PolyData (.vtk) du nuage de marqueurs
! Auteur        : Elhadji Ibrahima Thiam
! Standard      : Fortran 2018 / MPI-3
!===============================================================================
module module_markers
  use module_precision, only: dp
  use module_constants
  use module_mesh, only: mesh_2d_t, x_face, y_face, x_center, y_center
  use module_config, only: config_t
  use module_fields, only: flow_fields_t
  use module_mpi_comm, only: exchange_halos_scalar
  use module_interface_reconstruction, only: interface_segment_t, reconstruct_elvira, &
                                             reconstruct_interface, clip_segment_in_cell, &
                                             compute_alpha
  implicit none
  private
  include 'mpif.h'

  public :: marker_t, marker_cloud_t
  public :: init_markers_from_vof
  public :: interpolate_velocity_cvi
  public :: advect_markers_rk2
  public :: migrate_markers_mpi
  public :: prune_and_reseed_markers
  public :: compute_curvature_from_markers
  public :: write_markers_vtk
  public :: free_markers

  !-----------------------------------------------------------------------------
  ! Structure d'une particule-marqueur d'interface sans connectivité
  !-----------------------------------------------------------------------------
  type :: marker_t
     real(dp) :: x = 0.0_dp       ! Abscisse physique absolue [m]
     real(dp) :: y = 0.0_dp       ! Ordonnée physique absolue [m]
     real(dp) :: nx = 0.0_dp      ! Normale unitaire X orientée vers C=1
     real(dp) :: ny = 0.0_dp      ! Normale unitaire Y orientée vers C=1
     real(dp) :: kappa = 0.0_dp   ! Courbure locale portée par le marqueur [1/m]
     integer  :: i_cell = 1       ! Indice local X de la cellule hôte (1..nx_local)
     integer  :: j_cell = 1       ! Indice local Y de la cellule hôte (1..ny_local)
     logical  :: active = .true.  ! Indicateur de validité du marqueur
  end type marker_t

  !-----------------------------------------------------------------------------
  ! Conteneur dynamique de marqueurs pour le sous-domaine MPI local
  !-----------------------------------------------------------------------------
  type :: marker_cloud_t
     type(marker_t), allocatable :: list(:)
     integer :: count = 0
     integer :: capacity = 0
     integer, allocatable :: count_in_cell(:, :)
  end type marker_cloud_t

contains

  !=============================================================================
  ! Routine interne : ensure_capacity
  ! Rôle            : Redimensionne le tableau de marqueurs si nécessaire
  !=============================================================================
  subroutine ensure_capacity(cloud, min_cap)
    type(marker_cloud_t), intent(inout) :: cloud
    integer, intent(in)                 :: min_cap
    type(marker_t), allocatable :: tmp(:)
    integer :: new_cap

    if (.not. allocated(cloud%list)) then
       new_cap = max(min_cap, 1024)
       allocate(cloud%list(new_cap))
       cloud%capacity = new_cap
       cloud%count = 0
    else if (cloud%capacity < min_cap) then
       new_cap = max(min_cap, 2 * cloud%capacity)
       allocate(tmp(new_cap))
       if (cloud%count > 0) then
          tmp(1:cloud%count) = cloud%list(1:cloud%count)
       end if
       call move_alloc(tmp, cloud%list)
       cloud%capacity = new_cap
    end if
  end subroutine ensure_capacity

  !=============================================================================
  ! Routine : free_markers
  !=============================================================================
  subroutine free_markers(cloud)
    type(marker_cloud_t), intent(inout) :: cloud
    if (allocated(cloud%list)) deallocate(cloud%list)
    if (allocated(cloud%count_in_cell)) deallocate(cloud%count_in_cell)
    cloud%count = 0
    cloud%capacity = 0
  end subroutine free_markers

  !=============================================================================
  ! Routine : init_markers_from_vof
  ! Rôle    : Initialise le nuage de marqueurs à t = 0 sur les segments PLIC
  !           des cellules d'interface (0 < C < 1)
  !=============================================================================
  subroutine init_markers_from_vof(cloud, mesh, cfg, fields)
    type(marker_cloud_t), intent(inout) :: cloud
    type(mesh_2d_t),      intent(in)    :: mesh
    type(config_t),       intent(in)    :: cfg
    type(flow_fields_t),  intent(inout) :: fields

    type(interface_segment_t), allocatable :: segments(:, :)
    integer  :: i, j, k, nx, ny, n_per_cell
    real(dp) :: c_val, xmin, xmax, ymin, ymax, x1, y1, x2, y2, s
    real(dp) :: xp, yp, xp_proj, yp_proj, rx, ry, r_dist, sgn_c1
    logical  :: valid, is_circular_ic

    if (.not. cfg%use_markers) return

    nx = mesh%nx_local
    ny = mesh%ny_local
    n_per_cell = max(2, cfg%markers_per_cell)

    call ensure_capacity(cloud, max(1024, 4 * (nx + ny) * n_per_cell))
    if (allocated(cloud%count_in_cell)) deallocate(cloud%count_in_cell)
    allocate(cloud%count_in_cell(0:nx+1, 0:ny+1))
    cloud%count_in_cell = 0
    cloud%count = 0

    allocate(segments(0:nx+1, 0:ny+1))
    if (trim(cfg%reconstruction_method) == "NONE" .or. trim(cfg%reconstruction_method) == "SLIC") then
       call reconstruct_interface("ELVIRA", mesh, fields%C, segments)
    else
       call reconstruct_interface(cfg%reconstruction_method, mesh, fields%C, segments)
    end if

    is_circular_ic = (cfg%ic_type == IC_DROPLET .or. cfg%ic_type == IC_BUBBLE .or. &
                      cfg%ic_type == IC_CYLINDER .or. cfg%ic_type == IC_VORTEX_SHEAR) .and. &
                     (cfg%ic_radius > 0.0_dp) .and. &
                     (cfg%ic_axis_a <= 0.0_dp .or. cfg%ic_axis_b <= 0.0_dp .or. abs(cfg%ic_axis_a - cfg%ic_axis_b) < 1.0e-9_dp)

    do j = 1, ny
       ymin = y_face(mesh, j)
       ymax = y_face(mesh, j + 1)
       do i = 1, nx
          c_val = fields%C(i, j)
          if (c_val <= 1.0e-2_dp .or. c_val >= 1.0_dp - 1.0e-2_dp) cycle

          xmin = x_face(mesh, i)
          xmax = x_face(mesh, i + 1)
          call clip_segment_in_cell(xmin, xmax, ymin, ymax, &
                                    segments(i, j)%nx, segments(i, j)%ny, segments(i, j)%alpha, &
                                    x1, y1, x2, y2, valid)
          if (.not. valid) cycle

          call ensure_capacity(cloud, cloud%count + n_per_cell + 16)
          do k = 1, n_per_cell
             s = (real(k, dp) - 0.5_dp) / real(n_per_cell, dp)
             xp = x1 + s * (x2 - x1)
             yp = y1 + s * (y2 - y1)

             if (is_circular_ic) then
                rx = xp - cfg%ic_center_x
                ry = yp - cfg%ic_center_y
                r_dist = hypot(rx, ry)
                if (r_dist > 1.0e-12_dp .and. abs(r_dist - cfg%ic_radius) <= 0.75_dp * hypot(mesh%dx, mesh%dy)) then
                   xp_proj = cfg%ic_center_x + cfg%ic_radius * (rx / r_dist)
                   yp_proj = cfg%ic_center_y + cfg%ic_radius * (ry / r_dist)
                   if (xp_proj >= xmin - 1.0e-12_dp .and. xp_proj <= xmax + 1.0e-12_dp .and. &
                       yp_proj >= ymin - 1.0e-12_dp .and. yp_proj <= ymax + 1.0e-12_dp) then
                      xp = max(xmin, min(xmax, xp_proj))
                      yp = max(ymin, min(ymax, yp_proj))
                      cloud%count = cloud%count + 1
                      if (cfg%ic_type == IC_DROPLET .or. cfg%ic_type == IC_VORTEX_SHEAR) then
                         sgn_c1 = -1.0_dp
                      else
                         sgn_c1 = 1.0_dp
                      end if
                      cloud%list(cloud%count)%nx = sgn_c1 * (rx / r_dist)
                      cloud%list(cloud%count)%ny = sgn_c1 * (ry / r_dist)
                      cloud%list(cloud%count)%kappa = -sgn_c1 / cfg%ic_radius
                   else
                      cloud%count = cloud%count + 1
                      cloud%list(cloud%count)%nx = -segments(i, j)%nx
                      cloud%list(cloud%count)%ny = -segments(i, j)%ny
                      cloud%list(cloud%count)%kappa = 0.0_dp
                   end if
                else
                   cloud%count = cloud%count + 1
                   cloud%list(cloud%count)%nx = -segments(i, j)%nx
                   cloud%list(cloud%count)%ny = -segments(i, j)%ny
                   cloud%list(cloud%count)%kappa = 0.0_dp
                end if
             else
                cloud%count = cloud%count + 1
                cloud%list(cloud%count)%nx = -segments(i, j)%nx
                cloud%list(cloud%count)%ny = -segments(i, j)%ny
                cloud%list(cloud%count)%kappa = 0.0_dp
             end if

             cloud%list(cloud%count)%x = xp
             cloud%list(cloud%count)%y = yp
             cloud%list(cloud%count)%i_cell = i
             cloud%list(cloud%count)%j_cell = j
             cloud%list(cloud%count)%active = .true.
             cloud%count_in_cell(i, j) = cloud%count_in_cell(i, j) + 1
          end do
       end do
    end do

    deallocate(segments)
    call compute_curvature_from_markers(cloud, mesh, cfg, fields)
  end subroutine init_markers_from_vof

  !=============================================================================
  ! Routine pure : interpolate_velocity_cvi
  ! Rôle         : Interpolation Conservative de Vitesse (CVI) H(div)-conforme
  !                sur grille décalée MAC — Variante B exacte (Wang et al., 2015 ;
  !                Pusok et al., 2017) :
  !                - Interpole u sur sa boîte décalée [x_i^f, x_{i+1}^f] x [y_{j_u}^c, y_{j_u+1}^c]
  !                  avec j_u in {j-1, j} ;
  !                - Interpole v sur sa boîte décalée [x_{i_v}^c, x_{i_v+1}^c] x [y_j^f, y_{j+1}^f]
  !                  avec i_v in {i-1, i} ;
  !                - Garantit analytiquement dans les 4 sous-quadrants de (i, j) :
  !                  du_CVI/dx + dv_CVI/dy = (u(i+1,j)-u(i,j))/dx + (v(i,j+1)-v(i,j))/dy = D_{i,j}
  !                  en restant strictement dans le halo unique (0:nx+1, 0:ny+1).
  !=============================================================================
  pure subroutine interpolate_velocity_cvi(mesh, u, v, xp, yp, scheme, up, vp)
    type(mesh_2d_t),  intent(in)  :: mesh
    real(dp),         intent(in)  :: u(0:mesh%nx_local+1, 0:mesh%ny_local+1)
    real(dp),         intent(in)  :: v(0:mesh%nx_local+1, 0:mesh%ny_local+1)
    real(dp),         intent(in)  :: xp, yp
    character(len=*), intent(in)  :: scheme
    real(dp),         intent(out) :: up, vp

    integer  :: i, j, i_v, j_u, nx, ny
    real(dp) :: dx, dy, x_min_loc, y_min_loc
    real(dp) :: xi_u, eta_u, xi_v, eta_v
    real(dp) :: u_sw, u_se, u_nw, u_ne
    real(dp) :: v_sw, v_se, v_nw, v_ne
    real(dp) :: du_s, du_n, dv_w, dv_e
    real(dp) :: u_bilin, v_bilin, u_corr, v_corr

    if (.not. (xp == xp .and. abs(xp) < 1.0e30_dp .and. yp == yp .and. abs(yp) < 1.0e30_dp)) then
       up = 0.0_dp
       vp = 0.0_dp
       return
    end if

    nx = mesh%nx_local
    ny = mesh%ny_local
    dx = mesh%dx
    dy = mesh%dy

    x_min_loc = x_face(mesh, 1)
    y_min_loc = y_face(mesh, 1)

    i = floor((xp - x_min_loc) / dx) + 1
    j = floor((yp - y_min_loc) / dy) + 1
    i = max(1, min(nx, i))
    j = max(1, min(ny, j))

    ! Indices des boîtes décalées MAC (j_u in {j-1, j} subset [0, ny], i_v in {i-1, i} subset [0, nx])
    if (yp < y_center(mesh, j)) then
       j_u = j - 1
    else
       j_u = j
    end if
    if (xp < x_center(mesh, i)) then
       i_v = i - 1
    else
       i_v = i
    end if

    ! Coordonnées réduites locales dans chacune des boîtes décalées
    xi_u  = max(0.0_dp, min(1.0_dp, (xp - x_face(mesh, i))     / dx))
    eta_u = max(0.0_dp, min(1.0_dp, (yp - y_center(mesh, j_u)) / dy))
    xi_v  = max(0.0_dp, min(1.0_dp, (xp - x_center(mesh, i_v)) / dx))
    eta_v = max(0.0_dp, min(1.0_dp, (yp - y_face(mesh, j))     / dy))

    ! Sommets MAC de la boîte décalée u : [x_i^f, x_{i+1}^f] x [y_{j_u}^c, y_{j_u+1}^c]
    u_sw = u(i,   j_u)
    u_se = u(i+1, j_u)
    u_nw = u(i,   j_u+1)
    u_ne = u(i+1, j_u+1)

    ! Sommets MAC de la boîte décalée v : [x_{i_v}^c, x_{i_v+1}^c] x [y_j^f, y_{j+1}^f]
    v_sw = v(i_v,   j)
    v_se = v(i_v+1, j)
    v_nw = v(i_v,   j+1)
    v_ne = v(i_v+1, j+1)

    du_s = u_se - u_sw
    du_n = u_ne - u_nw
    dv_w = v_nw - v_sw
    dv_e = v_ne - v_se

    ! Interpolation bilinéaire sur chaque boîte décalée MAC
    u_bilin = (1.0_dp - eta_u) * ((1.0_dp - xi_u) * u_sw + xi_u * u_se) + &
              eta_u            * ((1.0_dp - xi_u) * u_nw + xi_u * u_ne)
    v_bilin = (1.0_dp - xi_v)  * ((1.0_dp - eta_v) * v_sw + eta_v * v_nw) + &
              xi_v             * ((1.0_dp - eta_v) * v_se + eta_v * v_ne)

    if (trim(scheme) == "BILINEAR_RK2" .or. trim(scheme) == "BILINEAR") then
       up = u_bilin
       vp = v_bilin
    else
       ! Correction quadratique CVI exacte (Variante B — Wang et al., 2015 ; Pusok et al., 2017) :
       ! d/dx[u_corr] = -(xi_u - 1/2)*(dv_e - dv_w)/dy
       ! d/dy[v_corr] = -(eta_v - 1/2)*(du_n - du_s)/dx
       u_corr = 0.5_dp * (dx / dy) * xi_u  * (1.0_dp - xi_u)  * (dv_e - dv_w)
       v_corr = 0.5_dp * (dy / dx) * eta_v * (1.0_dp - eta_v) * (du_n - du_s)
       up = u_bilin + u_corr
       vp = v_bilin + v_corr
    end if
  end subroutine interpolate_velocity_cvi

  !=============================================================================
  ! Routine : advect_markers_rk2
  ! Rôle    : Transporte les marqueurs par le schéma prédicteur-correcteur RK2
  !           (Heun d'ordre 2) avec la vitesse interpolée CVI
  !=============================================================================
  subroutine advect_markers_rk2(cloud, mesh, cfg, fields, dt)
    type(marker_cloud_t), intent(inout) :: cloud
    type(mesh_2d_t),      intent(in)    :: mesh
    type(config_t),       intent(in)    :: cfg
    type(flow_fields_t),  intent(in)    :: fields
    real(dp),             intent(in)    :: dt

    integer  :: p, k, nx, ny, n_rem, n_mid, n_eval, n_rep, n_fin, ierr
    integer  :: n_sL, n_sR, n_rL, n_rR, n_sB, n_sT, n_rB, n_rT
    integer  :: status(MPI_STATUS_SIZE)
    real(dp) :: x0, y0, x_star, y_star, u1, v1, u2, v2, xs, ys
    real(dp) :: x_min_loc, x_max_loc, y_min_loc, y_max_loc
    real(dp), parameter :: eps_b = 1.0e-12_dp
    logical  :: per_x, per_y, out_x, out_y

    integer,  allocatable :: rq_idx(:)
    real(dp), allocatable :: rq_x0(:), rq_y0(:), rq_u1(:), rq_v1(:)
    integer,  allocatable :: mid_qid(:), mid_dx(:), eval_qid(:), eval_dx(:), eval_dy(:)
    real(dp), allocatable :: mid_x(:), mid_y(:), eval_x(:), eval_y(:), eval_u2(:), eval_v2(:)
    integer,  allocatable :: rep_qid(:), rep_dx(:), fin_qid(:)
    real(dp), allocatable :: rep_u2(:), rep_v2(:), fin_u2(:), fin_v2(:)
    real(dp), allocatable :: send_w(:), send_e(:), recv_w(:), recv_e(:)
    real(dp), allocatable :: send_s(:), send_n(:), recv_s(:), recv_n(:)
    real(dp) :: cnt_s(1), cnt_r(1)

    if (.not. cfg%use_markers) return

    nx = mesh%nx_local
    ny = mesh%ny_local
    x_min_loc = x_face(mesh, 1)
    x_max_loc = x_face(mesh, nx + 1)
    y_min_loc = y_face(mesh, 1)
    y_max_loc = y_face(mesh, ny + 1)
    per_x = (cfg%bc_left%btype == BC_PERIODIC .or. cfg%bc_right%btype == BC_PERIODIC)
    per_y = (cfg%bc_bottom%btype == BC_PERIODIC .or. cfg%bc_top%btype == BC_PERIODIC)

    k = max(1, cloud%count)
    allocate(rq_idx(k), rq_x0(k), rq_y0(k), rq_u1(k), rq_v1(k))
    allocate(mid_qid(k), mid_dx(k), mid_x(k), mid_y(k))
    n_rem = 0
    n_mid = 0

    do p = 1, cloud%count
       if (.not. cloud%list(p)%active) cycle
       x0 = cloud%list(p)%x
       y0 = cloud%list(p)%y

       call interpolate_velocity_cvi(mesh, fields%u, fields%v, x0, y0, &
                                     cfg%marker_advection_scheme, u1, v1)
       if (.not. (u1 == u1 .and. abs(u1) < 1.0e30_dp .and. v1 == v1 .and. abs(v1) < 1.0e30_dp)) then
          u1 = 0.0_dp
          v1 = 0.0_dp
       end if

       x_star = x0 + dt * u1
       y_star = y0 + dt * v1

       if (per_x .and. mesh%px == 1) then
          x_star = x_star - mesh%lx * floor(x_star / mesh%lx)
          if (x_star >= mesh%lx) x_star = 0.0_dp
       else if (.not. per_x) then
          if (mesh%is_boundary_left  .and. x_star < 0.0_dp)   x_star = eps_b * mesh%dx
          if (mesh%is_boundary_right .and. x_star >= mesh%lx) x_star = mesh%lx - eps_b * mesh%dx
       end if

       if (per_y .and. mesh%py == 1) then
          y_star = y_star - mesh%ly * floor(y_star / mesh%ly)
          if (y_star >= mesh%ly) y_star = 0.0_dp
       else if (.not. per_y) then
          if (mesh%is_boundary_bottom .and. y_star < 0.0_dp)   y_star = eps_b * mesh%dy
          if (mesh%is_boundary_top    .and. y_star >= mesh%ly) y_star = mesh%ly - eps_b * mesh%dy
       end if

       out_x = (mesh%px > 1) .and. ((x_star < x_min_loc .and. mesh%neighbor_west /= MPI_PROC_NULL) .or. &
                                    (x_star >= x_max_loc .and. mesh%neighbor_east /= MPI_PROC_NULL))
       out_y = (mesh%py > 1) .and. ((y_star < y_min_loc .and. mesh%neighbor_south /= MPI_PROC_NULL) .or. &
                                    (y_star >= y_max_loc .and. mesh%neighbor_north /= MPI_PROC_NULL))

       if (mesh%size <= 1 .or. (.not. out_x .and. .not. out_y)) then
          call interpolate_velocity_cvi(mesh, fields%u, fields%v, x_star, y_star, &
                                        cfg%marker_advection_scheme, u2, v2)
          if (.not. (u2 == u2 .and. abs(u2) < 1.0e30_dp .and. v2 == v2 .and. abs(v2) < 1.0e30_dp)) then
             u2 = 0.0_dp
             v2 = 0.0_dp
          end if
          cloud%list(p)%x = x0 + 0.5_dp * dt * (u1 + u2)
          cloud%list(p)%y = y0 + 0.5_dp * dt * (v1 + v2)
       else
          n_rem = n_rem + 1
          rq_idx(n_rem) = p
          rq_x0(n_rem)  = x0
          rq_y0(n_rem)  = y0
          rq_u1(n_rem)  = u1
          rq_v1(n_rem)  = v1

          n_mid = n_mid + 1
          mid_qid(n_mid) = n_rem
          mid_dx(n_mid)  = 0
          mid_x(n_mid)   = x_star
          mid_y(n_mid)   = y_star
       end if
    end do

    if (mesh%size > 1) then
       ! --- PASSE X ALLER ---
       if (mesh%px > 1) then
          allocate(send_w(3 * max(1, n_mid)), send_e(3 * max(1, n_mid)))
          n_sL = 0; n_sR = 0
          k = n_mid
          n_mid = 0
          do p = 1, k
             if (mid_x(p) < x_min_loc .and. mesh%neighbor_west /= MPI_PROC_NULL) then
                xs = mid_x(p) + merge(mesh%lx, 0.0_dp, per_x .and. mesh%coord_x == 0)
                n_sL = n_sL + 1
                send_w(3*(n_sL-1) + 1) = real(mid_qid(p), dp)
                send_w(3*(n_sL-1) + 2) = xs
                send_w(3*(n_sL-1) + 3) = mid_y(p)
             else if (mid_x(p) >= x_max_loc .and. mesh%neighbor_east /= MPI_PROC_NULL) then
                xs = mid_x(p) - merge(mesh%lx, 0.0_dp, per_x .and. mesh%coord_x == mesh%px - 1)
                n_sR = n_sR + 1
                send_e(3*(n_sR-1) + 1) = real(mid_qid(p), dp)
                send_e(3*(n_sR-1) + 2) = xs
                send_e(3*(n_sR-1) + 3) = mid_y(p)
             else
                n_mid = n_mid + 1
                mid_qid(n_mid) = mid_qid(p)
                mid_dx(n_mid)  = mid_dx(p)
                mid_x(n_mid)   = mid_x(p)
                mid_y(n_mid)   = mid_y(p)
             end if
          end do

          cnt_s(1) = real(n_sL, dp); cnt_r(1) = 0.0_dp
          call MPI_Sendrecv(cnt_s, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 601, &
                            cnt_r, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 601, mesh%cart_comm, status, ierr)
          n_rR = nint(cnt_r(1))
          cnt_s(1) = real(n_sR, dp); cnt_r(1) = 0.0_dp
          call MPI_Sendrecv(cnt_s, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 602, &
                            cnt_r, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 602, mesh%cart_comm, status, ierr)
          n_rL = nint(cnt_r(1))

          allocate(recv_e(3 * max(1, n_rR)), recv_w(3 * max(1, n_rL)))
          call MPI_Sendrecv(send_w, 3 * n_sL, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 603, &
                            recv_e, 3 * n_rR, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 603, mesh%cart_comm, status, ierr)
          call MPI_Sendrecv(send_e, 3 * n_sR, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 604, &
                            recv_w, 3 * n_rL, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 604, mesh%cart_comm, status, ierr)

          if (n_mid + n_rR + n_rL > size(mid_qid)) then
             block
               integer,  allocatable :: tq(:), tdx(:)
               real(dp), allocatable :: tx(:), ty(:)
               integer :: new_sz
               new_sz = n_mid + n_rR + n_rL + 16
               allocate(tq(new_sz), tdx(new_sz), tx(new_sz), ty(new_sz))
               if (n_mid > 0) then
                  tq(1:n_mid) = mid_qid(1:n_mid); tdx(1:n_mid) = mid_dx(1:n_mid)
                  tx(1:n_mid) = mid_x(1:n_mid);   ty(1:n_mid)  = mid_y(1:n_mid)
               end if
               call move_alloc(tq, mid_qid); call move_alloc(tdx, mid_dx)
               call move_alloc(tx, mid_x);   call move_alloc(ty, mid_y)
             end block
          end if

          do p = 0, n_rR - 1
             n_mid = n_mid + 1
             mid_qid(n_mid) = nint(recv_e(3*p + 1))
             mid_dx(n_mid)  = -1
             mid_x(n_mid)   = recv_e(3*p + 2)
             mid_y(n_mid)   = recv_e(3*p + 3)
          end do
          do p = 0, n_rL - 1
             n_mid = n_mid + 1
             mid_qid(n_mid) = nint(recv_w(3*p + 1))
             mid_dx(n_mid)  = +1
             mid_x(n_mid)   = recv_w(3*p + 2)
             mid_y(n_mid)   = recv_w(3*p + 3)
          end do
          deallocate(send_w, send_e, recv_w, recv_e)
       end if

       ! --- PASSE Y ALLER ---
       allocate(eval_qid(max(1, n_mid + 16)), eval_dx(max(1, n_mid + 16)), eval_dy(max(1, n_mid + 16)))
       allocate(eval_x(max(1, n_mid + 16)), eval_y(max(1, n_mid + 16)), eval_u2(max(1, n_mid + 16)), eval_v2(max(1, n_mid + 16)))
       n_eval = 0

       if (mesh%py > 1) then
          allocate(send_s(4 * max(1, n_mid)), send_n(4 * max(1, n_mid)))
          n_sB = 0; n_sT = 0
          do p = 1, n_mid
             if (mid_y(p) < y_min_loc .and. mesh%neighbor_south /= MPI_PROC_NULL) then
                ys = mid_y(p) + merge(mesh%ly, 0.0_dp, per_y .and. mesh%coord_y == 0)
                n_sB = n_sB + 1
                send_s(4*(n_sB-1) + 1) = real(mid_qid(p), dp)
                send_s(4*(n_sB-1) + 2) = real(mid_dx(p), dp)
                send_s(4*(n_sB-1) + 3) = mid_x(p)
                send_s(4*(n_sB-1) + 4) = ys
             else if (mid_y(p) >= y_max_loc .and. mesh%neighbor_north /= MPI_PROC_NULL) then
                ys = mid_y(p) - merge(mesh%ly, 0.0_dp, per_y .and. mesh%coord_y == mesh%py - 1)
                n_sT = n_sT + 1
                send_n(4*(n_sT-1) + 1) = real(mid_qid(p), dp)
                send_n(4*(n_sT-1) + 2) = real(mid_dx(p), dp)
                send_n(4*(n_sT-1) + 3) = mid_x(p)
                send_n(4*(n_sT-1) + 4) = ys
             else
                n_eval = n_eval + 1
                eval_qid(n_eval) = mid_qid(p)
                eval_dx(n_eval)  = mid_dx(p)
                eval_dy(n_eval)  = 0
                eval_x(n_eval)   = mid_x(p)
                eval_y(n_eval)   = mid_y(p)
             end if
          end do

          cnt_s(1) = real(n_sB, dp); cnt_r(1) = 0.0_dp
          call MPI_Sendrecv(cnt_s, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 605, &
                            cnt_r, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 605, mesh%cart_comm, status, ierr)
          n_rT = nint(cnt_r(1))
          cnt_s(1) = real(n_sT, dp); cnt_r(1) = 0.0_dp
          call MPI_Sendrecv(cnt_s, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 606, &
                            cnt_r, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 606, mesh%cart_comm, status, ierr)
          n_rB = nint(cnt_r(1))

          allocate(recv_n(4 * max(1, n_rT)), recv_s(4 * max(1, n_rB)))
          call MPI_Sendrecv(send_s, 4 * n_sB, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 607, &
                            recv_n, 4 * n_rT, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 607, mesh%cart_comm, status, ierr)
          call MPI_Sendrecv(send_n, 4 * n_sT, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 608, &
                            recv_s, 4 * n_rB, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 608, mesh%cart_comm, status, ierr)

          if (n_eval + n_rT + n_rB > size(eval_qid)) then
             block
               integer,  allocatable :: tq(:), tdx(:), tdy(:)
               real(dp), allocatable :: tx(:), ty(:), tu(:), tv(:)
               integer :: new_sz
               new_sz = n_eval + n_rT + n_rB + 16
               allocate(tq(new_sz), tdx(new_sz), tdy(new_sz), tx(new_sz), ty(new_sz), tu(new_sz), tv(new_sz))
               if (n_eval > 0) then
                  tq(1:n_eval) = eval_qid(1:n_eval); tdx(1:n_eval) = eval_dx(1:n_eval); tdy(1:n_eval) = eval_dy(1:n_eval)
                  tx(1:n_eval) = eval_x(1:n_eval);   ty(1:n_eval)  = eval_y(1:n_eval)
               end if
               call move_alloc(tq, eval_qid); call move_alloc(tdx, eval_dx); call move_alloc(tdy, eval_dy)
               call move_alloc(tx, eval_x);   call move_alloc(ty, eval_y)
               call move_alloc(tu, eval_u2);  call move_alloc(tv, eval_v2)
             end block
          end if

          do p = 0, n_rT - 1
             n_eval = n_eval + 1
             eval_qid(n_eval) = nint(recv_n(4*p + 1))
             eval_dx(n_eval)  = nint(recv_n(4*p + 2))
             eval_dy(n_eval)  = -1
             eval_x(n_eval)   = recv_n(4*p + 3)
             eval_y(n_eval)   = recv_n(4*p + 4)
          end do
          do p = 0, n_rB - 1
             n_eval = n_eval + 1
             eval_qid(n_eval) = nint(recv_s(4*p + 1))
             eval_dx(n_eval)  = nint(recv_s(4*p + 2))
             eval_dy(n_eval)  = +1
             eval_x(n_eval)   = recv_s(4*p + 3)
             eval_y(n_eval)   = recv_s(4*p + 4)
          end do
          deallocate(send_s, send_n, recv_s, recv_n)
       else
          n_eval = n_mid
          if (n_eval > 0) then
             eval_qid(1:n_eval) = mid_qid(1:n_eval)
             eval_dx(1:n_eval)  = mid_dx(1:n_eval)
             eval_dy(1:n_eval)  = 0
             eval_x(1:n_eval)   = mid_x(1:n_eval)
             eval_y(1:n_eval)   = mid_y(1:n_eval)
          end if
       end if

       do p = 1, n_eval
          call interpolate_velocity_cvi(mesh, fields%u, fields%v, eval_x(p), eval_y(p), &
                                        cfg%marker_advection_scheme, eval_u2(p), eval_v2(p))
       end do

       ! --- PASSE Y RETOUR ---
       allocate(rep_qid(max(1, n_eval + 16)), rep_dx(max(1, n_eval + 16)), rep_u2(max(1, n_eval + 16)), rep_v2(max(1, n_eval + 16)))
       n_rep = 0
       if (mesh%py > 1) then
          allocate(send_s(4 * max(1, n_eval)), send_n(4 * max(1, n_eval)))
          n_sB = 0; n_sT = 0
          do p = 1, n_eval
             if (eval_dy(p) == +1) then
                n_sB = n_sB + 1
                send_s(4*(n_sB-1) + 1) = real(eval_qid(p), dp)
                send_s(4*(n_sB-1) + 2) = real(eval_dx(p), dp)
                send_s(4*(n_sB-1) + 3) = eval_u2(p)
                send_s(4*(n_sB-1) + 4) = eval_v2(p)
             else if (eval_dy(p) == -1) then
                n_sT = n_sT + 1
                send_n(4*(n_sT-1) + 1) = real(eval_qid(p), dp)
                send_n(4*(n_sT-1) + 2) = real(eval_dx(p), dp)
                send_n(4*(n_sT-1) + 3) = eval_u2(p)
                send_n(4*(n_sT-1) + 4) = eval_v2(p)
             else
                n_rep = n_rep + 1
                rep_qid(n_rep) = eval_qid(p); rep_dx(n_rep) = eval_dx(p)
                rep_u2(n_rep)  = eval_u2(p);  rep_v2(n_rep) = eval_v2(p)
             end if
          end do
          cnt_s(1) = real(n_sB, dp); cnt_r(1) = 0.0_dp
          call MPI_Sendrecv(cnt_s, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 615, &
                            cnt_r, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 615, mesh%cart_comm, status, ierr)
          n_rT = nint(cnt_r(1))
          cnt_s(1) = real(n_sT, dp); cnt_r(1) = 0.0_dp
          call MPI_Sendrecv(cnt_s, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 616, &
                            cnt_r, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 616, mesh%cart_comm, status, ierr)
          n_rB = nint(cnt_r(1))
          allocate(recv_n(4 * max(1, n_rT)), recv_s(4 * max(1, n_rB)))
          call MPI_Sendrecv(send_s, 4 * n_sB, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 617, &
                            recv_n, 4 * n_rT, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 617, mesh%cart_comm, status, ierr)
          call MPI_Sendrecv(send_n, 4 * n_sT, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 618, &
                            recv_s, 4 * n_rB, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 618, mesh%cart_comm, status, ierr)
          if (n_rep + n_rT + n_rB > size(rep_qid)) then
             block
               integer,  allocatable :: tq(:), tdx(:)
               real(dp), allocatable :: tu(:), tv(:)
               integer :: new_sz
               new_sz = n_rep + n_rT + n_rB + 16
               allocate(tq(new_sz), tdx(new_sz), tu(new_sz), tv(new_sz))
               if (n_rep > 0) then
                  tq(1:n_rep) = rep_qid(1:n_rep); tdx(1:n_rep) = rep_dx(1:n_rep)
                  tu(1:n_rep) = rep_u2(1:n_rep);  tv(1:n_rep)  = rep_v2(1:n_rep)
               end if
               call move_alloc(tq, rep_qid); call move_alloc(tdx, rep_dx)
               call move_alloc(tu, rep_u2);  call move_alloc(tv, rep_v2)
             end block
          end if
          do p = 0, n_rT - 1
             n_rep = n_rep + 1
             rep_qid(n_rep) = nint(recv_n(4*p + 1)); rep_dx(n_rep) = nint(recv_n(4*p + 2))
             rep_u2(n_rep)  = recv_n(4*p + 3);        rep_v2(n_rep) = recv_n(4*p + 4)
          end do
          do p = 0, n_rB - 1
             n_rep = n_rep + 1
             rep_qid(n_rep) = nint(recv_s(4*p + 1)); rep_dx(n_rep) = nint(recv_s(4*p + 2))
             rep_u2(n_rep)  = recv_s(4*p + 3);        rep_v2(n_rep) = recv_s(4*p + 4)
          end do
          deallocate(send_s, send_n, recv_s, recv_n)
       else
          n_rep = n_eval
          if (n_rep > 0) then
             rep_qid(1:n_rep) = eval_qid(1:n_rep); rep_dx(1:n_rep) = eval_dx(1:n_rep)
             rep_u2(1:n_rep)  = eval_u2(1:n_rep);  rep_v2(1:n_rep) = eval_v2(1:n_rep)
          end if
       end if

       ! --- PASSE X RETOUR ---
       allocate(fin_qid(max(1, n_rep + 16)), fin_u2(max(1, n_rep + 16)), fin_v2(max(1, n_rep + 16)))
       n_fin = 0
       if (mesh%px > 1) then
          allocate(send_w(3 * max(1, n_rep)), send_e(3 * max(1, n_rep)))
          n_sL = 0; n_sR = 0
          do p = 1, n_rep
             if (rep_dx(p) == +1) then
                n_sL = n_sL + 1
                send_w(3*(n_sL-1) + 1) = real(rep_qid(p), dp)
                send_w(3*(n_sL-1) + 2) = rep_u2(p)
                send_w(3*(n_sL-1) + 3) = rep_v2(p)
             else if (rep_dx(p) == -1) then
                n_sR = n_sR + 1
                send_e(3*(n_sR-1) + 1) = real(rep_qid(p), dp)
                send_e(3*(n_sR-1) + 2) = rep_u2(p)
                send_e(3*(n_sR-1) + 3) = rep_v2(p)
             else
                n_fin = n_fin + 1
                fin_qid(n_fin) = rep_qid(p); fin_u2(n_fin) = rep_u2(p); fin_v2(n_fin) = rep_v2(p)
             end if
          end do
          cnt_s(1) = real(n_sL, dp); cnt_r(1) = 0.0_dp
          call MPI_Sendrecv(cnt_s, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 611, &
                            cnt_r, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 611, mesh%cart_comm, status, ierr)
          n_rR = nint(cnt_r(1))
          cnt_s(1) = real(n_sR, dp); cnt_r(1) = 0.0_dp
          call MPI_Sendrecv(cnt_s, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 612, &
                            cnt_r, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 612, mesh%cart_comm, status, ierr)
          n_rL = nint(cnt_r(1))
          allocate(recv_e(3 * max(1, n_rR)), recv_w(3 * max(1, n_rL)))
          call MPI_Sendrecv(send_w, 3 * n_sL, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 613, &
                            recv_e, 3 * n_rR, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 613, mesh%cart_comm, status, ierr)
          call MPI_Sendrecv(send_e, 3 * n_sR, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 614, &
                            recv_w, 3 * n_rL, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 614, mesh%cart_comm, status, ierr)
          if (n_fin + n_rR + n_rL > size(fin_qid)) then
             block
               integer,  allocatable :: tq(:)
               real(dp), allocatable :: tu(:), tv(:)
               integer :: new_sz
               new_sz = n_fin + n_rR + n_rL + 16
               allocate(tq(new_sz), tu(new_sz), tv(new_sz))
               if (n_fin > 0) then
                  tq(1:n_fin) = fin_qid(1:n_fin); tu(1:n_fin) = fin_u2(1:n_fin); tv(1:n_fin) = fin_v2(1:n_fin)
               end if
               call move_alloc(tq, fin_qid); call move_alloc(tu, fin_u2); call move_alloc(tv, fin_v2)
             end block
          end if
          do p = 0, n_rR - 1
             n_fin = n_fin + 1
             fin_qid(n_fin) = nint(recv_e(3*p + 1)); fin_u2(n_fin) = recv_e(3*p + 2); fin_v2(n_fin) = recv_e(3*p + 3)
          end do
          do p = 0, n_rL - 1
             n_fin = n_fin + 1
             fin_qid(n_fin) = nint(recv_w(3*p + 1)); fin_u2(n_fin) = recv_w(3*p + 2); fin_v2(n_fin) = recv_w(3*p + 3)
          end do
          deallocate(send_w, send_e, recv_w, recv_e)
       else
          n_fin = n_rep
          if (n_fin > 0) then
             fin_qid(1:n_fin) = rep_qid(1:n_fin); fin_u2(1:n_fin) = rep_u2(1:n_fin); fin_v2(1:n_fin) = rep_v2(1:n_fin)
          end if
       end if

       do p = 1, n_fin
          k = fin_qid(p)
          if (k >= 1 .and. k <= n_rem) then
             u2 = fin_u2(p); v2 = fin_v2(p)
             if (.not. (u2 == u2 .and. abs(u2) < 1.0e30_dp .and. v2 == v2 .and. abs(v2) < 1.0e30_dp)) then
                u2 = 0.0_dp; v2 = 0.0_dp
             end if
             cloud%list(rq_idx(k))%x = rq_x0(k) + 0.5_dp * dt * (rq_u1(k) + u2)
             cloud%list(rq_idx(k))%y = rq_y0(k) + 0.5_dp * dt * (rq_v1(k) + v2)
          end if
       end do
       deallocate(eval_qid, eval_dx, eval_dy, eval_x, eval_y, eval_u2, eval_v2)
       deallocate(rep_qid, rep_dx, rep_u2, rep_v2, fin_qid, fin_u2, fin_v2)
    end if
    deallocate(rq_idx, rq_x0, rq_y0, rq_u1, rq_v1, mid_qid, mid_dx, mid_x, mid_y)
  end subroutine advect_markers_rk2

  !=============================================================================
  ! Routine : migrate_markers_mpi
  ! Rôle    : Gère les conditions aux limites périodiques/murs et migre les
  !           marqueurs sortants entre sous-domaines MPI en 2 passes (X puis Y).
  !=============================================================================
  subroutine migrate_markers_mpi(cloud, mesh, cfg)
    type(marker_cloud_t), intent(inout) :: cloud
    type(mesh_2d_t),      intent(in)    :: mesh
    type(config_t),       intent(in)    :: cfg

    integer  :: p, k, nx, ny, i_c, j_c, ierr
    real(dp) :: x_min_loc, x_max_loc, y_min_loc, y_max_loc, xp, yp, xs, ys
    real(dp), parameter :: eps_b = 1.0e-12_dp
    logical  :: per_x, per_y
    integer, parameter :: N_ATTR = 5 ! (x, y, nx, ny, kappa)

    real(dp), allocatable :: send_w(:), send_e(:), recv_w(:), recv_e(:)
    real(dp), allocatable :: send_s(:), send_n(:), recv_s(:), recv_n(:)
    real(dp) :: cnt_s(1), cnt_r(1)
    integer :: nw_out, ne_out, nw_in, ne_in
    integer :: ns_out, nn_out, ns_in, nn_in
    integer :: status(MPI_STATUS_SIZE)

    if (.not. cfg%use_markers) return

    nx = mesh%nx_local
    ny = mesh%ny_local
    x_min_loc = x_face(mesh, 1)
    x_max_loc = x_face(mesh, nx + 1)
    y_min_loc = y_face(mesh, 1)
    y_max_loc = y_face(mesh, ny + 1)

    per_x = (cfg%bc_left%btype == BC_PERIODIC .or. cfg%bc_right%btype == BC_PERIODIC)
    per_y = (cfg%bc_bottom%btype == BC_PERIODIC .or. cfg%bc_top%btype == BC_PERIODIC)

    ! 1. Bornage sur parois physiques non-périodiques, et repliement local uniquement si px/py == 1
    do p = 1, cloud%count
       if (.not. cloud%list(p)%active) cycle
       xp = cloud%list(p)%x
       yp = cloud%list(p)%y

       if (per_x) then
          if (mesh%px == 1) then
             xp = xp - mesh%lx * floor(xp / mesh%lx)
             if (xp >= mesh%lx) xp = 0.0_dp
          end if
       else
          if (mesh%is_boundary_left  .and. xp < 0.0_dp)   xp = eps_b * mesh%dx
          if (mesh%is_boundary_right .and. xp >= mesh%lx) xp = mesh%lx - eps_b * mesh%dx
       end if

       if (per_y) then
          if (mesh%py == 1) then
             yp = yp - mesh%ly * floor(yp / mesh%ly)
             if (yp >= mesh%ly) yp = 0.0_dp
          end if
       else
          if (mesh%is_boundary_bottom .and. yp < 0.0_dp)   yp = eps_b * mesh%dy
          if (mesh%is_boundary_top    .and. yp >= mesh%ly) yp = mesh%ly - eps_b * mesh%dy
       end if

       cloud%list(p)%x = xp
       cloud%list(p)%y = yp
    end do

    ! 2. Migration MPI inter-processus si mesh%size > 1
    if (mesh%size > 1) then
       ! --- PASSE 1 : Axe X (Ouest / Est) ---
       if (mesh%px > 1) then
          allocate(send_w(N_ATTR * max(1, cloud%count)))
          allocate(send_e(N_ATTR * max(1, cloud%count)))
          nw_out = 0; ne_out = 0

          do p = 1, cloud%count
             if (.not. cloud%list(p)%active) cycle
             xp = cloud%list(p)%x
             if (xp < x_min_loc .and. mesh%neighbor_west /= MPI_PROC_NULL) then
                xs = xp + merge(mesh%lx, 0.0_dp, per_x .and. mesh%coord_x == 0)
                send_w(nw_out*N_ATTR + 1) = xs
                send_w(nw_out*N_ATTR + 2) = cloud%list(p)%y
                send_w(nw_out*N_ATTR + 3) = cloud%list(p)%nx
                send_w(nw_out*N_ATTR + 4) = cloud%list(p)%ny
                send_w(nw_out*N_ATTR + 5) = cloud%list(p)%kappa
                nw_out = nw_out + 1
                cloud%list(p)%active = .false.
             else if (xp >= x_max_loc .and. mesh%neighbor_east /= MPI_PROC_NULL) then
                xs = xp - merge(mesh%lx, 0.0_dp, per_x .and. mesh%coord_x == mesh%px - 1)
                send_e(ne_out*N_ATTR + 1) = xs
                send_e(ne_out*N_ATTR + 2) = cloud%list(p)%y
                send_e(ne_out*N_ATTR + 3) = cloud%list(p)%nx
                send_e(ne_out*N_ATTR + 4) = cloud%list(p)%ny
                send_e(ne_out*N_ATTR + 5) = cloud%list(p)%kappa
                ne_out = ne_out + 1
                cloud%list(p)%active = .false.
             end if
          end do

          nw_in = 0; ne_in = 0
          cnt_s(1) = real(nw_out, dp); cnt_r(1) = 0.0_dp
          call MPI_Sendrecv(cnt_s, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 401, &
                            cnt_r, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 401, &
                            mesh%cart_comm, status, ierr)
          ne_in = nint(cnt_r(1))
          cnt_s(1) = real(ne_out, dp); cnt_r(1) = 0.0_dp
          call MPI_Sendrecv(cnt_s, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 402, &
                            cnt_r, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 402, &
                            mesh%cart_comm, status, ierr)
          nw_in = nint(cnt_r(1))

          allocate(recv_e(N_ATTR * max(1, ne_in)))
          allocate(recv_w(N_ATTR * max(1, nw_in)))

          call MPI_Sendrecv(send_w, nw_out * N_ATTR, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 403, &
                            recv_e, ne_in  * N_ATTR, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 403, &
                            mesh%cart_comm, status, ierr)
          call MPI_Sendrecv(send_e, ne_out * N_ATTR, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 404, &
                            recv_w, nw_in  * N_ATTR, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 404, &
                            mesh%cart_comm, status, ierr)

          call ensure_capacity(cloud, cloud%count + ne_in + nw_in + 16)
          do k = 0, ne_in - 1
             cloud%count = cloud%count + 1
             cloud%list(cloud%count)%x      = recv_e(k*N_ATTR + 1)
             cloud%list(cloud%count)%y      = recv_e(k*N_ATTR + 2)
             cloud%list(cloud%count)%nx     = recv_e(k*N_ATTR + 3)
             cloud%list(cloud%count)%ny     = recv_e(k*N_ATTR + 4)
             cloud%list(cloud%count)%kappa  = recv_e(k*N_ATTR + 5)
             cloud%list(cloud%count)%active = .true.
          end do
          do k = 0, nw_in - 1
             cloud%count = cloud%count + 1
             cloud%list(cloud%count)%x      = recv_w(k*N_ATTR + 1)
             cloud%list(cloud%count)%y      = recv_w(k*N_ATTR + 2)
             cloud%list(cloud%count)%nx     = recv_w(k*N_ATTR + 3)
             cloud%list(cloud%count)%ny     = recv_w(k*N_ATTR + 4)
             cloud%list(cloud%count)%kappa  = recv_w(k*N_ATTR + 5)
             cloud%list(cloud%count)%active = .true.
          end do

          deallocate(send_w, send_e, recv_w, recv_e)
       end if

       ! --- PASSE 2 : Axe Y (Sud / Nord) ---
       if (mesh%py > 1) then
          allocate(send_s(N_ATTR * max(1, cloud%count)))
          allocate(send_n(N_ATTR * max(1, cloud%count)))
          ns_out = 0; nn_out = 0

          do p = 1, cloud%count
             if (.not. cloud%list(p)%active) cycle
             yp = cloud%list(p)%y
             if (yp < y_min_loc .and. mesh%neighbor_south /= MPI_PROC_NULL) then
                ys = yp + merge(mesh%ly, 0.0_dp, per_y .and. mesh%coord_y == 0)
                send_s(ns_out*N_ATTR + 1) = cloud%list(p)%x
                send_s(ns_out*N_ATTR + 2) = ys
                send_s(ns_out*N_ATTR + 3) = cloud%list(p)%nx
                send_s(ns_out*N_ATTR + 4) = cloud%list(p)%ny
                send_s(ns_out*N_ATTR + 5) = cloud%list(p)%kappa
                ns_out = ns_out + 1
                cloud%list(p)%active = .false.
             else if (yp >= y_max_loc .and. mesh%neighbor_north /= MPI_PROC_NULL) then
                ys = yp - merge(mesh%ly, 0.0_dp, per_y .and. mesh%coord_y == mesh%py - 1)
                send_n(nn_out*N_ATTR + 1) = cloud%list(p)%x
                send_n(nn_out*N_ATTR + 2) = ys
                send_n(nn_out*N_ATTR + 3) = cloud%list(p)%nx
                send_n(nn_out*N_ATTR + 4) = cloud%list(p)%ny
                send_n(nn_out*N_ATTR + 5) = cloud%list(p)%kappa
                nn_out = nn_out + 1
                cloud%list(p)%active = .false.
             end if
          end do

          ns_in = 0; nn_in = 0
          cnt_s(1) = real(ns_out, dp); cnt_r(1) = 0.0_dp
          call MPI_Sendrecv(cnt_s, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 405, &
                            cnt_r, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 405, &
                            mesh%cart_comm, status, ierr)
          nn_in = nint(cnt_r(1))
          cnt_s(1) = real(nn_out, dp); cnt_r(1) = 0.0_dp
          call MPI_Sendrecv(cnt_s, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 406, &
                            cnt_r, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 406, &
                            mesh%cart_comm, status, ierr)
          ns_in = nint(cnt_r(1))

          allocate(recv_n(N_ATTR * max(1, nn_in)))
          allocate(recv_s(N_ATTR * max(1, ns_in)))

          call MPI_Sendrecv(send_s, ns_out * N_ATTR, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 407, &
                            recv_n, nn_in  * N_ATTR, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 407, &
                            mesh%cart_comm, status, ierr)
          call MPI_Sendrecv(send_n, nn_out * N_ATTR, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 408, &
                            recv_s, ns_in  * N_ATTR, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 408, &
                            mesh%cart_comm, status, ierr)

          call ensure_capacity(cloud, cloud%count + nn_in + ns_in + 16)
          do k = 0, nn_in - 1
             cloud%count = cloud%count + 1
             cloud%list(cloud%count)%x      = recv_n(k*N_ATTR + 1)
             cloud%list(cloud%count)%y      = recv_n(k*N_ATTR + 2)
             cloud%list(cloud%count)%nx     = recv_n(k*N_ATTR + 3)
             cloud%list(cloud%count)%ny     = recv_n(k*N_ATTR + 4)
             cloud%list(cloud%count)%kappa  = recv_n(k*N_ATTR + 5)
             cloud%list(cloud%count)%active = .true.
          end do
          do k = 0, ns_in - 1
             cloud%count = cloud%count + 1
             cloud%list(cloud%count)%x      = recv_s(k*N_ATTR + 1)
             cloud%list(cloud%count)%y      = recv_s(k*N_ATTR + 2)
             cloud%list(cloud%count)%nx     = recv_s(k*N_ATTR + 3)
             cloud%list(cloud%count)%ny     = recv_s(k*N_ATTR + 4)
             cloud%list(cloud%count)%kappa  = recv_s(k*N_ATTR + 5)
             cloud%list(cloud%count)%active = .true.
          end do

          deallocate(send_s, send_n, recv_s, recv_n)
       end if
    end if

    ! 3. Mise à jour des indices de cellule locale (i_cell, j_cell)
    do p = 1, cloud%count
       if (.not. cloud%list(p)%active) cycle
       i_c = floor((cloud%list(p)%x - x_min_loc) / mesh%dx) + 1
       j_c = floor((cloud%list(p)%y - y_min_loc) / mesh%dy) + 1
       cloud%list(p)%i_cell = max(1, min(nx, i_c))
       cloud%list(p)%j_cell = max(1, min(ny, j_c))
    end do
  end subroutine migrate_markers_mpi

  !=============================================================================
  ! Routine : prune_and_reseed_markers
  ! Rôle    : - Supprime les marqueurs tombés dans des cellules pures ou orphelins
  !           - Limite la sur-concentration (> marker_max_per_cell)
  !           - Ré-ensemence les mailles d'interface avec FILTRAGE GÉOMÉTRIQUE STRICT :
  !             Interdit l'injection si la cellule ne possède pas une normale PLIC univoque
  !             cohérente avec le front principal (continuité, présence biphasique,
  !             gradient seuil, cohérence angulaire et longueur minimale de segment).
  !=============================================================================
  subroutine prune_and_reseed_markers(cloud, mesh, cfg, fields)
    type(marker_cloud_t), intent(inout) :: cloud
    type(mesh_2d_t),      intent(in)    :: mesh
    type(config_t),       intent(in)    :: cfg
    type(flow_fields_t),  intent(in)    :: fields

    integer  :: p, i, j, k, q, s_idx, di, dj, nx, ny, n_per_cell, max_per_cell, new_count
    real(dp) :: c_val, xmin, xmax, ymin, ymax, x1, y1, x2, y2, s, xc_mid, yc_mid
    logical  :: valid, has_mixed_nb, has_marker_nb, normal_coherent, on_jet_nozzle
    real(dp) :: c_min_nb, c_max_nb, c_nb, xc_jet, half_w_jet, nbx, nby, nb_len
    real(dp) :: c_mm, c_0m, c_pm, c_m0, c_p0, c_mp, c_0p, c_pp
    real(dp) :: mx, my, nrm, snx, sny, alpha, h, min_grad, seg_len, dot_prod
    real(dp) :: q_x(4, 2), q_y(4, 2), q_nx(4, 2), q_ny(4, 2), q_kap(4, 2), ref_nx(4), ref_ny(4)
    integer  :: q_cnt(4, 2)
    real(dp), allocatable :: cnt_halo(:, :), nrm_x_halo(:, :), nrm_y_halo(:, :)

    if (.not. cfg%use_markers) return

    nx = mesh%nx_local
    ny = mesh%ny_local
    n_per_cell = max(2, cfg%markers_per_cell)
    max_per_cell = max(n_per_cell + 2, cfg%marker_max_per_cell)
    h = max(mesh%dx, mesh%dy)
    min_grad = 0.25_dp / h

    if (.not. allocated(cloud%count_in_cell)) then
       allocate(cloud%count_in_cell(0:nx+1, 0:ny+1))
    end if
    cloud%count_in_cell = 0

    ! 1. Élagage des marqueurs hors de la bande d'interface
    do p = 1, cloud%count
       if (.not. cloud%list(p)%active) cycle
       i = cloud%list(p)%i_cell
       j = cloud%list(p)%j_cell
       if (i < 1 .or. i > nx .or. j < 1 .or. j > ny) then
          cloud%list(p)%active = .false.
          cycle
       end if
       c_val = fields%C(i, j)

       if (c_val <= 0.02_dp .or. c_val >= 0.98_dp) then
          has_mixed_nb = .false.
          do dj = -1, 1
             do di = -1, 1
                if (di == 0 .and. dj == 0) cycle
                if (i+di >= 0 .and. i+di <= nx+1 .and. j+dj >= 0 .and. j+dj <= ny+1) then
                   if (fields%C(i+di, j+dj) > 0.02_dp .and. fields%C(i+di, j+dj) < 0.98_dp) then
                      has_mixed_nb = .true.
                      exit
                   end if
                end if
             end do
             if (has_mixed_nb) exit
          end do
          if (.not. has_mixed_nb) then
             cloud%list(p)%active = .false.
             cycle
          end if
       end if

       cloud%count_in_cell(i, j) = cloud%count_in_cell(i, j) + 1
    end do

    ! 1b. Fusion lagrangienne par sous-quadrants 2x2 avec séparation de normales opposées (films minces)
    do j = 1, ny
       ymin = y_face(mesh, j); ymax = y_face(mesh, j + 1)
       yc_mid = y_center(mesh, j)
       do i = 1, nx
          if (cloud%count_in_cell(i, j) <= max_per_cell) cycle
          xmin = x_face(mesh, i); xmax = x_face(mesh, i + 1)
          xc_mid = x_center(mesh, i)
          q_x = 0.0_dp; q_y = 0.0_dp; q_nx = 0.0_dp; q_ny = 0.0_dp; q_kap = 0.0_dp
          ref_nx = 0.0_dp; ref_ny = 0.0_dp; q_cnt = 0

          do p = 1, cloud%count
             if (.not. cloud%list(p)%active) cycle
             if (cloud%list(p)%i_cell == i .and. cloud%list(p)%j_cell == j) then
                q = 1
                if (cloud%list(p)%x >= xc_mid) q = q + 1
                if (cloud%list(p)%y >= yc_mid) q = q + 2
                s_idx = 1
                if (q_cnt(q, 1) == 0) then
                   ref_nx(q) = cloud%list(p)%nx
                   ref_ny(q) = cloud%list(p)%ny
                else if (cloud%list(p)%nx * ref_nx(q) + cloud%list(p)%ny * ref_ny(q) < 0.0_dp) then
                   s_idx = 2
                end if
                q_x(q, s_idx)   = q_x(q, s_idx)   + cloud%list(p)%x
                q_y(q, s_idx)   = q_y(q, s_idx)   + cloud%list(p)%y
                q_nx(q, s_idx)  = q_nx(q, s_idx)  + cloud%list(p)%nx
                q_ny(q, s_idx)  = q_ny(q, s_idx)  + cloud%list(p)%ny
                q_kap(q, s_idx) = q_kap(q, s_idx) + cloud%list(p)%kappa
                q_cnt(q, s_idx) = q_cnt(q, s_idx) + 1
                cloud%list(p)%active = .false.
             end if
          end do

          cloud%count_in_cell(i, j) = 0
          do q = 1, 4
             do s_idx = 1, 2
                if (q_cnt(q, s_idx) > 0) then
                   call ensure_capacity(cloud, cloud%count + 16)
                   cloud%count = cloud%count + 1
                   cloud%list(cloud%count)%x = max(xmin, min(xmax, q_x(q, s_idx) / real(q_cnt(q, s_idx), dp)))
                   cloud%list(cloud%count)%y = max(ymin, min(ymax, q_y(q, s_idx) / real(q_cnt(q, s_idx), dp)))
                   nrm = hypot(q_nx(q, s_idx), q_ny(q, s_idx))
                   if (nrm > 1.0e-12_dp) then
                      cloud%list(cloud%count)%nx = q_nx(q, s_idx) / nrm
                      cloud%list(cloud%count)%ny = q_ny(q, s_idx) / nrm
                   else
                      cloud%list(cloud%count)%nx = 0.0_dp
                      cloud%list(cloud%count)%ny = 1.0_dp
                   end if
                   cloud%list(cloud%count)%kappa  = q_kap(q, s_idx) / real(q_cnt(q, s_idx), dp)
                   cloud%list(cloud%count)%i_cell = i
                   cloud%list(cloud%count)%j_cell = j
                   cloud%list(cloud%count)%active = .true.
                   cloud%count_in_cell(i, j) = cloud%count_in_cell(i, j) + 1
                end if
             end do
          end do
       end do
    end do

    ! Synchronisation MPI / périodique du nombre de marqueurs et des normales moyennes dans les halos
    allocate(cnt_halo(0:nx+1, 0:ny+1), nrm_x_halo(0:nx+1, 0:ny+1), nrm_y_halo(0:nx+1, 0:ny+1))
    cnt_halo = 0.0_dp; nrm_x_halo = 0.0_dp; nrm_y_halo = 0.0_dp
    cnt_halo(1:nx, 1:ny) = real(cloud%count_in_cell(1:nx, 1:ny), dp)
    do p = 1, cloud%count
       if (.not. cloud%list(p)%active) cycle
       i = cloud%list(p)%i_cell; j = cloud%list(p)%j_cell
       if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
          nrm_x_halo(i, j) = nrm_x_halo(i, j) + cloud%list(p)%nx
          nrm_y_halo(i, j) = nrm_y_halo(i, j) + cloud%list(p)%ny
       end if
    end do
    call exchange_halos_scalar(cnt_halo, mesh)
    call exchange_halos_scalar(nrm_x_halo, mesh)
    call exchange_halos_scalar(nrm_y_halo, mesh)

    ! 2. Ré-ensemencement local avec FILTRAGE GÉOMÉTRIQUE STRICT
    xc_jet = merge(cfg%bc_top%jet_center, 0.5_dp * mesh%lx, cfg%bc_top%jet_center > 0.0_dp)
    half_w_jet = merge(0.5_dp * cfg%bc_top%jet_width, 0.04_dp, cfg%bc_top%jet_width > 0.0_dp)

    do j = 1, ny
       ymin = y_face(mesh, j)
       ymax = y_face(mesh, j + 1)
       do i = 1, nx
          c_val = fields%C(i, j)

          ! Condition 2.1 : Cellule mixte franche (seuil resserré 0.05 < C < 0.95) dépourvue de marqueurs
          if (c_val <= 0.05_dp .or. c_val >= 0.95_dp .or. cloud%count_in_cell(i, j) > 0) cycle

          on_jet_nozzle = (trim(cfg%bc_top%profile) == "JET_INLET" .and. mesh%is_boundary_top .and. &
                           j == ny .and. abs(x_center(mesh, i) - xc_jet) <= half_w_jet + mesh%dx)

          ! Condition 2.2 : Continuité avec le front principal existant (ou l'injecteur actif)
          has_marker_nb = on_jet_nozzle
          c_min_nb = 1.0_dp
          c_max_nb = 0.0_dp
          do dj = -1, 1
             do di = -1, 1
                if (di == 0 .and. dj == 0) cycle
                if (i+di >= 0 .and. i+di <= nx+1 .and. j+dj >= 0 .and. j+dj <= ny+1) then
                   if (cnt_halo(i+di, j+dj) > 0.5_dp) has_marker_nb = .true.
                   c_nb = fields%C(i+di, j+dj)
                   if (c_nb < c_min_nb) c_min_nb = c_nb
                   if (c_nb > c_max_nb) c_max_nb = c_nb
                end if
             end do
          end do
          if (.not. has_marker_nb) cycle ! Cellule isolée du front principal : INTERDIT

          ! Condition 2.3 : Présence obligatoire des deux phases dans le voisinage 3x3 (vrai saut d'interface)
          if (c_min_nb > 0.45_dp .or. c_max_nb < 0.55_dp) cycle

          ! Condition 2.4 : Normale univoque par gradient de Youngs avec magnitude suffisante
          c_mm = fields%C(i - 1, j - 1)
          c_0m = fields%C(i,     j - 1)
          c_pm = fields%C(i + 1, j - 1)
          c_m0 = fields%C(i - 1, j    )
          c_p0 = fields%C(i + 1, j    )
          c_mp = fields%C(i - 1, j + 1)
          c_0p = fields%C(i,     j + 1)
          c_pp = fields%C(i + 1, j + 1)

          mx = -((c_pp + 2.0_dp * c_p0 + c_pm) - (c_mp + 2.0_dp * c_m0 + c_mm)) / (8.0_dp * mesh%dx)
          my = -((c_pp + 2.0_dp * c_0p + c_mp) - (c_pm + 2.0_dp * c_0m + c_mm)) / (8.0_dp * mesh%dy)
          nrm = hypot(mx, my)

          if (nrm < min_grad) cycle ! Gradient diffus ou non univoque : INTERDIT

          snx = mx / nrm
          sny = my / nrm

          ! Condition 2.5 : Cohérence angulaire de la normale avec les marqueurs du front principal (incluant halos MPI)
          normal_coherent = on_jet_nozzle
          do dj = -1, 1
             do di = -1, 1
                if (di == 0 .and. dj == 0) cycle
                if (cnt_halo(i+di, j+dj) > 0.5_dp) then
                   nbx = nrm_x_halo(i+di, j+dj)
                   nby = nrm_y_halo(i+di, j+dj)
                   nb_len = hypot(nbx, nby)
                   if (nb_len > 1.0e-12_dp) then
                      dot_prod = (-snx) * (nbx / nb_len) + (-sny) * (nby / nb_len)
                      if (dot_prod >= 0.5_dp) then
                         normal_coherent = .true.
                         exit
                      end if
                   end if
                end if
             end do
             if (normal_coherent) exit
          end do
          if (.not. normal_coherent) then
             do p = 1, cloud%count
                if (.not. cloud%list(p)%active) cycle
                if (abs(cloud%list(p)%i_cell - i) <= 1 .and. abs(cloud%list(p)%j_cell - j) <= 1) then
                   dot_prod = (-snx) * cloud%list(p)%nx + (-sny) * cloud%list(p)%ny
                   if (dot_prod >= 0.5_dp) then
                      normal_coherent = .true.
                      exit
                   end if
                end if
             end do
          end if
          if (.not. normal_coherent) cycle ! Normale discordante avec le front principal : INTERDIT

          ! Condition 2.6 : Découpage PLIC et longueur minimale de segment
          xmin = x_face(mesh, i)
          xmax = x_face(mesh, i + 1)
          alpha = compute_alpha(xmin, xmax, ymin, ymax, snx, sny, c_val)

          call clip_segment_in_cell(xmin, xmax, ymin, ymax, &
                                    snx, sny, alpha, &
                                    x1, y1, x2, y2, valid)
          if (.not. valid) cycle

          seg_len = hypot(x2 - x1, y2 - y1)
          if (seg_len < 0.15_dp * h) cycle ! Segment dégénéré : INTERDIT

          call ensure_capacity(cloud, cloud%count + n_per_cell + 16)
          do k = 1, n_per_cell
             s = (real(k, dp) - 0.5_dp) / real(n_per_cell, dp)
             cloud%count = cloud%count + 1
             cloud%list(cloud%count)%x      = x1 + s * (x2 - x1)
             cloud%list(cloud%count)%y      = y1 + s * (y2 - y1)
             cloud%list(cloud%count)%nx     = -snx
             cloud%list(cloud%count)%ny     = -sny
             cloud%list(cloud%count)%kappa  = 0.0_dp
             cloud%list(cloud%count)%i_cell = i
             cloud%list(cloud%count)%j_cell = j
             cloud%list(cloud%count)%active = .true.
          end do
          cloud%count_in_cell(i, j) = n_per_cell
       end do
    end do

    deallocate(cnt_halo, nrm_x_halo, nrm_y_halo)

    ! 3. Compactage en place du tableau de marqueurs
    new_count = 0
    do p = 1, cloud%count
       if (cloud%list(p)%active) then
          new_count = new_count + 1
          if (new_count /= p) then
             cloud%list(new_count) = cloud%list(p)
          end if
       end if
    end do
    cloud%count = new_count
  end subroutine prune_and_reseed_markers

  !=============================================================================
  ! Routine : compute_curvature_from_markers
  ! Rôle    : Estime la normale unitaire (marker_nx, marker_ny) et la courbure
  !           locale (marker_kappa) dans chaque cellule de la bande interfaciale
  !           par régression parabolique pondérée dans le repère tourné (tau, eta)
  !           avec filtrage angulaire (n_p . n_ref >= 0.5).
  !=============================================================================
  subroutine compute_curvature_from_markers(cloud, mesh, cfg, fields)
    type(marker_cloud_t), intent(inout) :: cloud
    type(mesh_2d_t),      intent(in)    :: mesh
    type(config_t),       intent(in)    :: cfg
    type(flow_fields_t),  intent(inout) :: fields

    integer, parameter :: N_GHOST_ATTR = 4
    integer  :: i, j, p, k, nx, ny, di, dj, rx_rad, ry_rad, n_pts, n_local, n_tot, cnt_nb
    integer  :: nw_out, ne_out, ns_out, nn_out, nw_in, ne_in, ns_in, nn_in
    integer  :: ierr, status(MPI_STATUS_SIZE)
    logical  :: per_x, per_y
    real(dp) :: cnt_s(1), cnt_r(1)
    real(dp) :: dx, dy, h, r_fit, dc_dx, dc_dy, grad_mag, sum_nb
    real(dp) :: x_min_loc, x_max_loc, y_min_loc, y_max_loc, xp, yp
    real(dp) :: n_ref_x, n_ref_y, sum_nx, sum_ny, norm_avg
    real(dp) :: e_eta_x, e_eta_y, e_tau_x, e_tau_y
    real(dp) :: x_bar, y_bar, x_near, y_near, sum_w, d_min, d_cur, xc, yc
    real(dp) :: dx_p, dy_p, tau_h, eta_h, d2_h, w_p
    real(dp) :: S0, S1, S2, S3, S4, T0, T1, T2
    real(dp) :: M(3, 3), rhs(3), det_M, a_hat, b_hat, c_hat
    real(dp) :: kap_val, kap_max, m_nx, m_ny, n_len, slope_p

    ! Tableaux combinés (marqueurs locaux 1..n_local + fantômes MPI/périodiques)
    real(dp), allocatable :: pt_x(:), pt_y(:), pt_nx(:), pt_ny(:)
    integer,  allocatable :: orig_idx(:)
    real(dp), allocatable :: send_w(:), send_e(:), recv_w(:), recv_e(:)
    real(dp), allocatable :: send_s(:), send_n(:), recv_s(:), recv_n(:)
    real(dp), allocatable :: flag_kap(:, :), flag_nrm(:, :), kap_smooth(:, :)

    ! Indexation spatiale rapide par cellule étendue
    integer, allocatable :: head(:, :), next_pt(:)

    if (.not. cfg%use_markers) return

    nx = mesh%nx_local
    ny = mesh%ny_local
    dx = mesh%dx
    dy = mesh%dy
    h  = max(dx, dy)
    r_fit = 2.5_dp * h
    kap_max = 2.0_dp / min(dx, dy)
    rx_rad = max(2, ceiling(r_fit / dx))
    ry_rad = max(2, ceiling(r_fit / dy))

    x_min_loc = x_face(mesh, 1)
    x_max_loc = x_face(mesh, nx + 1)
    y_min_loc = y_face(mesh, 1)
    y_max_loc = y_face(mesh, ny + 1)

    per_x = (cfg%bc_left%btype == BC_PERIODIC .or. cfg%bc_right%btype == BC_PERIODIC)
    per_y = (cfg%bc_bottom%btype == BC_PERIODIC .or. cfg%bc_top%btype == BC_PERIODIC)

    fields%marker_nx         = 0.0_dp
    fields%marker_ny         = 0.0_dp
    fields%marker_kappa      = 0.0_dp
    fields%has_marker_normal = .false.
    fields%has_marker_kappa  = .false.

    ! 1. Collecte des marqueurs locaux actifs
    n_local = 0
    do p = 1, cloud%count
       if (cloud%list(p)%active) n_local = n_local + 1
    end do

    allocate(pt_x(max(16, 4 * n_local + 128)))
    allocate(pt_y(max(16, 4 * n_local + 128)))
    allocate(pt_nx(max(16, 4 * n_local + 128)))
    allocate(pt_ny(max(16, 4 * n_local + 128)))
    allocate(orig_idx(max(16, 4 * n_local + 128)))

    n_tot = 0
    do p = 1, cloud%count
       if (.not. cloud%list(p)%active) cycle
       n_tot = n_tot + 1
       pt_x(n_tot)     = cloud%list(p)%x
       pt_y(n_tot)     = cloud%list(p)%y
       pt_nx(n_tot)    = cloud%list(p)%nx
       pt_ny(n_tot)    = cloud%list(p)%ny
       orig_idx(n_tot) = p
    end do

    ! 2. Échange des marqueurs fantômes (bande de largeur r_fit) en X puis Y
    if (mesh%px > 1) then
       allocate(send_w(N_GHOST_ATTR * max(1, n_tot)))
       allocate(send_e(N_GHOST_ATTR * max(1, n_tot)))
       nw_out = 0; ne_out = 0
       do p = 1, n_tot
          xp = pt_x(p)
          if (xp < x_min_loc + r_fit .and. mesh%neighbor_west /= MPI_PROC_NULL) then
             send_w(nw_out*N_GHOST_ATTR + 1) = merge(xp + mesh%lx, xp, per_x .and. mesh%coord_x == 0)
             send_w(nw_out*N_GHOST_ATTR + 2) = pt_y(p)
             send_w(nw_out*N_GHOST_ATTR + 3) = pt_nx(p)
             send_w(nw_out*N_GHOST_ATTR + 4) = pt_ny(p)
             nw_out = nw_out + 1
          end if
          if (xp >= x_max_loc - r_fit .and. mesh%neighbor_east /= MPI_PROC_NULL) then
             send_e(ne_out*N_GHOST_ATTR + 1) = merge(xp - mesh%lx, xp, per_x .and. mesh%coord_x == mesh%px - 1)
             send_e(ne_out*N_GHOST_ATTR + 2) = pt_y(p)
             send_e(ne_out*N_GHOST_ATTR + 3) = pt_nx(p)
             send_e(ne_out*N_GHOST_ATTR + 4) = pt_ny(p)
             ne_out = ne_out + 1
          end if
       end do
       nw_in = 0; ne_in = 0
       cnt_s(1) = real(nw_out, dp); cnt_r(1) = 0.0_dp
       call MPI_Sendrecv(cnt_s, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 501, &
                         cnt_r, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 501, &
                         mesh%cart_comm, status, ierr)
       ne_in = nint(cnt_r(1))
       cnt_s(1) = real(ne_out, dp); cnt_r(1) = 0.0_dp
       call MPI_Sendrecv(cnt_s, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 502, &
                         cnt_r, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 502, &
                         mesh%cart_comm, status, ierr)
       nw_in = nint(cnt_r(1))
       allocate(recv_e(N_GHOST_ATTR * max(1, ne_in)))
       allocate(recv_w(N_GHOST_ATTR * max(1, nw_in)))
       call MPI_Sendrecv(send_w, nw_out * N_GHOST_ATTR, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 503, &
                         recv_e, ne_in  * N_GHOST_ATTR, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 503, &
                         mesh%cart_comm, status, ierr)
       call MPI_Sendrecv(send_e, ne_out * N_GHOST_ATTR, MPI_DOUBLE_PRECISION, mesh%neighbor_east, 504, &
                         recv_w, nw_in  * N_GHOST_ATTR, MPI_DOUBLE_PRECISION, mesh%neighbor_west, 504, &
                         mesh%cart_comm, status, ierr)
       call grow_ghost_arrays(pt_x, pt_y, pt_nx, pt_ny, orig_idx, n_tot + ne_in + nw_in + 16)
       do k = 0, ne_in - 1
          n_tot = n_tot + 1
          pt_x(n_tot)     = recv_e(k*N_GHOST_ATTR + 1)
          pt_y(n_tot)     = recv_e(k*N_GHOST_ATTR + 2)
          pt_nx(n_tot)    = recv_e(k*N_GHOST_ATTR + 3)
          pt_ny(n_tot)    = recv_e(k*N_GHOST_ATTR + 4)
          orig_idx(n_tot) = 0
       end do
       do k = 0, nw_in - 1
          n_tot = n_tot + 1
          pt_x(n_tot)     = recv_w(k*N_GHOST_ATTR + 1)
          pt_y(n_tot)     = recv_w(k*N_GHOST_ATTR + 2)
          pt_nx(n_tot)    = recv_w(k*N_GHOST_ATTR + 3)
          pt_ny(n_tot)    = recv_w(k*N_GHOST_ATTR + 4)
          orig_idx(n_tot) = 0
       end do
       deallocate(send_w, send_e, recv_w, recv_e)
    else if (per_x) then
       k = n_tot
       call grow_ghost_arrays(pt_x, pt_y, pt_nx, pt_ny, orig_idx, 3 * k + 16)
       do p = 1, k
          xp = pt_x(p)
          if (xp < x_min_loc + r_fit) then
             n_tot = n_tot + 1
             pt_x(n_tot) = xp + mesh%lx; pt_y(n_tot) = pt_y(p)
             pt_nx(n_tot) = pt_nx(p);    pt_ny(n_tot) = pt_ny(p)
             orig_idx(n_tot) = 0
          end if
          if (xp >= x_max_loc - r_fit) then
             n_tot = n_tot + 1
             pt_x(n_tot) = xp - mesh%lx; pt_y(n_tot) = pt_y(p)
             pt_nx(n_tot) = pt_nx(p);    pt_ny(n_tot) = pt_ny(p)
             orig_idx(n_tot) = 0
          end if
       end do
    end if

    if (mesh%py > 1) then
       allocate(send_s(N_GHOST_ATTR * max(1, n_tot)))
       allocate(send_n(N_GHOST_ATTR * max(1, n_tot)))
       ns_out = 0; nn_out = 0
       do p = 1, n_tot
          yp = pt_y(p)
          if (yp < y_min_loc + r_fit .and. mesh%neighbor_south /= MPI_PROC_NULL) then
             send_s(ns_out*N_GHOST_ATTR + 1) = pt_x(p)
             send_s(ns_out*N_GHOST_ATTR + 2) = merge(yp + mesh%ly, yp, per_y .and. mesh%coord_y == 0)
             send_s(ns_out*N_GHOST_ATTR + 3) = pt_nx(p)
             send_s(ns_out*N_GHOST_ATTR + 4) = pt_ny(p)
             ns_out = ns_out + 1
          end if
          if (yp >= y_max_loc - r_fit .and. mesh%neighbor_north /= MPI_PROC_NULL) then
             send_n(nn_out*N_GHOST_ATTR + 1) = pt_x(p)
             send_n(nn_out*N_GHOST_ATTR + 2) = merge(yp - mesh%ly, yp, per_y .and. mesh%coord_y == mesh%py - 1)
             send_n(nn_out*N_GHOST_ATTR + 3) = pt_nx(p)
             send_n(nn_out*N_GHOST_ATTR + 4) = pt_ny(p)
             nn_out = nn_out + 1
          end if
       end do
       ns_in = 0; nn_in = 0
       cnt_s(1) = real(ns_out, dp); cnt_r(1) = 0.0_dp
       call MPI_Sendrecv(cnt_s, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 505, &
                         cnt_r, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 505, &
                         mesh%cart_comm, status, ierr)
       nn_in = nint(cnt_r(1))
       cnt_s(1) = real(nn_out, dp); cnt_r(1) = 0.0_dp
       call MPI_Sendrecv(cnt_s, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 506, &
                         cnt_r, 1, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 506, &
                         mesh%cart_comm, status, ierr)
       ns_in = nint(cnt_r(1))
       allocate(recv_n(N_GHOST_ATTR * max(1, nn_in)))
       allocate(recv_s(N_GHOST_ATTR * max(1, ns_in)))
       call MPI_Sendrecv(send_s, ns_out * N_GHOST_ATTR, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 507, &
                         recv_n, nn_in  * N_GHOST_ATTR, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 507, &
                         mesh%cart_comm, status, ierr)
       call MPI_Sendrecv(send_n, nn_out * N_GHOST_ATTR, MPI_DOUBLE_PRECISION, mesh%neighbor_north, 508, &
                         recv_s, ns_in  * N_GHOST_ATTR, MPI_DOUBLE_PRECISION, mesh%neighbor_south, 508, &
                         mesh%cart_comm, status, ierr)
       call grow_ghost_arrays(pt_x, pt_y, pt_nx, pt_ny, orig_idx, n_tot + nn_in + ns_in + 16)
       do k = 0, nn_in - 1
          n_tot = n_tot + 1
          pt_x(n_tot)     = recv_n(k*N_GHOST_ATTR + 1)
          pt_y(n_tot)     = recv_n(k*N_GHOST_ATTR + 2)
          pt_nx(n_tot)    = recv_n(k*N_GHOST_ATTR + 3)
          pt_ny(n_tot)    = recv_n(k*N_GHOST_ATTR + 4)
          orig_idx(n_tot) = 0
       end do
       do k = 0, ns_in - 1
          n_tot = n_tot + 1
          pt_x(n_tot)     = recv_s(k*N_GHOST_ATTR + 1)
          pt_y(n_tot)     = recv_s(k*N_GHOST_ATTR + 2)
          pt_nx(n_tot)    = recv_s(k*N_GHOST_ATTR + 3)
          pt_ny(n_tot)    = recv_s(k*N_GHOST_ATTR + 4)
          orig_idx(n_tot) = 0
       end do
       deallocate(send_s, send_n, recv_s, recv_n)
    else if (per_y) then
       k = n_tot
       call grow_ghost_arrays(pt_x, pt_y, pt_nx, pt_ny, orig_idx, 3 * k + 16)
       do p = 1, k
          yp = pt_y(p)
          if (yp < y_min_loc + r_fit) then
             n_tot = n_tot + 1
             pt_x(n_tot) = pt_x(p);      pt_y(n_tot) = yp + mesh%ly
             pt_nx(n_tot) = pt_nx(p);    pt_ny(n_tot) = pt_ny(p)
             orig_idx(n_tot) = 0
          end if
          if (yp >= y_max_loc - r_fit) then
             n_tot = n_tot + 1
             pt_x(n_tot) = pt_x(p);      pt_y(n_tot) = yp - mesh%ly
             pt_nx(n_tot) = pt_nx(p);    pt_ny(n_tot) = pt_ny(p)
             orig_idx(n_tot) = 0
          end if
       end do
    end if

    allocate(head(1-rx_rad:nx+rx_rad, 1-ry_rad:ny+ry_rad))
    allocate(next_pt(max(1, n_tot)))
    head = 0
    next_pt = 0

    do p = 1, n_tot
       i = floor((pt_x(p) - x_min_loc) / dx) + 1
       j = floor((pt_y(p) - y_min_loc) / dy) + 1
       if (i >= 1 - rx_rad .and. i <= nx + rx_rad .and. j >= 1 - ry_rad .and. j <= ny + ry_rad) then
          next_pt(p) = head(i, j)
          head(i, j) = p
       end if
    end do

    do j = 1, ny
       do i = 1, nx
          ! Calcul du gradient local de C pour identifier la bande interfaciale
          dc_dx = ((fields%C(i+1, j+1) + 2.0_dp*fields%C(i+1, j) + fields%C(i+1, j-1)) - &
                   (fields%C(i-1, j+1) + 2.0_dp*fields%C(i-1, j) + fields%C(i-1, j-1))) / (8.0_dp * dx)
          dc_dy = ((fields%C(i+1, j+1) + 2.0_dp*fields%C(i, j+1) + fields%C(i-1, j+1)) - &
                   (fields%C(i+1, j-1) + 2.0_dp*fields%C(i, j-1) + fields%C(i-1, j-1))) / (8.0_dp * dy)
          grad_mag = hypot(dc_dx, dc_dy)

          if (grad_mag <= 1.0e-4_dp / h .and. head(i, j) == 0) cycle

          if (grad_mag > 1.0e-12_dp) then
             n_ref_x = dc_dx / grad_mag
             n_ref_y = dc_dy / grad_mag
          else
             p = head(i, j)
             if (p == 0) cycle
             n_ref_x = pt_nx(p)
             n_ref_y = pt_ny(p)
          end if

          ! 1. Recherche du centre local (x_bar, y_bar) et de la normale moyenne filtrée
          xc = x_center(mesh, i)
          yc = y_center(mesh, j)
          sum_nx = 0.0_dp; sum_ny = 0.0_dp
          x_bar = 0.0_dp;  y_bar = 0.0_dp; sum_w = 0.0_dp
          x_near = xc;     y_near = yc
          d_min = 1.0e30_dp
          n_pts = 0

          do dj = -ry_rad, ry_rad
             do di = -rx_rad, rx_rad
                p = head(i + di, j + dj)
                do while (p /= 0)
                   ! Filtrage angulaire pour protéger les films minces (cos(60 deg) = 0.5)
                   if (pt_nx(p) * n_ref_x + pt_ny(p) * n_ref_y >= 0.5_dp) then
                      d_cur = hypot(pt_x(p) - xc, pt_y(p) - yc)
                      if (d_cur <= r_fit) then
                         n_pts = n_pts + 1
                         sum_nx = sum_nx + pt_nx(p)
                         sum_ny = sum_ny + pt_ny(p)
                         if (di == 0 .and. dj == 0) then
                            x_bar = x_bar + pt_x(p)
                            y_bar = y_bar + pt_y(p)
                            sum_w = sum_w + 1.0_dp
                         end if
                         if (d_cur < d_min) then
                            d_min = d_cur
                            x_near = pt_x(p)
                            y_near = pt_y(p)
                         end if
                      end if
                   end if
                   p = next_pt(p)
                end do
             end do
          end do

          if (n_pts < 3) cycle
          if (sum_w > 0.0_dp) then
             x_bar = x_bar / sum_w
             y_bar = y_bar / sum_w
          else
             x_bar = x_near
             y_bar = y_near
          end if

          norm_avg = hypot(sum_nx, sum_ny)
          if (norm_avg > 1.0e-12_dp) then
             e_eta_x = sum_nx / norm_avg
             e_eta_y = sum_ny / norm_avg
          else
             e_eta_x = n_ref_x
             e_eta_y = n_ref_y
          end if
          e_tau_x = -e_eta_y
          e_tau_y =  e_eta_x

          ! 2. Assemblage du système normal 3x3 en coordonnées adimensionnées (tau/h, eta/h)
          S0 = 0.0_dp; S1 = 0.0_dp; S2 = 0.0_dp; S3 = 0.0_dp; S4 = 0.0_dp
          T0 = 0.0_dp; T1 = 0.0_dp; T2 = 0.0_dp

          do dj = -ry_rad, ry_rad
             do di = -rx_rad, rx_rad
                p = head(i + di, j + dj)
                do while (p /= 0)
                   if (pt_nx(p) * e_eta_x + pt_ny(p) * e_eta_y >= 0.5_dp) then
                      dx_p = pt_x(p) - x_bar
                      dy_p = pt_y(p) - y_bar
                      tau_h = (dx_p * e_tau_x + dy_p * e_tau_y) / h
                      eta_h = (dx_p * e_eta_x + dy_p * e_eta_y) / h
                      d2_h  = tau_h * tau_h + eta_h * eta_h
                      if (d2_h <= (2.5_dp)**2) then
                         w_p = exp(-d2_h / (1.5_dp**2))
                         S0 = S0 + w_p
                         S1 = S1 + w_p * tau_h
                         S2 = S2 + w_p * (tau_h**2)
                         S3 = S3 + w_p * (tau_h**3)
                         S4 = S4 + w_p * (tau_h**4)
                         T0 = T0 + w_p * eta_h
                         T1 = T1 + w_p * eta_h * tau_h
                         T2 = T2 + w_p * eta_h * (tau_h**2)
                      end if
                   end if
                   p = next_pt(p)
                end do
             end do
          end do

          ! Régularisation de Tikhonov légère sur la courbure pour conditionnement inconditionnel
          S4 = S4 + 1.0e-9_dp
          S2 = S2 + 1.0e-12_dp
          S0 = S0 + 1.0e-12_dp

          M(1, 1) = S4; M(1, 2) = S3; M(1, 3) = S2
          M(2, 1) = S3; M(2, 2) = S2; M(2, 3) = S1
          M(3, 1) = S2; M(3, 2) = S1; M(3, 3) = S0
          rhs(1)  = T2; rhs(2)  = T1; rhs(3)  = T0

          det_M = M(1,1)*(M(2,2)*M(3,3) - M(2,3)*M(3,2)) &
                - M(1,2)*(M(2,1)*M(3,3) - M(2,3)*M(3,1)) &
                + M(1,3)*(M(2,1)*M(3,2) - M(2,2)*M(3,1))

          if (abs(det_M) > 1.0e-14_dp) then
             a_hat = ( rhs(1)*(M(2,2)*M(3,3) - M(2,3)*M(3,2)) &
                     - M(1,2)*(rhs(2)*M(3,3) - M(2,3)*rhs(3)) &
                     + M(1,3)*(rhs(2)*M(3,2) - M(2,2)*rhs(3)) ) / det_M

             b_hat = ( M(1,1)*(rhs(2)*M(3,3) - M(2,3)*rhs(3)) &
                     - rhs(1)*(M(2,1)*M(3,3) - M(2,3)*M(3,1)) &
                     + M(1,3)*(M(2,1)*rhs(3) - rhs(2)*M(3,1)) ) / det_M

             c_hat = ( M(1,1)*(M(2,2)*rhs(3) - rhs(2)*M(3,2)) &
                     - M(1,2)*(M(2,1)*rhs(3) - rhs(2)*M(3,1)) &
                     + rhs(1)*(M(2,1)*M(3,2) - M(2,2)*M(3,1)) ) / det_M

             kap_val = (2.0_dp * (a_hat / h)) / ((1.0_dp + b_hat * b_hat)**1.5_dp)
             kap_val = max(-kap_max, min(kap_max, kap_val))

             m_nx = e_eta_x - b_hat * e_tau_x
             m_ny = e_eta_y - b_hat * e_tau_y
             n_len = hypot(m_nx, m_ny)
             if (n_len > 1.0e-12_dp) then
                m_nx = m_nx / n_len
                m_ny = m_ny / n_len
             else
                m_nx = e_eta_x
                m_ny = e_eta_y
             end if

             fields%marker_nx(i, j)         = m_nx
             fields%marker_ny(i, j)         = m_ny
             fields%marker_kappa(i, j)      = kap_val
             if (fields%C(i, j) > 0.01_dp .and. fields%C(i, j) < 1.0_dp - 0.01_dp .and. n_pts >= 3) then
                fields%has_marker_normal(i, j) = .true.
             end if
             fields%has_marker_kappa(i, j)  = .true.

             ! Mise à jour de la normale et de la courbure des marqueurs locaux de la cellule (i, j)
             p = head(i, j)
             do while (p /= 0)
                k = orig_idx(p)
                if (k > 0 .and. (pt_nx(p) * e_eta_x + pt_ny(p) * e_eta_y >= 0.5_dp)) then
                   dx_p = pt_x(p) - x_bar
                   dy_p = pt_y(p) - y_bar
                   tau_h = (dx_p * e_tau_x + dy_p * e_tau_y) / h
                   slope_p = 2.0_dp * a_hat * tau_h + b_hat
                   m_nx = e_eta_x - slope_p * e_tau_x
                   m_ny = e_eta_y - slope_p * e_tau_y
                   n_len = hypot(m_nx, m_ny)
                   if (n_len > 1.0e-12_dp) then
                      cloud%list(k)%nx = m_nx / n_len
                      cloud%list(k)%ny = m_ny / n_len
                   end if
                   cloud%list(k)%kappa = max(-kap_max, min(kap_max, &
                        (2.0_dp * (a_hat / h)) / ((1.0_dp + slope_p * slope_p)**1.5_dp)))
                end if
                p = next_pt(p)
             end do
          end if
       end do
    end do

    call exchange_halos_scalar(fields%marker_nx, mesh)
    call exchange_halos_scalar(fields%marker_ny, mesh)
    call exchange_halos_scalar(fields%marker_kappa, mesh)

    allocate(flag_kap(0:nx+1, 0:ny+1), flag_nrm(0:nx+1, 0:ny+1), kap_smooth(0:nx+1, 0:ny+1))
    flag_kap = merge(1.0_dp, 0.0_dp, fields%has_marker_kappa)
    flag_nrm = merge(1.0_dp, 0.0_dp, fields%has_marker_normal)
    call exchange_halos_scalar(flag_kap, mesh)
    call exchange_halos_scalar(flag_nrm, mesh)
    fields%has_marker_kappa  = (flag_kap > 0.5_dp)
    fields%has_marker_normal = (flag_nrm > 0.5_dp)

    ! Passe conservative de lissage tangentiel de courbure
    kap_smooth = fields%marker_kappa
    do j = 1, ny
       do i = 1, nx
          if (.not. fields%has_marker_kappa(i, j)) cycle
          m_nx = fields%marker_nx(i, j)
          m_ny = fields%marker_ny(i, j)
          sum_nb = 0.0_dp
          cnt_nb = 0
          do dj = -1, 1
             do di = -1, 1
                if (di == 0 .and. dj == 0) cycle
                if (fields%has_marker_kappa(i + di, j + dj) .and. &
                    (m_nx * fields%marker_nx(i + di, j + dj) + m_ny * fields%marker_ny(i + di, j + dj) >= 0.5_dp)) then
                   sum_nb = sum_nb + fields%marker_kappa(i + di, j + dj)
                   cnt_nb = cnt_nb + 1
                end if
             end do
          end do
          if (cnt_nb > 0) then
             kap_smooth(i, j) = max(-kap_max, min(kap_max, &
                  0.5_dp * fields%marker_kappa(i, j) + 0.5_dp * (sum_nb / real(cnt_nb, dp))))
          end if
       end do
    end do
    fields%marker_kappa = kap_smooth
    call exchange_halos_scalar(fields%marker_kappa, mesh)

    deallocate(pt_x, pt_y, pt_nx, pt_ny, orig_idx, head, next_pt, flag_kap, flag_nrm, kap_smooth)
  end subroutine compute_curvature_from_markers

  subroutine grow_ghost_arrays(px, py, pnx, pny, pidx, req_cap)
    real(dp), allocatable, intent(inout) :: px(:), py(:), pnx(:), pny(:)
    integer,  allocatable, intent(inout) :: pidx(:)
    integer, intent(in) :: req_cap
    real(dp), allocatable :: tx(:), ty(:), tnx(:), tny(:)
    integer,  allocatable :: tidx(:)
    integer :: cur_cap, new_cap

    cur_cap = size(px)
    if (req_cap <= cur_cap) return
    new_cap = max(req_cap, 2 * cur_cap)
    allocate(tx(new_cap), ty(new_cap), tnx(new_cap), tny(new_cap), tidx(new_cap))
    tx(1:cur_cap)   = px(1:cur_cap)
    ty(1:cur_cap)   = py(1:cur_cap)
    tnx(1:cur_cap)  = pnx(1:cur_cap)
    tny(1:cur_cap)  = pny(1:cur_cap)
    tidx(1:cur_cap) = pidx(1:cur_cap)
    call move_alloc(tx, px)
    call move_alloc(ty, py)
    call move_alloc(tnx, pnx)
    call move_alloc(tny, pny)
    call move_alloc(tidx, pidx)
  end subroutine grow_ghost_arrays

  !=============================================================================
  ! Routine : write_markers_vtk
  ! Rôle    : Exporte le nuage global de marqueurs (collecté sur le rang 0) au
  !           format VTK Legacy PolyData (.vtk) avec les vecteurs normaux et la
  !           courbure locale
  !=============================================================================
  subroutine write_markers_vtk(cloud, mesh, output_dir, case_name, step)
    type(marker_cloud_t), intent(in) :: cloud
    type(mesh_2d_t),      intent(in) :: mesh
    character(len=*),     intent(in) :: output_dir, case_name
    integer,              intent(in) :: step

    character(len=6)   :: step_str
    character(len=512) :: vtk_file
    integer :: vtk_unit, io_stat, p, n_loc, n_tot, r, idx, ierr
    real(dp), allocatable :: loc_buf(:), glob_buf(:)
    integer,  allocatable :: counts(:), counts5(:), displs5(:)

    n_loc = 0
    do p = 1, cloud%count
       if (cloud%list(p)%active) n_loc = n_loc + 1
    end do

    allocate(loc_buf(max(1, 5 * n_loc)))
    loc_buf = 0.0_dp
    idx = 0
    do p = 1, cloud%count
       if (cloud%list(p)%active) then
          loc_buf(idx + 1) = cloud%list(p)%x
          loc_buf(idx + 2) = cloud%list(p)%y
          loc_buf(idx + 3) = cloud%list(p)%nx
          loc_buf(idx + 4) = cloud%list(p)%ny
          loc_buf(idx + 5) = cloud%list(p)%kappa
          idx = idx + 5
       end if
    end do

    allocate(counts(mesh%size), counts5(mesh%size), displs5(mesh%size))
    counts = 0; counts5 = 0; displs5 = 0

    call MPI_Gather(n_loc, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, mesh%comm, ierr)

    n_tot = 0
    if (mesh%rank == 0) then
       do r = 1, mesh%size
          counts5(r) = 5 * counts(r)
          displs5(r) = 5 * n_tot
          n_tot = n_tot + counts(r)
       end do
    end if

    allocate(glob_buf(max(1, 5 * n_tot)))
    glob_buf = 0.0_dp

    call MPI_Gatherv(loc_buf, 5 * n_loc, MPI_DOUBLE_PRECISION, &
                     glob_buf, counts5, displs5, MPI_DOUBLE_PRECISION, &
                     0, mesh%comm, ierr)

    if (mesh%rank == 0) then
       write(step_str, '(I6.6)') step
       write(vtk_file, '(A,"/",A,"_markers_",A,".vtk")') &
         trim(output_dir), trim(case_name), step_str

       open(newunit=vtk_unit, file=trim(vtk_file), status="replace", action="write", iostat=io_stat)
       if (io_stat == 0) then
          write(vtk_unit, '(A)') "# vtk DataFile Version 3.0"
          write(vtk_unit, '(A)') "COVOF 2D Unconnected Interface Markers (CVI)"
          write(vtk_unit, '(A)') "ASCII"
          write(vtk_unit, '(A)') "DATASET POLYDATA"
          write(vtk_unit, '(A,I0,A)') "POINTS ", n_tot, " double"
          do p = 1, n_tot
             write(vtk_unit, '(ES16.8," ",ES16.8," 0.0")') glob_buf(5*(p-1)+1), glob_buf(5*(p-1)+2)
          end do

          write(vtk_unit, '(A,I0,A,I0)') "VERTICES ", n_tot, " ", 2 * n_tot
          do p = 0, n_tot - 1
             write(vtk_unit, '("1 ",I0)') p
          end do

          write(vtk_unit, '(A,I0)') "POINT_DATA ", n_tot
          write(vtk_unit, '(A)') "VECTORS MarkerNormal double"
          do p = 1, n_tot
             write(vtk_unit, '(ES16.8," ",ES16.8," 0.0")') glob_buf(5*(p-1)+3), glob_buf(5*(p-1)+4)
          end do

          write(vtk_unit, '(A)') "SCALARS MarkerCurvature double 1"
          write(vtk_unit, '(A)') "LOOKUP_TABLE default"
          do p = 1, n_tot
             write(vtk_unit, '(ES16.8)') glob_buf(5*(p-1)+5)
          end do

          close(vtk_unit)
       end if
    end if

    deallocate(loc_buf, glob_buf, counts, counts5, displs5)
  end subroutine write_markers_vtk

end module module_markers
