! [Procédure/Module] module module_interface_reconstruction
module module_interface_reconstruction
  use module_precision, only: dp
  use module_mesh, only: mesh_2d_t, x_face, y_face
  implicit none
  private

  ! Types publics
  public :: interface_segment_t

  ! Procédures publiques
  public :: cut_box_area
  public :: compute_alpha
  public :: reconstruct_slic
  public :: reconstruct_plic
  public :: reconstruct_elvira
  public :: reconstruct_lvira
  public :: reconstruct_interface
  public :: evaluate_limiter

  ! Type représentant un segment de droite d'interface affine : nx*x + ny*y + alpha = 0
  type :: interface_segment_t
     real(dp) :: nx = 0.0_dp
     real(dp) :: ny = 0.0_dp
     real(dp) :: alpha = 0.0_dp
  end type interface_segment_t

contains

  ! [Procédure/Module] function cut_box_area
  ! Calcule l'aire exacte de découpe d'un rectangle par un demi-plan (Sutherland-Hodgman + Shoelace)
  function cut_box_area(xmin, xmax, ymin, ymax, nx, ny, alpha) result(area)
    real(dp), intent(in) :: xmin, xmax, ymin, ymax, nx, ny, alpha
    real(dp) :: area

    real(dp) :: bx(4), by(4)
    real(dp) :: d(4)
    real(dp) :: cx(8), cy(8)
    integer  :: n_clipped, i, next_i, p, next_p
    real(dp) :: d1, d2, t
    real(dp) :: sum_cross

    if ((xmax - xmin) <= 1.0e-14_dp .or. (ymax - ymin) <= 1.0e-14_dp) then
       area = 0.0_dp
       return
    end if

    ! Définition des 4 sommets du rectangle en sens trigonométrique
    bx(1) = xmin; by(1) = ymin
    bx(2) = xmax; by(2) = ymin
    bx(3) = xmax; by(3) = ymax
    bx(4) = xmin; by(4) = ymax

    do i = 1, 4
       d(i) = nx * bx(i) + ny * by(i) + alpha
    end do

    n_clipped = 0
    do i = 1, 4
       next_i = mod(i, 4) + 1
       d1 = d(i)
       d2 = d(next_i)

       if (d1 <= 0.0_dp) then
          if (d2 <= 0.0_dp) then
             n_clipped = n_clipped + 1
             cx(n_clipped) = bx(next_i)
             cy(n_clipped) = by(next_i)
          else
             t = -d1 / (d2 - d1)
             n_clipped = n_clipped + 1
             cx(n_clipped) = bx(i) + t * (bx(next_i) - bx(i))
             cy(n_clipped) = by(i) + t * (by(next_i) - by(i))
          end if
       else
          if (d2 <= 0.0_dp) then
             t = -d1 / (d2 - d1)
             n_clipped = n_clipped + 1
             cx(n_clipped) = bx(i) + t * (bx(next_i) - bx(i))
             cy(n_clipped) = by(i) + t * (by(next_i) - by(i))
             n_clipped = n_clipped + 1
             cx(n_clipped) = bx(next_i)
             cy(n_clipped) = by(next_i)
          end if
       end if
    end do

    if (n_clipped < 3) then
       area = 0.0_dp
       return
    end if

    ! Formule de Gauss / Shoelace pour l'aire du polygone
    sum_cross = 0.0_dp
    do p = 1, n_clipped
       next_p = mod(p, n_clipped) + 1
       sum_cross = sum_cross + (cx(p) * cy(next_p) - cx(next_p) * cy(p))
    end do
    area = 0.5_dp * abs(sum_cross)
  end function cut_box_area

  ! [Procédure/Module] function compute_alpha
  ! Inversion géométrique par bisection : C_target -> paramètre alpha
  function compute_alpha(xmin, xmax, ymin, ymax, nx, ny, c_target) result(alpha)
    real(dp), intent(in) :: xmin, xmax, ymin, ymax, nx, ny, c_target
    real(dp) :: alpha

    real(dp) :: total_area, target_area
    real(dp) :: v1, v2, v3, v4, min_val, max_val
    real(dp) :: alpha_low, alpha_high, alpha_mid, a
    integer  :: iter

    total_area = (xmax - xmin) * (ymax - ymin)
    if (c_target <= 1.0e-7_dp) then
       alpha = 1.0e9_dp
       return
    end if
    if (c_target >= 1.0_dp - 1.0e-7_dp) then
       alpha = -1.0e9_dp
       return
    end if

    target_area = max(0.0_dp, min(1.0_dp, c_target)) * total_area

    v1 = nx * xmin + ny * ymin
    v2 = nx * xmax + ny * ymin
    v3 = nx * xmax + ny * ymax
    v4 = nx * xmin + ny * ymax

    min_val = min(v1, min(v2, min(v3, v4)))
    max_val = max(v1, max(v2, max(v3, v4)))

    alpha_low  = -max_val - 1.0e-12_dp
    alpha_high = -min_val + 1.0e-12_dp

    do iter = 1, 45
       alpha_mid = 0.5_dp * (alpha_low + alpha_high)
       a = cut_box_area(xmin, xmax, ymin, ymax, nx, ny, alpha_mid)
       if (a < target_area) then
          alpha_high = alpha_mid
       else
          alpha_low = alpha_mid
       end if
    end do
    alpha = 0.5_dp * (alpha_low + alpha_high)
  end function compute_alpha

  ! [Procédure/Module] subroutine reconstruct_slic
  subroutine reconstruct_slic(mesh, C, segments)
    type(mesh_2d_t), intent(in) :: mesh
    real(dp), intent(in) :: C(0:mesh%nx_local+1, 0:mesh%ny_local+1)
    type(interface_segment_t), intent(out) :: segments(0:mesh%nx_local+1, 0:mesh%ny_local+1)

    integer :: i, j, nx, ny
    real(dp) :: c_val, xmin, xmax, ymin, ymax, grad_x, grad_y, snx, sny, a

    nx = mesh%nx_local
    ny = mesh%ny_local

    do j = 0, ny + 1
       do i = 0, nx + 1
          c_val = C(i, j)
          if (c_val <= 1.0e-6_dp) then
             segments(i, j)%nx = 1.0_dp
             segments(i, j)%ny = 0.0_dp
             segments(i, j)%alpha = 1.0e9_dp
             cycle
          end if
          if (c_val >= 1.0_dp - 1.0e-6_dp) then
             segments(i, j)%nx = 1.0_dp
             segments(i, j)%ny = 0.0_dp
             segments(i, j)%alpha = -1.0e9_dp
             cycle
          end if

          xmin = x_face(mesh, i)
          xmax = x_face(mesh, i + 1)
          ymin = y_face(mesh, j)
          ymax = y_face(mesh, j + 1)

          if (i > 0 .and. i < nx + 1 .and. j > 0 .and. j < ny + 1) then
             grad_x = abs(C(i + 1, j) - C(i - 1, j))
             grad_y = abs(C(i, j + 1) - C(i, j - 1))
          else
             grad_x = 1.0_dp
             grad_y = 0.0_dp
          end if

          if (grad_x >= grad_y) then
             if (i > 0 .and. C(i - 1, j) >= C(i + 1, j)) then
                snx = 1.0_dp; sny = 0.0_dp
             else
                snx = -1.0_dp; sny = 0.0_dp
             end if
          else
             if (j > 0 .and. C(i, j - 1) >= C(i, j + 1)) then
                snx = 0.0_dp; sny = 1.0_dp
             else
                snx = 0.0_dp; sny = -1.0_dp
             end if
          end if

          a = compute_alpha(xmin, xmax, ymin, ymax, snx, sny, c_val)
          segments(i, j)%nx = snx
          segments(i, j)%ny = sny
          segments(i, j)%alpha = a
       end do
    end do
  end subroutine reconstruct_slic

  ! [Procédure/Module] subroutine reconstruct_plic
  subroutine reconstruct_plic(mesh, C, segments)
    type(mesh_2d_t), intent(in) :: mesh
    real(dp), intent(in) :: C(0:mesh%nx_local+1, 0:mesh%ny_local+1)
    type(interface_segment_t), intent(out) :: segments(0:mesh%nx_local+1, 0:mesh%ny_local+1)

    integer :: i, j, nx, ny
    real(dp) :: c_val, dx, dy, xmin, xmax, ymin, ymax
    real(dp) :: c_mm, c_0m, c_pm, c_m0, c_p0, c_mp, c_0p, c_pp
    real(dp) :: mx, my, norm, pnx, pny, a

    nx = mesh%nx_local
    ny = mesh%ny_local
    dx = mesh%dx
    dy = mesh%dy

    do j = 0, ny + 1
       do i = 0, nx + 1
          c_val = C(i, j)
          if (c_val <= 1.0e-6_dp) then
             segments(i, j)%nx = 1.0_dp
             segments(i, j)%ny = 0.0_dp
             segments(i, j)%alpha = 1.0e9_dp
             cycle
          end if
          if (c_val >= 1.0_dp - 1.0e-6_dp) then
             segments(i, j)%nx = 1.0_dp
             segments(i, j)%ny = 0.0_dp
             segments(i, j)%alpha = -1.0e9_dp
             cycle
          end if

          if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
             c_mm = C(i - 1, j - 1); c_0m = C(i,     j - 1); c_pm = C(i + 1, j - 1)
             c_m0 = C(i - 1, j    );                         c_p0 = C(i + 1, j    )
             c_mp = C(i - 1, j + 1); c_0p = C(i,     j + 1); c_pp = C(i + 1, j + 1)

             mx = -((c_pp + 2.0_dp*c_p0 + c_pm) - (c_mp + 2.0_dp*c_m0 + c_mm)) / (8.0_dp * dx)
             my = -((c_pp + 2.0_dp*c_0p + c_mp) - (c_pm + 2.0_dp*c_0m + c_mm)) / (8.0_dp * dy)
          else
             mx = 1.0_dp; my = 0.0_dp
          end if

          norm = hypot(mx, my)
          if (norm > 1.0e-12_dp) then
             pnx = mx / norm; pny = my / norm
          else
             pnx = 1.0_dp; pny = 0.0_dp
          end if

          xmin = x_face(mesh, i)
          xmax = x_face(mesh, i + 1)
          ymin = y_face(mesh, j)
          ymax = y_face(mesh, j + 1)

          a = compute_alpha(xmin, xmax, ymin, ymax, pnx, pny, c_val)
          segments(i, j)%nx = pnx
          segments(i, j)%ny = pny
          segments(i, j)%alpha = a
       end do
    end do
  end subroutine reconstruct_plic

  ! [Procédure/Module] subroutine reconstruct_elvira
  subroutine reconstruct_elvira(mesh, C, segments)
    type(mesh_2d_t), intent(in) :: mesh
    real(dp), intent(in) :: C(0:mesh%nx_local+1, 0:mesh%ny_local+1)
    type(interface_segment_t), intent(out) :: segments(0:mesh%nx_local+1, 0:mesh%ny_local+1)

    integer :: i, j, nx, ny, di, dj, k
    real(dp) :: c_val, dx, dy, xmin, xmax, ymin, ymax
    real(dp) :: H_m1, H_0, H_p1, W_m1, W_0, W_p1
    real(dp) :: sgn_x, sgn_y, mx, my
    real(dp) :: cand_nx(7), cand_ny(7)
    real(dp) :: best_err, err, cnorm, cnx, cny, calpha, diff
    real(dp) :: best_nx, best_ny, best_alpha
    real(dp) :: v_xmin, v_xmax, v_ymin, v_ymax, v_target, v_recon

    nx = mesh%nx_local
    ny = mesh%ny_local
    dx = mesh%dx
    dy = mesh%dy

    do j = 0, ny + 1
       do i = 0, nx + 1
          c_val = C(i, j)
          if (c_val <= 1.0e-6_dp) then
             segments(i, j)%nx = 1.0_dp
             segments(i, j)%ny = 0.0_dp
             segments(i, j)%alpha = 1.0e9_dp
             cycle
          end if
          if (c_val >= 1.0_dp - 1.0e-6_dp) then
             segments(i, j)%nx = 1.0_dp
             segments(i, j)%ny = 0.0_dp
             segments(i, j)%alpha = -1.0e9_dp
             cycle
          end if

          xmin = x_face(mesh, i)
          xmax = x_face(mesh, i + 1)
          ymin = y_face(mesh, j)
          ymax = y_face(mesh, j + 1)

          if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
             H_m1 = (C(i - 1, j - 1) + C(i - 1, j) + C(i - 1, j + 1)) * dy
             H_0  = (C(i,     j - 1) + C(i,     j) + C(i,     j + 1)) * dy
             H_p1 = (C(i + 1, j - 1) + C(i + 1, j) + C(i + 1, j + 1)) * dy

             W_m1 = (C(i - 1, j - 1) + C(i, j - 1) + C(i + 1, j - 1)) * dx
             W_0  = (C(i - 1, j    ) + C(i, j    ) + C(i + 1, j    )) * dx
             W_p1 = (C(i - 1, j + 1) + C(i, j + 1) + C(i + 1, j + 1)) * dx

             if (C(i, j + 1) < C(i, j - 1)) then; sgn_y = 1.0_dp; else; sgn_y = -1.0_dp; end if
             if (C(i + 1, j) < C(i - 1, j)) then; sgn_x = 1.0_dp; else; sgn_x = -1.0_dp; end if

             ! 6 candidats Height Functions
             cand_nx(1) = -sgn_y * (H_p1 - H_m1) / (2.0_dp * dx); cand_ny(1) = sgn_y * 1.0_dp
             cand_nx(2) = -sgn_y * (H_p1 - H_0)  / (dx);          cand_ny(2) = sgn_y * 1.0_dp
             cand_nx(3) = -sgn_y * (H_0  - H_m1) / (dx);          cand_ny(3) = sgn_y * 1.0_dp

             cand_nx(4) = sgn_x * 1.0_dp; cand_ny(4) = -sgn_x * (W_p1 - W_m1) / (2.0_dp * dy)
             cand_nx(5) = sgn_x * 1.0_dp; cand_ny(5) = -sgn_x * (W_p1 - W_0)  / (dy)
             cand_nx(6) = sgn_x * 1.0_dp; cand_ny(6) = -sgn_x * (W_0  - W_m1) / (dy)

             ! Candidat 7 : Youngs
             mx = -((C(i+1,j+1) + 2.0_dp*C(i+1,j) + C(i+1,j-1)) - (C(i-1,j+1) + 2.0_dp*C(i-1,j) + C(i-1,j-1))) / (8.0_dp * dx)
             my = -((C(i+1,j+1) + 2.0_dp*C(i,j+1) + C(i-1,j+1)) - (C(i+1,j-1) + 2.0_dp*C(i,j-1) + C(i-1,j-1))) / (8.0_dp * dy)
             cand_nx(7) = mx; cand_ny(7) = my

             best_err = 1.0e30_dp
             best_nx = 1.0_dp; best_ny = 0.0_dp; best_alpha = 0.0_dp

             do k = 1, 7
                cnorm = hypot(cand_nx(k), cand_ny(k))
                if (cnorm < 1.0e-12_dp) cycle
                cnx = cand_nx(k) / cnorm
                cny = cand_ny(k) / cnorm

                calpha = compute_alpha(xmin, xmax, ymin, ymax, cnx, cny, c_val)

                err = 0.0_dp
                do dj = -1, 1
                   do di = -1, 1
                      if (di == 0 .and. dj == 0) cycle
                      v_xmin = x_face(mesh, i + di)
                      v_xmax = x_face(mesh, i + di + 1)
                      v_ymin = y_face(mesh, j + dj)
                      v_ymax = y_face(mesh, j + dj + 1)

                      v_target = C(i + di, j + dj) * dx * dy
                      v_recon  = cut_box_area(v_xmin, v_xmax, v_ymin, v_ymax, cnx, cny, calpha)
                      diff = v_recon - v_target
                      err = err + diff * diff
                   end do
                end do

                if (err < best_err) then
                   best_err = err
                   best_nx = cnx
                   best_ny = cny
                   best_alpha = calpha
                end if
             end do

             segments(i, j)%nx = best_nx
             segments(i, j)%ny = best_ny
             segments(i, j)%alpha = best_alpha
          else
             segments(i, j)%nx = 1.0_dp
             segments(i, j)%ny = 0.0_dp
             segments(i, j)%alpha = compute_alpha(xmin, xmax, ymin, ymax, 1.0_dp, 0.0_dp, c_val)
          end if
       end do
    end do
  end subroutine reconstruct_elvira

  ! [Procédure/Module] subroutine reconstruct_lvira
  subroutine reconstruct_lvira(mesh, C, segments)
    type(mesh_2d_t), intent(in) :: mesh
    real(dp), intent(in) :: C(0:mesh%nx_local+1, 0:mesh%ny_local+1)
    type(interface_segment_t), intent(out) :: segments(0:mesh%nx_local+1, 0:mesh%ny_local+1)

    integer :: i, j, nx, ny, iter, di, dj
    real(dp) :: c_val, dx, dy, xmin, xmax, ymin, ymax, theta0
    real(dp) :: invphi, invphi2, a_pt, b_pt, h_pt, c_pt, d_pt, fc, fd
    real(dp) :: best_theta, best_nx, best_ny, best_alpha
    real(dp) :: v_xmin, v_xmax, v_ymin, v_ymax, v_target, v_recon, diff

    ! Étape 1 : Initialisation par ELVIRA
    call reconstruct_elvira(mesh, C, segments)

    nx = mesh%nx_local
    ny = mesh%ny_local
    dx = mesh%dx
    dy = mesh%dy
    invphi = (sqrt(5.0_dp) - 1.0_dp) / 2.0_dp
    invphi2 = (3.0_dp - sqrt(5.0_dp)) / 2.0_dp

    do j = 1, ny
       do i = 1, nx
          c_val = C(i, j)
          if (c_val <= 1.0e-6_dp .or. c_val >= 1.0_dp - 1.0e-6_dp) cycle

          xmin = x_face(mesh, i)
          xmax = x_face(mesh, i + 1)
          ymin = y_face(mesh, j)
          ymax = y_face(mesh, j + 1)

          theta0 = atan2(segments(i, j)%ny, segments(i, j)%nx)
          a_pt = theta0 - 0.5_dp
          b_pt = theta0 + 0.5_dp
          h_pt = b_pt - a_pt

          c_pt = a_pt + invphi2 * h_pt
          d_pt = a_pt + invphi * h_pt

          ! Évaluation fc
          fc = 0.0_dp
          do dj = -1, 1
             do di = -1, 1
                if (di == 0 .and. dj == 0) cycle
                v_xmin = x_face(mesh, i + di); v_xmax = x_face(mesh, i + di + 1)
                v_ymin = y_face(mesh, j + dj); v_ymax = y_face(mesh, j + dj + 1)
                v_target = C(i + di, j + dj) * dx * dy
                v_recon  = cut_box_area(v_xmin, v_xmax, v_ymin, v_ymax, cos(c_pt), sin(c_pt), &
                                       compute_alpha(xmin, xmax, ymin, ymax, cos(c_pt), sin(c_pt), c_val))
                diff = v_recon - v_target
                fc = fc + diff * diff
             end do
          end do

          ! Évaluation fd
          fd = 0.0_dp
          do dj = -1, 1
             do di = -1, 1
                if (di == 0 .and. dj == 0) cycle
                v_xmin = x_face(mesh, i + di); v_xmax = x_face(mesh, i + di + 1)
                v_ymin = y_face(mesh, j + dj); v_ymax = y_face(mesh, j + dj + 1)
                v_target = C(i + di, j + dj) * dx * dy
                v_recon  = cut_box_area(v_xmin, v_xmax, v_ymin, v_ymax, cos(d_pt), sin(d_pt), &
                                       compute_alpha(xmin, xmax, ymin, ymax, cos(d_pt), sin(d_pt), c_val))
                diff = v_recon - v_target
                fd = fd + diff * diff
             end do
          end do

          do iter = 1, 10
             if (fc < fd) then
                b_pt = d_pt; d_pt = c_pt; fd = fc
                h_pt = invphi * h_pt
                c_pt = a_pt + invphi2 * h_pt
                fc = 0.0_dp
                do dj = -1, 1
                   do di = -1, 1
                      if (di == 0 .and. dj == 0) cycle
                      v_xmin = x_face(mesh, i + di); v_xmax = x_face(mesh, i + di + 1)
                      v_ymin = y_face(mesh, j + dj); v_ymax = y_face(mesh, j + dj + 1)
                      v_target = C(i + di, j + dj) * dx * dy
                      v_recon  = cut_box_area(v_xmin, v_xmax, v_ymin, v_ymax, cos(c_pt), sin(c_pt), &
                                             compute_alpha(xmin, xmax, ymin, ymax, cos(c_pt), sin(c_pt), c_val))
                      diff = v_recon - v_target
                      fc = fc + diff * diff
                   end do
                end do
             else
                a_pt = c_pt; c_pt = d_pt; fc = fd
                h_pt = invphi * h_pt
                d_pt = a_pt + invphi * h_pt
                fd = 0.0_dp
                do dj = -1, 1
                   do di = -1, 1
                      if (di == 0 .and. dj == 0) cycle
                      v_xmin = x_face(mesh, i + di); v_xmax = x_face(mesh, i + di + 1)
                      v_ymin = y_face(mesh, j + dj); v_ymax = y_face(mesh, j + dj + 1)
                      v_target = C(i + di, j + dj) * dx * dy
                      v_recon  = cut_box_area(v_xmin, v_xmax, v_ymin, v_ymax, cos(d_pt), sin(d_pt), &
                                             compute_alpha(xmin, xmax, ymin, ymax, cos(d_pt), sin(d_pt), c_val))
                      diff = v_recon - v_target
                      fd = fd + diff * diff
                   end do
                end do
             end if
          end do

          best_theta = 0.5_dp * (a_pt + b_pt)
          best_nx = cos(best_theta)
          best_ny = sin(best_theta)
          best_alpha = compute_alpha(xmin, xmax, ymin, ymax, best_nx, best_ny, c_val)

          segments(i, j)%nx = best_nx
          segments(i, j)%ny = best_ny
          segments(i, j)%alpha = best_alpha
       end do
    end do
  end subroutine reconstruct_lvira

  ! [Procédure/Module] subroutine reconstruct_interface
  subroutine reconstruct_interface(method_name, mesh, C, segments)
    character(len=*), intent(in) :: method_name
    type(mesh_2d_t), intent(in) :: mesh
    real(dp), intent(in) :: C(0:mesh%nx_local+1, 0:mesh%ny_local+1)
    type(interface_segment_t), intent(out) :: segments(0:mesh%nx_local+1, 0:mesh%ny_local+1)

    integer :: i, j, nx, ny, i_left_src, i_right_src, j_bot_src, j_top_src
    real(dp) :: c_val

    select case (trim(method_name))
    case ("SLIC")
       call reconstruct_slic(mesh, C, segments)
    case ("LVIRA")
       call reconstruct_lvira(mesh, C, segments)
    case ("ELVIRA")
       call reconstruct_elvira(mesh, C, segments)
    case default
       call reconstruct_plic(mesh, C, segments)
    end select

    nx = mesh%nx_local
    ny = mesh%ny_local
    if (.not. mesh%is_boundary_left .and. mesh%size <= 1) then
       i_left_src = nx
    else
       i_left_src = 1
    end if
    if (.not. mesh%is_boundary_right .and. mesh%size <= 1) then
       i_right_src = 1
    else
       i_right_src = nx
    end if
    if (.not. mesh%is_boundary_bottom .and. mesh%size <= 1) then
       j_bot_src = ny
    else
       j_bot_src = 1
    end if
    if (.not. mesh%is_boundary_top .and. mesh%size <= 1) then
       j_top_src = 1
    else
       j_top_src = ny
    end if

    do j = 1, ny
       c_val = C(0, j)
       if (c_val > 1.0e-6_dp .and. c_val < 1.0_dp - 1.0e-6_dp) then
          segments(0, j)%nx = segments(i_left_src, j)%nx
          segments(0, j)%ny = segments(i_left_src, j)%ny
          segments(0, j)%alpha = compute_alpha(x_face(mesh, 0), x_face(mesh, 1), &
                                               y_face(mesh, j), y_face(mesh, j + 1), &
                                               segments(0, j)%nx, segments(0, j)%ny, c_val)
       end if
       c_val = C(nx + 1, j)
       if (c_val > 1.0e-6_dp .and. c_val < 1.0_dp - 1.0e-6_dp) then
          segments(nx + 1, j)%nx = segments(i_right_src, j)%nx
          segments(nx + 1, j)%ny = segments(i_right_src, j)%ny
          segments(nx + 1, j)%alpha = compute_alpha(x_face(mesh, nx + 1), x_face(mesh, nx + 2), &
                                                    y_face(mesh, j), y_face(mesh, j + 1), &
                                                    segments(nx + 1, j)%nx, segments(nx + 1, j)%ny, c_val)
       end if
    end do

    do i = 0, nx + 1
       c_val = C(i, 0)
       if (c_val > 1.0e-6_dp .and. c_val < 1.0_dp - 1.0e-6_dp) then
          segments(i, 0)%nx = segments(i, j_bot_src)%nx
          segments(i, 0)%ny = segments(i, j_bot_src)%ny
          segments(i, 0)%alpha = compute_alpha(x_face(mesh, i), x_face(mesh, i + 1), &
                                               y_face(mesh, 0), y_face(mesh, 1), &
                                               segments(i, 0)%nx, segments(i, 0)%ny, c_val)
       end if
       c_val = C(i, ny + 1)
       if (c_val > 1.0e-6_dp .and. c_val < 1.0_dp - 1.0e-6_dp) then
          segments(i, ny + 1)%nx = segments(i, j_top_src)%nx
          segments(i, ny + 1)%ny = segments(i, j_top_src)%ny
          segments(i, ny + 1)%alpha = compute_alpha(x_face(mesh, i), x_face(mesh, i + 1), &
                                                    y_face(mesh, ny + 1), y_face(mesh, ny + 2), &
                                                    segments(i, ny + 1)%nx, segments(i, ny + 1)%ny, c_val)
       end if
    end do
  end subroutine reconstruct_interface

  ! [Procédure/Module] function evaluate_limiter
  function evaluate_limiter(limiter_name, r, c_donor, c_upwind, c_acceptor) result(phi)
    character(len=*), intent(in) :: limiter_name
    real(dp), intent(in) :: r, c_donor, c_upwind, c_acceptor
    real(dp) :: phi
    real(dp) :: denom, c_tilde_D, c_tilde_f

    select case (trim(limiter_name))
    case ("VAN_LEER")
       if (r > 0.0_dp) then
          phi = 2.0_dp * r / (1.0_dp + r)
       else
          phi = 0.0_dp
       end if
    case ("MUSCL")
       phi = max(0.0_dp, min(2.0_dp * r, min(0.5_dp * (1.0_dp + r), 2.0_dp)))
    case ("HRIC")
       denom = c_acceptor - c_upwind
       if (abs(denom) < 1.0e-12_dp) then
          phi = 0.0_dp
          return
       end if
       c_tilde_D = (c_donor - c_upwind) / denom
       c_tilde_f = c_tilde_D
       if (c_tilde_D >= 0.0_dp .and. c_tilde_D <= 1.0_dp) then
          if (c_tilde_D < 0.5_dp) then
             c_tilde_f = 2.0_dp * c_tilde_D
          else
             c_tilde_f = 1.0_dp
          end if
       end if
       if (abs(1.0_dp - c_tilde_D) > 1.0e-12_dp) then
          phi = 2.0_dp * (c_tilde_f - c_tilde_D) / (1.0_dp - c_tilde_D)
          phi = max(0.0_dp, min(2.0_dp, phi))
       else
          phi = 1.0_dp
       end if
    case default ! "SUPERBEE"
       phi = max(0.0_dp, max(min(2.0_dp * r, 1.0_dp), min(r, 2.0_dp)))
    end select
  end function evaluate_limiter

end module module_interface_reconstruction
