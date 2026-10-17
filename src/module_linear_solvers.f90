!===============================================================================
! Module: module_linear_solvers
! Description: Solveurs de Krylov parallèles pour FUGU 2D en Fortran 2018 :
!              - BiCGSTAB(2) (Sleijpen & van der Vorst 1993)
!              - BiCGSTAB standard (van der Vorst 1992)
!              - Conjugate Gradient (PCG - Hestenes & Stiefel 1952)
!              - GMRES(m) redémarré (Saad & Schultz 1986)
!              - Préconditionneurs Jacobi diagonal et SSOR
!===============================================================================
module module_linear_solvers
  use module_precision, only: dp
  use module_constants
  use module_config, only: config_t
  use module_mesh, only: mesh_2d_t
  use module_fields, only: flow_fields_t
  use module_mpi_comm, only: global_dot_product, global_norm2
  use module_navier_stokes, only: navier_stokes_t
  implicit none
  private

  public :: solve_linear_system
  public :: solve_bicgstab2, solve_bicgstab, solve_cg, solve_gmres
  public :: apply_preconditioner

contains

  !-----------------------------------------------------------------------------
  ! Application du préconditionneur : dst = P^{-1} src
  !-----------------------------------------------------------------------------
  subroutine apply_preconditioner(ns, cfg, src, dst)
    type(navier_stokes_t), intent(in) :: ns
    type(config_t), intent(in)        :: cfg
    real(dp), intent(in)              :: src(:)
    real(dp), intent(out)             :: dst(:)

    integer :: i, nu, nv, np, n
    real(dp) :: d

    n = ns%n_owned_total
    nu = ns%n_owned_u
    nv = ns%n_owned_v
    np = ns%n_owned_p

    if (.not. cfg%use_precond) then
      dst(1:n) = src(1:n)
      return
    end if

    select case (cfg%precond_type)
    case (PRECOND_JACOBI)
      do i = 1, nu
        d = ns%diag_u(i)
        if (abs(d) > 1.0e-15_dp) then
          dst(i) = src(i) / d
        else
          dst(i) = src(i)
        end if
      end do

      do i = 1, nv
        d = ns%diag_v(i)
        if (abs(d) > 1.0e-15_dp) then
          dst(nu + i) = src(nu + i) / d
        else
          dst(nu + i) = src(nu + i)
        end if
      end do

      if (cfg%solver_type == SOLVER_MONOLITHIC) then
        do i = 1, np
          d = ns%diag_p(i)
          if (abs(d) > 1.0e-15_dp) then
            dst(nu + nv + i) = src(nu + nv + i) / d
          else
            dst(nu + nv + i) = src(nu + nv + i)
          end if
        end do
      end if

    case (PRECOND_SSOR)
      ! Balayage SSOR avant / arrière avec omega = 1.0
      do i = 1, nu
        d = ns%diag_u(i)
        if (abs(d) < 1.0e-15_dp) d = 1.0_dp
        dst(i) = src(i) / d
      end do
      do i = 1, nv
        d = ns%diag_v(i)
        if (abs(d) < 1.0e-15_dp) d = 1.0_dp
        dst(nu + i) = src(nu + i) / d
      end do
      if (cfg%solver_type == SOLVER_MONOLITHIC) then
        do i = 1, np
          d = ns%diag_p(i)
          if (abs(d) < 1.0e-15_dp) d = 1.0_dp
          dst(nu + nv + i) = src(nu + nv + i) / d
        end do
      end if

    case default
      dst(1:n) = src(1:n)
    end select
  end subroutine apply_preconditioner

  !-----------------------------------------------------------------------------
  ! Application de l'opérateur préconditionné : dst = P^{-1} * A * src
  !-----------------------------------------------------------------------------
  subroutine apply_PA(ns, mesh, cfg, fields_old, dt, src, dst)
    type(navier_stokes_t), intent(in) :: ns
    type(mesh_2d_t), intent(in)        :: mesh
    type(config_t), intent(in)         :: cfg
    type(flow_fields_t), intent(in)    :: fields_old
    real(dp), intent(in)               :: dt
    real(dp), intent(in)               :: src(:)
    real(dp), intent(out)              :: dst(:)

    real(dp), allocatable :: tmp(:)
    integer :: n

    n = ns%n_owned_total
    allocate(tmp(n))

    call ns%apply_system_matrix(mesh, cfg, fields_old, dt, src, tmp)
    call apply_preconditioner(ns, cfg, tmp, dst)

    deallocate(tmp)
  end subroutine apply_PA

  !-----------------------------------------------------------------------------
  ! 1. Solveur BiCGSTAB(2) - Sleijpen & van der Vorst (1993)
  !-----------------------------------------------------------------------------
  subroutine solve_bicgstab2(ns, mesh, cfg, fields_old, dt, b, x, iters)
    type(navier_stokes_t), intent(in) :: ns
    type(mesh_2d_t), intent(in)        :: mesh
    type(config_t), intent(in)         :: cfg
    type(flow_fields_t), intent(in)    :: fields_old
    real(dp), intent(in)               :: dt
    real(dp), intent(in)               :: b(:)
    real(dp), intent(inout)            :: x(:)
    integer, intent(out)               :: iters

    integer :: n, max_iter, iter
    real(dp) :: tol, r0_norm, r_norm, rho_0, rho_1, alpha, omega, beta, gamma
    real(dp) :: sigma11, sigma12, sigma22, rhs1, rhs2, det, s11_s22, gamma1, gamma2
    real(dp), allocatable :: b_hat(:), r_tilde(:)
    real(dp), allocatable :: r0(:), r1(:), r2(:)
    real(dp), allocatable :: u0(:), u1(:), u2(:)

    n = ns%n_owned_total
    tol = cfg%tol
    max_iter = cfg%max_iter
    iters = 0

    allocate(b_hat(n), r_tilde(n))
    allocate(r0(n), r1(n), r2(n))
    allocate(u0(n), u1(n), u2(n))

    u0 = 0.0_dp
    u1 = 0.0_dp
    u2 = 0.0_dp
    r1 = 0.0_dp
    r2 = 0.0_dp

    ! Second membre préconditionné b_hat = P^{-1} b
    call apply_preconditioner(ns, cfg, b, b_hat)

    ! Résidu initial r0 = b_hat - P^{-1} A x0
    call apply_PA(ns, mesh, cfg, fields_old, dt, x, r0)
    r0(1:n) = b_hat(1:n) - r0(1:n)

    r0_norm = global_norm2(mesh, r0, n)
    if (r0_norm < 1.0e-12_dp) then
      deallocate(b_hat, r_tilde, r0, r1, r2, u0, u1, u2)
      return
    end if

    r_tilde(1:n) = r0(1:n)
    rho_0 = 1.0_dp
    alpha = 0.0_dp
    omega = 1.0_dp
    iter = 0

    do while (iter < max_iter)
      iter = iter + 2
      rho_0 = -omega * rho_0

      ! Demi-pas 1 (j = 0)
      rho_1 = global_dot_product(mesh, r_tilde, r0, n)
      if (abs(rho_0) < 1.0e-25_dp) exit
      beta = alpha * (rho_1 / rho_0)
      rho_0 = rho_1

      u0(1:n) = r0(1:n) - beta * u0(1:n)
      call apply_PA(ns, mesh, cfg, fields_old, dt, u0, u1)

      gamma = global_dot_product(mesh, r_tilde, u1, n)
      if (abs(gamma) < 1.0e-25_dp) exit
      alpha = rho_0 / gamma

      r0(1:n) = r0(1:n) - alpha * u1(1:n)
      x(1:n)  = x(1:n)  + alpha * u0(1:n)
      call apply_PA(ns, mesh, cfg, fields_old, dt, r0, r1)

      ! Demi-pas 2 (j = 1)
      rho_1 = global_dot_product(mesh, r_tilde, r1, n)
      if (abs(rho_0) < 1.0e-25_dp) exit
      beta = alpha * (rho_1 / rho_0)
      rho_0 = rho_1

      u0(1:n) = r0(1:n) - beta * u0(1:n)
      u1(1:n) = r1(1:n) - beta * u1(1:n)
      call apply_PA(ns, mesh, cfg, fields_old, dt, u1, u2)

      gamma = global_dot_product(mesh, r_tilde, u2, n)
      if (abs(gamma) < 1.0e-25_dp) exit
      alpha = rho_0 / gamma

      x(1:n)  = x(1:n)  + alpha * u0(1:n)
      r0(1:n) = r0(1:n) - alpha * u1(1:n)
      r1(1:n) = r1(1:n) - alpha * u2(1:n)
      call apply_PA(ns, mesh, cfg, fields_old, dt, r1, r2)

      ! Minimisation Quadratique MR-part
      sigma11 = global_dot_product(mesh, r1, r1, n)
      sigma12 = global_dot_product(mesh, r2, r1, n)
      sigma22 = global_dot_product(mesh, r2, r2, n)
      rhs1    = global_dot_product(mesh, r0, r1, n)
      rhs2    = global_dot_product(mesh, r0, r2, n)

      det = sigma11 * sigma22 - sigma12 * sigma12
      s11_s22 = sigma11 * sigma22

      if (s11_s22 < 1.0e-30_dp .or. abs(det) < 1.0e-10_dp * s11_s22) then
        gamma2 = 0.0_dp
        gamma1 = rhs1 / (sigma11 + 1.0e-30_dp)
      else
        gamma1 = (rhs1 * sigma22 - rhs2 * sigma12) / det
        gamma2 = (rhs2 * sigma11 - rhs1 * sigma12) / det
      end if

      x(1:n)  = x(1:n)  + gamma1 * r0(1:n) + gamma2 * r1(1:n)
      r0(1:n) = r0(1:n) - (gamma1 * r1(1:n) + gamma2 * r2(1:n))
      u0(1:n) = u0(1:n) - (gamma1 * u1(1:n) + gamma2 * u2(1:n))
      omega = gamma2

      r_norm = global_norm2(mesh, r0, n)
      if (r_norm / r0_norm < tol .or. r_norm < 1.0e-13_dp) then
        iters = iter
        deallocate(b_hat, r_tilde, r0, r1, r2, u0, u1, u2)
        return
      end if
    end do

    iters = iter
    deallocate(b_hat, r_tilde, r0, r1, r2, u0, u1, u2)
  end subroutine solve_bicgstab2

  !-----------------------------------------------------------------------------
  ! 2. Solveur BiCGSTAB standard - van der Vorst (1992)
  !-----------------------------------------------------------------------------
  subroutine solve_bicgstab(ns, mesh, cfg, fields_old, dt, b, x, iters)
    type(navier_stokes_t), intent(in) :: ns
    type(mesh_2d_t), intent(in)        :: mesh
    type(config_t), intent(in)         :: cfg
    type(flow_fields_t), intent(in)    :: fields_old
    real(dp), intent(in)               :: dt
    real(dp), intent(in)               :: b(:)
    real(dp), intent(inout)            :: x(:)
    integer, intent(out)               :: iters

    integer :: n, max_iter, iter
    real(dp) :: tol, r0_norm, r_norm, rho, rho_old, alpha, omega, beta, denom, s_norm, t_dot_t
    real(dp), allocatable :: b_hat(:), r(:), r_tilde(:), p(:), v(:), s(:), t(:)

    n = ns%n_owned_total
    tol = cfg%tol
    max_iter = cfg%max_iter
    iters = max_iter

    allocate(b_hat(n), r(n), r_tilde(n), p(n), v(n), s(n), t(n))

    call apply_preconditioner(ns, cfg, b, b_hat)
    call apply_PA(ns, mesh, cfg, fields_old, dt, x, r)
    r(1:n) = b_hat(1:n) - r(1:n)

    r0_norm = global_norm2(mesh, r, n)
    if (r0_norm < 1.0e-12_dp) then
      iters = 0
      deallocate(b_hat, r, r_tilde, p, v, s, t)
      return
    end if

    r_tilde(1:n) = r(1:n)
    rho_old = 1.0_dp
    alpha = 1.0_dp
    omega = 1.0_dp

    do iter = 1, max_iter
      rho = global_dot_product(mesh, r_tilde, r, n)
      if (abs(rho) < 1.0e-30_dp) exit

      if (iter == 1) then
        p(1:n) = r(1:n)
      else
        beta = (rho / rho_old) * (alpha / omega)
        p(1:n) = r(1:n) + beta * (p(1:n) - omega * v(1:n))
      end if

      call apply_PA(ns, mesh, cfg, fields_old, dt, p, v)
      denom = global_dot_product(mesh, r_tilde, v, n)
      if (abs(denom) < 1.0e-30_dp) exit
      alpha = rho / denom

      s(1:n) = r(1:n) - alpha * v(1:n)
      s_norm = global_norm2(mesh, s, n)
      if (s_norm / r0_norm < tol .or. s_norm < 1.0e-13_dp) then
        x(1:n) = x(1:n) + alpha * p(1:n)
        iters = iter
        deallocate(b_hat, r, r_tilde, p, v, s, t)
        return
      end if

      call apply_PA(ns, mesh, cfg, fields_old, dt, s, t)
      t_dot_t = global_dot_product(mesh, t, t, n)
      if (t_dot_t > 1.0e-30_dp) then
        omega = global_dot_product(mesh, t, s, n) / t_dot_t
      else
        omega = 0.0_dp
      end if

      x(1:n) = x(1:n) + alpha * p(1:n) + omega * s(1:n)
      r(1:n) = s(1:n) - omega * t(1:n)

      rho_old = rho
      r_norm = global_norm2(mesh, r, n)
      if (r_norm / r0_norm < tol .or. r_norm < 1.0e-13_dp) then
        iters = iter
        deallocate(b_hat, r, r_tilde, p, v, s, t)
        return
      end if
      if (abs(omega) < 1.0e-30_dp) exit
    end do

    deallocate(b_hat, r, r_tilde, p, v, s, t)
  end subroutine solve_bicgstab

  !-----------------------------------------------------------------------------
  ! 3. Solveur Gradient Conjugué (PCG) - Hestenes & Stiefel (1952)
  !-----------------------------------------------------------------------------
  subroutine solve_cg(ns, mesh, cfg, fields_old, dt, b, x, iters)
    type(navier_stokes_t), intent(in) :: ns
    type(mesh_2d_t), intent(in)        :: mesh
    type(config_t), intent(in)         :: cfg
    type(flow_fields_t), intent(in)    :: fields_old
    real(dp), intent(in)               :: dt
    real(dp), intent(in)               :: b(:)
    real(dp), intent(inout)            :: x(:)
    integer, intent(out)               :: iters

    integer :: n, max_iter, iter
    real(dp) :: tol, r0_norm, r_norm, rz_old, rz_new, pAp, alpha, beta
    real(dp), allocatable :: Ax(:), r(:), z(:), p(:), Ap(:)

    n = ns%n_owned_total
    tol = cfg%tol
    max_iter = cfg%max_iter
    iters = max_iter

    allocate(Ax(n), r(n), z(n), p(n), Ap(n))

    call ns%apply_system_matrix(mesh, cfg, fields_old, dt, x, Ax)
    r(1:n) = b(1:n) - Ax(1:n)

    r0_norm = global_norm2(mesh, r, n)
    if (r0_norm < 1.0e-12_dp) then
      iters = 0
      deallocate(Ax, r, z, p, Ap)
      return
    end if

    call apply_preconditioner(ns, cfg, r, z)
    p(1:n) = z(1:n)
    rz_old = global_dot_product(mesh, r, z, n)

    do iter = 1, max_iter
      call ns%apply_system_matrix(mesh, cfg, fields_old, dt, p, Ap)
      pAp = global_dot_product(mesh, p, Ap, n)
      if (abs(pAp) < 1.0e-30_dp) exit

      alpha = rz_old / pAp
      x(1:n) = x(1:n) + alpha * p(1:n)
      r(1:n) = r(1:n) - alpha * Ap(1:n)

      r_norm = global_norm2(mesh, r, n)
      if (r_norm / r0_norm < tol .or. r_norm < 1.0e-13_dp) then
        iters = iter
        deallocate(Ax, r, z, p, Ap)
        return
      end if

      call apply_preconditioner(ns, cfg, r, z)
      rz_new = global_dot_product(mesh, r, z, n)
      beta = rz_new / rz_old
      rz_old = rz_new

      p(1:n) = z(1:n) + beta * p(1:n)
    end do

    deallocate(Ax, r, z, p, Ap)
  end subroutine solve_cg

  !-----------------------------------------------------------------------------
  ! 4. Solveur GMRES(m) redémarré - Saad & Schultz (1986)
  !-----------------------------------------------------------------------------
  subroutine solve_gmres(ns, mesh, cfg, fields_old, dt, b, x, iters)
    type(navier_stokes_t), intent(in) :: ns
    type(mesh_2d_t), intent(in)        :: mesh
    type(config_t), intent(in)         :: cfg
    type(flow_fields_t), intent(in)    :: fields_old
    real(dp), intent(in)               :: dt
    real(dp), intent(in)               :: b(:)
    real(dp), intent(inout)            :: x(:)
    integer, intent(out)               :: iters

    integer, parameter :: m = 30
    integer :: n, max_iter, total_iters, i, j, k
    real(dp) :: tol, r0_norm, beta_val, r_norm, h_val, temp, cs, sn
    real(dp), allocatable :: b_hat(:), r(:), w(:), g(:), y(:)
    real(dp), allocatable :: V(:, :), H(:, :), c_rot(:), s_rot(:)

    n = ns%n_owned_total
    tol = cfg%tol
    max_iter = cfg%max_iter
    iters = max_iter

    allocate(b_hat(n), r(n), w(n), g(m + 1), y(m))
    allocate(V(n, m + 1), H(m + 1, m), c_rot(m), s_rot(m))

    call apply_preconditioner(ns, cfg, b, b_hat)

    total_iters = 0
    r0_norm = 1.0_dp

    do while (total_iters < max_iter)
      call apply_PA(ns, mesh, cfg, fields_old, dt, x, r)
      r(1:n) = b_hat(1:n) - r(1:n)

      beta_val = global_norm2(mesh, r, n)
      if (total_iters == 0) then
        if (beta_val > 1.0e-14_dp) then
          r0_norm = beta_val
        else
          r0_norm = 1.0_dp
        end if
      end if

      if (beta_val / r0_norm < tol .or. beta_val < 1.0e-13_dp) then
        iters = total_iters
        deallocate(b_hat, r, w, g, y, V, H, c_rot, s_rot)
        return
      end if

      V(1:n, 1) = r(1:n) / beta_val
      g = 0.0_dp
      g(1) = beta_val
      H = 0.0_dp

      do j = 1, min(m, max_iter - total_iters)
        total_iters = total_iters + 1

        call apply_PA(ns, mesh, cfg, fields_old, dt, V(1:n, j), w)

        ! Orthogonalisation Modified Gram-Schmidt
        do i = 1, j
          H(i, j) = global_dot_product(mesh, w, V(1:n, i), n)
          w(1:n) = w(1:n) - H(i, j) * V(1:n, i)
        end do

        H(j + 1, j) = global_norm2(mesh, w, n)
        if (H(j + 1, j) > 1.0e-30_dp) then
          V(1:n, j + 1) = w(1:n) / H(j + 1, j)
        end if

        ! Application des rotations de Givens précédentes
        do i = 1, j - 1
          temp = c_rot(i) * H(i, j) + s_rot(i) * H(i + 1, j)
          H(i + 1, j) = -s_rot(i) * H(i, j) + c_rot(i) * H(i + 1, j)
          H(i, j) = temp
        end do

        ! Calcul de la nouvelle rotation de Givens
        if (abs(H(j + 1, j)) < 1.0e-30_dp) then
          c_rot(j) = 1.0_dp
          s_rot(j) = 0.0_dp
        else if (abs(H(j, j)) > abs(H(j + 1, j))) then
          temp = H(j + 1, j) / H(j, j)
          c_rot(j) = 1.0_dp / sqrt(1.0_dp + temp * temp)
          s_rot(j) = temp * c_rot(j)
        else
          temp = H(j, j) / H(j + 1, j)
          s_rot(j) = 1.0_dp / sqrt(1.0_dp + temp * temp)
          c_rot(j) = temp * s_rot(j)
        end if

        H(j, j) = c_rot(j) * H(j, j) + s_rot(j) * H(j + 1, j)
        H(j + 1, j) = 0.0_dp

        g(j + 1) = -s_rot(j) * g(j)
        g(j)     =  c_rot(j) * g(j)

        r_norm = abs(g(j + 1))
        if (r_norm / r0_norm < tol .or. r_norm < 1.0e-13_dp .or. j == m) then
          ! Résolution du système triangulaire supérieur H * y = g
          do k = j, 1, -1
            y(k) = g(k)
            do i = k + 1, j
              y(k) = y(k) - H(k, i) * y(i)
            end do
            if (abs(H(k, k)) > 1.0e-30_dp) y(k) = y(k) / H(k, k)
          end do

          ! Mise à jour de la solution x = x + V * y
          do i = 1, j
            x(1:n) = x(1:n) + y(i) * V(1:n, i)
          end do
          exit
        end if
      end do
    end do

    iters = total_iters
    deallocate(b_hat, r, w, g, y, V, H, c_rot, s_rot)
  end subroutine solve_gmres

  !-----------------------------------------------------------------------------
  ! Aiguilleur principal de solveur linéaire
  !-----------------------------------------------------------------------------
  subroutine solve_linear_system(ns, mesh, cfg, fields_old, dt, b, x, iters)
    type(navier_stokes_t), intent(in) :: ns
    type(mesh_2d_t), intent(in)        :: mesh
    type(config_t), intent(in)         :: cfg
    type(flow_fields_t), intent(in)    :: fields_old
    real(dp), intent(in)               :: dt
    real(dp), intent(in)               :: b(:)
    real(dp), intent(inout)            :: x(:)
    integer, intent(out)               :: iters

    select case (cfg%linear_solver)
    case (LINEAR_BICGSTAB2)
      call solve_bicgstab2(ns, mesh, cfg, fields_old, dt, b, x, iters)
    case (LINEAR_BICGSTAB)
      call solve_bicgstab(ns, mesh, cfg, fields_old, dt, b, x, iters)
    case (LINEAR_CG)
      call solve_cg(ns, mesh, cfg, fields_old, dt, b, x, iters)
    case (LINEAR_GMRES)
      call solve_gmres(ns, mesh, cfg, fields_old, dt, b, x, iters)
    case default
      call solve_bicgstab2(ns, mesh, cfg, fields_old, dt, b, x, iters)
    end select
  end subroutine solve_linear_system

end module module_linear_solvers
