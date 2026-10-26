!===============================================================================
! Module: module_linear_solvers
! Description: Solveurs de Krylov parallèles pour COVOF 2D en Fortran 2018 :
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
  use module_navier_stokes, only: navier_stokes_t, is_u_face_in_obstacle, is_v_face_in_obstacle
  implicit none
  private

  public :: solve_linear_system
  public :: solve_bicgstab2, solve_bicgstab, solve_cg, solve_gmres
  public :: apply_preconditioner

contains

  !-----------------------------------------------------------------------------
  ! Application du préconditionneur : dst = P^{-1} src
  !-----------------------------------------------------------------------------
  subroutine apply_preconditioner(ns, mesh, cfg, src, dst)
    ! Déclaration de variable locale / paramètre de calcul
    type(navier_stokes_t), intent(in) :: ns
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)        :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)        :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)              :: src(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(out)             :: dst(:)

    ! Déclaration de variable locale / paramètre de calcul
    integer :: i, j, dof, nu, nv, np, n, u_nx, v_ny, sweep
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: d, val, grad_p
    ! Déclaration de variable locale / paramètre de calcul
    logical :: is_wall
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: rp_2d(:, :), zp_2d(:, :)
    real(dp), allocatable :: ru(:), rv(:), rhs_u(:), rhs_v(:)

    ! Affectation et calcul numérique
    n = ns%n_owned_total
    ! Affectation et calcul numérique
    nu = ns%n_owned_u
    ! Affectation et calcul numérique
    nv = ns%n_owned_v
    ! Affectation et calcul numérique
    np = ns%n_owned_p

    ! Évaluation de la condition logique
    if (.not. cfg%use_precond) then
      ! Affectation et calcul numérique
      dst(1:n) = src(1:n)
      ! Retour à la procédure appelante
      return
    end if

    select case (cfg%precond_type)
    case (PRECOND_JACOBI)
      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 1, nu
        ! Affectation et calcul numérique
        d = ns%diag_u(i)
        ! Évaluation de la condition logique
        if (abs(d) > 1.0e-15_dp) then
          ! Affectation et calcul numérique
          dst(i) = src(i) / d
        ! Branche par défaut si condition non satisfaite
        else
          ! Affectation et calcul numérique
          dst(i) = src(i)
        end if
      end do

      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 1, nv
        ! Affectation et calcul numérique
        d = ns%diag_v(i)
        ! Évaluation de la condition logique
        if (abs(d) > 1.0e-15_dp) then
          ! Affectation et calcul numérique
          dst(nu + i) = src(nu + i) / d
        ! Branche par défaut si condition non satisfaite
        else
          ! Affectation et calcul numérique
          dst(nu + i) = src(nu + i)
        end if
      end do

      ! Évaluation de la condition logique
      if (cfg%solver_type == SOLVER_MONOLITHIC) then
        ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
        do i = 1, np
          ! Affectation et calcul numérique
          d = ns%diag_p(i)
          ! Évaluation de la condition logique
          if (abs(d) > 1.0e-15_dp) then
            ! Affectation et calcul numérique
            dst(nu + nv + i) = src(nu + nv + i) / d
          ! Branche par défaut si condition non satisfaite
          else
            ! Affectation et calcul numérique
            dst(nu + nv + i) = src(nu + nv + i)
          end if
        end do
      end if

    case (PRECOND_SSOR)
      ! Balayage SSOR avant / arrière avec omega = 1.0
      do i = 1, nu
        ! Affectation et calcul numérique
        d = ns%diag_u(i)
        ! Affectation et calcul numérique
        if (abs(d) < 1.0e-15_dp) d = 1.0_dp
        ! Affectation et calcul numérique
        dst(i) = src(i) / d
      end do
      ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
      do i = 1, nv
        ! Affectation et calcul numérique
        d = ns%diag_v(i)
        ! Affectation et calcul numérique
        if (abs(d) < 1.0e-15_dp) d = 1.0_dp
        ! Affectation et calcul numérique
        dst(nu + i) = src(nu + i) / d
      end do
      ! Évaluation de la condition logique
      if (cfg%solver_type == SOLVER_MONOLITHIC) then
        ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
        do i = 1, np
          ! Affectation et calcul numérique
          d = ns%diag_p(i)
          ! Affectation et calcul numérique
          if (abs(d) < 1.0e-15_dp) d = 1.0_dp
          ! Affectation et calcul numérique
          dst(nu + nv + i) = src(nu + nv + i) / d
        end do
      end if

    case (PRECOND_BLOCK_TRIANGULAR)
      if (cfg%solver_type == SOLVER_MONOLITHIC) then
        u_nx = mesh%nx_local
        if (mesh%is_boundary_right) u_nx = u_nx + 1
        v_ny = mesh%ny_local
        if (mesh%is_boundary_top) v_ny = v_ny + 1

        allocate(rp_2d(mesh%nx_local, mesh%ny_local))
        allocate(zp_2d(mesh%nx_local, mesh%ny_local))

        do j = 1, mesh%ny_local
          do i = 1, mesh%nx_local
            zp_2d(i, j) = src(nu + nv + (j - 1) * mesh%nx_local + i)
          end do
        end do

        rp_2d = 0.0_dp

        ! 1. 8 balayages de Gauss-Seidel symétrique sur le complément de Schur S
        do sweep = 1, 8
          ! Balayage avant
          do j = 1, mesh%ny_local
            do i = 1, mesh%nx_local
              val = zp_2d(i, j)
              if (i > 1) then
                val = val + ns%schur_aW(i, j) * rp_2d(i - 1, j)
              else if (cfg%bc_left%btype == BC_PERIODIC) then
                val = val + ns%schur_aW(i, j) * rp_2d(mesh%nx_local, j)
              end if

              if (i < mesh%nx_local) then
                val = val + ns%schur_aE(i, j) * rp_2d(i + 1, j)
              else if (cfg%bc_right%btype == BC_PERIODIC) then
                val = val + ns%schur_aE(i, j) * rp_2d(1, j)
              end if

              if (j > 1) then
                val = val + ns%schur_aS(i, j) * rp_2d(i, j - 1)
              else if (cfg%bc_bottom%btype == BC_PERIODIC) then
                val = val + ns%schur_aS(i, j) * rp_2d(i, mesh%ny_local)
              end if

              if (j < mesh%ny_local) then
                val = val + ns%schur_aN(i, j) * rp_2d(i, j + 1)
              else if (cfg%bc_top%btype == BC_PERIODIC) then
                val = val + ns%schur_aN(i, j) * rp_2d(i, 1)
              end if

              rp_2d(i, j) = val / ns%schur_aP(i, j)
            end do
          end do

          ! Balayage arrière
          do j = mesh%ny_local, 1, -1
            do i = mesh%nx_local, 1, -1
              val = zp_2d(i, j)
              if (i > 1) then
                val = val + ns%schur_aW(i, j) * rp_2d(i - 1, j)
              else if (cfg%bc_left%btype == BC_PERIODIC) then
                val = val + ns%schur_aW(i, j) * rp_2d(mesh%nx_local, j)
              end if

              if (i < mesh%nx_local) then
                val = val + ns%schur_aE(i, j) * rp_2d(i + 1, j)
              else if (cfg%bc_right%btype == BC_PERIODIC) then
                val = val + ns%schur_aE(i, j) * rp_2d(1, j)
              end if

              if (j > 1) then
                val = val + ns%schur_aS(i, j) * rp_2d(i, j - 1)
              else if (cfg%bc_bottom%btype == BC_PERIODIC) then
                val = val + ns%schur_aS(i, j) * rp_2d(i, mesh%ny_local)
              end if

              if (j < mesh%ny_local) then
                val = val + ns%schur_aN(i, j) * rp_2d(i, j + 1)
              else if (cfg%bc_top%btype == BC_PERIODIC) then
                val = val + ns%schur_aN(i, j) * rp_2d(i, 1)
              end if

              rp_2d(i, j) = val / ns%schur_aP(i, j)
            end do
          end do
        end do

        ! 2. Second membre pour vitesse u : rhs_u = zu - G_u * rp
        allocate(rhs_u(nu), ru(nu))
        rhs_u = 0.0_dp
        do j = 1, mesh%ny_local
          do i = 1, u_nx
            dof = (j - 1) * u_nx + i
            is_wall = (i == 1 .and. mesh%is_boundary_left) .or. &
                      (i == u_nx .and. mesh%is_boundary_right) .or. &
                      is_u_face_in_obstacle(mesh, cfg, i, j)
            if (is_wall) then
              rhs_u(dof) = src(dof)
            else
              grad_p = (rp_2d(i, j) - rp_2d(i - 1, j)) * mesh%dy
              rhs_u(dof) = src(dof) - grad_p
            end if
          end do
        end do

        ! 4 balayages de Gauss-Seidel symétrique sur Fu * ru = rhs_u
        ru = 0.0_dp
        do sweep = 1, 4
          ! Forward
          do j = 1, mesh%ny_local
            do i = 1, u_nx
              dof = (j - 1) * u_nx + i
              is_wall = (i == 1 .and. mesh%is_boundary_left) .or. &
                        (i == u_nx .and. mesh%is_boundary_right) .or. &
                        is_u_face_in_obstacle(mesh, cfg, i, j)
              if (is_wall) then
                ru(dof) = rhs_u(dof) / ns%u_aP(dof)
                cycle
              end if
              val = rhs_u(dof)
              if (i > 1) val = val + ns%u_aW(dof) * ru(dof - 1)
              if (i < u_nx) val = val + ns%u_aE(dof) * ru(dof + 1)
              if (j > 1) val = val + ns%u_aS(dof) * ru(dof - u_nx)
              if (j < mesh%ny_local) val = val + ns%u_aN(dof) * ru(dof + u_nx)
              ru(dof) = val / ns%u_aP(dof)
            end do
          end do
          ! Backward
          do j = mesh%ny_local, 1, -1
            do i = u_nx, 1, -1
              dof = (j - 1) * u_nx + i
              is_wall = (i == 1 .and. mesh%is_boundary_left) .or. &
                        (i == u_nx .and. mesh%is_boundary_right) .or. &
                        is_u_face_in_obstacle(mesh, cfg, i, j)
              if (is_wall) then
                ru(dof) = rhs_u(dof) / ns%u_aP(dof)
                cycle
              end if
              val = rhs_u(dof)
              if (i > 1) val = val + ns%u_aW(dof) * ru(dof - 1)
              if (i < u_nx) val = val + ns%u_aE(dof) * ru(dof + 1)
              if (j > 1) val = val + ns%u_aS(dof) * ru(dof - u_nx)
              if (j < mesh%ny_local) val = val + ns%u_aN(dof) * ru(dof + u_nx)
              ru(dof) = val / ns%u_aP(dof)
            end do
          end do
        end do
        dst(1:nu) = ru(1:nu)
        deallocate(rhs_u, ru)

        ! 3. Second membre pour vitesse v : rhs_v = zv - G_v * rp
        allocate(rhs_v(nv), rv(nv))
        rhs_v = 0.0_dp
        do j = 1, v_ny
          do i = 1, mesh%nx_local
            dof = (j - 1) * mesh%nx_local + i
            is_wall = (j == 1 .and. mesh%is_boundary_bottom .and. &
                       (cfg%bc_bottom%btype == BC_DIRICHLET .or. cfg%bc_bottom%btype == BC_FREESLIP)) .or. &
                      (j == v_ny .and. mesh%is_boundary_top .and. &
                       (cfg%bc_top%btype == BC_DIRICHLET .or. cfg%bc_top%btype == BC_FREESLIP)) .or. &
                      is_v_face_in_obstacle(mesh, cfg, i, j)
            if (is_wall) then
              rhs_v(dof) = src(nu + dof)
            else
              grad_p = (rp_2d(i, j) - rp_2d(i, j - 1)) * mesh%dx
              rhs_v(dof) = src(nu + dof) - grad_p
            end if
          end do
        end do

        ! 4 balayages de Gauss-Seidel symétrique sur Fv * rv = rhs_v
        rv = 0.0_dp
        do sweep = 1, 4
          ! Forward
          do j = 1, v_ny
            do i = 1, mesh%nx_local
              dof = (j - 1) * mesh%nx_local + i
              is_wall = (j == 1 .and. mesh%is_boundary_bottom .and. &
                         (cfg%bc_bottom%btype == BC_DIRICHLET .or. cfg%bc_bottom%btype == BC_FREESLIP)) .or. &
                        (j == v_ny .and. mesh%is_boundary_top .and. &
                         (cfg%bc_top%btype == BC_DIRICHLET .or. cfg%bc_top%btype == BC_FREESLIP)) .or. &
                        is_v_face_in_obstacle(mesh, cfg, i, j)
              if (is_wall) then
                rv(dof) = rhs_v(dof) / ns%v_aP(dof)
                cycle
              end if
              val = rhs_v(dof)
              if (i > 1) val = val + ns%v_aW(dof) * rv(dof - 1)
              if (i < mesh%nx_local) val = val + ns%v_aE(dof) * rv(dof + 1)
              if (j > 1) val = val + ns%v_aS(dof) * rv(dof - mesh%nx_local)
              if (j < v_ny) val = val + ns%v_aN(dof) * rv(dof + mesh%nx_local)
              rv(dof) = val / ns%v_aP(dof)
            end do
          end do
          ! Backward
          do j = v_ny, 1, -1
            do i = mesh%nx_local, 1, -1
              dof = (j - 1) * mesh%nx_local + i
              is_wall = (j == 1 .and. mesh%is_boundary_bottom .and. &
                         (cfg%bc_bottom%btype == BC_DIRICHLET .or. cfg%bc_bottom%btype == BC_FREESLIP)) .or. &
                        (j == v_ny .and. mesh%is_boundary_top .and. &
                         (cfg%bc_top%btype == BC_DIRICHLET .or. cfg%bc_top%btype == BC_FREESLIP)) .or. &
                        is_v_face_in_obstacle(mesh, cfg, i, j)
              if (is_wall) then
                rv(dof) = rhs_v(dof) / ns%v_aP(dof)
                cycle
              end if
              val = rhs_v(dof)
              if (i > 1) val = val + ns%v_aW(dof) * rv(dof - 1)
              if (i < mesh%nx_local) val = val + ns%v_aE(dof) * rv(dof + 1)
              if (j > 1) val = val + ns%v_aS(dof) * rv(dof - mesh%nx_local)
              if (j < v_ny) val = val + ns%v_aN(dof) * rv(dof + mesh%nx_local)
              rv(dof) = val / ns%v_aP(dof)
            end do
          end do
        end do
        dst(nu + 1 : nu + nv) = rv(1:nv)
        deallocate(rhs_v, rv)

        ! 4. Pression
        do j = 1, mesh%ny_local
          do i = 1, mesh%nx_local
            dof = nu + nv + (j - 1) * mesh%nx_local + i
            dst(dof) = rp_2d(i, j)
          end do
        end do

        deallocate(rp_2d, zp_2d)

      else
        ! Mode Lagrangien augmenté : diagonale Jacobi
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
      end if

    case default
      ! Affectation et calcul numérique
      dst(1:n) = src(1:n)
    end select
  end subroutine apply_preconditioner

  !-----------------------------------------------------------------------------
  ! Application de l'opérateur préconditionné : dst = P^{-1} * A * src
  !-----------------------------------------------------------------------------
  subroutine apply_PA(ns, mesh, cfg, fields_old, dt, src, dst)
    ! Déclaration de variable locale / paramètre de calcul
    type(navier_stokes_t), intent(in) :: ns
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)        :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)         :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(in)    :: fields_old
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)               :: dt
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)               :: src(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(out)              :: dst(:)

    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: tmp(:)
    ! Déclaration de variable locale / paramètre de calcul
    integer :: n

    ! Affectation et calcul numérique
    n = ns%n_owned_total
    ! Allocation dynamique de tableau contigu en mémoire
    allocate(tmp(n))

    ! Appel de sous-programme modulaire
    call ns%apply_system_matrix(mesh, cfg, fields_old, dt, src, tmp)
    ! Appel de sous-programme modulaire
    call apply_preconditioner(ns, mesh, cfg, tmp, dst)

    ! Libération de la mémoire allouée
    deallocate(tmp)
  end subroutine apply_PA

  !-----------------------------------------------------------------------------
  ! 1. Solveur BiCGSTAB(2) - Sleijpen & van der Vorst (1993)
  !-----------------------------------------------------------------------------
  subroutine solve_bicgstab2(ns, mesh, cfg, fields_old, dt, b, x, iters)
    ! Déclaration de variable locale / paramètre de calcul
    type(navier_stokes_t), intent(in) :: ns
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)        :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)         :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(in)    :: fields_old
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)               :: dt
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)               :: b(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(inout)            :: x(:)
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(out)               :: iters

    ! Déclaration de variable locale / paramètre de calcul
    integer :: n, max_iter, iter
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: tol, r0_norm, r_norm, rho_0, rho_1, alpha, omega, beta, gamma
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: sigma11, sigma12, sigma22, rhs1, rhs2, det, s11_s22, gamma1, gamma2
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: b_hat(:), r_tilde(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: r0(:), r1(:), r2(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: u0(:), u1(:), u2(:)

    ! Affectation et calcul numérique
    n = ns%n_owned_total
    ! Affectation et calcul numérique
    tol = cfg%tol
    ! Affectation et calcul numérique
    max_iter = cfg%max_iter
    ! Affectation et calcul numérique
    iters = 0

    ! Allocation dynamique de tableau contigu en mémoire
    allocate(b_hat(n), r_tilde(n))
    ! Allocation dynamique de tableau contigu en mémoire
    allocate(r0(n), r1(n), r2(n))
    ! Allocation dynamique de tableau contigu en mémoire
    allocate(u0(n), u1(n), u2(n))

    ! Affectation et calcul numérique
    u0 = 0.0_dp
    ! Affectation et calcul numérique
    u1 = 0.0_dp
    ! Affectation et calcul numérique
    u2 = 0.0_dp
    ! Affectation et calcul numérique
    r1 = 0.0_dp
    ! Affectation et calcul numérique
    r2 = 0.0_dp

    ! Second membre préconditionné b_hat = P^{-1} b
    call apply_preconditioner(ns, mesh, cfg, b, b_hat)

    ! Résidu initial r0 = b_hat - P^{-1} A x0
    call apply_PA(ns, mesh, cfg, fields_old, dt, x, r0)
    ! Affectation et calcul numérique
    r0(1:n) = b_hat(1:n) - r0(1:n)

    ! Évaluation de la norme du résidu ou de l'erreur numérique
    r0_norm = global_norm2(mesh, r0, n)
    ! Évaluation de la condition logique
    if (r0_norm < 1.0e-12_dp) then
      ! Libération de la mémoire allouée
      deallocate(b_hat, r_tilde, r0, r1, r2, u0, u1, u2)
      ! Retour à la procédure appelante
      return
    end if

    ! Affectation et calcul numérique
    r_tilde(1:n) = r0(1:n)
    ! Calcul de la masse volumique locale ou interpolée aux faces
    rho_0 = 1.0_dp
    ! Affectation et calcul numérique
    alpha = 0.0_dp
    ! Affectation et calcul numérique
    omega = 1.0_dp
    ! Affectation et calcul numérique
    iter = 0

    ! Boucle d'itération numérique
    do while (iter < max_iter)
      ! Affectation et calcul numérique
      iter = iter + 2
      ! Calcul de la masse volumique locale ou interpolée aux faces
      rho_0 = -omega * rho_0

      ! Demi-pas 1 (j = 0)
      rho_1 = global_dot_product(mesh, r_tilde, r0, n)
      if (abs(rho_0) < 1.0e-25_dp) exit
      ! Affectation et calcul numérique
      beta = alpha * (rho_1 / rho_0)
      ! Calcul de la masse volumique locale ou interpolée aux faces
      rho_0 = rho_1

      ! Affectation et calcul numérique
      u0(1:n) = r0(1:n) - beta * u0(1:n)
      ! Appel de sous-programme modulaire
      call apply_PA(ns, mesh, cfg, fields_old, dt, u0, u1)

      ! Affectation et calcul numérique
      gamma = global_dot_product(mesh, r_tilde, u1, n)
      if (abs(gamma) < 1.0e-25_dp) exit
      ! Affectation et calcul numérique
      alpha = rho_0 / gamma

      ! Affectation et calcul numérique
      r0(1:n) = r0(1:n) - alpha * u1(1:n)
      ! Affectation et calcul numérique
      x(1:n)  = x(1:n)  + alpha * u0(1:n)
      ! Appel de sous-programme modulaire
      call apply_PA(ns, mesh, cfg, fields_old, dt, r0, r1)

      ! Demi-pas 2 (j = 1)
      rho_1 = global_dot_product(mesh, r_tilde, r1, n)
      if (abs(rho_0) < 1.0e-25_dp) exit
      ! Affectation et calcul numérique
      beta = alpha * (rho_1 / rho_0)
      ! Calcul de la masse volumique locale ou interpolée aux faces
      rho_0 = rho_1

      ! Affectation et calcul numérique
      u0(1:n) = r0(1:n) - beta * u0(1:n)
      ! Affectation et calcul numérique
      u1(1:n) = r1(1:n) - beta * u1(1:n)
      ! Appel de sous-programme modulaire
      call apply_PA(ns, mesh, cfg, fields_old, dt, u1, u2)

      ! Affectation et calcul numérique
      gamma = global_dot_product(mesh, r_tilde, u2, n)
      if (abs(gamma) < 1.0e-25_dp) exit
      ! Affectation et calcul numérique
      alpha = rho_0 / gamma

      ! Affectation et calcul numérique
      x(1:n)  = x(1:n)  + alpha * u0(1:n)
      ! Affectation et calcul numérique
      r0(1:n) = r0(1:n) - alpha * u1(1:n)
      ! Affectation et calcul numérique
      r1(1:n) = r1(1:n) - alpha * u2(1:n)
      ! Appel de sous-programme modulaire
      call apply_PA(ns, mesh, cfg, fields_old, dt, r1, r2)

      ! Minimisation Quadratique MR-part
      sigma11 = global_dot_product(mesh, r1, r1, n)
      ! Affectation et calcul numérique
      sigma12 = global_dot_product(mesh, r2, r1, n)
      ! Affectation et calcul numérique
      sigma22 = global_dot_product(mesh, r2, r2, n)
      ! Affectation et calcul numérique
      rhs1    = global_dot_product(mesh, r0, r1, n)
      ! Affectation et calcul numérique
      rhs2    = global_dot_product(mesh, r0, r2, n)

      ! Affectation et calcul numérique
      det = sigma11 * sigma22 - sigma12 * sigma12
      ! Affectation et calcul numérique
      s11_s22 = sigma11 * sigma22

      ! Évaluation de la condition logique
      if (s11_s22 < 1.0e-30_dp .or. abs(det) < 1.0e-10_dp * s11_s22) then
        ! Affectation et calcul numérique
        gamma2 = 0.0_dp
        ! Affectation et calcul numérique
        gamma1 = rhs1 / (sigma11 + 1.0e-30_dp)
      ! Branche par défaut si condition non satisfaite
      else
        ! Affectation et calcul numérique
        gamma1 = (rhs1 * sigma22 - rhs2 * sigma12) / det
        ! Affectation et calcul numérique
        gamma2 = (rhs2 * sigma11 - rhs1 * sigma12) / det
      end if

      ! Affectation et calcul numérique
      x(1:n)  = x(1:n)  + gamma1 * r0(1:n) + gamma2 * r1(1:n)
      ! Affectation et calcul numérique
      r0(1:n) = r0(1:n) - (gamma1 * r1(1:n) + gamma2 * r2(1:n))
      ! Affectation et calcul numérique
      u0(1:n) = u0(1:n) - (gamma1 * u1(1:n) + gamma2 * u2(1:n))
      ! Affectation et calcul numérique
      omega = gamma2

      ! Évaluation de la norme du résidu ou de l'erreur numérique
      r_norm = global_norm2(mesh, r0, n)
      ! Évaluation de la condition logique
      if (r_norm / r0_norm < tol .or. r_norm < 1.0e-13_dp) then
        ! Affectation et calcul numérique
        iters = iter
        ! Libération de la mémoire allouée
        deallocate(b_hat, r_tilde, r0, r1, r2, u0, u1, u2)
        ! Retour à la procédure appelante
        return
      end if
    end do

    ! Affectation et calcul numérique
    iters = iter
    ! Libération de la mémoire allouée
    deallocate(b_hat, r_tilde, r0, r1, r2, u0, u1, u2)
  end subroutine solve_bicgstab2

  !-----------------------------------------------------------------------------
  ! 2. Solveur BiCGSTAB standard - van der Vorst (1992)
  !-----------------------------------------------------------------------------
  subroutine solve_bicgstab(ns, mesh, cfg, fields_old, dt, b, x, iters)
    ! Déclaration de variable locale / paramètre de calcul
    type(navier_stokes_t), intent(in) :: ns
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)        :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)         :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(in)    :: fields_old
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)               :: dt
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)               :: b(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(inout)            :: x(:)
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(out)               :: iters

    ! Déclaration de variable locale / paramètre de calcul
    integer :: n, max_iter, iter
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: tol, r0_norm, r_norm, rho, rho_old, alpha, omega, beta, denom, s_norm, t_dot_t
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: b_hat(:), r(:), r_tilde(:), p(:), v(:), s(:), t(:)

    ! Affectation et calcul numérique
    n = ns%n_owned_total
    ! Affectation et calcul numérique
    tol = cfg%tol
    ! Affectation et calcul numérique
    max_iter = cfg%max_iter
    ! Affectation et calcul numérique
    iters = max_iter

    ! Allocation dynamique de tableau contigu en mémoire
    allocate(b_hat(n), r(n), r_tilde(n), p(n), v(n), s(n), t(n))

    ! Appel de sous-programme modulaire
    call apply_preconditioner(ns, mesh, cfg, b, b_hat)
    ! Appel de sous-programme modulaire
    call apply_PA(ns, mesh, cfg, fields_old, dt, x, r)
    ! Affectation et calcul numérique
    r(1:n) = b_hat(1:n) - r(1:n)

    ! Évaluation de la norme du résidu ou de l'erreur numérique
    r0_norm = global_norm2(mesh, r, n)
    ! Évaluation de la condition logique
    if (r0_norm < 1.0e-12_dp) then
      ! Affectation et calcul numérique
      iters = 0
      ! Libération de la mémoire allouée
      deallocate(b_hat, r, r_tilde, p, v, s, t)
      ! Retour à la procédure appelante
      return
    end if

    ! Affectation et calcul numérique
    r_tilde(1:n) = r(1:n)
    ! Calcul de la masse volumique locale ou interpolée aux faces
    rho_old = 1.0_dp
    ! Affectation et calcul numérique
    alpha = 1.0_dp
    ! Affectation et calcul numérique
    omega = 1.0_dp

    ! Boucle d'itération numérique
    do iter = 1, max_iter
      ! Calcul de la masse volumique locale ou interpolée aux faces
      rho = global_dot_product(mesh, r_tilde, r, n)
      if (abs(rho) < 1.0e-30_dp) exit

      ! Évaluation de la condition logique
      if (iter == 1) then
        ! Actualisation de la pression statique au centre de la maille
        p(1:n) = r(1:n)
      ! Branche par défaut si condition non satisfaite
      else
        ! Affectation et calcul numérique
        beta = (rho / rho_old) * (alpha / omega)
        ! Actualisation de la pression statique au centre de la maille
        p(1:n) = r(1:n) + beta * (p(1:n) - omega * v(1:n))
      end if

      ! Appel de sous-programme modulaire
      call apply_PA(ns, mesh, cfg, fields_old, dt, p, v)
      ! Affectation et calcul numérique
      denom = global_dot_product(mesh, r_tilde, v, n)
      if (abs(denom) < 1.0e-30_dp) exit
      ! Affectation et calcul numérique
      alpha = rho / denom

      ! Affectation et calcul numérique
      s(1:n) = r(1:n) - alpha * v(1:n)
      ! Évaluation de la norme du résidu ou de l'erreur numérique
      s_norm = global_norm2(mesh, s, n)
      ! Évaluation de la condition logique
      if (s_norm / r0_norm < tol .or. s_norm < 1.0e-13_dp) then
        ! Affectation et calcul numérique
        x(1:n) = x(1:n) + alpha * p(1:n)
        ! Affectation et calcul numérique
        iters = iter
        ! Libération de la mémoire allouée
        deallocate(b_hat, r, r_tilde, p, v, s, t)
        ! Retour à la procédure appelante
        return
      end if

      ! Appel de sous-programme modulaire
      call apply_PA(ns, mesh, cfg, fields_old, dt, s, t)
      ! Affectation et calcul numérique
      t_dot_t = global_dot_product(mesh, t, t, n)
      ! Évaluation de la condition logique
      if (t_dot_t > 1.0e-30_dp) then
        ! Affectation et calcul numérique
        omega = global_dot_product(mesh, t, s, n) / t_dot_t
      ! Branche par défaut si condition non satisfaite
      else
        ! Affectation et calcul numérique
        omega = 0.0_dp
      end if

      ! Affectation et calcul numérique
      x(1:n) = x(1:n) + alpha * p(1:n) + omega * s(1:n)
      ! Affectation et calcul numérique
      r(1:n) = s(1:n) - omega * t(1:n)

      ! Calcul de la masse volumique locale ou interpolée aux faces
      rho_old = rho
      ! Évaluation de la norme du résidu ou de l'erreur numérique
      r_norm = global_norm2(mesh, r, n)
      ! Évaluation de la condition logique
      if (r_norm / r0_norm < tol .or. r_norm < 1.0e-13_dp) then
        ! Affectation et calcul numérique
        iters = iter
        ! Libération de la mémoire allouée
        deallocate(b_hat, r, r_tilde, p, v, s, t)
        ! Retour à la procédure appelante
        return
      end if
      if (abs(omega) < 1.0e-30_dp) exit
    end do

    ! Libération de la mémoire allouée
    deallocate(b_hat, r, r_tilde, p, v, s, t)
  end subroutine solve_bicgstab

  !-----------------------------------------------------------------------------
  ! 3. Solveur Gradient Conjugué (PCG) - Hestenes & Stiefel (1952)
  !-----------------------------------------------------------------------------
  subroutine solve_cg(ns, mesh, cfg, fields_old, dt, b, x, iters)
    ! Déclaration de variable locale / paramètre de calcul
    type(navier_stokes_t), intent(in) :: ns
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)        :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)         :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(in)    :: fields_old
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)               :: dt
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)               :: b(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(inout)            :: x(:)
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(out)               :: iters

    ! Déclaration de variable locale / paramètre de calcul
    integer :: n, max_iter, iter
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: tol, r0_norm, r_norm, rz_old, rz_new, pAp, alpha, beta
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: Ax(:), r(:), z(:), p(:), Ap(:)

    ! Affectation et calcul numérique
    n = ns%n_owned_total
    ! Affectation et calcul numérique
    tol = cfg%tol
    ! Affectation et calcul numérique
    max_iter = cfg%max_iter
    ! Affectation et calcul numérique
    iters = max_iter

    ! Allocation dynamique de tableau contigu en mémoire
    allocate(Ax(n), r(n), z(n), p(n), Ap(n))

    ! Appel de sous-programme modulaire
    call ns%apply_system_matrix(mesh, cfg, fields_old, dt, x, Ax)
    ! Affectation et calcul numérique
    r(1:n) = b(1:n) - Ax(1:n)

    ! Évaluation de la norme du résidu ou de l'erreur numérique
    r0_norm = global_norm2(mesh, r, n)
    ! Évaluation de la condition logique
    if (r0_norm < 1.0e-12_dp) then
      ! Affectation et calcul numérique
      iters = 0
      ! Libération de la mémoire allouée
      deallocate(Ax, r, z, p, Ap)
      ! Retour à la procédure appelante
      return
    end if

    ! Appel de sous-programme modulaire
    call apply_preconditioner(ns, mesh, cfg, r, z)
    ! Actualisation de la pression statique au centre de la maille
    p(1:n) = z(1:n)
    ! Affectation et calcul numérique
    rz_old = global_dot_product(mesh, r, z, n)

    ! Boucle d'itération numérique
    do iter = 1, max_iter
      ! Appel de sous-programme modulaire
      call ns%apply_system_matrix(mesh, cfg, fields_old, dt, p, Ap)
      ! Affectation et calcul numérique
      pAp = global_dot_product(mesh, p, Ap, n)
      if (abs(pAp) < 1.0e-30_dp) exit

      ! Affectation et calcul numérique
      alpha = rz_old / pAp
      ! Affectation et calcul numérique
      x(1:n) = x(1:n) + alpha * p(1:n)
      ! Affectation et calcul numérique
      r(1:n) = r(1:n) - alpha * Ap(1:n)

      ! Évaluation de la norme du résidu ou de l'erreur numérique
      r_norm = global_norm2(mesh, r, n)
      ! Évaluation de la condition logique
      if (r_norm / r0_norm < tol .or. r_norm < 1.0e-13_dp) then
        ! Affectation et calcul numérique
        iters = iter
        ! Libération de la mémoire allouée
        deallocate(Ax, r, z, p, Ap)
        ! Retour à la procédure appelante
        return
      end if

      ! Appel de sous-programme modulaire
      call apply_preconditioner(ns, mesh, cfg, r, z)
      ! Affectation et calcul numérique
      rz_new = global_dot_product(mesh, r, z, n)
      ! Affectation et calcul numérique
      beta = rz_new / rz_old
      ! Affectation et calcul numérique
      rz_old = rz_new

      ! Actualisation de la pression statique au centre de la maille
      p(1:n) = z(1:n) + beta * p(1:n)
    end do

    ! Libération de la mémoire allouée
    deallocate(Ax, r, z, p, Ap)
  end subroutine solve_cg

  !-----------------------------------------------------------------------------
  ! 4. Solveur GMRES(m) redémarré - Saad & Schultz (1986)
  !-----------------------------------------------------------------------------
  subroutine solve_gmres(ns, mesh, cfg, fields_old, dt, b, x, iters)
    ! Déclaration de variable locale / paramètre de calcul
    type(navier_stokes_t), intent(in) :: ns
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)        :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)         :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(in)    :: fields_old
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)               :: dt
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)               :: b(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(inout)            :: x(:)
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(out)               :: iters

    ! Déclaration de variable locale / paramètre de calcul
    integer, parameter :: m = 30
    ! Déclaration de variable locale / paramètre de calcul
    integer :: n, max_iter, total_iters, i, j, k
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: tol, r0_norm, beta_val, r_norm, h_val, temp, cs, sn
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: b_hat(:), r(:), w(:), g(:), y(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: V(:, :), H(:, :), c_rot(:), s_rot(:)

    ! Affectation et calcul numérique
    n = ns%n_owned_total
    ! Affectation et calcul numérique
    tol = cfg%tol
    ! Affectation et calcul numérique
    max_iter = cfg%max_iter
    ! Affectation et calcul numérique
    iters = max_iter

    ! Allocation dynamique de tableau contigu en mémoire
    allocate(b_hat(n), r(n), w(n), g(m + 1), y(m))
    ! Allocation dynamique de tableau contigu en mémoire
    allocate(V(n, m + 1), H(m + 1, m), c_rot(m), s_rot(m))

    ! Appel de sous-programme modulaire
    call apply_preconditioner(ns, mesh, cfg, b, b_hat)

    ! Affectation et calcul numérique
    total_iters = 0
    ! Évaluation de la norme du résidu ou de l'erreur numérique
    r0_norm = 1.0_dp

    ! Boucle d'itération numérique
    do while (total_iters < max_iter)
      ! Appel de sous-programme modulaire
      call apply_PA(ns, mesh, cfg, fields_old, dt, x, r)
      ! Affectation et calcul numérique
      r(1:n) = b_hat(1:n) - r(1:n)

      ! Affectation et calcul numérique
      beta_val = global_norm2(mesh, r, n)
      ! Évaluation de la condition logique
      if (total_iters == 0) then
        ! Évaluation de la condition logique
        if (beta_val > 1.0e-14_dp) then
          ! Évaluation de la norme du résidu ou de l'erreur numérique
          r0_norm = beta_val
        ! Branche par défaut si condition non satisfaite
        else
          ! Évaluation de la norme du résidu ou de l'erreur numérique
          r0_norm = 1.0_dp
        end if
      end if

      ! Évaluation de la condition logique
      if (beta_val / r0_norm < tol .or. beta_val < 1.0e-13_dp) then
        ! Affectation et calcul numérique
        iters = total_iters
        ! Libération de la mémoire allouée
        deallocate(b_hat, r, w, g, y, V, H, c_rot, s_rot)
        ! Retour à la procédure appelante
        return
      end if

      ! Affectation et calcul numérique
      V(1:n, 1) = r(1:n) / beta_val
      ! Affectation et calcul numérique
      g = 0.0_dp
      ! Affectation et calcul numérique
      g(1) = beta_val
      ! Affectation et calcul numérique
      H = 0.0_dp

      ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
      do j = 1, min(m, max_iter - total_iters)
        ! Affectation et calcul numérique
        total_iters = total_iters + 1

        ! Appel de sous-programme modulaire
        call apply_PA(ns, mesh, cfg, fields_old, dt, V(1:n, j), w)

        ! Orthogonalisation Modified Gram-Schmidt
        do i = 1, j
          ! Affectation et calcul numérique
          H(i, j) = global_dot_product(mesh, w, V(1:n, i), n)
          ! Affectation et calcul numérique
          w(1:n) = w(1:n) - H(i, j) * V(1:n, i)
        end do

        ! Affectation et calcul numérique
        H(j + 1, j) = global_norm2(mesh, w, n)
        ! Évaluation de la condition logique
        if (H(j + 1, j) > 1.0e-30_dp) then
          ! Affectation et calcul numérique
          V(1:n, j + 1) = w(1:n) / H(j + 1, j)
        end if

        ! Application des rotations de Givens précédentes
        do i = 1, j - 1
          ! Affectation et calcul numérique
          temp = c_rot(i) * H(i, j) + s_rot(i) * H(i + 1, j)
          ! Affectation et calcul numérique
          H(i + 1, j) = -s_rot(i) * H(i, j) + c_rot(i) * H(i + 1, j)
          ! Affectation et calcul numérique
          H(i, j) = temp
        end do

        ! Calcul de la nouvelle rotation de Givens
        if (abs(H(j + 1, j)) < 1.0e-30_dp) then
          ! Affectation et calcul numérique
          c_rot(j) = 1.0_dp
          ! Affectation et calcul numérique
          s_rot(j) = 0.0_dp
        ! Branche conditionnelle alternative
        else if (abs(H(j, j)) > abs(H(j + 1, j))) then
          ! Affectation et calcul numérique
          temp = H(j + 1, j) / H(j, j)
          ! Affectation et calcul numérique
          c_rot(j) = 1.0_dp / sqrt(1.0_dp + temp * temp)
          ! Affectation et calcul numérique
          s_rot(j) = temp * c_rot(j)
        ! Branche par défaut si condition non satisfaite
        else
          ! Affectation et calcul numérique
          temp = H(j, j) / H(j + 1, j)
          ! Affectation et calcul numérique
          s_rot(j) = 1.0_dp / sqrt(1.0_dp + temp * temp)
          ! Affectation et calcul numérique
          c_rot(j) = temp * s_rot(j)
        end if

        ! Affectation et calcul numérique
        H(j, j) = c_rot(j) * H(j, j) + s_rot(j) * H(j + 1, j)
        ! Affectation et calcul numérique
        H(j + 1, j) = 0.0_dp

        ! Affectation et calcul numérique
        g(j + 1) = -s_rot(j) * g(j)
        ! Affectation et calcul numérique
        g(j)     =  c_rot(j) * g(j)

        ! Évaluation de la norme du résidu ou de l'erreur numérique
        r_norm = abs(g(j + 1))
        ! Évaluation de la condition logique
        if (r_norm / r0_norm < tol .or. r_norm < 1.0e-13_dp .or. j == m) then
          ! Résolution du système triangulaire supérieur H * y = g
          do k = j, 1, -1
            ! Affectation et calcul numérique
            y(k) = g(k)
            ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
            do i = k + 1, j
              ! Affectation et calcul numérique
              y(k) = y(k) - H(k, i) * y(i)
            end do
            ! Affectation et calcul numérique
            if (abs(H(k, k)) > 1.0e-30_dp) y(k) = y(k) / H(k, k)
          end do

          ! Mise à jour de la solution x = x + V * y
          do i = 1, j
            ! Affectation et calcul numérique
            x(1:n) = x(1:n) + y(i) * V(1:n, i)
          end do
          exit
        end if
      end do
    end do

    ! Affectation et calcul numérique
    iters = total_iters
    ! Libération de la mémoire allouée
    deallocate(b_hat, r, w, g, y, V, H, c_rot, s_rot)
  end subroutine solve_gmres

  !-----------------------------------------------------------------------------
  ! Aiguilleur principal de solveur linéaire
  !-----------------------------------------------------------------------------
  subroutine solve_linear_system(ns, mesh, cfg, fields_old, dt, b, x, iters)
    ! Déclaration de variable locale / paramètre de calcul
    type(navier_stokes_t), intent(in) :: ns
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)        :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in)         :: cfg
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(in)    :: fields_old
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)               :: dt
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)               :: b(:)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(inout)            :: x(:)
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(out)               :: iters

    select case (cfg%linear_solver)
    case (LINEAR_BICGSTAB2)
      ! Appel de sous-programme modulaire
      call solve_bicgstab2(ns, mesh, cfg, fields_old, dt, b, x, iters)
    case (LINEAR_BICGSTAB)
      ! Appel de sous-programme modulaire
      call solve_bicgstab(ns, mesh, cfg, fields_old, dt, b, x, iters)
    case (LINEAR_CG)
      ! Appel de sous-programme modulaire
      call solve_cg(ns, mesh, cfg, fields_old, dt, b, x, iters)
    case (LINEAR_GMRES)
      ! Appel de sous-programme modulaire
      call solve_gmres(ns, mesh, cfg, fields_old, dt, b, x, iters)
    case default
      ! Appel de sous-programme modulaire
      call solve_bicgstab2(ns, mesh, cfg, fields_old, dt, b, x, iters)
    end select
  end subroutine solve_linear_system

end module module_linear_solvers
