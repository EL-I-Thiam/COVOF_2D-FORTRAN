! [Procédure/Module] module module_interface_vof
module module_interface_vof
  use module_precision, only: dp
  use module_constants
  use module_mesh, only: mesh_2d_t, x_center, y_center, x_face, y_face, is_in_obstacle
  use module_fields, only: flow_fields_t
  use module_config, only: config_t
  use module_mpi_comm, only: exchange_halos_scalar, apply_vof_bc
  use module_interface_reconstruction, only: interface_segment_t, cut_box_area, compute_alpha, &
                                            reconstruct_interface, evaluate_limiter
  implicit none
  private

  public :: interface_vof_t
  public :: set_initial_conditions
  public :: advect_vof
  public :: update_physical_properties
  public :: compute_curvature_and_csf

  ! Compteur de pas pour Strang splitting alterné
  integer, save :: vof_step_count = 0

  type :: interface_vof_t
     ! Déclaration de variable locale / paramètre de calcul
     type(config_t) :: cfg
  contains
     ! Affectation et calcul numérique
     procedure :: init => vof_init
     ! Affectation et calcul numérique
     procedure :: set_initial_condition => vof_set_ic
     ! Affectation et calcul numérique
     procedure :: advect => vof_advect_step
  end type interface_vof_t

contains

  ! [Procédure/Module] subroutine set_initial_conditions
  subroutine set_initial_conditions(fields, mesh, cfg)
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(inout) :: fields
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in) :: cfg

    ! Déclaration de variable locale / paramètre de calcul
    integer :: i, j, nx, ny
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: x, y, dist, eps, int_y, x_c_jet, half_w_jet, stub_h_jet

    ! Affectation et calcul numérique
    nx = mesh%nx_local
    ! Affectation et calcul numérique
    ny = mesh%ny_local
    ! Affectation et calcul numérique
    eps = 1.5_dp * max(mesh%dx, mesh%dy)

    ! 1. Initialisation des vitesses (par défaut 0 ou imposé)
    fields%u = cfg%ic_u
    ! Affectation et calcul numérique
    fields%v = cfg%ic_v

    ! 2. Initialisation de la fraction volumique C
    do j = 0, ny + 1
       ! Affectation et calcul numérique
       y = y_center(mesh, j)
       ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
       do i = 0, nx + 1
          ! Affectation et calcul numérique
          x = x_center(mesh, i)

          select case (cfg%ic_type)
          case (IC_UNIFORM)
             ! Transport / actualisation de la fraction volumique VOF
             fields%C(i, j) = 1.0_dp

          case (IC_DROPLET)
             ! Affectation et calcul numérique
             dist = sqrt((x - cfg%ic_center_x)**2 + (y - cfg%ic_center_y)**2) - cfg%ic_radius
             ! Évaluation de la condition logique
             if (dist < -eps) then
                ! Transport / actualisation de la fraction volumique VOF
                fields%C(i, j) = 1.0_dp
             ! Branche conditionnelle alternative
             else if (dist > eps) then
                ! Transport / actualisation de la fraction volumique VOF
                fields%C(i, j) = 0.0_dp
             ! Branche par défaut si condition non satisfaite
             else
                ! Transport / actualisation de la fraction volumique VOF
                fields%C(i, j) = 0.5_dp - dist / (2.0_dp * eps) - (1.0_dp / (2.0_dp * PI)) * sin(PI * dist / eps)
             end if

          case (IC_BUBBLE)
             ! Affectation et calcul numérique
             dist = sqrt((x - cfg%ic_center_x)**2 + (y - cfg%ic_center_y)**2) - cfg%ic_radius
             ! Évaluation de la condition logique
             if (dist < -eps) then
                ! Transport / actualisation de la fraction volumique VOF
                fields%C(i, j) = 0.0_dp ! Intérieur bulle (fluide 2)
             ! Branche conditionnelle alternative
             else if (dist > eps) then
                ! Transport / actualisation de la fraction volumique VOF
                fields%C(i, j) = 1.0_dp ! Extérieur (fluide 1)
             ! Branche par défaut si condition non satisfaite
             else
                ! Transport / actualisation de la fraction volumique VOF
                fields%C(i, j) = 0.5_dp + dist / (2.0_dp * eps) + (1.0_dp / (2.0_dp * PI)) * sin(PI * dist / eps)
             end if

          case (IC_DAM_BREAK)
             ! Affectation et calcul numérique
             if (x >= cfg%ic_box_x0 .and. x <= cfg%ic_box_x1 .and. &
                 ! Affectation et calcul numérique
                 y >= cfg%ic_box_y0 .and. y <= cfg%ic_box_y1) then
                ! Transport / actualisation de la fraction volumique VOF
                fields%C(i, j) = 1.0_dp
             ! Branche par défaut si condition non satisfaite
             else
                ! Transport / actualisation de la fraction volumique VOF
                fields%C(i, j) = 0.0_dp
             end if

          case (IC_LIQUID_SHEET)
             ! Nappe liquide cisaillée (Cas 3.2 de l'article Fluids 2021)
             if (cfg%ic_center_y > 0.0_dp) then
                ! Calcul de distance à l'axe de la nappe
                dist = abs(y - cfg%ic_center_y) - 0.5_dp * cfg%ic_thickness
             else
                ! Calcul de distance au plan médian
                dist = abs(y - 0.5_dp * mesh%ly) - 0.5_dp * cfg%ic_thickness
             end if
             ! Évaluation de l'indicatrice C
             if (dist < -eps) then
                fields%C(i, j) = 0.0_dp ! Intérieur nappe liquide (fluide 2)
             else if (dist > eps) then
                fields%C(i, j) = 1.0_dp ! Gaz ambiant (fluide 1)
             else
                fields%C(i, j) = 0.5_dp + dist / (2.0_dp * eps) + &
                     (1.0_dp / (2.0_dp * PI)) * sin(PI * dist / eps)
             end if

          case (IC_VISCOUS_JET)
             ! Jet liquide visqueux injecté (Cas 3.3 de l'article Fluids 2021)
             if (cfg%ic_center_x > 0.0_dp) then
                x_c_jet = cfg%ic_center_x
             else
                x_c_jet = 0.5_dp * mesh%lx
             end if
             if (cfg%ic_jet_width > 0.0_dp) then
                half_w_jet = 0.5_dp * cfg%ic_jet_width
             else
                half_w_jet = 0.04_dp
             end if
             stub_h_jet = 2.0_dp * mesh%dy
             if (y >= mesh%ly - stub_h_jet .and. abs(x - x_c_jet) <= half_w_jet) then
                fields%C(i, j) = 0.0_dp ! Amorce liquide (fluide 2)
             else
                fields%C(i, j) = 1.0_dp ! Cavité de gaz (fluide 1)
             end if

          case (IC_VORTEX_SHEAR)
             ! Affectation et calcul numérique
             dist = sqrt((x - cfg%ic_center_x)**2 + (y - cfg%ic_center_y)**2) - cfg%ic_radius
             ! Évaluation de la condition logique
             if (dist < -eps) then
                ! Transport / actualisation de la fraction volumique VOF
                fields%C(i, j) = 1.0_dp
             ! Branche conditionnelle alternative
             else if (dist > eps) then
                ! Transport / actualisation de la fraction volumique VOF
                fields%C(i, j) = 0.0_dp
             ! Branche par défaut si condition non satisfaite
             else
                ! Transport / actualisation de la fraction volumique VOF
                fields%C(i, j) = 0.5_dp - dist / (2.0_dp * eps) - (1.0_dp / (2.0_dp * PI)) * sin(PI * dist / eps)
             end if

          case (IC_RAYLEIGH_TAYLOR)
             ! Affectation et calcul numérique
             int_y = cfg%ic_center_y + cfg%ic_pert_amp * cos(2.0_dp * PI * cfg%ic_pert_k * x / mesh%lx)
             ! Évaluation de la condition logique
             if (y >= int_y) then
                ! Transport / actualisation de la fraction volumique VOF
                fields%C(i, j) = 1.0_dp ! Fluide lourd au-dessus
             ! Branche par défaut si condition non satisfaite
             else
                ! Transport / actualisation de la fraction volumique VOF
                fields%C(i, j) = 0.0_dp ! Fluide léger en dessous
             end if

          case (IC_LAPLACIAN_TEST)
             ! Transport / actualisation de la fraction volumique VOF
             fields%C(i, j) = 1.0_dp

          case default
             ! Transport / actualisation de la fraction volumique VOF
             fields%C(i, j) = 1.0_dp
          end select
       end do
    end do

    ! 3. Détection des obstacles solides pour la pénalisation de Brinkman
    fields%solid_mask = 0
    ! Évaluation de la condition logique
    if (cfg%obs_enabled) then
       ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
       do j = 0, ny + 1
          ! Affectation et calcul numérique
          y = y_center(mesh, j)
          ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
          do i = 0, nx + 1
             ! Affectation et calcul numérique
             x = x_center(mesh, i)
             if (is_in_obstacle(x, y, cfg%obs_type, cfg%obs_x0, cfg%obs_y0, cfg%obs_x1, cfg%obs_y1, &
                                cfg%obs_xc, cfg%obs_yc, cfg%obs_radius)) then
                ! Affectation et calcul numérique
                fields%solid_mask(i, j) = 1
             end if
          end do
       end do
    end if

    ! Appel de sous-programme modulaire
    call exchange_halos_scalar(fields%C, mesh)
    ! Appel de sous-programme modulaire
    call apply_vof_bc(fields, mesh, cfg)
    ! Appel de sous-programme modulaire
    call update_physical_properties(fields, mesh, cfg)
    ! Appel de sous-programme modulaire
    call compute_curvature_and_csf(fields, mesh, cfg)
  end subroutine set_initial_conditions

  ! [Procédure/Module] subroutine update_physical_properties
  subroutine update_physical_properties(fields, mesh, cfg)
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(inout) :: fields
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in) :: cfg

    ! Déclaration de variable locale / paramètre de calcul
    integer :: i, j, nx, ny
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: c, y, z, mu_val, m_sw, m_se, m_nw, m_ne

    ! Affectation et calcul numérique
    nx = mesh%nx_local
    ! Affectation et calcul numérique
    ny = mesh%ny_local

    ! Évaluation de la condition logique
    if (.not. cfg%enable_two_phase) then
       ! Calcul de la masse volumique locale ou interpolée aux faces
       fields%rho = cfg%rho1
       ! Calcul de la masse volumique locale ou interpolée aux faces
       fields%rho_u = cfg%rho1
       ! Calcul de la masse volumique locale ou interpolée aux faces
       fields%rho_v = cfg%rho1
       ! Évaluation de la condition logique
       if (cfg%viscosity_model == VISC_EXPONENTIAL) then
          ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
          do j = 0, ny + 1
             ! Affectation et calcul numérique
             y = y_center(mesh, j)
             ! Affectation et calcul numérique
             z = y - cfg%visc_e
             ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
             mu_val = cfg%visc_mu0 * exp(-cfg%visc_a * z / cfg%visc_e)
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             fields%mu(:, j) = mu_val
          end do
          ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
          do j = 0, ny + 1
             ! Affectation et calcul numérique
             y = y_face(mesh, j)
             ! Affectation et calcul numérique
             z = y - cfg%visc_e
             ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
             mu_val = cfg%visc_mu0 * exp(-cfg%visc_a * z / cfg%visc_e)
             ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
             fields%mu_node(:, j) = mu_val
          end do
       ! Branche par défaut si condition non satisfaite
       else
          ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
          fields%mu = cfg%mu1
          ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
          fields%mu_node = cfg%mu1
       end if
       ! Retour à la procédure appelante
       return
    end if

    ! Coins de C
    fields%C(0, 0)       = fields%C(1, 1)
    ! Transport / actualisation de la fraction volumique VOF
    fields%C(nx+1, 0)    = fields%C(nx, 1)
    ! Transport / actualisation de la fraction volumique VOF
    fields%C(0, ny+1)    = fields%C(1, ny)
    ! Transport / actualisation de la fraction volumique VOF
    fields%C(nx+1, ny+1) = fields%C(nx, ny)

    ! 1. Masse volumique arithmétique et viscosité au centre
    do j = 0, ny + 1
       ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
       do i = 0, nx + 1
          ! Affectation et calcul numérique
          c = max(0.0_dp, min(1.0_dp, fields%C(i, j)))
          ! Calcul de la masse volumique locale ou interpolée aux faces
          fields%rho(i, j) = c * cfg%rho1 + (1.0_dp - c) * cfg%rho2
          ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
          fields%mu(i, j)  = c * cfg%mu1  + (1.0_dp - c) * cfg%mu2
       end do
    end do

    ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
    fields%mu(0, 0)       = fields%mu(1, 1)
    ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
    fields%mu(nx+1, 0)    = fields%mu(nx, 1)
    ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
    fields%mu(0, ny+1)    = fields%mu(1, ny)
    ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
    fields%mu(nx+1, ny+1) = fields%mu(nx, ny)

    ! 2. Masse volumique sur les faces u
    do j = 1, ny
       ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
       do i = 1, nx + 1
          ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
          fields%rho_u(i, j) = 0.5_dp * (fields%rho(i - 1, j) + fields%rho(i, j))
       end do
    end do

    ! 3. Masse volumique sur les faces v
    do j = 1, ny + 1
       ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
       do i = 1, nx
          ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
          fields%rho_v(i, j) = 0.5_dp * (fields%rho(i, j - 1) + fields%rho(i, j))
       end do
    end do

    ! 4. Viscosité harmonique aux nœuds (coins de mailles)
    do j = 1, ny + 1
       ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
       do i = 1, nx + 1
          ! Affectation et calcul numérique
          m_sw = fields%mu(i-1, j-1)
          ! Affectation et calcul numérique
          m_se = fields%mu(i, j-1)
          ! Affectation et calcul numérique
          m_nw = fields%mu(i-1, j)
          ! Affectation et calcul numérique
          m_ne = fields%mu(i, j)
          ! Calcul de la viscosité dynamique locale (harmonique aux nœuds)
          fields%mu_node(i, j) = 4.0_dp / (1.0_dp / m_sw + 1.0_dp / m_se + 1.0_dp / m_nw + 1.0_dp / m_ne)
       end do
    end do
  end subroutine update_physical_properties

  ! [Procédure/Module] pure real
  pure real(dp) function superbee(r)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in) :: r
    ! Affectation et calcul numérique
    superbee = max(0.0_dp, max(min(2.0_dp * r, 1.0_dp), min(r, 2.0_dp)))
  end function superbee

  ! [Procédure/Module] subroutine compute_curvature_and_csf
  subroutine compute_curvature_and_csf(fields, mesh, cfg)
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(inout) :: fields
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in) :: cfg

    ! Déclaration de variable locale / paramètre de calcul
    integer :: i, j, nx, ny
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: dx, dy, dc_dx, dc_dy, grad_mag, kap_avg, grad_c
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: norm_x(:, :), norm_y(:, :)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: kappa(:, :)

    ! Évaluation de la condition logique
    if (.not. cfg%enable_surface_tension .or. cfg%sigma <= 0.0_dp) then
       ! Affectation et calcul numérique
       fields%fx_sigma = 0.0_dp
       ! Affectation et calcul numérique
       fields%fy_sigma = 0.0_dp
       ! Affectation et calcul numérique
       fields%curvature = 0.0_dp
       ! Retour à la procédure appelante
       return
    end if

    ! Affectation et calcul numérique
    nx = mesh%nx_local
    ! Affectation et calcul numérique
    ny = mesh%ny_local
    ! Calcul de la métrique spatiale élémentaire
    dx = mesh%dx
    ! Calcul de la métrique spatiale élémentaire
    dy = mesh%dy

    ! Allocation dynamique de tableau contigu en mémoire
    allocate(norm_x(0:nx+1, 0:ny+1), norm_y(0:nx+1, 0:ny+1))
    ! Allocation dynamique de tableau contigu en mémoire
    allocate(kappa(0:nx+1, 0:ny+1))
    ! Évaluation de la norme du résidu ou de l'erreur numérique
    norm_x = 0.0_dp
    ! Évaluation de la norme du résidu ou de l'erreur numérique
    norm_y = 0.0_dp
    ! Affectation et calcul numérique
    kappa = 0.0_dp

    ! 1. Gradient et normale unitaire aux centres
    do j = 1, ny
       ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
       do i = 1, nx
          ! Calcul de la métrique spatiale élémentaire
          dc_dx = (fields%C(i+1, j) - fields%C(i-1, j)) / (2.0_dp * dx)
          ! Calcul de la métrique spatiale élémentaire
          dc_dy = (fields%C(i, j+1) - fields%C(i, j-1)) / (2.0_dp * dy)
          ! Affectation et calcul numérique
          grad_mag = sqrt(dc_dx * dc_dx + dc_dy * dc_dy)
          ! Évaluation de la condition logique
          if (grad_mag > 1.0e-3_dp / dx) then
             ! Évaluation de la norme du résidu ou de l'erreur numérique
             norm_x(i, j) = dc_dx / grad_mag
             ! Évaluation de la norme du résidu ou de l'erreur numérique
             norm_y(i, j) = dc_dy / grad_mag
          end if
       end do
    end do

    ! Appel de sous-programme modulaire
    call exchange_halos_scalar(norm_x, mesh)
    ! Appel de sous-programme modulaire
    call exchange_halos_scalar(norm_y, mesh)

    ! 2. Courbure kappa = -div(n)
    do j = 1, ny
       ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
       do i = 1, nx
          ! Évaluation de la condition logique
          if (norm_x(i, j) /= 0.0_dp .or. norm_y(i, j) /= 0.0_dp) then
             ! Affectation et calcul numérique
             kappa(i, j) = - ((norm_x(i+1, j) - norm_x(i-1, j)) / (2.0_dp * dx) + &
                              (norm_y(i, j+1) - norm_y(i, j-1)) / (2.0_dp * dy))
          end if
       end do
    end do

    ! Appel de sous-programme modulaire
    call exchange_halos_scalar(kappa, mesh)
    ! Affectation et calcul numérique
    fields%curvature = kappa

    ! 3. Forces de tension superficielle CSF (Balanced-Force)
    do j = 1, ny
       ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
       do i = 1, nx + 1
          ! Affectation et calcul numérique
          kap_avg = 0.5_dp * (kappa(i-1, j) + kappa(i, j))
          ! Affectation et calcul numérique
          grad_c  = (fields%C(i, j) - fields%C(i-1, j)) / dx
          ! Affectation et calcul numérique
          fields%fx_sigma(i, j) = cfg%sigma * kap_avg * grad_c
       end do
    end do

    ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
    do j = 1, ny + 1
       ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
       do i = 1, nx
          ! Affectation et calcul numérique
          kap_avg = 0.5_dp * (kappa(i, j-1) + kappa(i, j))
          ! Affectation et calcul numérique
          grad_c  = (fields%C(i, j) - fields%C(i, j-1)) / dy
          ! Affectation et calcul numérique
          fields%fy_sigma(i, j) = cfg%sigma * kap_avg * grad_c
       end do
    end do

    ! Libération de la mémoire allouée
    deallocate(norm_x, norm_y, kappa)
  end subroutine compute_curvature_and_csf

  ! [Procédure/Module] subroutine advect_vof
  subroutine advect_vof(fields, mesh, cfg, dt)
    type(flow_fields_t), intent(inout) :: fields
    type(mesh_2d_t), intent(in) :: mesh
    type(config_t), intent(in) :: cfg
    real(dp), intent(in) :: dt

    select case (trim(cfg%advection_method))
    case ("CIAM")
       call advect_vof_ciam(fields, mesh, cfg, dt)
    case ("WEYMOUTH_YUE")
       call advect_vof_weymouth_yue(fields, mesh, cfg, dt)
    case default
       call advect_vof_tvd(fields, mesh, cfg, dt)
    end select

    vof_step_count = vof_step_count + 1

    call exchange_halos_scalar(fields%C, mesh)
    call apply_vof_bc(fields, mesh, cfg)
    call update_physical_properties(fields, mesh, cfg)
    call compute_curvature_and_csf(fields, mesh, cfg)
  end subroutine advect_vof

  ! [Procédure/Module] subroutine advect_vof_tvd
  subroutine advect_vof_tvd(fields, mesh, cfg, dt)
    type(flow_fields_t), intent(inout) :: fields
    type(mesh_2d_t), intent(in) :: mesh
    type(config_t), intent(in) :: cfg
    real(dp), intent(in) :: dt

    integer :: i, j, nx, ny, i_up, i_down, j_up, j_down
    real(dp) :: dx, dy, u_face, v_face, c_l, c_r, c_b, c_t, c_u
    real(dp) :: dc_up, dc_loc, r, phi, fx_in, fx_out, fy_in, fy_out
    real(dp), allocatable :: flux_x(:, :), flux_y(:, :), c_new(:, :)

    nx = mesh%nx_local
    ny = mesh%ny_local
    dx = mesh%dx
    dy = mesh%dy

    allocate(flux_x(0:nx+1, 0:ny+1), flux_y(0:nx+1, 0:ny+1))
    allocate(c_new(0:nx+1, 0:ny+1))
    flux_x = 0.0_dp
    flux_y = 0.0_dp
    c_new = fields%C

    ! Flux selon X
    do j = 1, ny
       do i = 1, nx + 1
          u_face = fields%u(i, j)
          c_l = fields%C(i-1, j)
          c_r = fields%C(i, j)

          if (u_face >= 0.0_dp) then
             i_up = max(0, i - 2)
             c_u = fields%C(i_up, j)
             dc_up = c_l - c_u
             dc_loc = c_r - c_l
             if (abs(dc_loc) > 1.0e-12_dp) then
                r = dc_up / dc_loc
             else
                r = 0.0_dp
             end if
             phi = evaluate_limiter(cfg%tvd_limiter, r, c_l, c_u, c_r)
             flux_x(i, j) = u_face * (c_l + 0.5_dp * phi * (1.0_dp - min(1.0_dp, abs(u_face) * dt / dx)) * dc_loc)
          else
             i_down = min(nx + 1, i + 1)
             c_u = fields%C(i_down, j)
             dc_up = c_u - c_r
             dc_loc = c_r - c_l
             if (abs(dc_loc) > 1.0e-12_dp) then
                r = dc_up / dc_loc
             else
                r = 0.0_dp
             end if
             phi = evaluate_limiter(cfg%tvd_limiter, r, c_r, c_u, c_l)
             flux_x(i, j) = u_face * (c_r - 0.5_dp * phi * (1.0_dp - min(1.0_dp, abs(u_face) * dt / dx)) * dc_loc)
          end if
       end do
    end do

    ! Flux selon Y
    do j = 1, ny + 1
       do i = 1, nx
          v_face = fields%v(i, j)
          c_b = fields%C(i, j-1)
          c_t = fields%C(i, j)

          if (v_face >= 0.0_dp) then
             j_up = max(0, j - 2)
             c_u = fields%C(i, j_up)
             dc_up = c_b - c_u
             dc_loc = c_t - c_b
             if (abs(dc_loc) > 1.0e-12_dp) then
                r = dc_up / dc_loc
             else
                r = 0.0_dp
             end if
             phi = evaluate_limiter(cfg%tvd_limiter, r, c_b, c_u, c_t)
             flux_y(i, j) = v_face * (c_b + 0.5_dp * phi * (1.0_dp - min(1.0_dp, abs(v_face) * dt / dy)) * dc_loc)
          else
             j_down = min(ny + 1, j + 1)
             c_u = fields%C(i, j_down)
             dc_up = c_u - c_t
             dc_loc = c_t - c_b
             if (abs(dc_loc) > 1.0e-12_dp) then
                r = dc_up / dc_loc
             else
                r = 0.0_dp
             end if
             phi = evaluate_limiter(cfg%tvd_limiter, r, c_t, c_u, c_b)
             flux_y(i, j) = v_face * (c_t - 0.5_dp * phi * (1.0_dp - min(1.0_dp, abs(v_face) * dt / dy)) * dc_loc)
          end if
       end do
    end do

    ! Mise à jour conservative
    do j = 1, ny
       do i = 1, nx
          fx_in  = flux_x(i, j)
          fx_out = flux_x(i+1, j)
          fy_in  = flux_y(i, j)
          fy_out = flux_y(i, j+1)

          c_new(i, j) = fields%C(i, j) - (dt / dx) * (fx_out - fx_in) - (dt / dy) * (fy_out - fy_in)
          c_new(i, j) = max(0.0_dp, min(1.0_dp, c_new(i, j)))
       end do
    end do

    fields%C = c_new
    deallocate(flux_x, flux_y, c_new)
  end subroutine advect_vof_tvd

  ! [Procédure/Module] subroutine advect_vof_ciam
  subroutine advect_vof_ciam(fields, mesh, cfg, dt)
    type(flow_fields_t), intent(inout) :: fields
    type(mesh_2d_t), intent(in) :: mesh
    type(config_t), intent(in) :: cfg
    real(dp), intent(in) :: dt

    integer :: i, j, nx, ny, donor_idx
    real(dp) :: dx, dy, delta_x, delta_y, u_face, v_face
    real(dp) :: x_pos, y_pos, b_xmin, b_xmax, b_ymin, b_ymax, denom
    type(interface_segment_t), allocatable :: segments(:, :)
    real(dp), allocatable :: flux_vol(:, :), c_star(:, :), c_final(:, :)
    logical :: x_first

    nx = mesh%nx_local
    ny = mesh%ny_local
    dx = mesh%dx
    dy = mesh%dy

    allocate(segments(0:nx+1, 0:ny+1))
    allocate(flux_vol(0:nx+1, 0:ny+1))
    allocate(c_star(0:nx+1, 0:ny+1))
    allocate(c_final(0:nx+1, 0:ny+1))

    x_first = (mod(vof_step_count, 2) == 0)

    if (x_first) then
       ! 1. Passe X
       call reconstruct_interface(cfg%reconstruction_method, mesh, fields%C, segments)
       flux_vol = 0.0_dp
       do j = 1, ny
          do i = 1, nx + 1
             u_face = fields%u(i, j)
             delta_x = abs(u_face) * dt
             if (delta_x <= 1.0e-14_dp) cycle
             x_pos = x_face(mesh, i)
             if (u_face > 0.0_dp) then
                b_xmin = x_pos - delta_x; b_xmax = x_pos
                b_ymin = y_face(mesh, j); b_ymax = y_face(mesh, j + 1)
                flux_vol(i, j) = cut_box_area(b_xmin, b_xmax, b_ymin, b_ymax, &
                                             segments(i-1, j)%nx, segments(i-1, j)%ny, segments(i-1, j)%alpha)
             else
                b_xmin = x_pos; b_xmax = x_pos + delta_x
                b_ymin = y_face(mesh, j); b_ymax = y_face(mesh, j + 1)
                flux_vol(i, j) = -cut_box_area(b_xmin, b_xmax, b_ymin, b_ymax, &
                                              segments(i, j)%nx, segments(i, j)%ny, segments(i, j)%alpha)
             end if
          end do
       end do

       c_star = fields%C
       do j = 1, ny
          do i = 1, nx
             denom = dx - dt * (fields%u(i + 1, j) - fields%u(i, j))
             if (denom > 1.0e-12_dp) then
                c_star(i, j) = (fields%C(i, j) * dx - (flux_vol(i + 1, j) - flux_vol(i, j)) / dy) / denom
             else
                c_star(i, j) = fields%C(i, j) - (flux_vol(i + 1, j) - flux_vol(i, j)) / (dx * dy)
             end if
             c_star(i, j) = max(0.0_dp, min(1.0_dp, c_star(i, j)))
          end do
       end do
       call exchange_halos_scalar(c_star, mesh)

       ! 2. Passe Y
       call reconstruct_interface(cfg%reconstruction_method, mesh, c_star, segments)
       flux_vol = 0.0_dp
       do j = 1, ny + 1
          do i = 1, nx
             v_face = fields%v(i, j)
             delta_y = abs(v_face) * dt
             if (delta_y <= 1.0e-14_dp) cycle
             y_pos = y_face(mesh, j)
             if (v_face > 0.0_dp) then
                b_xmin = x_face(mesh, i); b_xmax = x_face(mesh, i + 1)
                b_ymin = y_pos - delta_y; b_ymax = y_pos
                flux_vol(i, j) = cut_box_area(b_xmin, b_xmax, b_ymin, b_ymax, &
                                             segments(i, j-1)%nx, segments(i, j-1)%ny, segments(i, j-1)%alpha)
             else
                b_xmin = x_face(mesh, i); b_xmax = x_face(mesh, i + 1)
                b_ymin = y_pos; b_ymax = y_pos + delta_y
                flux_vol(i, j) = -cut_box_area(b_xmin, b_xmax, b_ymin, b_ymax, &
                                              segments(i, j)%nx, segments(i, j)%ny, segments(i, j)%alpha)
             end if
          end do
       end do

       c_final = c_star
       do j = 1, ny
          do i = 1, nx
             denom = dy - dt * (fields%v(i, j + 1) - fields%v(i, j))
             if (denom > 1.0e-12_dp) then
                c_final(i, j) = (c_star(i, j) * dy - (flux_vol(i, j + 1) - flux_vol(i, j)) / dx) / denom
             else
                c_final(i, j) = c_star(i, j) - (flux_vol(i, j + 1) - flux_vol(i, j)) / (dx * dy)
             end if
             c_final(i, j) = max(0.0_dp, min(1.0_dp, c_final(i, j)))
          end do
       end do
    else
       ! 1. Passe Y
       call reconstruct_interface(cfg%reconstruction_method, mesh, fields%C, segments)
       flux_vol = 0.0_dp
       do j = 1, ny + 1
          do i = 1, nx
             v_face = fields%v(i, j)
             delta_y = abs(v_face) * dt
             if (delta_y <= 1.0e-14_dp) cycle
             y_pos = y_face(mesh, j)
             if (v_face > 0.0_dp) then
                b_xmin = x_face(mesh, i); b_xmax = x_face(mesh, i + 1)
                b_ymin = y_pos - delta_y; b_ymax = y_pos
                flux_vol(i, j) = cut_box_area(b_xmin, b_xmax, b_ymin, b_ymax, &
                                             segments(i, j-1)%nx, segments(i, j-1)%ny, segments(i, j-1)%alpha)
             else
                b_xmin = x_face(mesh, i); b_xmax = x_face(mesh, i + 1)
                b_ymin = y_pos; b_ymax = y_pos + delta_y
                flux_vol(i, j) = -cut_box_area(b_xmin, b_xmax, b_ymin, b_ymax, &
                                              segments(i, j)%nx, segments(i, j)%ny, segments(i, j)%alpha)
             end if
          end do
       end do

       c_star = fields%C
       do j = 1, ny
          do i = 1, nx
             denom = dy - dt * (fields%v(i, j + 1) - fields%v(i, j))
             if (denom > 1.0e-12_dp) then
                c_star(i, j) = (fields%C(i, j) * dy - (flux_vol(i, j + 1) - flux_vol(i, j)) / dx) / denom
             else
                c_star(i, j) = fields%C(i, j) - (flux_vol(i, j + 1) - flux_vol(i, j)) / (dx * dy)
             end if
             c_star(i, j) = max(0.0_dp, min(1.0_dp, c_star(i, j)))
          end do
       end do
       call exchange_halos_scalar(c_star, mesh)

       ! 2. Passe X
       call reconstruct_interface(cfg%reconstruction_method, mesh, c_star, segments)
       flux_vol = 0.0_dp
       do j = 1, ny
          do i = 1, nx + 1
             u_face = fields%u(i, j)
             delta_x = abs(u_face) * dt
             if (delta_x <= 1.0e-14_dp) cycle
             x_pos = x_face(mesh, i)
             if (u_face > 0.0_dp) then
                b_xmin = x_pos - delta_x; b_xmax = x_pos
                b_ymin = y_face(mesh, j); b_ymax = y_face(mesh, j + 1)
                flux_vol(i, j) = cut_box_area(b_xmin, b_xmax, b_ymin, b_ymax, &
                                             segments(i-1, j)%nx, segments(i-1, j)%ny, segments(i-1, j)%alpha)
             else
                b_xmin = x_pos; b_xmax = x_pos + delta_x
                b_ymin = y_face(mesh, j); b_ymax = y_face(mesh, j + 1)
                flux_vol(i, j) = -cut_box_area(b_xmin, b_xmax, b_ymin, b_ymax, &
                                              segments(i, j)%nx, segments(i, j)%ny, segments(i, j)%alpha)
             end if
          end do
       end do

       c_final = c_star
       do j = 1, ny
          do i = 1, nx
             denom = dx - dt * (fields%u(i + 1, j) - fields%u(i, j))
             if (denom > 1.0e-12_dp) then
                c_final(i, j) = (c_star(i, j) * dx - (flux_vol(i + 1, j) - flux_vol(i, j)) / dy) / denom
             else
                c_final(i, j) = c_star(i, j) - (flux_vol(i + 1, j) - flux_vol(i, j)) / (dx * dy)
             end if
             c_final(i, j) = max(0.0_dp, min(1.0_dp, c_final(i, j)))
          end do
       end do
    end if

    fields%C = c_final
    deallocate(segments, flux_vol, c_star, c_final)
  end subroutine advect_vof_ciam

  ! [Procédure/Module] subroutine advect_vof_weymouth_yue
  subroutine advect_vof_weymouth_yue(fields, mesh, cfg, dt)
    type(flow_fields_t), intent(inout) :: fields
    type(mesh_2d_t), intent(in) :: mesh
    type(config_t), intent(in) :: cfg
    real(dp), intent(in) :: dt

    integer :: i, j, nx, ny
    real(dp) :: dx, dy, cell_area, u_face, v_face, delta_x, delta_y
    real(dp) :: x_pos, y_pos, b_xmin, b_xmax, b_ymin, b_ymax
    logical :: x_first
    type(interface_segment_t), allocatable :: segments(:, :)
    real(dp), allocatable :: flux_vol(:, :), c_star(:, :), c_final(:, :), c_d(:, :)

    nx = mesh%nx_local
    ny = mesh%ny_local
    dx = mesh%dx
    dy = mesh%dy
    cell_area = dx * dy
    x_first = (mod(vof_step_count, 2) == 0)

    allocate(segments(0:nx+1, 0:ny+1))
    allocate(flux_vol(0:nx+1, 0:ny+1))
    allocate(c_star(0:nx+1, 0:ny+1))
    allocate(c_final(0:nx+1, 0:ny+1))
    allocate(c_d(0:nx+1, 0:ny+1))

    ! Fonction indicatrice binaire de dilatation figée au temps n (Weymouth & Yue, JCP 2010, Eq. 14)
    do j = 0, ny + 1
       do i = 0, nx + 1
          if (fields%C(i, j) >= 0.5_dp) then
             c_d(i, j) = 1.0_dp
          else
             c_d(i, j) = 0.0_dp
          end if
       end do
    end do

    if (x_first) then
       ! 1. Passe X (Eq. 15 Weymouth & Yue 2010)
       call reconstruct_interface(cfg%reconstruction_method, mesh, fields%C, segments)
       flux_vol = 0.0_dp
       do j = 1, ny
          do i = 1, nx + 1
             u_face = fields%u(i, j)
             delta_x = abs(u_face) * dt
             if (delta_x <= 1.0e-14_dp) cycle
             x_pos = x_face(mesh, i)
             if (u_face > 0.0_dp) then
                b_xmin = x_pos - delta_x; b_xmax = x_pos
                b_ymin = y_face(mesh, j); b_ymax = y_face(mesh, j + 1)
                flux_vol(i, j) = cut_box_area(b_xmin, b_xmax, b_ymin, b_ymax, &
                                             segments(i-1, j)%nx, segments(i-1, j)%ny, segments(i-1, j)%alpha)
             else
                b_xmin = x_pos; b_xmax = x_pos + delta_x
                b_ymin = y_face(mesh, j); b_ymax = y_face(mesh, j + 1)
                flux_vol(i, j) = -cut_box_area(b_xmin, b_xmax, b_ymin, b_ymax, &
                                              segments(i, j)%nx, segments(i, j)%ny, segments(i, j)%alpha)
             end if
          end do
       end do

       c_star = fields%C
       do j = 1, ny
          do i = 1, nx
             c_star(i, j) = fields%C(i, j) - (flux_vol(i + 1, j) - flux_vol(i, j)) / cell_area + &
                            c_d(i, j) * (dt / dx) * (fields%u(i + 1, j) - fields%u(i, j))
          end do
       end do
       call exchange_halos_scalar(c_star, mesh)

       ! 2. Passe Y (Eq. 16 Weymouth & Yue 2010)
       call reconstruct_interface(cfg%reconstruction_method, mesh, c_star, segments)
       flux_vol = 0.0_dp
       do j = 1, ny + 1
          do i = 1, nx
             v_face = fields%v(i, j)
             delta_y = abs(v_face) * dt
             if (delta_y <= 1.0e-14_dp) cycle
             y_pos = y_face(mesh, j)
             if (v_face > 0.0_dp) then
                b_xmin = x_face(mesh, i); b_xmax = x_face(mesh, i + 1)
                b_ymin = y_pos - delta_y; b_ymax = y_pos
                flux_vol(i, j) = cut_box_area(b_xmin, b_xmax, b_ymin, b_ymax, &
                                             segments(i, j-1)%nx, segments(i, j-1)%ny, segments(i, j-1)%alpha)
             else
                b_xmin = x_face(mesh, i); b_xmax = x_face(mesh, i + 1)
                b_ymin = y_pos; b_ymax = y_pos + delta_y
                flux_vol(i, j) = -cut_box_area(b_xmin, b_xmax, b_ymin, b_ymax, &
                                              segments(i, j)%nx, segments(i, j)%ny, segments(i, j)%alpha)
             end if
          end do
       end do

       c_final = c_star
       do j = 1, ny
          do i = 1, nx
             c_final(i, j) = c_star(i, j) - (flux_vol(i, j + 1) - flux_vol(i, j)) / cell_area + &
                             c_d(i, j) * (dt / dy) * (fields%v(i, j + 1) - fields%v(i, j))
             c_final(i, j) = max(0.0_dp, min(1.0_dp, c_final(i, j)))
          end do
       end do
    else
       ! 1. Passe Y (Eq. 16 Weymouth & Yue 2010)
       call reconstruct_interface(cfg%reconstruction_method, mesh, fields%C, segments)
       flux_vol = 0.0_dp
       do j = 1, ny + 1
          do i = 1, nx
             v_face = fields%v(i, j)
             delta_y = abs(v_face) * dt
             if (delta_y <= 1.0e-14_dp) cycle
             y_pos = y_face(mesh, j)
             if (v_face > 0.0_dp) then
                b_xmin = x_face(mesh, i); b_xmax = x_face(mesh, i + 1)
                b_ymin = y_pos - delta_y; b_ymax = y_pos
                flux_vol(i, j) = cut_box_area(b_xmin, b_xmax, b_ymin, b_ymax, &
                                             segments(i, j-1)%nx, segments(i, j-1)%ny, segments(i, j-1)%alpha)
             else
                b_xmin = x_face(mesh, i); b_xmax = x_face(mesh, i + 1)
                b_ymin = y_pos; b_ymax = y_pos + delta_y
                flux_vol(i, j) = -cut_box_area(b_xmin, b_xmax, b_ymin, b_ymax, &
                                              segments(i, j)%nx, segments(i, j)%ny, segments(i, j)%alpha)
             end if
          end do
       end do

       c_star = fields%C
       do j = 1, ny
          do i = 1, nx
             c_star(i, j) = fields%C(i, j) - (flux_vol(i, j + 1) - flux_vol(i, j)) / cell_area + &
                            c_d(i, j) * (dt / dy) * (fields%v(i, j + 1) - fields%v(i, j))
          end do
       end do
       call exchange_halos_scalar(c_star, mesh)

       ! 2. Passe X (Eq. 15 Weymouth & Yue 2010)
       call reconstruct_interface(cfg%reconstruction_method, mesh, c_star, segments)
       flux_vol = 0.0_dp
       do j = 1, ny
          do i = 1, nx + 1
             u_face = fields%u(i, j)
             delta_x = abs(u_face) * dt
             if (delta_x <= 1.0e-14_dp) cycle
             x_pos = x_face(mesh, i)
             if (u_face > 0.0_dp) then
                b_xmin = x_pos - delta_x; b_xmax = x_pos
                b_ymin = y_face(mesh, j); b_ymax = y_face(mesh, j + 1)
                flux_vol(i, j) = cut_box_area(b_xmin, b_xmax, b_ymin, b_ymax, &
                                             segments(i-1, j)%nx, segments(i-1, j)%ny, segments(i-1, j)%alpha)
             else
                b_xmin = x_pos; b_xmax = x_pos + delta_x
                b_ymin = y_face(mesh, j); b_ymax = y_face(mesh, j + 1)
                flux_vol(i, j) = -cut_box_area(b_xmin, b_xmax, b_ymin, b_ymax, &
                                              segments(i, j)%nx, segments(i, j)%ny, segments(i, j)%alpha)
             end if
          end do
       end do

       c_final = c_star
       do j = 1, ny
          do i = 1, nx
             c_final(i, j) = c_star(i, j) - (flux_vol(i + 1, j) - flux_vol(i, j)) / cell_area + &
                             c_d(i, j) * (dt / dx) * (fields%u(i + 1, j) - fields%u(i, j))
             c_final(i, j) = max(0.0_dp, min(1.0_dp, c_final(i, j)))
          end do
       end do
    end if

    fields%C = c_final
    deallocate(segments, flux_vol, c_star, c_final, c_d)
  end subroutine advect_vof_weymouth_yue

  ! [Procédure/Module] subroutine vof_init
  subroutine vof_init(this, cfg)
    ! Déclaration de variable locale / paramètre de calcul
    class(interface_vof_t), intent(inout) :: this
    ! Déclaration de variable locale / paramètre de calcul
    type(config_t), intent(in) :: cfg
    ! Affectation et calcul numérique
    this%cfg = cfg
  end subroutine vof_init

  ! [Procédure/Module] subroutine vof_set_ic
  subroutine vof_set_ic(this, mesh, fields)
    ! Déclaration de variable locale / paramètre de calcul
    class(interface_vof_t), intent(in) :: this
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(inout) :: fields
    ! Appel de sous-programme modulaire
    call set_initial_conditions(fields, mesh, this%cfg)
    ! Appel de sous-programme modulaire
    call update_physical_properties(fields, mesh, this%cfg)
    ! Appel de sous-programme modulaire
    call compute_curvature_and_csf(fields, mesh, this%cfg)
  end subroutine vof_set_ic

  ! [Procédure/Module] subroutine vof_advect_step
  subroutine vof_advect_step(this, mesh, fields, dt)
    ! Déclaration de variable locale / paramètre de calcul
    class(interface_vof_t), intent(in) :: this
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in) :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(inout) :: fields
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in) :: dt
    ! Appel de sous-programme modulaire
    call advect_vof(fields, mesh, this%cfg, dt)
    ! Appel de sous-programme modulaire
    call update_physical_properties(fields, mesh, this%cfg)
    ! Appel de sous-programme modulaire
    call compute_curvature_and_csf(fields, mesh, this%cfg)
  end subroutine vof_advect_step

end module module_interface_vof
