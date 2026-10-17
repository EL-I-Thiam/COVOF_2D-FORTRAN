module module_interface_vof
  use module_precision, only: dp
  use module_constants
  use module_mesh, only: mesh_2d_t, x_center, y_center, x_face, y_face, is_in_obstacle
  use module_fields, only: flow_fields_t
  use module_config, only: config_t
  use module_mpi_comm, only: exchange_halos_scalar, apply_vof_bc
  implicit none
  private

  public :: interface_vof_t
  public :: set_initial_conditions
  public :: advect_vof
  public :: update_physical_properties
  public :: compute_curvature_and_csf

  type :: interface_vof_t
     type(config_t) :: cfg
  contains
     procedure :: init => vof_init
     procedure :: set_initial_condition => vof_set_ic
     procedure :: advect => vof_advect_step
  end type interface_vof_t

contains

  subroutine set_initial_conditions(fields, mesh, cfg)
    type(flow_fields_t), intent(inout) :: fields
    type(mesh_2d_t), intent(in) :: mesh
    type(config_t), intent(in) :: cfg

    integer :: i, j, nx, ny
    real(dp) :: x, y, dist, eps, int_y

    nx = mesh%nx_local
    ny = mesh%ny_local
    eps = 1.5_dp * max(mesh%dx, mesh%dy)

    ! 1. Initialisation des vitesses (par défaut 0 ou imposé)
    fields%u = cfg%ic_u
    fields%v = cfg%ic_v

    ! 2. Initialisation de la fraction volumique C
    do j = 0, ny + 1
       y = y_center(mesh, j)
       do i = 0, nx + 1
          x = x_center(mesh, i)

          select case (cfg%ic_type)
          case (IC_UNIFORM)
             fields%C(i, j) = 1.0_dp

          case (IC_DROPLET)
             dist = sqrt((x - cfg%ic_center_x)**2 + (y - cfg%ic_center_y)**2) - cfg%ic_radius
             if (dist < -eps) then
                fields%C(i, j) = 1.0_dp
             else if (dist > eps) then
                fields%C(i, j) = 0.0_dp
             else
                fields%C(i, j) = 0.5_dp - dist / (2.0_dp * eps) - (1.0_dp / (2.0_dp * PI)) * sin(PI * dist / eps)
             end if

          case (IC_BUBBLE)
             dist = sqrt((x - cfg%ic_center_x)**2 + (y - cfg%ic_center_y)**2) - cfg%ic_radius
             if (dist < -eps) then
                fields%C(i, j) = 0.0_dp ! Intérieur bulle (fluide 2)
             else if (dist > eps) then
                fields%C(i, j) = 1.0_dp ! Extérieur (fluide 1)
             else
                fields%C(i, j) = 0.5_dp + dist / (2.0_dp * eps) + (1.0_dp / (2.0_dp * PI)) * sin(PI * dist / eps)
             end if

          case (IC_DAM_BREAK)
             if (x >= cfg%ic_box_x0 .and. x <= cfg%ic_box_x1 .and. &
                 y >= cfg%ic_box_y0 .and. y <= cfg%ic_box_y1) then
                fields%C(i, j) = 1.0_dp
             else
                fields%C(i, j) = 0.0_dp
             end if

          case (IC_VORTEX_SHEAR)
             dist = sqrt((x - cfg%ic_center_x)**2 + (y - cfg%ic_center_y)**2) - cfg%ic_radius
             if (dist < -eps) then
                fields%C(i, j) = 1.0_dp
             else if (dist > eps) then
                fields%C(i, j) = 0.0_dp
             else
                fields%C(i, j) = 0.5_dp - dist / (2.0_dp * eps) - (1.0_dp / (2.0_dp * PI)) * sin(PI * dist / eps)
             end if

          case (IC_RAYLEIGH_TAYLOR)
             int_y = cfg%ic_center_y + cfg%ic_pert_amp * cos(2.0_dp * PI * cfg%ic_pert_k * x / mesh%lx)
             if (y >= int_y) then
                fields%C(i, j) = 1.0_dp ! Fluide lourd au-dessus
             else
                fields%C(i, j) = 0.0_dp ! Fluide léger en dessous
             end if

          case (IC_LAPLACIAN_TEST)
             fields%C(i, j) = 1.0_dp

          case default
             fields%C(i, j) = 1.0_dp
          end select
       end do
    end do

    ! 3. Détection des obstacles solides pour la pénalisation de Brinkman
    fields%solid_mask = 0
    if (cfg%obs_enabled) then
       do j = 0, ny + 1
          y = y_center(mesh, j)
          do i = 0, nx + 1
             x = x_center(mesh, i)
             if (is_in_obstacle(x, y, cfg%obs_type, cfg%obs_x0, cfg%obs_y0, cfg%obs_x1, cfg%obs_y1, &
                                cfg%obs_xc, cfg%obs_yc, cfg%obs_radius)) then
                fields%solid_mask(i, j) = 1
             end if
          end do
       end do
    end if

    call exchange_halos_scalar(fields%C, mesh)
    call apply_vof_bc(fields, mesh, cfg)
    call update_physical_properties(fields, mesh, cfg)
    call compute_curvature_and_csf(fields, mesh, cfg)
  end subroutine set_initial_conditions

  subroutine update_physical_properties(fields, mesh, cfg)
    type(flow_fields_t), intent(inout) :: fields
    type(mesh_2d_t), intent(in) :: mesh
    type(config_t), intent(in) :: cfg

    integer :: i, j, nx, ny
    real(dp) :: c, y, z, mu_val, m_sw, m_se, m_nw, m_ne

    nx = mesh%nx_local
    ny = mesh%ny_local

    if (.not. cfg%enable_two_phase) then
       fields%rho = cfg%rho1
       fields%rho_u = cfg%rho1
       fields%rho_v = cfg%rho1
       if (cfg%viscosity_model == VISC_EXPONENTIAL) then
          do j = 0, ny + 1
             y = y_center(mesh, j)
             z = y - cfg%visc_e
             mu_val = cfg%visc_mu0 * exp(cfg%visc_a * (z * z) / (cfg%visc_e * cfg%visc_e))
             fields%mu(:, j) = mu_val
          end do
          do j = 0, ny + 1
             y = y_face(mesh, j)
             z = y - cfg%visc_e
             mu_val = cfg%visc_mu0 * exp(cfg%visc_a * (z * z) / (cfg%visc_e * cfg%visc_e))
             fields%mu_node(:, j) = mu_val
          end do
       else
          fields%mu = cfg%mu1
          fields%mu_node = cfg%mu1
       end if
       return
    end if

    ! Coins de C
    fields%C(0, 0)       = fields%C(1, 1)
    fields%C(nx+1, 0)    = fields%C(nx, 1)
    fields%C(0, ny+1)    = fields%C(1, ny)
    fields%C(nx+1, ny+1) = fields%C(nx, ny)

    ! 1. Masse volumique arithmétique et viscosité au centre
    do j = 0, ny + 1
       do i = 0, nx + 1
          c = max(0.0_dp, min(1.0_dp, fields%C(i, j)))
          fields%rho(i, j) = c * cfg%rho1 + (1.0_dp - c) * cfg%rho2
          fields%mu(i, j)  = c * cfg%mu1  + (1.0_dp - c) * cfg%mu2
       end do
    end do

    fields%mu(0, 0)       = fields%mu(1, 1)
    fields%mu(nx+1, 0)    = fields%mu(nx, 1)
    fields%mu(0, ny+1)    = fields%mu(1, ny)
    fields%mu(nx+1, ny+1) = fields%mu(nx, ny)

    ! 2. Masse volumique sur les faces u
    do j = 1, ny
       do i = 1, nx + 1
          fields%rho_u(i, j) = 0.5_dp * (fields%rho(i - 1, j) + fields%rho(i, j))
       end do
    end do

    ! 3. Masse volumique sur les faces v
    do j = 1, ny + 1
       do i = 1, nx
          fields%rho_v(i, j) = 0.5_dp * (fields%rho(i, j - 1) + fields%rho(i, j))
       end do
    end do

    ! 4. Viscosité harmonique aux nœuds (coins de mailles)
    do j = 1, ny + 1
       do i = 1, nx + 1
          m_sw = fields%mu(i-1, j-1)
          m_se = fields%mu(i, j-1)
          m_nw = fields%mu(i-1, j)
          m_ne = fields%mu(i, j)
          fields%mu_node(i, j) = 4.0_dp / (1.0_dp / m_sw + 1.0_dp / m_se + 1.0_dp / m_nw + 1.0_dp / m_ne)
       end do
    end do
  end subroutine update_physical_properties

  pure real(dp) function superbee(r)
    real(dp), intent(in) :: r
    superbee = max(0.0_dp, max(min(2.0_dp * r, 1.0_dp), min(r, 2.0_dp)))
  end function superbee

  subroutine compute_curvature_and_csf(fields, mesh, cfg)
    type(flow_fields_t), intent(inout) :: fields
    type(mesh_2d_t), intent(in) :: mesh
    type(config_t), intent(in) :: cfg

    integer :: i, j, nx, ny
    real(dp) :: dx, dy, dc_dx, dc_dy, grad_mag, kap_avg, grad_c
    real(dp), allocatable :: norm_x(:, :), norm_y(:, :)
    real(dp), allocatable :: kappa(:, :)

    if (.not. cfg%enable_surface_tension .or. cfg%sigma <= 0.0_dp) then
       fields%fx_sigma = 0.0_dp
       fields%fy_sigma = 0.0_dp
       fields%curvature = 0.0_dp
       return
    end if

    nx = mesh%nx_local
    ny = mesh%ny_local
    dx = mesh%dx
    dy = mesh%dy

    allocate(norm_x(0:nx+1, 0:ny+1), norm_y(0:nx+1, 0:ny+1))
    allocate(kappa(0:nx+1, 0:ny+1))
    norm_x = 0.0_dp
    norm_y = 0.0_dp
    kappa = 0.0_dp

    ! 1. Gradient et normale unitaire aux centres
    do j = 1, ny
       do i = 1, nx
          dc_dx = (fields%C(i+1, j) - fields%C(i-1, j)) / (2.0_dp * dx)
          dc_dy = (fields%C(i, j+1) - fields%C(i, j-1)) / (2.0_dp * dy)
          grad_mag = sqrt(dc_dx * dc_dx + dc_dy * dc_dy)
          if (grad_mag > 1.0e-3_dp / dx) then
             norm_x(i, j) = dc_dx / grad_mag
             norm_y(i, j) = dc_dy / grad_mag
          end if
       end do
    end do

    call exchange_halos_scalar(norm_x, mesh)
    call exchange_halos_scalar(norm_y, mesh)

    ! 2. Courbure kappa = -div(n)
    do j = 1, ny
       do i = 1, nx
          if (norm_x(i, j) /= 0.0_dp .or. norm_y(i, j) /= 0.0_dp) then
             kappa(i, j) = - ((norm_x(i+1, j) - norm_x(i-1, j)) / (2.0_dp * dx) + &
                              (norm_y(i, j+1) - norm_y(i, j-1)) / (2.0_dp * dy))
          end if
       end do
    end do

    call exchange_halos_scalar(kappa, mesh)
    fields%curvature = kappa

    ! 3. Forces de tension superficielle CSF (Balanced-Force)
    do j = 1, ny
       do i = 1, nx + 1
          kap_avg = 0.5_dp * (kappa(i-1, j) + kappa(i, j))
          grad_c  = (fields%C(i, j) - fields%C(i-1, j)) / dx
          fields%fx_sigma(i, j) = cfg%sigma * kap_avg * grad_c
       end do
    end do

    do j = 1, ny + 1
       do i = 1, nx
          kap_avg = 0.5_dp * (kappa(i, j-1) + kappa(i, j))
          grad_c  = (fields%C(i, j) - fields%C(i, j-1)) / dy
          fields%fy_sigma(i, j) = cfg%sigma * kap_avg * grad_c
       end do
    end do

    deallocate(norm_x, norm_y, kappa)
  end subroutine compute_curvature_and_csf

  subroutine advect_vof(fields, mesh, cfg, dt)
    type(flow_fields_t), intent(inout) :: fields
    type(mesh_2d_t), intent(in) :: mesh
    type(config_t), intent(in) :: cfg
    real(dp), intent(in) :: dt

    integer :: i, j, nx, ny, i_up, i_down, j_up, j_down
    real(dp) :: dx, dy, u_face, v_face, c_l, c_r, c_b, c_t
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

    ! Flux selon X avec limiteur Superbee
    do j = 1, ny
       do i = 1, nx + 1
          u_face = fields%u(i, j)
          c_l = fields%C(i-1, j)
          c_r = fields%C(i, j)

          if (u_face >= 0.0_dp) then
             i_up = max(0, i - 2)
             dc_up = c_l - fields%C(i_up, j)
             dc_loc = c_r - c_l
             if (abs(dc_loc) > 1.0e-12_dp) then
                r = dc_up / dc_loc
             else
                r = 0.0_dp
             end if
             phi = superbee(r)
             flux_x(i, j) = u_face * (c_l + 0.5_dp * phi * (1.0_dp - u_face * dt / dx) * dc_loc)
          else
             i_down = min(nx + 1, i + 1)
             dc_up = fields%C(i_down, j) - c_r
             dc_loc = c_r - c_l
             if (abs(dc_loc) > 1.0e-12_dp) then
                r = dc_up / dc_loc
             else
                r = 0.0_dp
             end if
             phi = superbee(r)
             flux_x(i, j) = u_face * (c_r - 0.5_dp * phi * (1.0_dp + u_face * dt / dx) * dc_loc)
          end if
       end do
    end do

    ! Flux selon Y avec limiteur Superbee
    do j = 1, ny + 1
       do i = 1, nx
          v_face = fields%v(i, j)
          c_b = fields%C(i, j-1)
          c_t = fields%C(i, j)

          if (v_face >= 0.0_dp) then
             j_up = max(0, j - 2)
             dc_up = c_b - fields%C(i, j_up)
             dc_loc = c_t - c_b
             if (abs(dc_loc) > 1.0e-12_dp) then
                r = dc_up / dc_loc
             else
                r = 0.0_dp
             end if
             phi = superbee(r)
             flux_y(i, j) = v_face * (c_b + 0.5_dp * phi * (1.0_dp - v_face * dt / dy) * dc_loc)
          else
             j_down = min(ny + 1, j + 1)
             dc_up = fields%C(i, j_down) - c_t
             dc_loc = c_t - c_b
             if (abs(dc_loc) > 1.0e-12_dp) then
                r = dc_up / dc_loc
             else
                r = 0.0_dp
             end if
             phi = superbee(r)
             flux_y(i, j) = v_face * (c_t - 0.5_dp * phi * (1.0_dp + v_face * dt / dy) * dc_loc)
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
    call exchange_halos_scalar(fields%C, mesh)
    call apply_vof_bc(fields, mesh, cfg)

    call update_physical_properties(fields, mesh, cfg)
    call compute_curvature_and_csf(fields, mesh, cfg)

    deallocate(flux_x, flux_y, c_new)
  end subroutine advect_vof

  subroutine vof_init(this, cfg)
    class(interface_vof_t), intent(inout) :: this
    type(config_t), intent(in) :: cfg
    this%cfg = cfg
  end subroutine vof_init

  subroutine vof_set_ic(this, mesh, fields)
    class(interface_vof_t), intent(in) :: this
    type(mesh_2d_t), intent(in) :: mesh
    type(flow_fields_t), intent(inout) :: fields
    call set_initial_conditions(fields, mesh, this%cfg)
    call update_physical_properties(fields, mesh, this%cfg)
    call compute_curvature_and_csf(fields, mesh, this%cfg)
  end subroutine vof_set_ic

  subroutine vof_advect_step(this, mesh, fields, dt)
    class(interface_vof_t), intent(in) :: this
    type(mesh_2d_t), intent(in) :: mesh
    type(flow_fields_t), intent(inout) :: fields
    real(dp), intent(in) :: dt
    call advect_vof(fields, mesh, this%cfg, dt)
    call update_physical_properties(fields, mesh, this%cfg)
    call compute_curvature_and_csf(fields, mesh, this%cfg)
  end subroutine vof_advect_step

end module module_interface_vof
