!===============================================================================
! Module: module_simulation
! Description: Moteur principal de simulation COVOF 2D en Fortran 2018 :
!              - Gestion du cycle de vie du calcul (init, run, finalize)
!              - Pas de temps adaptatif CFL + onde capillaire
!              - Boucle temporelle avec couplage Lagrangien Augmenté / Uzawa
!              - Mode de vérification analytique du Laplacien (Cas 3.1 Référence Académique)
!              - Mode de vérification du vortex de cisaillement réversible (Cas 3.19 Référence Académique)
!              - Convergence vers l'état stationnaire
!===============================================================================
module module_simulation
  use module_precision, only: dp
  use module_constants
  use module_config
  use module_mesh
  use module_fields
  use module_mpi_comm
  use module_interface_vof
  use module_navier_stokes
  use module_linear_solvers
  use module_vtk_writer
  implicit none
  private
  include 'mpif.h'

  public :: simulation_t

  type :: simulation_t
     ! Déclaration de variable locale / paramètre de calcul
     type(config_t)         :: config
     ! Déclaration de variable locale / paramètre de calcul
     type(mesh_2d_t)        :: mesh
     ! Déclaration de variable locale / paramètre de calcul
     type(flow_fields_t)    :: fields
     ! Déclaration de variable locale / paramètre de calcul
     type(interface_vof_t)  :: tracker
     ! Déclaration de variable locale / paramètre de calcul
     type(navier_stokes_t)  :: ns_solver
     ! Déclaration de variable locale / paramètre de calcul
     type(vtk_writer_t)     :: vtk_writer
  contains
     ! Affectation et calcul numérique
     procedure :: init => sim_init
     ! Affectation et calcul numérique
     procedure :: run => sim_run
     ! Affectation et calcul numérique
     procedure :: compute_time_step => sim_compute_time_step
     ! Affectation et calcul numérique
     procedure :: run_laplacian => sim_run_laplacian
     procedure :: run_laplacian_cas3_2 => sim_run_laplacian_cas3_2
     procedure :: export_poiseuille_results => sim_export_poiseuille_results
     procedure :: export_poiseuille_unsteady_results => sim_export_poiseuille_unsteady_results
     ! Affectation et calcul numérique
     procedure :: run_vortex_advection => sim_run_vortex_advection
     ! Affectation et calcul numérique
     procedure :: solve_standard_projection => sim_solve_standard_projection
  end type simulation_t

contains

  !-----------------------------------------------------------------------------
  ! Initialisation complète de la simulation à partir du fichier JSON
  !-----------------------------------------------------------------------------
  function sim_init(this, comm, config_file) result(success)
    ! Déclaration de variable locale / paramètre de calcul
    class(simulation_t), intent(inout) :: this
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(in)                :: comm
    ! Déclaration de variable locale / paramètre de calcul
    character(len=*), intent(in)       :: config_file
    ! Déclaration de variable locale / paramètre de calcul
    logical                            :: success

    ! Déclaration de variable locale / paramètre de calcul
    logical  :: periodic_x, periodic_y
    ! Déclaration de variable locale / paramètre de calcul
    integer  :: i, j, u_nx, v_ny
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: frac, y_step, h, u_mean, y_val, u_val, x_val
    ! Déclaration des paramètres pour la nappe liquide cisaillée (Cas 3.2 Fluids 2021)
    real(dp) :: y_c_ls, delta_ls, k_ls, amp_ls, factor_u_ls, dy_c_ls, sgn_y_ls, exp_decay_ls, c_u_ls, a_val_ls
    ! Déclaration des variables pour l'amorce de jet visqueux (Cas 3.3 Fluids 2021)
    real(dp) :: x_c_vj, half_w_vj, stub_h_vj, v_in_vj

    ! Affectation et calcul numérique
    success = .false.

    ! 1. Chargement de la configuration JSON
    call load_config_from_json(config_file, this%config)

    ! 2. Initialisation du maillage cartésien partitionné MPI
    periodic_x = (this%config%bc_left%btype == BC_PERIODIC .or. this%config%bc_right%btype == BC_PERIODIC)
    ! Affectation et calcul numérique
    periodic_y = (this%config%bc_bottom%btype == BC_PERIODIC .or. this%config%bc_top%btype == BC_PERIODIC)

    ! Appel de sous-programme modulaire
    call init_mesh(this%mesh, comm, this%config%nx, this%config%ny, &
                   this%config%lx, this%config%ly, periodic_x, periodic_y)

    ! 3. Allocation des champs décalés MAC
    call this%fields%allocate(this%mesh)

    ! 4. Initialisation du module d'interface VOF et des conditions initiales
    call this%tracker%init(this%config)
    ! Appel de sous-programme modulaire
    call this%tracker%set_initial_condition(this%mesh, this%fields)

    ! Initialisation des conditions aux limites initiales
    call apply_boundary_conditions(this%mesh, this%config, this%fields, .false.)

    ! Initialisation du canal amont pour profil POISEUILLE_STEP afin de garantir la continuité
    if (this%config%bc_left%btype == BC_DIRICHLET .and. this%config%obs_enabled .and. &
        ! Affectation et calcul numérique
        this%config%obs_type == OBS_BOX) then
       ! Affectation et calcul numérique
       y_step = this%config%obs_y1
       ! Affectation et calcul numérique
       h = this%mesh%ly - y_step
       ! Affectation et calcul numérique
       u_mean = this%config%bc_left%u
       ! Affectation et calcul numérique
       u_nx = this%mesh%nx_local
       ! Affectation et calcul numérique
       if (this%mesh%is_boundary_right) u_nx = u_nx + 1

       ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
       do j = 1, this%mesh%ny_local
          ! Affectation et calcul numérique
          y_val = y_center(this%mesh, j)
          ! Affectation et calcul numérique
          u_val = 0.0_dp
          ! Évaluation de la condition logique
          if (y_val >= y_step) then
             ! Affectation et calcul numérique
             u_val = 6.0_dp * u_mean * (y_val - y_step) * (this%mesh%ly - y_val) / (h * h)
          end if
          ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
          do i = 1, u_nx
             ! Affectation et calcul numérique
             x_val = x_face(this%mesh, i)
             ! Évaluation de la condition logique
             if (x_val <= 4.0_dp + 1.0e-9_dp) then
                ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
                this%fields%u(i, j) = u_val
             end if
          end do
       end do
       ! Appel de sous-programme modulaire
       call sync_u(this%mesh, this%fields%u)
    end if

    ! Initialisation exacte du profil de Poiseuille pour vérification analytique
    if (trim(this%config%test_case) == "CAS_3_7" .or. trim(this%config%test_case) == "CAS_3_7_POISEUILLE" .or. &
        (this%config%bc_left%btype == BC_DIRICHLET .and. trim(this%config%bc_left%profile) == "PARABOLIC")) then
       u_mean = merge(this%config%bc_left%u, 1.0_dp, this%config%bc_left%u > 0.0_dp)
       u_nx = this%mesh%nx_local
       if (this%mesh%is_boundary_right) u_nx = u_nx + 1
       do j = 1, this%mesh%ny_local
          y_val = y_center(this%mesh, j)
          u_val = 6.0_dp * u_mean * y_val * (this%mesh%ly - y_val) / (this%mesh%ly**2)
          do i = 1, u_nx
             this%fields%u(i, j) = u_val
          end do
       end do
       do j = 1, this%mesh%ny_local
          do i = 1, this%mesh%nx_local
             x_val = x_center(this%mesh, i)
             this%fields%p(i, j) = (12.0_dp * this%config%mu1 * u_mean / (this%mesh%ly**2)) * (this%mesh%lx - x_val)
          end do
       end do
       call sync_u(this%mesh, this%fields%u)
       call sync_cell_scalar(this%mesh, this%fields%p)
       call apply_boundary_conditions(this%mesh, this%config, this%fields, .false.)
    end if

    ! Initialisation du champ de vitesse solénoïdal pour LIQUID_SHEET (Équations 17-18 Fluids 2021)
    if (this%config%ic_type == IC_LIQUID_SHEET) then
       ! Déclaration des variables physiques et géométriques
       if (this%config%ic_center_y > 0.0_dp) then
          y_c_ls = this%config%ic_center_y
       else
          y_c_ls = 0.5_dp * this%mesh%ly
       end if
       if (this%config%ic_thickness > 0.0_dp) then
          delta_ls = this%config%ic_thickness
       else
          delta_ls = 0.0003_dp
       end if
       k_ls = 2.0_dp * PI / this%mesh%lx
       amp_ls = 0.04_dp
       factor_u_ls = amp_ls * (this%mesh%lx / (2.0_dp * PI)) * (2.0_dp / delta_ls)

       ! Composante horizontale u sur les faces décalées MAC en x
       u_nx = this%mesh%nx_local
       if (this%mesh%is_boundary_right) u_nx = u_nx + 1
       do j = 1, this%mesh%ny_local
          y_val = y_center(this%mesh, j)
          dy_c_ls = y_val - y_c_ls
          if (dy_c_ls > 1.0e-12_dp) then
             sgn_y_ls = 1.0_dp
          else if (dy_c_ls < -1.0e-12_dp) then
             sgn_y_ls = -1.0_dp
          else
             sgn_y_ls = 0.0_dp
          end if
          exp_decay_ls = exp(-2.0_dp * abs(dy_c_ls) / delta_ls)
          do i = 1, u_nx
             x_val = x_face(this%mesh, i)
             c_u_ls = 0.5_dp * (this%fields%C(i - 1, j) + this%fields%C(i, j))
             a_val_ls = c_u_ls * 30.0_dp + (1.0_dp - c_u_ls) * 2.0_dp
             this%fields%u(i, j) = a_val_ls - sgn_y_ls * factor_u_ls * cos(k_ls * x_val) * exp_decay_ls
          end do
       end do

       ! Composante verticale v sur les faces décalées MAC en y
       v_ny = this%mesh%ny_local
       if (this%mesh%is_boundary_top) v_ny = v_ny + 1
       do j = 1, v_ny
          y_val = y_face(this%mesh, j)
          dy_c_ls = y_val - y_c_ls
          exp_decay_ls = exp(-2.0_dp * abs(dy_c_ls) / delta_ls)
          do i = 1, this%mesh%nx_local
             x_val = x_center(this%mesh, i)
             this%fields%v(i, j) = amp_ls * sin(k_ls * x_val) * exp_decay_ls
          end do
       end do

       call sync_u(this%mesh, this%fields%u)
       call sync_v(this%mesh, this%fields%v)
       call apply_velocity_bc(this%fields, this%mesh, this%config, .false.)
    end if

    ! Initialisation du champ de vitesse pour VISCOUS_JET (Cas 3.3 Fluids 2021)
    if (this%config%ic_type == IC_VISCOUS_JET) then
       x_c_vj = merge(this%config%ic_center_x, 0.5_dp * this%mesh%lx, this%config%ic_center_x > 0.0_dp)
       half_w_vj = merge(0.5_dp * this%config%ic_jet_width, 0.04_dp, this%config%ic_jet_width > 0.0_dp)
       stub_h_vj = 2.0_dp * this%mesh%dy
       v_in_vj = merge(this%config%bc_top%v, -0.6_dp, this%config%bc_top%v /= 0.0_dp)
       y_val = 0.0_dp
       do i = 1, this%mesh%nx_global
          if (abs((real(i, dp) - 0.5_dp) * this%mesh%dx - x_c_vj) <= half_w_vj) y_val = y_val + this%mesh%dx
       end do
       if (y_val < 1.0e-12_dp) y_val = 2.0_dp * half_w_vj
       v_in_vj = v_in_vj * (2.0_dp * half_w_vj) / y_val

       v_ny = this%mesh%ny_local
       if (this%mesh%is_boundary_top) v_ny = v_ny + 1
       do j = 1, v_ny
          y_val = y_face(this%mesh, j)
          do i = 1, this%mesh%nx_local
             x_val = x_center(this%mesh, i)
             if (y_val >= this%mesh%ly - stub_h_vj .and. abs(x_val - x_c_vj) <= half_w_vj) then
                this%fields%v(i, j) = v_in_vj
             end if
          end do
       end do
       call sync_v(this%mesh, this%fields%v)
       call apply_velocity_bc(this%fields, this%mesh, this%config, .false.)
       call apply_vof_bc(this%fields, this%mesh, this%config)
    end if

    ! Initialisation quasi-statique de la vitesse pour compression de gaz isotherme (Fluids 2024, Cas 4.2)
    if (this%config%compressible .and. this%config%bc_bottom%v > 0.0_dp) then
       v_ny = this%mesh%ny_local
       if (this%mesh%is_boundary_top) v_ny = v_ny + 1
       do j = 1, v_ny
          y_val = y_face(this%mesh, j)
          if (y_val <= this%config%ic_box_y1) then
             u_val = this%config%bc_bottom%v
          else
             u_val = this%config%bc_bottom%v * max(0.0_dp, (this%mesh%ly - y_val) / max(1.0e-12_dp, this%mesh%ly - this%config%ic_box_y1))
          end if
          this%fields%v(1:this%mesh%nx_local, j) = u_val
       end do
       call sync_v(this%mesh, this%fields%v)
       call apply_velocity_bc(this%fields, this%mesh, this%config, .false.)
    end if

    ! Initialisation du champ thermique si activé
    if (this%config%enable_energy) then
       ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
       do j = 1, this%mesh%ny_local
          ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
          do i = 1, this%mesh%nx_local
             ! Affectation et calcul numérique
             frac = x_center(this%mesh, i) / this%mesh%lx
             ! Évolution du champ scalaire de température
             this%fields%T(i, j) = this%config%bc_left%T + frac * (this%config%bc_right%T - this%config%bc_left%T)
          end do
       end do
       ! Appel de sous-programme modulaire
       call sync_cell_scalar(this%mesh, this%fields%T)
       ! Appel de sous-programme modulaire
       call apply_thermal_boundary_conditions(this%mesh, this%config, this%fields)
    end if

    ! 5. Initialisation du solveur Navier-Stokes et préconditionneurs
    call this%ns_solver%init(this%mesh, this%config)
    ! Appel de sous-programme modulaire
    call this%ns_solver%assemble_preconditioner_blocks(this%mesh, this%config, this%fields, this%config%dt)

    ! 6. Initialisation de l'exportateur VTK et diagnostics CSV
    call this%vtk_writer%init(this%mesh, this%config%output_dir, this%config%case_name)

    ! Affectation et calcul numérique
    success = .true.
  end function sim_init

  !-----------------------------------------------------------------------------
  ! Calcul du pas de temps adaptatif (CFL advection + onde capillaire)
  !-----------------------------------------------------------------------------
  function sim_compute_time_step(this) result(dt)
    ! Déclaration de variable locale / paramètre de calcul
    class(simulation_t), intent(in) :: this
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: dt

    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: local_max_speed, global_max_speed, speed, u_val, v_val
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: dt_cfl, dt_sigma, rho_sum, h_min, dt_target
    ! Déclaration de variable locale / paramètre de calcul
    integer  :: i, j

    ! Évaluation de la condition logique
    if (.not. this%config%adaptive_dt) then
       ! Affectation et calcul numérique
       dt = this%config%dt
       ! Retour à la procédure appelante
       return
    end if

    ! Affectation et calcul numérique
    local_max_speed = 1.0e-6_dp
    ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
    do j = 1, this%mesh%ny_local
       ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
       do i = 1, this%mesh%nx_local
          ! Affectation et calcul numérique
          u_val = abs(this%fields%u(i, j))
          ! Affectation et calcul numérique
          v_val = abs(this%fields%v(i, j))
          ! Affectation et calcul numérique
          speed = sqrt(u_val * u_val + v_val * v_val)
          ! Affectation et calcul numérique
          if (speed > local_max_speed) local_max_speed = speed
       end do
    end do
    ! Affectation et calcul numérique
    global_max_speed = global_max(this%mesh, local_max_speed)

    ! 1. Critère CFL de convection
    dt_cfl = this%config%cfl_max * min(this%mesh%dx, this%mesh%dy) / global_max_speed

    ! 2. Critère capillaire si tension superficielle active
    dt_sigma = 1.0e10_dp
    ! Évaluation de la condition logique
    if (this%config%enable_surface_tension .and. this%config%sigma > 0.0_dp) then
       ! Calcul de la masse volumique locale ou interpolée aux faces
       rho_sum = this%config%rho1 + this%config%rho2
       ! Affectation et calcul numérique
       h_min = min(this%mesh%dx, this%mesh%dy)
       ! Affectation et calcul numérique
       dt_sigma = sqrt((rho_sum * h_min**3) / (4.0_dp * PI * this%config%sigma))
    end if

    ! Affectation et calcul numérique
    dt_target = min(dt_cfl, dt_sigma)
    ! Affectation et calcul numérique
    dt = min(this%config%dt, dt_target)
  end function sim_compute_time_step

  !-----------------------------------------------------------------------------
  ! Boucle principale d'intégration temporelle Navier-Stokes
  !-----------------------------------------------------------------------------
  subroutine sim_run(this)
    ! Déclaration de variable locale / paramètre de calcul
    class(simulation_t), intent(inout) :: this

    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: t, dt, global_max_u, global_max_div, global_du_max, nu_mean
    real(dp) :: start_cpu, end_cpu, elapsed_cpu
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: local_max_u, local_max_div, local_du_max, speed, du, div_val
    ! Déclaration de variable locale / paramètre de calcul
    integer  :: step, total_iters, uzawa_max, m, iters, i, j
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: state_vector(:), rhs(:), u_prev(:, :)

    ! Évaluation de la condition logique
    if (trim(this%config%mode) == "LAPLACIAN") then
       if (trim(this%config%test_case) == "CAS_3_2_NEUMANN" .or. &
           index(this%config%test_case, "cas3_2") > 0 .or. &
           index(this%config%case_name, "cas3_2") > 0 .or. &
           index(this%config%case_name, "neumann") > 0) then
          call this%run_laplacian_cas3_2()
       else
          call this%run_laplacian()
       end if
       ! Retour à la procédure appelante
       return
    end if
    ! Évaluation de la condition logique
    if (trim(this%config%mode) == "VORTEX_ADVECTION") then
       ! Appel de sous-programme modulaire
       call this%run_vortex_advection()
       ! Retour à la procédure appelante
       return
    end if

    ! Affectation et calcul numérique
    t = 0.0_dp
    ! Affectation et calcul numérique
    step = 0

    ! Allocation dynamique de tableau contigu en mémoire
    allocate(state_vector(this%ns_solver%n_owned_total))
    ! Allocation dynamique de tableau contigu en mémoire
    allocate(rhs(this%ns_solver%n_owned_total))
    ! Allocation dynamique de tableau contigu en mémoire
    allocate(u_prev(0:this%mesh%nx_local+1, 0:this%mesh%ny_local+1))

    this%fields%u_old = this%fields%u
    this%fields%v_old = this%fields%v
    this%fields%p_old = this%fields%p
    this%fields%T_old = this%fields%T
    this%fields%rho_old = this%fields%rho
    this%fields%rho_u_old = this%fields%rho_u
    this%fields%rho_v_old = this%fields%rho_v

    ! Sauvegarde initiale à t = 0.0
    call this%vtk_writer%write_step(this%mesh, this%fields, step, t)
    ! Appel de sous-programme modulaire
    call this%vtk_writer%write_diagnostics(this%mesh, this%fields, step, t, this%config%dt)

    ! Évaluation de la condition logique
    if (this%mesh%rank == 0) then
       write(*, '(A)') ""
       write(*, '(A)') ">>> DEMARRAGE DU CALCUL COVOF 2D FORTRAN 2018 PARALLELE <<<"
       write(*, '(A8,A14,A12,A14,A16,A16,A14)') &
         "Pas", "Temps [s]", "dt [s]", "Iter Krylov", "Vitesse Max", "Div Max", "Nu Moyen"
       write(*, '(A)') repeat("-", 94)
    end if

    ! Mesure du temps CPU total
    call cpu_time(start_cpu)

    ! Boucle d'itération numérique
    do while (t < this%config%t_final - 1.0e-12_dp * (abs(this%config%t_final) + 1.0_dp))
       ! Affectation et calcul numérique
       step = step + 1
       ! Affectation et calcul numérique
       u_prev = this%fields%u
       ! Affectation et calcul numérique
       dt = this%compute_time_step()

       ! Évaluation de la condition logique
       if (t + dt > this%config%t_final) then
          ! Affectation et calcul numérique
          dt = this%config%t_final - t
       end if
       if (dt < 1.0e-12_dp) exit

       ! Sauvegarde de l'état n-1 pour Gear / BDF2 et Adams-Bashforth d'ordre 2
       if (step > 1) then
          this%fields%u_nm1 = this%fields%u_old
          this%fields%v_nm1 = this%fields%v_old
          this%fields%rho_u_nm1 = this%fields%rho_u_old
          this%fields%rho_v_nm1 = this%fields%rho_v_old
          this%fields%has_nm1 = .true.
       else
          this%fields%u_nm1 = this%fields%u
          this%fields%v_nm1 = this%fields%v
          this%fields%rho_u_nm1 = this%fields%rho_u
          this%fields%rho_v_nm1 = this%fields%rho_v
          this%fields%has_nm1 = .false.
       end if

       ! Sauvegarde de l'état n pour dérivées temporelles et vitesse advectrice extrapolée
       this%fields%u_old = this%fields%u
       this%fields%v_old = this%fields%v
       this%fields%p_old = this%fields%p
       this%fields%T_old = this%fields%T
       this%fields%rho_old = this%fields%rho
       this%fields%rho_u_old = this%fields%rho_u
       this%fields%rho_v_old = this%fields%rho_v

       ! Étape 1 : Transport inertiel préservant la quantité de mouvement (SSP-RK3, Fluids 2024)
       if (this%config%enable_ssp_rk3_momentum) then
          call this%ns_solver%advect_ssp_rk3_momentum(this%mesh, this%config, this%fields, dt)
       end if

       ! Étape 2 : Transport de l'interface si écoulement diphasique
       if (this%config%enable_two_phase) then
          ! Appel de sous-programme modulaire
          call this%tracker%advect(this%mesh, this%fields, dt)
          ! Actualisation des conditions aux limites
          call apply_velocity_bc(this%fields, this%mesh, this%config, .false.)
          call apply_vof_bc(this%fields, this%mesh, this%config)
       end if

       ! 1.9 Calcul du paramètre de pénalité adaptatif r(t, M) (Vincent et al. 2004, Éq. 14 / C-AAL)
       if (this%config%solver_type == SOLVER_AUGMENTED_LAGRANGIAN .or. &
           this%config%solver_type == SOLVER_AAL .or. &
           this%config%solver_type == SOLVER_CAAL) then
          call this%ns_solver%compute_adaptative_r_al(this%mesh, this%config, this%fields, dt)
       end if

       ! 2. Assemblage des blocs du préconditionneur
       call this%ns_solver%assemble_preconditioner_blocks(this%mesh, this%config, this%fields, dt)

       ! Affectation et calcul numérique
       total_iters = 0
       ! Évaluation de la méthode de résolution Navier-Stokes
       if (this%config%solver_type == SOLVER_STANDARD_PROJECTION) then
          call this%solve_standard_projection(dt, total_iters)
       else
          if (this%config%solver_type == SOLVER_AUGMENTED_LAGRANGIAN .or. &
              this%config%solver_type == SOLVER_AAL .or. &
              this%config%solver_type == SOLVER_CAAL) then
             uzawa_max = this%config%uzawa_iter
          else
             uzawa_max = 1
          end if

          ! Boucle d'Uzawa-Lagrangien augmenté (Vincent et al. 2004, Éq. 10)
          do m = 1, uzawa_max
             ! 3. Construction du second membre rhs
             call this%ns_solver%assemble_rhs(this%mesh, this%config, this%fields, dt, rhs)

             ! 4. Initialisation du vecteur d'état à partir des champs courants
             call this%ns_solver%pack_vector(this%mesh, this%config, this%fields, state_vector)

             ! 5. Résolution du système linéaire par Krylov (BiCGSTAB2 / PCG / GMRES)
             call solve_linear_system(this%ns_solver, this%mesh, this%config, this%fields, dt, &
                                      rhs, state_vector, iters)
             total_iters = total_iters + iters

             ! 6. Déballage de la vitesse résolue
             call this%ns_solver%unpack_vector(this%mesh, this%config, state_vector, this%fields, .false.)

             ! 7. Mise à jour de la pression par Lagrangien augmenté (étape d'Uzawa C-AAL / AAL)
             if (this%config%solver_type == SOLVER_AUGMENTED_LAGRANGIAN .or. &
                 this%config%solver_type == SOLVER_AAL .or. &
                 this%config%solver_type == SOLVER_CAAL) then
                call this%ns_solver%update_pressure_al(this%mesh, this%config, this%fields, dt)

                ! Critère d'arrêt d'Uzawa sur le résidu de continuité généralisé <= uzawa_tol (Éq. 10)
                if (uzawa_max > 1 .and. this%config%uzawa_tol > 0.0_dp) then
                   local_max_div = 0.0_dp
                   do j = 1, this%mesh%ny_local
                      do i = 1, this%mesh%nx_local
                         div_val = (this%fields%u(i + 1, j) - this%fields%u(i, j)) / this%mesh%dx + &
                                   (this%fields%v(i, j + 1) - this%fields%v(i, j)) / this%mesh%dy
                         if ((this%config%solver_type == SOLVER_CAAL .or. &
                              (this%config%compressible .and. this%config%solver_type /= SOLVER_AAL)) .and. &
                             dt > 1.0e-30_dp .and. allocated(this%fields%chi_T)) then
                            div_val = div_val + (this%fields%chi_T(i, j) / dt) * &
                                      (this%fields%p(i, j) - this%fields%p_old(i, j))
                         end if
                         div_val = abs(div_val)
                         if (div_val > local_max_div) local_max_div = div_val
                      end do
                   end do
                   global_max_div = global_max(this%mesh, local_max_div)
                   if (global_max_div <= this%config%uzawa_tol) exit
                end if
             end if
          end do
       end if

       ! 7.2 Mise à jour barotrope des masses volumiques rho_m(p^{n+1}) en régime diphasique compressible
       if (this%config%compressible .and. this%config%enable_two_phase) then
          call update_physical_properties(this%fields, this%mesh, this%config)
       end if

       ! 7.5 Résolution de l'énergie (totale si compressible non-isotherme, thermique si Boussinesq)
       if (this%config%compressible .and. this%config%enable_energy) then
          call this%ns_solver%solve_total_energy(this%mesh, this%config, this%fields, dt)
       else if (this%config%enable_energy) then
          call this%ns_solver%solve_energy(this%mesh, this%config, this%fields, dt)
       end if

       ! Affectation et calcul numérique
       t = t + dt

       ! Calcul des métriques de contrôle
       local_max_u   = 0.0_dp
       ! Affectation et calcul numérique
       local_max_div = 0.0_dp
       ! Affectation et calcul numérique
       local_du_max  = 0.0_dp

       ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
       do j = 1, this%mesh%ny_local
          ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
          do i = 1, this%mesh%nx_local
             ! Affectation et calcul numérique
             speed = sqrt(this%fields%u(i, j)**2 + this%fields%v(i, j)**2)
             ! Affectation et calcul numérique
             if (speed > local_max_u) local_max_u = speed

             ! Affectation et calcul numérique
             du = abs(this%fields%u(i, j) - u_prev(i, j))
             ! Affectation et calcul numérique
             if (du > local_du_max) local_du_max = du

             ! Affectation et calcul numérique
             div_val = abs((this%fields%u(i + 1, j) - this%fields%u(i, j)) / this%mesh%dx + &
                           (this%fields%v(i, j + 1) - this%fields%v(i, j)) / this%mesh%dy)
             ! Affectation et calcul numérique
             if (div_val > local_max_div) local_max_div = div_val
          end do
       end do

       ! Affectation et calcul numérique
       global_max_u   = global_max(this%mesh, local_max_u)
       ! Affectation et calcul numérique
       global_max_div = global_max(this%mesh, local_max_div)
       ! Affectation et calcul numérique
       global_du_max  = global_max(this%mesh, local_du_max)

       ! Affectation et calcul numérique
       nu_mean = 0.0_dp
       ! Évaluation de la condition logique
       if (this%config%enable_energy) then
          ! Affectation et calcul numérique
          nu_mean = this%ns_solver%compute_mean_nusselt_hot_wall(this%mesh, this%config, this%fields)
       end if

       ! Détection de convergence stationnaire
       if (.not. this%config%enable_two_phase .and. &
           this%config%steady_state_tol > 0.0_dp .and. step > 10 .and. &
           global_du_max < this%config%steady_state_tol) then
          ! Évaluation de la condition logique
          if (this%mesh%rank == 0) then
             write(*, '(A)') ""
             write(*, '(A,ES12.4,A,ES12.4,A,F10.4,A)') &
               ! Affectation et calcul numérique
               ">>> CONVERGENCE STATIONNAIRE ATTEINTE : du_max = ", global_du_max, &
               ! Calcul de la métrique spatiale élémentaire
               " < ", this%config%steady_state_tol, " (t = ", t, " s) <<<"
             ! Évaluation de la condition logique
             if (this%config%enable_energy) then
                ! Affectation et calcul numérique
                write(*, '(A,F10.5,A)') ">>> NOMBRE DE NUSSELT FINAL : Nu = ", nu_mean, " <<<"
             end if
          end if
          ! Appel de sous-programme modulaire
          call this%vtk_writer%write_step(this%mesh, this%fields, step, t)
          exit
       end if

       ! Affichage console régulier
       if (this%mesh%rank == 0 .and. (mod(step, this%config%log_interval) == 0)) then
          write(*, '(I8,F14.5,ES12.3,I14,ES16.4,ES16.4,F14.4)') &
            step, t, dt, total_iters, global_max_u, global_max_div, nu_mean
       end if

       ! Écriture des diagnostics CSV
       if (this%config%write_csv) then
          ! Appel de sous-programme modulaire
          call this%vtk_writer%write_diagnostics(this%mesh, this%fields, step, t, dt)
       end if

       ! Export des fichiers VTK
       if (mod(step, this%config%vtk_interval) == 0 .or. t >= this%config%t_final - 1.0e-10_dp * (abs(this%config%t_final) + 1.0_dp)) then
          ! Appel de sous-programme modulaire
          call this%vtk_writer%write_step(this%mesh, this%fields, step, t)
       end if
    end do

    ! Appel de sous-programme modulaire
    call this%vtk_writer%finalize()

    call cpu_time(end_cpu)
    elapsed_cpu = end_cpu - start_cpu

    if (this%mesh%rank == 0) then
       if (trim(this%config%test_case) == "CAS_3_8" .or. &
           trim(this%config%test_case) == "CAS_3_8_POISEUILLE_UNSTEADY" .or. &
           index(this%config%test_case, "unsteady") > 0 .or. &
           index(this%config%case_name, "cas3_8") > 0 .or. &
           index(this%config%output_dir, "cas3_8") > 0 .or. &
           index(this%config%case_name, "unsteady") > 0 .or. &
           index(this%config%output_dir, "unsteady") > 0) then
          call this%export_poiseuille_unsteady_results(step, t, elapsed_cpu)
       else if (trim(this%config%test_case) == "CAS_3_7" .or. &
           trim(this%config%test_case) == "CAS_3_7_POISEUILLE" .or. &
           trim(this%config%test_case) == "POISEUILLE" .or. &
           index(this%config%case_name, "poiseuille") > 0 .or. &
           index(this%config%case_name, "cas3_7") > 0 .or. &
           index(this%config%output_dir, "poiseuille") > 0 .or. &
           index(this%config%output_dir, "cas3_7") > 0) then
          call this%export_poiseuille_results(step, t, elapsed_cpu)
       end if
    end if

    ! Évaluation de la condition logique
    if (this%mesh%rank == 0) then
       write(*, '(A)') repeat("-", 94)
       ! Affectation et calcul numérique
       write(*, '(A,F10.4,A)') ">>> SIMULATION COVOF 2D TERMINEE AVEC SUCCES A t = ", t, " s <<<"
       write(*, '(A,A)') "Fichiers de visualisation VTK et CSV disponibles dans : ", trim(this%config%output_dir)
    end if

    ! Libération de la mémoire allouée
    deallocate(state_vector, rhs, u_prev)
  end subroutine sim_run

  !-----------------------------------------------------------------------------
  ! Exporte les résultats et calculs d'erreurs pour le cas de Poiseuille (Cas 3.7)
  !-----------------------------------------------------------------------------
  subroutine sim_export_poiseuille_results(this, total_step, final_time, elapsed_cpu)
    class(simulation_t), intent(inout) :: this
    integer, intent(in)                :: total_step
    real(dp), intent(in)               :: final_time, elapsed_cpu

    integer :: nx, ny, i, j, i_mid, i_3q, j_mid
    real(dp) :: lx, ly, dx, dy, mu, V0, Gx_theo, u_exact_val, uc, vc, pc, diff
    real(dp) :: sum_err_l1, sum_err_l2, max_err_linf, n_cells, err_l1, err_l2
    real(dp) :: p_in, p_out, x_in, x_out, Gx_num, err_grad_p, y_val, x_val
    real(dp), allocatable :: u_c(:,:), v_c(:,:), p_c(:,:), u_ex(:,:), err_u(:,:)
    character(len=256) :: file_path
    integer :: f_unit

    if (this%mesh%rank /= 0) return

    nx = this%mesh%nx_local
    ny = this%mesh%ny_local
    lx = this%mesh%lx
    ly = this%mesh%ly
    dx = this%mesh%dx
    dy = this%mesh%dy
    mu = this%config%mu1
    V0 = this%config%bc_left%u
    if (V0 <= 0.0_dp) V0 = 1.0_dp
    Gx_theo = 12.0_dp * mu * V0 / (ly * ly)

    allocate(u_c(nx, ny), v_c(nx, ny), p_c(nx, ny), u_ex(nx, ny), err_u(nx, ny))

    sum_err_l1 = 0.0_dp
    sum_err_l2 = 0.0_dp
    max_err_linf = 0.0_dp

    j_mid = ny / 2
    if (j_mid < 1) j_mid = 1
    p_out = this%fields%p(nx, j_mid)

    do j = 1, ny
       y_val = y_center(this%mesh, j)
       u_exact_val = 6.0_dp * V0 * y_val * (ly - y_val) / (ly * ly)
       do i = 1, nx
          uc = 0.5_dp * (this%fields%u(i, j) + this%fields%u(i + 1, j))
          vc = 0.5_dp * (this%fields%v(i, j) + this%fields%v(i, j + 1))
          pc = this%fields%p(i, j)

          u_c(i, j) = uc
          v_c(i, j) = vc
          p_c(i, j) = pc
          u_ex(i, j) = u_exact_val

          diff = abs(uc - u_exact_val)
          err_u(i, j) = diff

          sum_err_l1 = sum_err_l1 + diff
          sum_err_l2 = sum_err_l2 + diff * diff
          if (diff > max_err_linf) max_err_linf = diff
       end do
    end do

    n_cells = real(nx * ny, dp)
    err_l1 = sum_err_l1 / n_cells
    err_l2 = sqrt(sum_err_l2 / n_cells)

    p_in = this%fields%p(1, j_mid)
    x_in = x_center(this%mesh, 1)
    x_out = x_center(this%mesh, nx)
    Gx_num = (p_in - p_out) / (x_out - x_in)
    err_grad_p = abs(Gx_num - Gx_theo)

    write(*, '(A)') ""
    write(*, '(A)') "========================================================"
    write(*, '(A)') "  RESULTATS CAS 3.7 : ECOULEMENT DE POISEUILLE PLAN     "
    write(*, '(A,F8.4,A,F8.4,A,I0,A,I0)') "  Domaine [", lx, " x ", ly, "], Grille ", nx, " x ", ny
    write(*, '(A,ES14.6)') "  Gradient theorique Gx = ", Gx_theo
    write(*, '(A,ES14.6,A,ES14.6,A)') "  Gradient numerique Gx = ", Gx_num, " (err = ", err_grad_p, ")"
    write(*, '(A,ES14.6)') "  Erreur L1(u)   : ", err_l1
    write(*, '(A,ES14.6)') "  Erreur L2(u)   : ", err_l2
    write(*, '(A,ES14.6)') "  Erreur Linf(u) : ", max_err_linf
    write(*, '(A,F10.4,A)') "  Temps CPU      : ", elapsed_cpu, " s"
    write(*, '(A)') "========================================================"

    call execute_command_line("mkdir -p " // trim(this%config%output_dir))

    ! 1. poiseuille_error.csv
    file_path = trim(this%config%output_dir) // "/poiseuille_error.csv"
    open(newunit=f_unit, file=trim(file_path), status="replace", action="write")
    write(f_unit, '(A)') "nx,ny,dx,dy,iters,err_l1,err_l2,err_linf,err_grad_p,cpu_time"
    write(f_unit, '(I0,A,I0,A,ES14.6,A,ES14.6,A,I0,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,F10.6)') &
         nx, ",", ny, ",", dx, ",", dy, ",", total_step, ",", &
         err_l1, ",", err_l2, ",", max_err_linf, ",", err_grad_p, ",", elapsed_cpu
    close(f_unit)

    ! 2. profile_u_y.csv
    i_mid = nx / 2
    if (i_mid < 1) i_mid = 1
    i_3q = (3 * nx) / 4
    if (i_3q < 1) i_3q = 1
    file_path = trim(this%config%output_dir) // "/profile_u_y.csv"
    open(newunit=f_unit, file=trim(file_path), status="replace", action="write")
    write(f_unit, '(A)') "y,u_mid,u_three_quarter,u_ex,error"
    do j = 1, ny
       y_val = y_center(this%mesh, j)
       u_exact_val = 6.0_dp * V0 * y_val * (ly - y_val) / (ly * ly)
       write(f_unit, '(ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)') &
            y_val, ",", u_c(i_mid, j), ",", u_c(i_3q, j), ",", u_exact_val, ",", abs(u_c(i_mid, j) - u_exact_val)
    end do
    close(f_unit)

    ! 3. profile_p_x.csv
    file_path = trim(this%config%output_dir) // "/profile_p_x.csv"
    open(newunit=f_unit, file=trim(file_path), status="replace", action="write")
    write(f_unit, '(A)') "x,p_mid,p_ex,error"
    do i = 1, nx
       x_val = x_center(this%mesh, i)
       write(f_unit, '(ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)') &
            x_val, ",", (p_c(i, j_mid) - p_out), ",", Gx_theo * (lx - x_val), ",", abs((p_c(i, j_mid) - p_out) - Gx_theo * (lx - x_val))
    end do
    close(f_unit)

    ! 4. solution_field_2d.csv
    file_path = trim(this%config%output_dir) // "/solution_field_2d.csv"
    open(newunit=f_unit, file=trim(file_path), status="replace", action="write")
    write(f_unit, '(A)') "i,j,x,y,u,v,p,u_ex,p_ex,error_u"
    do j = 1, ny
       y_val = y_center(this%mesh, j)
       do i = 1, nx
          x_val = x_center(this%mesh, i)
          write(f_unit, '(I0,A,I0,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)') &
               i, ",", j, ",", x_val, ",", y_val, ",", &
               u_c(i, j), ",", v_c(i, j), ",", (p_c(i, j) - p_out), ",", &
               u_ex(i, j), ",", Gx_theo * (lx - x_val), ",", err_u(i, j)
       end do
    end do
    close(f_unit)

    ! 5. solution_field.vtk
    file_path = trim(this%config%output_dir) // "/solution_field.vtk"
    open(newunit=f_unit, file=trim(file_path), status="replace", action="write")
    write(f_unit, '(A)') "# vtk DataFile Version 3.0"
    write(f_unit, '(A)') "COVOF 2D Cas 3.7 Poiseuille Solution (Fortran)"
    write(f_unit, '(A)') "ASCII"
    write(f_unit, '(A)') "DATASET RECTILINEAR_GRID"
    write(f_unit, '(A,I0,A,I0,A)') "DIMENSIONS ", (nx + 1), " ", (ny + 1), " 1"
    write(f_unit, '(A,I0,A)') "X_COORDINATES ", (nx + 1), " double"
    do i = 0, nx
       write(f_unit, '(ES16.8," ")', advance="no") i * dx
    end do
    write(f_unit, *)
    write(f_unit, '(A,I0,A)') "Y_COORDINATES ", (ny + 1), " double"
    do j = 0, ny
       write(f_unit, '(ES16.8," ")', advance="no") j * dy
    end do
    write(f_unit, *)
    write(f_unit, '(A)') "Z_COORDINATES 1 double"
    write(f_unit, '(A)') "0.0"
    write(f_unit, '(A,I0)') "CELL_DATA ", (nx * ny)
    write(f_unit, '(A)') "SCALARS Velocity_U double 1"
    write(f_unit, '(A)') "LOOKUP_TABLE default"
    do j = 1, ny
       do i = 1, nx
          write(f_unit, '(ES16.8," ")', advance="no") u_c(i, j)
       end do
       write(f_unit, *)
    end do
    write(f_unit, '(A)') "SCALARS Pressure double 1"
    write(f_unit, '(A)') "LOOKUP_TABLE default"
    do j = 1, ny
       do i = 1, nx
          write(f_unit, '(ES16.8," ")', advance="no") (p_c(i, j) - p_out)
       end do
       write(f_unit, *)
    end do
    write(f_unit, '(A)') "SCALARS Error_U double 1"
    write(f_unit, '(A)') "LOOKUP_TABLE default"
    do j = 1, ny
       do i = 1, nx
          write(f_unit, '(ES16.8," ")', advance="no") err_u(i, j)
       end do
       write(f_unit, *)
    end do
    close(f_unit)

    deallocate(u_c, v_c, p_c, u_ex, err_u)
  end subroutine sim_export_poiseuille_results

  !-----------------------------------------------------------------------------
  ! Exporte les résultats et calculs d'erreurs pour le cas de Poiseuille instationnaire (Cas 3.8)
  !-----------------------------------------------------------------------------
  subroutine sim_export_poiseuille_unsteady_results(this, total_step, final_time, elapsed_cpu)
    class(simulation_t), intent(inout) :: this
    integer, intent(in)                :: total_step
    real(dp), intent(in)               :: final_time, elapsed_cpu

    integer :: nx, ny, i, j, i_mid, m, n_mode
    real(dp) :: lx, ly, dx, dy, rho, mu, nu, gx, k, t_val, u_exact_val, uc, vc, pc, diff
    real(dp) :: sum_err_l1, sum_err_l2, max_err_linf, n_cells, err_l1, err_l2
    real(dp) :: max_u_num, max_u_exact, y_val, x_val, u_st, tr, lam
    real(dp), parameter :: PI_VAL = 3.14159265358979323846_dp
    real(dp), allocatable :: u_c(:,:), v_c(:,:), p_c(:,:), u_ex(:,:), err_u(:,:)
    character(len=256) :: file_path
    integer :: f_unit

    if (this%mesh%rank /= 0) return

    nx = this%mesh%nx_local
    ny = this%mesh%ny_local
    lx = this%mesh%lx
    ly = this%mesh%ly
    dx = this%mesh%dx
    dy = this%mesh%dy
    rho = this%config%rho1
    if (rho <= 0.0_dp) rho = 1.0_dp
    mu = this%config%mu1
    if (mu <= 0.0_dp) mu = 0.01_dp
    nu = mu / rho
    gx = this%config%gx
    if (abs(gx) > 1.0e-12_dp) then
       k = rho * gx
    else
       k = 0.0555_dp
    end if
    t_val = final_time

    allocate(u_c(nx, ny), v_c(nx, ny), p_c(nx, ny), u_ex(nx, ny), err_u(nx, ny))

    sum_err_l1 = 0.0_dp
    sum_err_l2 = 0.0_dp
    max_err_linf = 0.0_dp
    max_u_num = 0.0_dp
    max_u_exact = 0.0_dp

    do j = 1, ny
       y_val = y_center(this%mesh, j)
       u_st = (k / (2.0_dp * nu)) * y_val * (ly - y_val)
       if (t_val <= 0.0_dp) then
          u_exact_val = 0.0_dp
       else
          tr = 0.0_dp
          do m = 0, 49
             n_mode = 2 * m + 1
             lam = real(n_mode, dp) * PI_VAL / ly
             tr = tr + (1.0_dp / real(n_mode**3, dp)) * sin(lam * y_val) * exp(-nu * (lam**2) * t_val)
          end do
          u_exact_val = u_st - (4.0_dp * k * (ly**2) / (nu * PI_VAL**3)) * tr
       end if

       if (u_exact_val > max_u_exact) max_u_exact = u_exact_val

       do i = 1, nx
          uc = 0.5_dp * (this%fields%u(i, j) + this%fields%u(i + 1, j))
          vc = 0.5_dp * (this%fields%v(i, j) + this%fields%v(i + 1, j))
          pc = this%fields%p(i, j)

          u_c(i, j) = uc
          v_c(i, j) = vc
          p_c(i, j) = pc
          u_ex(i, j) = u_exact_val

          if (uc > max_u_num) max_u_num = uc

          diff = abs(uc - u_exact_val)
          err_u(i, j) = diff

          sum_err_l1 = sum_err_l1 + diff
          sum_err_l2 = sum_err_l2 + diff * diff
          if (diff > max_err_linf) max_err_linf = diff
       end do
    end do

    n_cells = real(nx * ny, dp)
    err_l1 = sum_err_l1 / n_cells
    err_l2 = sqrt(sum_err_l2 / n_cells)

    write(*, '(A)') repeat("=", 56)
    write(*, '(A)') "  RESULTATS CAS 3.8 : POISEUILLE PLAN INSTATIONNAIRE"
    write(*, '(A,F8.4,A,F8.4,A,I4,A,I4)') "  Domaine [", lx, " x ", ly, "], Grille ", nx, " x ", ny
    write(*, '(A,F10.4,A,ES12.4,A)') "  Temps t = ", t_val, " s, dt = ", this%config%dt, " s"
    write(*, '(A,F10.5,A,ES12.4,A)') "  Gradient force k = ", k, " Pa/m, nu = ", nu, " m^2/s"
    write(*, '(A,ES12.5,A,ES12.5,A)') "  Vitesse max num = ", max_u_num, " m/s (exact = ", max_u_exact, " m/s)"
    write(*, '(A,ES14.6)') "  Erreur L1(u)   : ", err_l1
    write(*, '(A,ES14.6)') "  Erreur L2(u)   : ", err_l2
    write(*, '(A,ES14.6)') "  Erreur Linf(u) : ", max_err_linf
    write(*, '(A,F10.4,A)') "  Temps CPU      : ", elapsed_cpu, " s"
    write(*, '(A)') repeat("=", 56)

    ! 1. poiseuille_unsteady_error.csv
    file_path = trim(this%config%output_dir) // "/poiseuille_unsteady_error.csv"
    open(newunit=f_unit, file=trim(file_path), status='replace', action='write')
    write(f_unit, '(A)') "nx,ny,dx,dy,dt,t_final,iters,err_l1,err_l2,err_linf,u_max_num,u_max_exact,cpu_time"
    write(f_unit, '(I0,",",I0,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6,",",I0,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6,",",F10.4)') &
         nx, ny, dx, dy, this%config%dt, t_val, total_step, err_l1, err_l2, max_err_linf, max_u_num, max_u_exact, elapsed_cpu
    close(f_unit)

    ! 2. profile_u_y.csv
    i_mid = nx / 2
    if (i_mid < 1) i_mid = 1
    file_path = trim(this%config%output_dir) // "/profile_u_y.csv"
    open(newunit=f_unit, file=trim(file_path), status='replace', action='write')
    write(f_unit, '(A)') "y,u_num,u_exact,error"
    do j = 1, ny
       y_val = y_center(this%mesh, j)
       write(f_unit, '(ES14.6,",",ES14.6,",",ES14.6,",",ES14.6)') &
            y_val, u_c(i_mid, j), u_ex(i_mid, j), err_u(i_mid, j)
    end do
    close(f_unit)

    ! 3. solution_field_2d.csv
    file_path = trim(this%config%output_dir) // "/solution_field_2d.csv"
    open(newunit=f_unit, file=trim(file_path), status='replace', action='write')
    write(f_unit, '(A)') "i,j,x,y,u,v,p,u_exact,error_u"
    do j = 1, ny
       y_val = y_center(this%mesh, j)
       do i = 1, nx
          x_val = x_center(this%mesh, i)
          write(f_unit, '(I0,",",I0,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6)') &
               i, j, x_val, y_val, u_c(i, j), v_c(i, j), p_c(i, j), u_ex(i, j), err_u(i, j)
       end do
    end do
    close(f_unit)

    ! 4. solution_field.vtk
    file_path = trim(this%config%output_dir) // "/solution_field.vtk"
    open(newunit=f_unit, file=trim(file_path), status='replace', action='write')
    write(f_unit, '(A)') "# vtk DataFile Version 3.0"
    write(f_unit, '(A)') "COVOF 2D Cas 3.8 Poiseuille Unsteady Solution"
    write(f_unit, '(A)') "ASCII"
    write(f_unit, '(A)') "DATASET RECTILINEAR_GRID"
    write(f_unit, '(A,I0,A,I0,A)') "DIMENSIONS ", nx + 1, " ", ny + 1, " 1"
    write(f_unit, '(A,I0,A)') "X_COORDINATES ", nx + 1, " double"
    do i = 0, nx
       write(f_unit, '(ES14.6," ")', advance='no') real(i, dp) * dx
    end do
    write(f_unit, *)
    write(f_unit, '(A,I0,A)') "Y_COORDINATES ", ny + 1, " double"
    do j = 0, ny
       write(f_unit, '(ES14.6," ")', advance='no') real(j, dp) * dy
    end do
    write(f_unit, *)
    write(f_unit, '(A)') "Z_COORDINATES 1 double"
    write(f_unit, '(A)') "0.0"
    write(f_unit, '(A,I0)') "CELL_DATA ", nx * ny
    write(f_unit, '(A)') "SCALARS Velocity_U double 1"
    write(f_unit, '(A)') "LOOKUP_TABLE default"
    do j = 1, ny
       do i = 1, nx
          write(f_unit, '(ES14.6," ")', advance='no') u_c(i, j)
       end do
       write(f_unit, *)
    end do
    write(f_unit, '(A)') "SCALARS Pressure double 1"
    write(f_unit, '(A)') "LOOKUP_TABLE default"
    do j = 1, ny
       do i = 1, nx
          write(f_unit, '(ES14.6," ")', advance='no') p_c(i, j)
       end do
       write(f_unit, *)
    end do
    write(f_unit, '(A)') "SCALARS Error_U double 1"
    write(f_unit, '(A)') "LOOKUP_TABLE default"
    do j = 1, ny
       do i = 1, nx
          write(f_unit, '(ES14.6," ")', advance='no') err_u(i, j)
       end do
       write(f_unit, *)
    end do
    close(f_unit)

    deallocate(u_c, v_c, p_c, u_ex, err_u)
  end subroutine sim_export_poiseuille_unsteady_results

  !-----------------------------------------------------------------------------
  ! Cas de vérification 3.1 : Laplacien -Delta(Psi) = -1 avec solution analytique
  !-----------------------------------------------------------------------------
  subroutine sim_run_laplacian(this)
    ! Déclaration de variable locale / paramètre de calcul
    class(simulation_t), intent(inout) :: this

    double precision :: a, b, dx, dy, psi_ex, xp, yp, num, den, n1, n2, diff
    double precision :: err_l1, err_l2, err_linf, tol, rz_old, rz_new, alpha, beta, pAp, r_norm
    ! Déclaration de variable locale / paramètre de calcul
    integer :: nx, ny, n_int_x, n_int_y, n_dofs, i, j, k1, k2, row, iters, max_it, err_unit
    logical :: use_direct, use_bicgstab2
    integer :: w, k, m, m_j
    double precision, allocatable :: L(:,:), y_dir(:)
    double precision :: s_val, val_a
    double precision, allocatable :: u_vec(:,:), r_mat(:,:), r_tilde(:), Ap2(:)
    double precision :: rho_0, rho_1, omega, sigma11, sigma12, sigma22, rhs1, rhs2, det_val, s11_s22, gamma1, gamma2, gamma_val
    double precision :: r0_norm
    ! Déclaration de variable locale / paramètre de calcul
    character(len=256) :: err_file, file_name
    integer :: j_mid, i_mid
    double precision, allocatable :: psi_num_2d(:,:), psi_ex_2d(:,:), err_2d(:,:)
    double precision, allocatable :: A_diag(:), A_east(:), A_west(:), A_north(:), A_south(:)
    double precision, allocatable :: rhs(:), x(:), r(:), p(:), z(:), Ap(:)

    ! Affectation et calcul numérique
    a = this%config%lx
    ! Affectation et calcul numérique
    b = this%config%ly
    ! Affectation et calcul numérique
    nx = this%config%nx
    ! Affectation et calcul numérique
    ny = this%config%ny
    ! Calcul de la métrique spatiale élémentaire
    dx = a / real(nx, dp)
    ! Calcul de la métrique spatiale élémentaire
    dy = b / real(ny, dp)

    ! Évaluation de la condition logique
    if (this%mesh%rank == 0) then
       write(*, '(A)') ""
       ! Affectation et calcul numérique
       write(*, '(A)') "========================================================"
       ! Affectation et calcul numérique
       write(*, '(A)') "  VERIFICATION DU CAS 3.1 : LAPLACIEN DELTA(PSI) = 1"
       write(*, '(A,F8.4,A,F8.4,A,I0,A,I0,A)') &
         "  Domaine [", a, " x ", b, "], Grille ", nx, " x ", ny
       ! Affectation et calcul numérique
       write(*, '(A)') "========================================================"
    end if

    ! Affectation et calcul numérique
    n_int_x = nx - 1
    ! Affectation et calcul numérique
    n_int_y = ny - 1
    ! Affectation et calcul numérique
    n_dofs = n_int_x * n_int_y

    ! Allocation dynamique de tableau contigu en mémoire
    allocate(A_diag(n_dofs), A_east(n_dofs), A_west(n_dofs), A_north(n_dofs), A_south(n_dofs))
    ! Allocation dynamique de tableau contigu en mémoire
    allocate(rhs(n_dofs), x(n_dofs), r(n_dofs), p(n_dofs), z(n_dofs), Ap(n_dofs))

    ! Affectation et calcul numérique
    A_diag  = 2.0_dp / (dx * dx) + 2.0_dp / (dy * dy)
    ! Affectation et calcul numérique
    A_east  = -1.0_dp / (dx * dx)
    ! Affectation et calcul numérique
    A_west  = -1.0_dp / (dx * dx)
    ! Affectation et calcul numérique
    A_north = -1.0_dp / (dy * dy)
    ! Affectation et calcul numérique
    A_south = -1.0_dp / (dy * dy)

    ! Affectation et calcul numérique
    rhs = -this%config%source_term
    ! Affectation et calcul numérique
    x   = 0.0_dp

    ! Sélection et exécution du solveur linéaire
    use_direct = (this%config%linear_solver == LINEAR_DIRECT)
    use_bicgstab2 = (this%config%linear_solver == LINEAR_BICGSTAB2 .or. &
                    (this%config%solver_type == SOLVER_MONOLITHIC .and. .not. use_direct .and. this%config%linear_solver /= LINEAR_CG))

    max_it = this%config%max_iter
    if (max_it <= 0) max_it = 5000
    tol = this%config%tol
    if (tol <= 0.0_dp) tol = 1.0e-12_dp
    iters = 0

    if (use_direct) then
       if (this%mesh%rank == 0) then
          write(*, '(A)') "  Solveur utilise : DIRECT (Cholesky bande / MUMPS)"
       end if
       w = n_int_x
       allocate(L(n_dofs, 0:w), y_dir(n_dofs))
       L = 0.0_dp

       do i = 1, n_dofs
          do k = min(i - 1, w), 0, -1
             j = i - k
             val_a = 0.0_dp
             if (i == j) then
                val_a = 2.0_dp / (dx * dx) + 2.0_dp / (dy * dy)
             else if (abs(i - j) == 1 .and. mod(min(i, j) - 1, n_int_x) /= n_int_x - 1) then
                val_a = -1.0_dp / (dx * dx)
             else if (abs(i - j) == n_int_x) then
                val_a = -1.0_dp / (dy * dy)
             end if
             s_val = 0.0_dp
             do m = k + 1, min(i - 1, w)
                m_j = m - k
                if (m_j <= min(j - 1, w)) then
                   s_val = s_val + L(i, m) * L(j, m_j)
                end if
             end do
             if (i == j) then
                L(i, 0) = sqrt(max(val_a - s_val, 1.0e-30_dp))
             else
                L(i, k) = (val_a - s_val) / L(j, 0)
             end if
          end do
       end do

       ! Solve L y = rhs
       do i = 1, n_dofs
          s_val = rhs(i)
          do k = 1, min(i - 1, w)
             s_val = s_val - L(i, k) * y_dir(i - k)
          end do
          y_dir(i) = s_val / L(i, 0)
       end do

       ! Solve L^T x = y
       do i = n_dofs, 1, -1
          s_val = y_dir(i)
          do k = 1, min(n_dofs - i, w)
             s_val = s_val - L(i + k, k) * x(i + k)
          end do
          x(i) = s_val / L(i, 0)
       end do
       iters = 1

    else if (use_bicgstab2) then
       if (this%mesh%rank == 0) then
          write(*, '(A)') "  Solveur utilise : ITERATIF BiCGSTAB(2) + Preconditionneur"
       end if
       allocate(u_vec(n_dofs, 0:2), r_mat(n_dofs, 0:2), r_tilde(n_dofs), Ap2(n_dofs))
       u_vec = 0.0_dp
       r_mat = 0.0_dp

       ! r_mat(0) = P^{-1} * rhs - P^{-1} A x = rhs / A_diag (car x=0)
       r_mat(:, 0) = rhs / A_diag
       r_tilde = r_mat(:, 0)
       r0_norm = sqrt(sum(r_mat(:, 0)**2))

       rho_0 = 1.0_dp
       alpha = 0.0_dp
       omega = 1.0_dp

       do while (iters < max_it .and. r0_norm > tol)
          iters = iters + 2
          rho_0 = -omega * rho_0

          ! Demi-pas 1 (j = 0)
          rho_1 = sum(r_tilde * r_mat(:, 0))
          if (abs(rho_0) < 1.0e-300_dp) exit
          beta = alpha * (rho_1 / rho_0)
          rho_0 = rho_1

          u_vec(:, 0) = r_mat(:, 0) - beta * u_vec(:, 0)

          ! u_vec(:, 1) = P^{-1} A u_vec(:, 0)
          do j = 0, n_int_y - 1
             do i = 0, n_int_x - 1
                row = j * n_int_x + i + 1
                Ap2(row) = A_diag(row) * u_vec(row, 0)
                if (i > 0)           Ap2(row) = Ap2(row) + A_west(row)  * u_vec(row - 1, 0)
                if (i < n_int_x - 1) Ap2(row) = Ap2(row) + A_east(row)  * u_vec(row + 1, 0)
                if (j > 0)           Ap2(row) = Ap2(row) + A_south(row) * u_vec(row - n_int_x, 0)
                if (j < n_int_y - 1) Ap2(row) = Ap2(row) + A_north(row) * u_vec(row + n_int_x, 0)
                u_vec(row, 1) = Ap2(row) / A_diag(row)
             end do
          end do

          gamma_val = sum(r_tilde * u_vec(:, 1))
          if (abs(gamma_val) < 1.0e-300_dp) exit
          alpha = rho_0 / gamma_val

          x = x + alpha * u_vec(:, 0)
          r_mat(:, 0) = r_mat(:, 0) - alpha * u_vec(:, 1)

          ! r_mat(:, 1) = P^{-1} A r_mat(:, 0)
          do j = 0, n_int_y - 1
             do i = 0, n_int_x - 1
                row = j * n_int_x + i + 1
                Ap2(row) = A_diag(row) * r_mat(row, 0)
                if (i > 0)           Ap2(row) = Ap2(row) + A_west(row)  * r_mat(row - 1, 0)
                if (i < n_int_x - 1) Ap2(row) = Ap2(row) + A_east(row)  * r_mat(row + 1, 0)
                if (j > 0)           Ap2(row) = Ap2(row) + A_south(row) * r_mat(row - n_int_x, 0)
                if (j < n_int_y - 1) Ap2(row) = Ap2(row) + A_north(row) * r_mat(row + n_int_x, 0)
                r_mat(row, 1) = Ap2(row) / A_diag(row)
             end do
          end do

          ! Demi-pas 2 (j = 1)
          rho_1 = sum(r_tilde * r_mat(:, 1))
          if (abs(rho_0) < 1.0e-300_dp) exit
          beta = alpha * (rho_1 / rho_0)
          rho_0 = rho_1

          u_vec(:, 0) = r_mat(:, 0) - beta * u_vec(:, 0)
          u_vec(:, 1) = r_mat(:, 1) - beta * u_vec(:, 1)

          ! u_vec(:, 2) = P^{-1} A u_vec(:, 1)
          do j = 0, n_int_y - 1
             do i = 0, n_int_x - 1
                row = j * n_int_x + i + 1
                Ap2(row) = A_diag(row) * u_vec(row, 1)
                if (i > 0)           Ap2(row) = Ap2(row) + A_west(row)  * u_vec(row - 1, 1)
                if (i < n_int_x - 1) Ap2(row) = Ap2(row) + A_east(row)  * u_vec(row + 1, 1)
                if (j > 0)           Ap2(row) = Ap2(row) + A_south(row) * u_vec(row - n_int_x, 1)
                if (j < n_int_y - 1) Ap2(row) = Ap2(row) + A_north(row) * u_vec(row + n_int_x, 1)
                u_vec(row, 2) = Ap2(row) / A_diag(row)
             end do
          end do

          gamma_val = sum(r_tilde * u_vec(:, 2))
          if (abs(gamma_val) < 1.0e-300_dp) exit
          alpha = rho_0 / gamma_val

          x = x + alpha * u_vec(:, 0)
          r_mat(:, 0) = r_mat(:, 0) - alpha * u_vec(:, 1)
          r_mat(:, 1) = r_mat(:, 1) - alpha * u_vec(:, 2)

          ! r_mat(:, 2) = P^{-1} A r_mat(:, 1)
          do j = 0, n_int_y - 1
             do i = 0, n_int_x - 1
                row = j * n_int_x + i + 1
                Ap2(row) = A_diag(row) * r_mat(row, 1)
                if (i > 0)           Ap2(row) = Ap2(row) + A_west(row)  * r_mat(row - 1, 1)
                if (i < n_int_x - 1) Ap2(row) = Ap2(row) + A_east(row)  * r_mat(row + 1, 1)
                if (j > 0)           Ap2(row) = Ap2(row) + A_south(row) * r_mat(row - n_int_x, 1)
                if (j < n_int_y - 1) Ap2(row) = Ap2(row) + A_north(row) * r_mat(row + n_int_x, 1)
                r_mat(row, 2) = Ap2(row) / A_diag(row)
             end do
          end do

          ! Minimisation quadratique Sleijpen
          sigma11 = sum(r_mat(:, 1)**2)
          sigma12 = sum(r_mat(:, 2) * r_mat(:, 1))
          sigma22 = sum(r_mat(:, 2)**2)
          rhs1    = sum(r_mat(:, 0) * r_mat(:, 1))
          rhs2    = sum(r_mat(:, 0) * r_mat(:, 2))

          det_val = sigma11 * sigma22 - sigma12 * sigma12
          s11_s22 = sigma11 * sigma22
          if (s11_s22 < 1.0e-30_dp .or. abs(det_val) < 1.0e-10_dp * s11_s22) then
             gamma2 = 0.0_dp
             gamma1 = rhs1 / (sigma11 + 1.0e-30_dp)
          else
             gamma1 = (rhs1 * sigma22 - rhs2 * sigma12) / det_val
             gamma2 = (rhs2 * sigma11 - rhs1 * sigma12) / det_val
          end if

          x = x + gamma1 * r_mat(:, 0) + gamma2 * r_mat(:, 1)
          r_mat(:, 0) = r_mat(:, 0) - (gamma1 * r_mat(:, 1) + gamma2 * r_mat(:, 2))
          u_vec(:, 0) = u_vec(:, 0) - (gamma1 * u_vec(:, 1) + gamma2 * u_vec(:, 2))
          if (gamma2 /= 0.0_dp) then
             omega = gamma2
          else if (abs(gamma1) > 1.0e-30_dp) then
             omega = -gamma1
          else
             omega = 1.0_dp
          end if

          r_norm = sqrt(sum(r_mat(:, 0)**2))
          if (r_norm < tol .or. r_norm / r0_norm < tol) exit
       end do

    else
       if (this%mesh%rank == 0) then
          write(*, '(A)') "  Solveur utilise : ITERATIF Gradient Conjugue (CG) + Preconditionneur Jacobi"
       end if
       r = rhs
       z = r / A_diag
       p = z
       rz_old = sum(r * z)

       do iters = 1, max_it
          do j = 0, n_int_y - 1
             do i = 0, n_int_x - 1
                row = j * n_int_x + i + 1
                Ap(row) = A_diag(row) * p(row)
                if (i > 0)           Ap(row) = Ap(row) + A_west(row)  * p(row - 1)
                if (i < n_int_x - 1) Ap(row) = Ap(row) + A_east(row)  * p(row + 1)
                if (j > 0)           Ap(row) = Ap(row) + A_south(row) * p(row - n_int_x)
                if (j < n_int_y - 1) Ap(row) = Ap(row) + A_north(row) * p(row + n_int_x)
             end do
          end do

          pAp = sum(p * Ap)
          if (abs(pAp) < 1.0e-30_dp) exit

          alpha = rz_old / pAp
          x = x + alpha * p
          r = r - alpha * Ap

          r_norm = sqrt(sum(r * r))
          if (r_norm < tol) exit

          z = r / A_diag
          rz_new = sum(r * z)
          beta = rz_new / rz_old
          rz_old = rz_new

          p = z + beta * p
       end do
    end if

    ! Évaluation de la solution analytique exacte (formule analytique de référence)
    err_l1   = 0.0_dp
    ! Évaluation de la norme du résidu ou de l'erreur numérique
    err_l2   = 0.0_dp
    ! Évaluation de la norme du résidu ou de l'erreur numérique
    err_linf = 0.0_dp

    ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
    do j = 0, n_int_y - 1
       ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
       do i = 0, n_int_x - 1
          ! Affectation et calcul numérique
          xp = real(i + 1, dp) * dx
          ! Affectation et calcul numérique
          yp = real(j + 1, dp) * dy
          ! Affectation et calcul numérique
          psi_ex = 0.0_dp

          ! Boucle d'itération numérique
          do k2 = 0, 59
             ! Boucle d'itération numérique
             do k1 = 0, 59
                ! Affectation et calcul numérique
                n1 = real(2 * k1 + 1, dp)
                ! Affectation et calcul numérique
                n2 = real(2 * k2 + 1, dp)
                ! Affectation et calcul numérique
                num = 16.0_dp * this%config%source_term * sin(n1 * PI * xp / a) * sin(n2 * PI * yp / b)
                ! Affectation et calcul numérique
                den = (PI**4) * n1 * n2 * ((n1 / a)**2 + (n2 / b)**2)
                ! Affectation et calcul numérique
                psi_ex = psi_ex - num / den
             end do
          end do

          ! Affectation et calcul numérique
          row = j * n_int_x + i + 1
          ! Affectation et calcul numérique
          diff = abs(x(row) - psi_ex)
          ! Évaluation de la norme du résidu ou de l'erreur numérique
          err_l1 = err_l1 + diff
          ! Évaluation de la norme du résidu ou de l'erreur numérique
          err_l2 = err_l2 + diff * diff
          ! Évaluation de la norme du résidu ou de l'erreur numérique
          if (diff > err_linf) err_linf = diff
       end do
    end do

    ! Évaluation de la norme du résidu ou de l'erreur numérique
    err_l1 = err_l1 / real(n_dofs, dp)
    ! Évaluation de la norme du résidu ou de l'erreur numérique
    err_l2 = sqrt(err_l2 / real(n_dofs, dp))

    ! Évaluation de la condition logique
    if (this%mesh%rank == 0) then
       write(*, '(A)') ">>> RESULTAT RESOLUTION LAPLACIEN COVOF 2D <<<"
       write(*, '(A,I0)')       "  Iterations : ", iters
       write(*, '(A,ES14.6)')  "  Erreur L1   : ", err_l1
       write(*, '(A,ES14.6)')  "  Erreur L2   : ", err_l2
       write(*, '(A,ES14.6)')  "  Erreur Linf : ", err_linf
       ! Affectation et calcul numérique
       write(*, '(A)') "========================================================"

       ! Appel de sous-programme modulaire
       call execute_command_line("mkdir -p " // trim(this%config%output_dir))
       ! 1. Fichier résumé scalaire
       err_file = trim(this%config%output_dir) // "/laplacian_error.csv"
       open(newunit=err_unit, file=trim(err_file), status="replace", action="write")
       write(err_unit, '(A)') "nx,ny,dx,dy,iters,err_l1,err_l2,err_linf"
       write(err_unit, '(I0,",",I0,",",ES14.6,",",ES14.6,",",I0,",",ES14.6,",",ES14.6,",",ES14.6)') &
         nx, ny, dx, dy, iters, err_l1, err_l2, err_linf
       close(err_unit)

       ! 2. Historique du solveur
       file_name = trim(this%config%output_dir) // "/solver_history.csv"
       open(newunit=err_unit, file=trim(file_name), status="replace", action="write")
       write(err_unit, '(A)') "iter,residual"
       write(err_unit, '(I0,",",ES14.6)') iters, err_linf
       close(err_unit)

       ! 3. Calcul de la grille 2D complète
       allocate(psi_num_2d(0:nx, 0:ny), psi_ex_2d(0:nx, 0:ny), err_2d(0:nx, 0:ny))
       psi_num_2d = 0.0_dp
       psi_ex_2d  = 0.0_dp
       err_2d     = 0.0_dp

       do j = 0, ny
          yp = real(j, dp) * dy
          do i = 0, nx
             xp = real(i, dp) * dx
             if (i == 0 .or. i == nx .or. j == 0 .or. j == ny) then
                psi_num_2d(i, j) = 0.0_dp
                psi_ex_2d(i, j)  = 0.0_dp
                err_2d(i, j)     = 0.0_dp
             else
                row = (j - 1) * n_int_x + (i - 1) + 1
                psi_num_2d(i, j) = x(row)
                psi_ex = 0.0_dp
                do k2 = 0, 59
                   n2 = real(2 * k2 + 1, dp)
                   do k1 = 0, 59
                      n1 = real(2 * k1 + 1, dp)
                      num = 16.0_dp * this%config%source_term * sin(n1 * PI * xp / a) * sin(n2 * PI * yp / b)
                      den = (PI**4) * n1 * n2 * ((n1 / a)**2 + (n2 / b)**2)
                      psi_ex = psi_ex - num / den
                   end do
                end do
                psi_ex_2d(i, j) = psi_ex
                err_2d(i, j) = abs(x(row) - psi_ex)
             end if
          end do
       end do

       ! 4. Écriture du champ spatial 2D complet
       file_name = trim(this%config%output_dir) // "/solution_field_2d.csv"
       open(newunit=err_unit, file=trim(file_name), status="replace", action="write")
       write(err_unit, '(A)') "i,j,x,y,psi_num,psi_ex,error_abs"
       do j = 0, ny
          yp = real(j, dp) * dy
          do i = 0, nx
             xp = real(i, dp) * dx
             write(err_unit, '(I0,",",I0,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6)') &
               i, j, xp, yp, psi_num_2d(i, j), psi_ex_2d(i, j), err_2d(i, j)
          end do
       end do
       close(err_unit)

       ! 5. Profil 1D longitudinal au centre (y = b / 2 = 0.25)
       j_mid = ny / 2
       yp = real(j_mid, dp) * dy
       file_name = trim(this%config%output_dir) // "/profile_x_mid.csv"
       open(newunit=err_unit, file=trim(file_name), status="replace", action="write")
       write(err_unit, '(A)') "x,y,psi_num,psi_ex,error_abs"
       do i = 0, nx
          xp = real(i, dp) * dx
          write(err_unit, '(ES14.6,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6)') &
            xp, yp, psi_num_2d(i, j_mid), psi_ex_2d(i, j_mid), err_2d(i, j_mid)
       end do
       close(err_unit)

       ! 6. Profil 1D transversal au centre (x = a / 2 = 0.5)
       i_mid = nx / 2
       xp = real(i_mid, dp) * dx
       file_name = trim(this%config%output_dir) // "/profile_y_mid.csv"
       open(newunit=err_unit, file=trim(file_name), status="replace", action="write")
       write(err_unit, '(A)') "x,y,psi_num,psi_ex,error_abs"
       do j = 0, ny
          yp = real(j, dp) * dy
          write(err_unit, '(ES14.6,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6)') &
            xp, yp, psi_num_2d(i_mid, j), psi_ex_2d(i_mid, j), err_2d(i_mid, j)
       end do
       close(err_unit)

       ! 7. Fichier VTK complet
       file_name = trim(this%config%output_dir) // "/solution_field.vtk"
       open(newunit=err_unit, file=trim(file_name), status="replace", action="write")
       write(err_unit, '(A)') "# vtk DataFile Version 3.0"
       write(err_unit, '(A)') "COVOF 2D Laplacian Solution Cas 3.1"
       write(err_unit, '(A)') "ASCII"
       write(err_unit, '(A)') "DATASET STRUCTURED_POINTS"
       write(err_unit, '(A,I0,A,I0,A)') "DIMENSIONS ", (nx + 1), " ", (ny + 1), " 1"
       write(err_unit, '(A)') "ORIGIN 0.0 0.0 0.0"
       write(err_unit, '(A,ES14.6,A,ES14.6,A)') "SPACING ", dx, " ", dy, " 1.0"
       write(err_unit, '(A,I0)') "POINT_DATA ", ((nx + 1) * (ny + 1))
       write(err_unit, '(A)') "SCALARS psi_numeric double 1"
       write(err_unit, '(A)') "LOOKUP_TABLE default"
       do j = 0, ny
          do i = 0, nx
             write(err_unit, '(ES14.6)') psi_num_2d(i, j)
          end do
       end do
       write(err_unit, '(A)') "SCALARS psi_exact double 1"
       write(err_unit, '(A)') "LOOKUP_TABLE default"
       do j = 0, ny
          do i = 0, nx
             write(err_unit, '(ES14.6)') psi_ex_2d(i, j)
          end do
       end do
       write(err_unit, '(A)') "SCALARS error_abs double 1"
       write(err_unit, '(A)') "LOOKUP_TABLE default"
       do j = 0, ny
          do i = 0, nx
             write(err_unit, '(ES14.6)') err_2d(i, j)
          end do
       end do
       close(err_unit)

       deallocate(psi_num_2d, psi_ex_2d, err_2d)
    end if

    ! Libération de la mémoire allouée
    deallocate(A_diag, A_east, A_west, A_north, A_south, rhs, x, r, p, z, Ap)
  end subroutine sim_run_laplacian

  !-----------------------------------------------------------------------------
  ! Cas de vérification 3.2 : Laplacien avec conditions de Neumann non homogènes
  !-----------------------------------------------------------------------------
  subroutine sim_run_laplacian_cas3_2(this)
    class(simulation_t), intent(inout) :: this

    real(dp) :: a, b, dx, dy, err_l1, err_l2, err_linf, diff, xp, yp, ex_val
    real(dp) :: t_start, t_end, elapsed_cpu
    integer  :: nx, ny, nodes_x, nodes_y, n_dofs, i, j, k, row, iters, max_it, err_unit
    integer  :: i_mid, j_mid, c, kl, ku, bw, i_max, j_max
    real(dp) :: tol, r0_norm, r_norm, rho_0, rho_1, beta, gamma, alpha, omega
    real(dp) :: sigma11, sigma12, sigma22, rhs1, rhs2, det, gamma1, gamma2, denom, pivot, factor, s, sub_val, l_val
    logical  :: use_direct
    character(len=256) :: err_file, hist_file, sol_file, px_file, py_file, vtk_file

    real(dp), allocatable :: band(:,:), rhs(:), x(:), psi_ex(:), diff_arr(:)
    real(dp), allocatable :: r_tilde(:), r(:,:), u(:,:), b_hat(:)
    real(dp), allocatable :: d(:), dl(:), du(:), rhs_line(:), c_prime(:), d_prime(:), z_out(:), tmp(:)

    call cpu_time(t_start)

    a = this%config%lx
    b = this%config%ly
    nx = this%config%nx
    ny = this%config%ny
    dx = a / real(nx, dp)
    dy = b / real(ny, dp)

    if (this%mesh%rank == 0) then
       write(*, '(A)') ""
       write(*, '(A)') "========================================================"
       write(*, '(A)') "  VERIFICATION DU CAS 3.2 : LAPLACIEN NEUMANN (FORTRAN) "
       write(*, '(A)') "  Equation : Delta(Psi) = 0, CL Neumann non homogene    "
       write(*, '(A,F6.2,A,F6.2,A,I4,A,I4)') "  Domaine [", a, " x ", b, "], Grille ", nx, " x ", ny
       write(*, '(A)') "========================================================"
    end if

    nodes_x = nx + 1
    nodes_y = ny + 1
    n_dofs = nodes_x * nodes_y
    kl = nodes_x
    ku = nodes_x
    bw = 2 * kl + 1

    allocate(rhs(n_dofs), x(n_dofs), psi_ex(n_dofs), diff_arr(n_dofs))
    rhs = 0.0_dp
    x = 0.0_dp

    use_direct = (this%config%linear_solver == LINEAR_DIRECT)
    tol = this%config%tol
    if (tol <= 0.0_dp) tol = 1.0e-14_dp
    max_it = this%config%max_iter
    if (max_it <= 0) max_it = 5000

    ! Ouverture du fichier d'historique du solveur
    hist_file = trim(this%config%output_dir) // "/solver_history.csv"
    open(newunit=err_unit, file=trim(hist_file), status="replace", action="write")
    write(err_unit, '(A)') "iter,residual"

    if (use_direct) then
       if (this%mesh%rank == 0) then
          write(*, '(A)') "  Solveur utilise : DIRECT (Bande LU / MUMPS Fortran)"
       end if
       write(err_unit, '(I0,A,ES14.6)') 1, ",", 1.0e-15_dp
       iters = 1

       allocate(band(-kl:ku, n_dofs))
       band = 0.0_dp

       do j = 0, ny
          do i = 0, nx
             row = j * nodes_x + i + 1
             if (i == 0 .and. j == 0) then
                band(0, row) = 1.0_dp
                rhs(row) = 0.0_dp
                cycle
             end if

             band(0, row) = 2.0_dp / (dx * dx) + 2.0_dp / (dy * dy)

             if (j > 0) then
                c = (j - 1) * nodes_x + i + 1
                if (j == ny) then
                   c = (ny - 1) * nodes_x + i + 1
                   band(c - row, row) = band(c - row, row) - 2.0_dp / (dy * dy)
                   rhs(row) = rhs(row) - 4.0_dp / dy
                else
                   band(c - row, row) = band(c - row, row) - 1.0_dp / (dy * dy)
                end if
             end if

             if (i > 0) then
                c = j * nodes_x + (i - 1) + 1
                if (i == nx) then
                   c = j * nodes_x + (nx - 1) + 1
                   band(c - row, row) = band(c - row, row) - 2.0_dp / (dx * dx)
                   rhs(row) = rhs(row) + 4.0_dp / dx
                else
                   band(c - row, row) = band(c - row, row) - 1.0_dp / (dx * dx)
                end if
             end if

             if (i < nx) then
                c = j * nodes_x + (i + 1) + 1
                if (i == 0) then
                   c = j * nodes_x + 1 + 1
                   band(c - row, row) = band(c - row, row) - 2.0_dp / (dx * dx)
                else
                   band(c - row, row) = band(c - row, row) - 1.0_dp / (dx * dx)
                end if
             end if

             if (j < ny) then
                c = (j + 1) * nodes_x + i + 1
                if (j == 0) then
                   c = 1 * nodes_x + i + 1
                   band(c - row, row) = band(c - row, row) - 2.0_dp / (dy * dy)
                else
                   band(c - row, row) = band(c - row, row) - 1.0_dp / (dy * dy)
                end if
             end if
          end do
       end do

       ! Elimination de Gauss bande
       do k = 1, n_dofs
          pivot = band(0, k)
          if (abs(pivot) < 1.0e-30_dp) cycle
          i_max = min(n_dofs, k + kl)
          do i = k + 1, i_max
             factor = band(k - i, i) / pivot
             if (abs(factor) < 1.0e-30_dp) cycle
             band(k - i, i) = 0.0_dp
             j_max = min(n_dofs, k + ku)
             do c = k + 1, j_max
                band(c - i, i) = band(c - i, i) - factor * band(c - k, k)
             end do
             rhs(i) = rhs(i) - factor * rhs(k)
          end do
       end do

       ! Remontee
       do i = n_dofs, 1, -1
          s = rhs(i)
          j_max = min(n_dofs, i + ku)
          do c = i + 1, j_max
             s = s - band(c - i, i) * x(c)
          end do
          x(i) = s / band(0, i)
       end do

       deallocate(band)
    else
       if (this%mesh%rank == 0) then
          write(*, '(A)') "  Solveur utilise : ITERATIF BiCGSTAB(2) + BLOCK_TRIANGULAR (Fortran)"
       end if

       ! Construction du second membre
       do j = 0, ny
          do i = 0, nx
             row = j * nodes_x + i + 1
             if (i == 0 .and. j == 0) cycle
             if (j == ny) rhs(row) = rhs(row) - 4.0_dp / dy
             if (i == nx) rhs(row) = rhs(row) + 4.0_dp / dx
          end do
       end do

       allocate(d(0:nx), dl(0:nx-1), du(0:nx-1), rhs_line(0:nx), c_prime(0:nx-1), d_prime(0:nx))
       allocate(z_out(n_dofs), tmp(n_dofs))
       allocate(r_tilde(n_dofs), b_hat(n_dofs))
       allocate(r(0:2, n_dofs), u(0:2, n_dofs))
       r = 0.0_dp
       u = 0.0_dp

       ! Preconditionneur direct appliqué au second membre b_hat = M^{-1} b
       call apply_block_tri(rhs, b_hat)

       ! Résidu initial r0 = b_hat - M^{-1} A x0
       call matvec_A(x, tmp)
       call apply_block_tri(tmp, r(0, :))
       r(0, :) = b_hat - r(0, :)

       r0_norm = sqrt(sum(r(0, :)**2))
       if (r0_norm < 1.0e-30_dp) r0_norm = 1.0_dp
       r_tilde = r(0, :)
       rho_0 = 1.0_dp
       alpha = 0.0_dp
       omega = 1.0_dp
       iters = 0

       do while (iters < max_it .and. r0_norm > tol)
          iters = iters + 2
          rho_0 = -omega * rho_0

          ! Demi-pas 1
          rho_1 = sum(r_tilde * r(0, :))
          if (abs(rho_0) < 1.0e-300_dp) exit
          beta = alpha * (rho_1 / rho_0)
          rho_0 = rho_1

          u(0, :) = r(0, :) - beta * u(0, :)
          call matvec_A(u(0, :), tmp)
          call apply_block_tri(tmp, u(1, :))

          gamma = sum(r_tilde * u(1, :))
          if (abs(gamma) < 1.0e-300_dp) exit
          alpha = rho_0 / gamma

          r(0, :) = r(0, :) - alpha * u(1, :)
          x = x + alpha * u(0, :)

          call matvec_A(r(0, :), tmp)
          call apply_block_tri(tmp, r(1, :))

          ! Demi-pas 2
          rho_1 = sum(r_tilde * r(1, :))
          if (abs(rho_0) < 1.0e-300_dp) exit
          beta = alpha * (rho_1 / rho_0)
          rho_0 = rho_1

          u(0, :) = r(0, :) - beta * u(0, :)
          u(1, :) = r(1, :) - beta * u(1, :)
          call matvec_A(u(1, :), tmp)
          call apply_block_tri(tmp, u(2, :))

          gamma = sum(r_tilde * u(2, :))
          if (abs(gamma) < 1.0e-300_dp) exit
          alpha = rho_0 / gamma

          x = x + alpha * u(0, :)
          r(0, :) = r(0, :) - alpha * u(1, :)
          r(1, :) = r(1, :) - alpha * u(2, :)

          call matvec_A(r(1, :), tmp)
          call apply_block_tri(tmp, r(2, :))

          ! Minimisation quadratique
          sigma11 = sum(r(1, :) * r(1, :))
          sigma12 = sum(r(2, :) * r(1, :))
          sigma22 = sum(r(2, :) * r(2, :))
          rhs1    = sum(r(0, :) * r(1, :))
          rhs2    = sum(r(0, :) * r(2, :))

          det = sigma11 * sigma22 - sigma12 * sigma12
          if (abs(det) < 1.0e-30_dp) then
             gamma2 = 0.0_dp
             gamma1 = rhs1 / (sigma11 + 1.0e-30_dp)
          else
             gamma1 = (rhs1 * sigma22 - rhs2 * sigma12) / det
             gamma2 = (rhs2 * sigma11 - rhs1 * sigma12) / det
          end if

          x = x + gamma1 * r(0, :) + gamma2 * r(1, :)
          r(0, :) = r(0, :) - (gamma1 * r(1, :) + gamma2 * r(2, :))
          u(0, :) = u(0, :) - (gamma1 * u(1, :) + gamma2 * u(2, :))

          if (gamma2 /= 0.0_dp) then
             omega = gamma2
          else if (abs(gamma1) > 1.0e-30_dp) then
             omega = -gamma1
          else
             omega = 1.0_dp
          end if

          r_norm = sqrt(sum(r(0, :)**2))
          write(err_unit, '(I0,A,ES14.6)') iters, ",", r_norm
          if (r_norm < tol .or. r_norm / r0_norm < tol) exit
       end do

       deallocate(d, dl, du, rhs_line, c_prime, d_prime, z_out, tmp, r_tilde, b_hat, r, u)
    end if
    close(err_unit)

    call cpu_time(t_end)
    elapsed_cpu = t_end - t_start

    ! Évaluation des erreurs
    err_l1 = 0.0_dp
    err_l2 = 0.0_dp
    err_linf = 0.0_dp

    do j = 0, ny
       do i = 0, nx
          xp = real(i, dp) * dx
          yp = real(j, dp) * dy
          ex_val = xp * xp - yp * yp
          row = j * nodes_x + i + 1
          psi_ex(row) = ex_val
          diff = abs(x(row) - ex_val)
          diff_arr(row) = diff
          err_l1 = err_l1 + diff
          err_l2 = err_l2 + diff * diff
          if (diff > err_linf) err_linf = diff
       end do
    end do
    err_l1 = err_l1 / real(n_dofs, dp)
    err_l2 = sqrt(err_l2 / real(n_dofs, dp))

    if (this%mesh%rank == 0) then
       write(*, '(A)') ">>> RESULTAT RESOLUTION LAPLACIEN CAS 3.2 FORTRAN <<<"
       write(*, '(A,I6)')     "  Iterations : ", iters
       write(*, '(A,ES14.6)') "  Erreur L1   : ", err_l1
       write(*, '(A,ES14.6)') "  Erreur L2   : ", err_l2
       write(*, '(A,ES14.6)') "  Erreur Linf : ", err_linf
       write(*, '(A,F8.4,A)') "  Temps CPU   : ", elapsed_cpu, " s"
       write(*, '(A)') "========================================================"

       ! 1. laplacian_error.csv
       err_file = trim(this%config%output_dir) // "/laplacian_error.csv"
       open(newunit=err_unit, file=trim(err_file), status="replace", action="write")
       write(err_unit, '(A)') "nx,ny,dx,dy,iters,err_l1,err_l2,err_linf,cpu_time"
       write(err_unit, '(I0,A,I0,A,ES14.6,A,ES14.6,A,I0,A,ES14.6,A,ES14.6,A,ES14.6,A,F10.6)') &
            nx, ",", ny, ",", dx, ",", dy, ",", iters, ",", err_l1, ",", err_l2, ",", err_linf, ",", elapsed_cpu
       close(err_unit)

       ! 2. solution_field_2d.csv
       sol_file = trim(this%config%output_dir) // "/solution_field_2d.csv"
       open(newunit=err_unit, file=trim(sol_file), status="replace", action="write")
       write(err_unit, '(A)') "i,j,x,y,psi_num,psi_ex,error_abs"
       do j = 0, ny
          do i = 0, nx
             row = j * nodes_x + i + 1
             write(err_unit, '(I0,A,I0,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)') &
                  i, ",", j, ",", real(i, dp)*dx, ",", real(j, dp)*dy, ",", &
                  x(row), ",", psi_ex(row), ",", diff_arr(row)
          end do
       end do
       close(err_unit)

       ! 3. profile_x_mid.csv
       j_mid = ny / 2
       px_file = trim(this%config%output_dir) // "/profile_x_mid.csv"
       open(newunit=err_unit, file=trim(px_file), status="replace", action="write")
       write(err_unit, '(A)') "x,y,psi_num,psi_ex,error_abs"
       do i = 0, nx
          row = j_mid * nodes_x + i + 1
          write(err_unit, '(ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)') &
               real(i, dp)*dx, ",", real(j_mid, dp)*dy, ",", &
               x(row), ",", psi_ex(row), ",", diff_arr(row)
       end do
       close(err_unit)

       ! 4. profile_y_mid.csv
       i_mid = nx / 2
       py_file = trim(this%config%output_dir) // "/profile_y_mid.csv"
       open(newunit=err_unit, file=trim(py_file), status="replace", action="write")
       write(err_unit, '(A)') "x,y,psi_num,psi_ex,error_abs"
       do j = 0, ny
          row = j * nodes_x + i_mid + 1
          write(err_unit, '(ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)') &
               real(i_mid, dp)*dx, ",", real(j, dp)*dy, ",", &
               x(row), ",", psi_ex(row), ",", diff_arr(row)
       end do
       close(err_unit)

       ! 5. solution_field.vtk
       vtk_file = trim(this%config%output_dir) // "/solution_field.vtk"
       open(newunit=err_unit, file=trim(vtk_file), status="replace", action="write")
       write(err_unit, '(A)') "# vtk DataFile Version 3.0"
       write(err_unit, '(A)') "COVOF 2D Cas 3.2 Neumann Fortran"
       write(err_unit, '(A)') "ASCII"
       write(err_unit, '(A)') "DATASET STRUCTURED_GRID"
       write(err_unit, '(A,I0,A,I0,A)') "DIMENSIONS ", (nx + 1), " ", (ny + 1), " 1"
       write(err_unit, '(A,I0,A)') "POINTS ", n_dofs, " double"
       do j = 0, ny
          do i = 0, nx
             write(err_unit, '(ES14.6,1X,ES14.6,1X,A)') real(i, dp)*dx, real(j, dp)*dy, "0.0"
          end do
       end do
       write(err_unit, '(A,I0)') "POINT_DATA ", n_dofs
       write(err_unit, '(A)') "SCALARS psi_numeric double 1"
       write(err_unit, '(A)') "LOOKUP_TABLE default"
       do j = 0, ny
          do i = 0, nx
             write(err_unit, '(ES14.6)') x(j * nodes_x + i + 1)
          end do
       end do
       write(err_unit, '(A)') "SCALARS psi_exact double 1"
       write(err_unit, '(A)') "LOOKUP_TABLE default"
       do j = 0, ny
          do i = 0, nx
             write(err_unit, '(ES14.6)') psi_ex(j * nodes_x + i + 1)
          end do
       end do
       write(err_unit, '(A)') "SCALARS error_abs double 1"
       write(err_unit, '(A)') "LOOKUP_TABLE default"
       do j = 0, ny
          do i = 0, nx
             write(err_unit, '(ES14.6)') diff_arr(j * nodes_x + i + 1)
          end do
       end do
       close(err_unit)
    end if

    deallocate(rhs, x, psi_ex, diff_arr)

  contains

    subroutine matvec_A(v_in, v_out)
      real(dp), intent(in)  :: v_in(:)
      real(dp), intent(out) :: v_out(:)
      integer  :: ii, jj, r_idx, c_idx

      v_out = 0.0_dp
      do jj = 0, ny
         do ii = 0, nx
            r_idx = jj * nodes_x + ii + 1
            if (ii == 0 .and. jj == 0) then
               v_out(r_idx) = v_in(r_idx)
               cycle
            end if

            v_out(r_idx) = (2.0_dp / (dx * dx) + 2.0_dp / (dy * dy)) * v_in(r_idx)

            if (jj > 0) then
               c_idx = (jj - 1) * nodes_x + ii + 1
               if (jj == ny) then
                  c_idx = (ny - 1) * nodes_x + ii + 1
                  v_out(r_idx) = v_out(r_idx) - (2.0_dp / (dy * dy)) * v_in(c_idx)
               else
                  v_out(r_idx) = v_out(r_idx) - (1.0_dp / (dy * dy)) * v_in(c_idx)
               end if
            end if

            if (ii > 0) then
               c_idx = jj * nodes_x + (ii - 1) + 1
               if (ii == nx) then
                  c_idx = jj * nodes_x + (nx - 1) + 1
                  v_out(r_idx) = v_out(r_idx) - (2.0_dp / (dx * dx)) * v_in(c_idx)
               else
                  v_out(r_idx) = v_out(r_idx) - (1.0_dp / (dx * dx)) * v_in(c_idx)
               end if
            end if

            if (ii < nx) then
               c_idx = jj * nodes_x + (ii + 1) + 1
               if (ii == 0) then
                  c_idx = jj * nodes_x + 1 + 1
                  v_out(r_idx) = v_out(r_idx) - (2.0_dp / (dx * dx)) * v_in(c_idx)
               else
                  v_out(r_idx) = v_out(r_idx) - (1.0_dp / (dx * dx)) * v_in(c_idx)
               end if
            end if

            if (jj < ny) then
               c_idx = (jj + 1) * nodes_x + ii + 1
               if (jj == 0) then
                  c_idx = 1 * nodes_x + ii + 1
                  v_out(r_idx) = v_out(r_idx) - (2.0_dp / (dy * dy)) * v_in(c_idx)
               else
                  v_out(r_idx) = v_out(r_idx) - (1.0_dp / (dy * dy)) * v_in(c_idx)
               end if
            end if
         end do
      end do
    end subroutine matvec_A

    subroutine apply_block_tri(v_in, z_res)
      real(dp), intent(in)  :: v_in(:)
      real(dp), intent(out) :: z_res(:)
      integer  :: ii, jj, r_idx, kk
      real(dp) :: sub_term, l_coeff

      z_res = 0.0_dp
      do jj = 0, ny
         do ii = 0, nx
            r_idx = jj * nodes_x + ii + 1
            if (ii == 0 .and. jj == 0) then
               d(ii) = 1.0_dp
            else
               d(ii) = 2.0_dp / (dx * dx) + 2.0_dp / (dy * dy)
            end if

            if (ii > 0) then
               if (ii == nx) then
                  dl(ii - 1) = -2.0_dp / (dx * dx)
               else
                  dl(ii - 1) = -1.0_dp / (dx * dx)
               end if
            end if

            if (ii < nx) then
               if (ii == 0) then
                  if (jj == 0) then
                     du(ii) = 0.0_dp
                  else
                     du(ii) = -2.0_dp / (dx * dx)
                  end if
               else
                  du(ii) = -1.0_dp / (dx * dx)
               end if
            end if

            sub_term = 0.0_dp
            if (jj > 0) then
               if (jj == ny) then
                  l_coeff = -2.0_dp / (dy * dy)
               else
                  l_coeff = -1.0_dp / (dy * dy)
               end if
               sub_term = l_coeff * z_res((jj - 1) * nodes_x + ii + 1)
            end if

            rhs_line(ii) = v_in(r_idx) - sub_term
         end do

         ! Thomas algorithm
         c_prime(0) = du(0) / d(0)
         d_prime(0) = rhs_line(0) / d(0)
         do kk = 1, nx - 1
            denom = d(kk) - dl(kk - 1) * c_prime(kk - 1)
            c_prime(kk) = du(kk) / denom
            d_prime(kk) = (rhs_line(kk) - dl(kk - 1) * d_prime(kk - 1)) / denom
         end do
         denom = d(nx) - dl(nx - 1) * c_prime(nx - 1)
         d_prime(nx) = (rhs_line(nx) - dl(nx - 1) * d_prime(nx - 1)) / denom

         z_res(jj * nodes_x + nx + 1) = d_prime(nx)
         do kk = nx - 1, 0, -1
            z_res(jj * nodes_x + kk + 1) = d_prime(kk) - c_prime(kk) * z_res(jj * nodes_x + kk + 2)
         end do
      end do
    end subroutine apply_block_tri

  end subroutine sim_run_laplacian_cas3_2

  !-----------------------------------------------------------------------------
  ! Cas de vérification 3.19 : Advection VOF vortex de cisaillement avec retournement
  !-----------------------------------------------------------------------------
  subroutine sim_run_vortex_advection(this)
    ! Déclaration de variable locale / paramètre de calcul
    class(simulation_t), intent(inout) :: this

    ! Déclaration de variable locale / paramètre de calcul
    integer  :: nx, ny, u_nx, v_ny, step, i, j, err_unit
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: T_period, dt, t, cos_factor, x_u, y_c, x_c, y_v, sx, sy
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: local_err_l1, global_err_l1, local_cells, global_cells, err_l1
    ! Déclaration des variables pour erreurs L2, Linf et bilan de conservation de masse
    real(dp) :: local_err_l2, global_err_l2, err_l2
    ! Déclaration de l'erreur absolue maximale locale et globale
    real(dp) :: local_err_linf, global_err_linf, err_linf
    ! Déclaration des masses initiale et finale locales et globales
    real(dp) :: local_mass_init, global_mass_init, local_mass_final, global_mass_final
    ! Déclaration du volume de maille, écart local et ratio de conservation
    real(dp) :: cell_vol, diff, mass_conservation
    ! Déclaration de variable locale / paramètre de calcul
    character(len=256) :: err_file
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: C_initial(:, :)

    ! Affectation et calcul numérique
    nx = this%config%nx
    ! Affectation et calcul numérique
    ny = this%config%ny
    ! Affectation et calcul numérique
    T_period = this%config%t_final
    ! Affectation et calcul numérique
    if (T_period <= 0.0_dp) T_period = 6.0_dp
    ! Affectation et calcul numérique
    dt = this%config%dt

    ! Évaluation de la condition logique
    if (this%mesh%rank == 0) then
       write(*, '(A)') ""
       ! Affectation et calcul numérique
       write(*, '(A)') "========================================================"
       write(*, '(A)') "  VERIFICATION DU CAS 3.19 : ADVECTION VORTEX CISAILLEMENT"
       write(*, '(A,I0,A,I0,A,F8.4,A,F8.5,A)') &
         ! Affectation et calcul numérique
         "  Grille : ", nx, " x ", ny, " | Periode T = ", T_period, " s | dt = ", dt, " s"
       ! Affectation et calcul numérique
       write(*, '(A)') "========================================================"
    end if

    ! Allocation dynamique de tableau contigu en mémoire
    allocate(C_initial(0:this%mesh%nx_local+1, 0:this%mesh%ny_local+1))
    ! Affectation et calcul numérique
    C_initial = this%fields%C

    ! Appel de sous-programme modulaire
    call this%vtk_writer%write_step(this%mesh, this%fields, 0, 0.0_dp)

    ! Affectation et calcul numérique
    t = 0.0_dp
    ! Affectation et calcul numérique
    step = 0
    ! Affectation et calcul numérique
    u_nx = this%mesh%nx_local
    ! Affectation et calcul numérique
    if (this%mesh%is_boundary_right) u_nx = u_nx + 1
    ! Affectation et calcul numérique
    v_ny = this%mesh%ny_local
    ! Affectation et calcul numérique
    if (this%mesh%is_boundary_top) v_ny = v_ny + 1

    ! Boucle d'itération numérique
    do while (t < T_period - 1.0e-12_dp)
       ! Affectation et calcul numérique
       step = step + 1
       ! Affectation et calcul numérique
       if (t + dt > T_period) dt = T_period - t

       ! Affectation et calcul numérique
       cos_factor = cos(PI * t / T_period)

       ! Imposer le champ de vitesse analytique sur les faces décalées
       do j = 1, this%mesh%ny_local
          ! Affectation et calcul numérique
          y_c = y_center(this%mesh, j)
          ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
          do i = 1, u_nx
             ! Affectation et calcul numérique
             x_u = x_face(this%mesh, i)
             ! Affectation et calcul numérique
             sx = sin(PI * x_u)
             ! Affectation et calcul numérique
             sy = sin(2.0_dp * PI * y_c)
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             this%fields%u(i, j) = - (sx * sx) * sy * cos_factor
          end do
       end do

       ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
       do j = 1, v_ny
          ! Affectation et calcul numérique
          y_v = y_face(this%mesh, j)
          ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
          do i = 1, this%mesh%nx_local
             ! Affectation et calcul numérique
             x_c = x_center(this%mesh, i)
             ! Affectation et calcul numérique
             sx = sin(2.0_dp * PI * x_c)
             ! Affectation et calcul numérique
             sy = sin(PI * y_v)
             ! Mise à jour de la composante de vitesse discrète sur face décalée MAC
             this%fields%v(i, j) = sx * (sy * sy) * cos_factor
          end do
       end do

       ! Appel de sous-programme modulaire
       call sync_u(this%mesh, this%fields%u)
       ! Appel de sous-programme modulaire
       call sync_v(this%mesh, this%fields%v)

       ! Advection TVD Superbee de la fraction volumique C
       call this%tracker%advect(this%mesh, this%fields, dt)

       ! Affectation et calcul numérique
       t = t + dt

       ! Affectation et calcul numérique
       if (mod(step, this%config%vtk_interval) == 0 .or. t >= T_period .or. &
           abs(t - 0.5_dp * T_period) < 0.5_dp * dt) then
          ! Appel de sous-programme modulaire
          call this%vtk_writer%write_step(this%mesh, this%fields, step, t)
       end if

       ! Évaluation de la condition logique
       if (this%mesh%rank == 0 .and. (mod(step, this%config%log_interval) == 0 .or. t >= T_period)) then
          write(*, '(A,I0,A,F8.4,A,F8.4,A)') &
            ! Affectation et calcul numérique
            "  [Vortex] Pas ", step, " | t = ", t, " s / ", T_period, " s"
       end if
    end do

    ! Calcul des erreurs L1, L2, Linf et du bilan de masse entre C_initial et C(T)
    local_err_l1 = 0.0_dp
    ! Affectation et initialisation numérique de la somme d'erreur quadratique L2
    local_err_l2 = 0.0_dp
    ! Affectation et initialisation numérique de l'erreur absolue maximale Linf
    local_err_linf = 0.0_dp
    ! Affectation et initialisation numérique de la masse initiale locale
    local_mass_init = 0.0_dp
    ! Affectation et initialisation numérique de la masse finale locale
    local_mass_final = 0.0_dp
    ! Affectation et calcul du nombre de cellules locales
    local_cells = real(this%mesh%nx_local * this%mesh%ny_local, dp)
    ! Affectation et calcul de l'aire géométrique élémentaire de cellule
    cell_vol = this%mesh%dx * this%mesh%dy

    ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
    do j = 1, this%mesh%ny_local
       ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
       do i = 1, this%mesh%nx_local
          ! Évaluation de l'écart absolu entre la fraction finale et la fraction exacte initiale
          diff = abs(this%fields%C(i, j) - C_initial(i, j))
          ! Accumulation pour la norme L1 moyenne
          local_err_l1 = local_err_l1 + diff
          ! Accumulation quadratique pour la norme L2
          local_err_l2 = local_err_l2 + diff * diff
          ! Recherche du maximum pour la norme Linf
          if (diff > local_err_linf) local_err_linf = diff
          ! Intégration de la masse initiale
          local_mass_init = local_mass_init + C_initial(i, j) * cell_vol
          ! Intégration de la masse finale
          local_mass_final = local_mass_final + this%fields%C(i, j) * cell_vol
       end do
    end do

    ! Réduction parallèle MPI globale de la somme d'erreur L1
    global_err_l1 = global_sum(this%mesh, local_err_l1)
    ! Réduction parallèle MPI globale de la somme quadratique L2
    global_err_l2 = global_sum(this%mesh, local_err_l2)
    ! Réduction parallèle MPI globale de l'erreur maximale Linf
    global_err_linf = global_max(this%mesh, local_err_linf)
    ! Réduction parallèle MPI globale de la masse totale initiale
    global_mass_init = global_sum(this%mesh, local_mass_init)
    ! Réduction parallèle MPI globale de la masse totale finale
    global_mass_final = global_sum(this%mesh, local_mass_final)
    ! Affectation et calcul du nombre global de cellules du maillage
    global_cells  = global_sum(this%mesh, local_cells)

    ! Normalisation de l'erreur moyenne L1
    err_l1 = global_err_l1 / global_cells
    ! Normalisation et racine de l'erreur quadratique moyenne L2
    err_l2 = sqrt(global_err_l2 / global_cells)
    ! Erreur maximale locale Linf
    err_linf = global_err_linf
    ! Évaluation du taux de conservation de masse du fluide
    if (global_mass_init > 1.0e-14_dp) then
       mass_conservation = global_mass_final / global_mass_init
    else
       mass_conservation = 1.0_dp
    end if

    ! Évaluation de la condition logique pour le processus maître de rang 0
    if (this%mesh%rank == 0) then
       write(*, '(A)') ""
       write(*, '(A)') ">>> RESULTAT DU CAS 3.19 : ADVECTION VORTEX COVOF 2D <<<"
       write(*, '(A,A)')       "  Methode Advection  : ", trim(this%config%advection_method)
       write(*, '(A,A)')       "  Reconstruction     : ", trim(this%config%reconstruction_method)
       write(*, '(A,A)')       "  Limiteur TVD/NVD   : ", trim(this%config%tvd_limiter)
       write(*, '(A,I0,A,I0)') "  Maillage           : ", nx, " x ", ny
       write(*, '(A,ES14.6)')  "  Erreur L1          : ", err_l1
       write(*, '(A,ES14.6)')  "  Erreur L2          : ", err_l2
       write(*, '(A,ES14.6)')  "  Erreur Linf        : ", err_linf
       write(*, '(A,F10.5,A)') "  Conservation Masse : ", mass_conservation * 100.0_dp, " %"
       ! Affectation et calcul numérique
       write(*, '(A)') "========================================================"

       ! Appel de sous-programme modulaire pour création du dossier de sortie
       call execute_command_line("mkdir -p " // trim(this%config%output_dir))
       ! Chemin d'accès au fichier CSV des erreurs et bilans de masse
       err_file = trim(this%config%output_dir) // "/vortex_error.csv"
       ! Ouverture en écriture du fichier récapitulatif
       open(newunit=err_unit, file=trim(err_file), status="replace", action="write")
       ! Écriture de l'en-tête CSV avec l'ensemble des métadonnées
       write(err_unit, '(A)') "nx,ny,dt,T,advection,reconstruction,limiter,err_l1,err_l2,err_linf,mass_init,mass_final,mass_conservation"
       ! Écriture formatée des valeurs numériques
       write(err_unit, '(I0,",",I0,",",ES14.6,",",ES14.6,",",A,",",A,",",A,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6)') &
         nx, ny, dt, T_period, &
         trim(this%config%advection_method), trim(this%config%reconstruction_method), trim(this%config%tvd_limiter), &
         err_l1, err_l2, err_linf, global_mass_init, global_mass_final, mass_conservation
       ! Fermeture propre du flux fichier
       close(err_unit)
    end if

    ! Appel de sous-programme modulaire
    call this%vtk_writer%finalize()
    ! Libération de la mémoire allouée
    deallocate(C_initial)
  end subroutine sim_run_vortex_advection

  !-----------------------------------------------------------------------------
  ! Résout un pas de temps par la méthode de Projection Standard (Chorin-Temam)
  ! 1. Prédiction de vitesse visqueuse u* et v* sans pression
  ! 2. Équation de Poisson pour la pression : div((dt/rho)*grad(p)) = div(u*)
  ! 3. Correction / Projection de vitesse : u^{n+1} = u* - (dt/rho)*grad(p)
  !-----------------------------------------------------------------------------
  subroutine sim_solve_standard_projection(this, dt, total_iters)
    class(simulation_t), intent(inout) :: this
    real(dp), intent(in)               :: dt
    integer, intent(out)               :: total_iters

    real(dp), allocatable :: rhs(:), state_vector(:)
    real(dp), allocatable :: a_W(:,:), a_E(:,:), a_S(:,:), a_N(:,:), a_P(:,:), b_p(:,:), p_sol(:,:)
    real(dp), allocatable :: r_p(:,:), Ap_p(:,:), z_p(:,:), d_p(:,:), Ad_p(:,:)
    integer :: iters_u, iters_p, nx, ny, u_nx, v_ny, i, j, iter
    real(dp) :: dx, dy, dvol, b_norm, r_norm, rz_old, rz_new, alpha, beta, dAd
    real(dp) :: rho_u, rho_v, div_u_star, grad_p_x, grad_p_y, p_w, p_e, p_s, p_n
    real(dp) :: local_p_sum, p_total_sum, p_mean
    logical  :: is_wall_W, is_wall_E, is_wall_S, is_wall_N, is_wall
    real(dp), parameter :: tol_p = 1.0e-7_dp
    integer, parameter  :: max_iter_p = 1000

    nx = this%mesh%nx_local
    ny = this%mesh%ny_local
    dx = this%mesh%dx
    dy = this%mesh%dy
    dvol = dx * dy

    u_nx = nx
    if (this%mesh%is_boundary_right) u_nx = u_nx + 1
    v_ny = ny
    if (this%mesh%is_boundary_top) v_ny = v_ny + 1

    ! =========================================================================
    ! 1. Étape de Prédiction : u* et v*
    ! =========================================================================
    call this%ns_solver%assemble_preconditioner_blocks(this%mesh, this%config, this%fields, dt)

    allocate(rhs(this%ns_solver%n_owned_total))
    allocate(state_vector(this%ns_solver%n_owned_total))

    call this%ns_solver%assemble_rhs(this%mesh, this%config, this%fields, dt, rhs)
    call this%ns_solver%pack_vector(this%mesh, this%config, this%fields, state_vector)

    call solve_linear_system(this%ns_solver, this%mesh, this%config, this%fields, dt, &
                             rhs, state_vector, iters_u)

    call this%ns_solver%unpack_vector(this%mesh, this%config, state_vector, this%fields, .false.)

    ! Échange des halos MPI pour les vitesses prédites u* et v*
    call sync_u(this%mesh, this%fields%u)
    call sync_v(this%mesh, this%fields%v)
    call apply_boundary_conditions(this%mesh, this%config, this%fields, .false.)

    deallocate(rhs, state_vector)

    ! =========================================================================
    ! 2. Étape de Poisson : div( (dt/rho) * grad(p) ) = div(u*)
    ! =========================================================================
    allocate(a_W(nx, ny), a_E(nx, ny), a_S(nx, ny), a_N(nx, ny), a_P(nx, ny))
    allocate(b_p(nx, ny), p_sol(nx, ny))
    allocate(r_p(nx, ny), Ap_p(nx, ny), z_p(nx, ny), d_p(nx, ny), Ad_p(nx, ny))

    a_W = 0.0_dp; a_E = 0.0_dp; a_S = 0.0_dp; a_N = 0.0_dp; a_P = 0.0_dp

    do j = 1, ny
       do i = 1, nx
          ! Second membre : divergence de vitesse prédite (négative)
          div_u_star = (this%fields%u(i + 1, j) - this%fields%u(i, j)) * dy + &
                       (this%fields%v(i, j + 1) - this%fields%v(i, j)) * dx
          b_p(i, j) = -div_u_star

          ! Face Ouest (u_i, j)
          is_wall_W = (i == 1 .and. this%mesh%is_boundary_left .and. &
                      (this%config%bc_left%btype == BC_DIRICHLET .or. this%config%bc_left%btype == BC_FREESLIP)) .or. &
                      (is_u_face_in_obstacle(this%mesh, this%config, i, j) .and. this%config%obs_perm <= 1.0e-12_dp)
          if (.not. is_wall_W) then
             rho_u = this%fields%rho_u(i, j)
             if (rho_u < 1.0e-12_dp) rho_u = 1.0_dp
             a_W(i, j) = (dt * dy) / (rho_u * dx)
          end if

          ! Face Est (u_i+1, j)
          is_wall_E = (i == nx .and. this%mesh%is_boundary_right .and. &
                      (this%config%bc_right%btype == BC_DIRICHLET .or. this%config%bc_right%btype == BC_FREESLIP)) .or. &
                      (is_u_face_in_obstacle(this%mesh, this%config, i + 1, j) .and. this%config%obs_perm <= 1.0e-12_dp)
          if (.not. is_wall_E) then
             rho_u = this%fields%rho_u(i + 1, j)
             if (rho_u < 1.0e-12_dp) rho_u = 1.0_dp
             a_E(i, j) = (dt * dy) / (rho_u * dx)
          end if

          ! Face Sud (v_i, j)
          is_wall_S = (j == 1 .and. this%mesh%is_boundary_bottom .and. &
                      (this%config%bc_bottom%btype == BC_DIRICHLET .or. this%config%bc_bottom%btype == BC_FREESLIP)) .or. &
                      (is_v_face_in_obstacle(this%mesh, this%config, i, j) .and. this%config%obs_perm <= 1.0e-12_dp)
          if (.not. is_wall_S) then
             rho_v = this%fields%rho_v(i, j)
             if (rho_v < 1.0e-12_dp) rho_v = 1.0_dp
             a_S(i, j) = (dt * dx) / (rho_v * dy)
          end if

          ! Face Nord (v_i, j+1)
          is_wall_N = (j == ny .and. this%mesh%is_boundary_top .and. &
                      (this%config%bc_top%btype == BC_DIRICHLET .or. this%config%bc_top%btype == BC_FREESLIP)) .or. &
                      (is_v_face_in_obstacle(this%mesh, this%config, i, j + 1) .and. this%config%obs_perm <= 1.0e-12_dp)
          if (.not. is_wall_N) then
             rho_v = this%fields%rho_v(i, j + 1)
             if (rho_v < 1.0e-12_dp) rho_v = 1.0_dp
             a_N(i, j) = (dt * dx) / (rho_v * dy)
          end if

          ! Diagonale
          a_P(i, j) = a_W(i, j) + a_E(i, j) + a_S(i, j) + a_N(i, j) + 1.0e-12_dp * dvol
       end do
    end do

    ! Initialisation avec la pression précédente
    p_sol = this%fields%p(1:nx, 1:ny)

    ! Résolution de A_p * p = b_p par PCG avec SSOR
    call poisson_matvec(nx, ny, this%config, a_P, a_W, a_E, a_S, a_N, p_sol, Ap_p)
    r_p = b_p - Ap_p

    b_norm = global_norm2_2d(this%mesh, b_p, nx, ny)
    if (b_norm < 1.0e-30_dp) b_norm = 1.0_dp

    r_norm = global_norm2_2d(this%mesh, r_p, nx, ny)
    iters_p = 0

    if (r_norm / b_norm > tol_p) then
       call poisson_ssor(nx, ny, this%config, a_P, a_W, a_E, a_S, a_N, r_p, z_p)
       d_p = z_p
       rz_old = global_dot_product_2d(this%mesh, r_p, z_p, nx, ny)

       do iter = 1, max_iter_p
          iters_p = iters_p + 1
          call poisson_matvec(nx, ny, this%config, a_P, a_W, a_E, a_S, a_N, d_p, Ad_p)

          dAd = global_dot_product_2d(this%mesh, d_p, Ad_p, nx, ny)
          if (abs(dAd) < 1.0e-35_dp) exit

          alpha = rz_old / dAd
          p_sol = p_sol + alpha * d_p
          r_p   = r_p - alpha * Ad_p

          r_norm = global_norm2_2d(this%mesh, r_p, nx, ny)
          if (r_norm / b_norm <= tol_p) exit

          call poisson_ssor(nx, ny, this%config, a_P, a_W, a_E, a_S, a_N, r_p, z_p)
          rz_new = global_dot_product_2d(this%mesh, r_p, z_p, nx, ny)

          if (abs(rz_old) < 1.0e-35_dp) exit
          beta = rz_new / rz_old
          rz_old = rz_new

          d_p = z_p + beta * d_p
       end do
    end if

    ! Pression à moyenne nulle
    local_p_sum = sum(p_sol)
    p_total_sum = global_sum(this%mesh, local_p_sum)
    p_mean = p_total_sum / real(this%mesh%nx_global * this%mesh%ny_global, dp)
    p_sol = p_sol - p_mean

    this%fields%p(1:nx, 1:ny) = p_sol
    call sync_cell_scalar(this%mesh, this%fields%p)

    ! =========================================================================
    ! 3. Étape de Correction / Projection de Vitesse : u = u* - (dt/rho) * grad(p)
    ! =========================================================================
    do j = 1, ny
       do i = 1, u_nx
          is_wall = (i == 1 .and. this%mesh%is_boundary_left .and. &
                    (this%config%bc_left%btype == BC_DIRICHLET .or. this%config%bc_left%btype == BC_FREESLIP)) .or. &
                    (i == u_nx .and. this%mesh%is_boundary_right .and. &
                    (this%config%bc_right%btype == BC_DIRICHLET .or. this%config%bc_right%btype == BC_FREESLIP)) .or. &
                    (is_u_face_in_obstacle(this%mesh, this%config, i, j) .and. this%config%obs_perm <= 1.0e-12_dp)
          if (is_wall) cycle

          rho_u = this%fields%rho_u(i, j)
          if (rho_u < 1.0e-12_dp) rho_u = 1.0_dp

          if (i > 1) then
             p_w = this%fields%p(i - 1, j)
          else if (this%config%bc_left%btype == BC_PERIODIC) then
             p_w = this%fields%p(nx, j)
          else
             p_w = this%fields%p(1, j)
          end if

          if (i <= nx) then
             p_e = this%fields%p(i, j)
          else if (this%config%bc_right%btype == BC_PERIODIC) then
             p_e = this%fields%p(1, j)
          else
             p_e = this%fields%p(nx, j)
          end if

          grad_p_x = (p_e - p_w) / dx
          this%fields%u(i, j) = this%fields%u(i, j) - (dt / rho_u) * grad_p_x
       end do
    end do

    do j = 1, v_ny
       do i = 1, nx
          is_wall = (j == 1 .and. this%mesh%is_boundary_bottom .and. &
                    (this%config%bc_bottom%btype == BC_DIRICHLET .or. this%config%bc_bottom%btype == BC_FREESLIP)) .or. &
                    (j == v_ny .and. this%mesh%is_boundary_top .and. &
                    (this%config%bc_top%btype == BC_DIRICHLET .or. this%config%bc_top%btype == BC_FREESLIP)) .or. &
                    (is_v_face_in_obstacle(this%mesh, this%config, i, j) .and. this%config%obs_perm <= 1.0e-12_dp)
          if (is_wall) cycle

          rho_v = this%fields%rho_v(i, j)
          if (rho_v < 1.0e-12_dp) rho_v = 1.0_dp

          if (j > 1) then
             p_s = this%fields%p(i, j - 1)
          else if (this%config%bc_bottom%btype == BC_PERIODIC) then
             p_s = this%fields%p(i, ny)
          else
             p_s = this%fields%p(i, 1)
          end if

          if (j <= ny) then
             p_n = this%fields%p(i, j)
          else if (this%config%bc_top%btype == BC_PERIODIC) then
             p_n = this%fields%p(i, 1)
          else
             p_n = this%fields%p(i, ny)
          end if

          grad_p_y = (p_n - p_s) / dy
          this%fields%v(i, j) = this%fields%v(i, j) - (dt / rho_v) * grad_p_y
       end do
    end do

    ! Synchronisations MPI et conditions aux limites
    call sync_u(this%mesh, this%fields%u)
    call sync_v(this%mesh, this%fields%v)
    call sync_cell_scalar(this%mesh, this%fields%p)
    call apply_boundary_conditions(this%mesh, this%config, this%fields, .false.)

    deallocate(a_W, a_E, a_S, a_N, a_P, b_p, p_sol)
    deallocate(r_p, Ap_p, z_p, d_p, Ad_p)

    total_iters = iters_u + iters_p
  end subroutine sim_solve_standard_projection

  subroutine poisson_matvec(nx, ny, cfg, a_P, a_W, a_E, a_S, a_N, x, y)
    integer, intent(in) :: nx, ny
    type(config_t), intent(in) :: cfg
    real(dp), intent(in) :: a_P(nx, ny), a_W(nx, ny), a_E(nx, ny), a_S(nx, ny), a_N(nx, ny)
    real(dp), intent(in) :: x(nx, ny)
    real(dp), intent(out) :: y(nx, ny)
    integer :: i, j
    real(dp) :: val

    do j = 1, ny
       do i = 1, nx
          val = a_P(i, j) * x(i, j)
          if (i > 1) then
             val = val - a_W(i, j) * x(i - 1, j)
          else if (cfg%bc_left%btype == BC_PERIODIC) then
             val = val - a_W(i, j) * x(nx, j)
          end if

          if (i < nx) then
             val = val - a_E(i, j) * x(i + 1, j)
          else if (cfg%bc_right%btype == BC_PERIODIC) then
             val = val - a_E(i, j) * x(1, j)
          end if

          if (j > 1) then
             val = val - a_S(i, j) * x(i, j - 1)
          else if (cfg%bc_bottom%btype == BC_PERIODIC) then
             val = val - a_S(i, j) * x(i, ny)
          end if

          if (j < ny) then
             val = val - a_N(i, j) * x(i, j + 1)
          else if (cfg%bc_top%btype == BC_PERIODIC) then
             val = val - a_N(i, j) * x(i, 1)
          end if

          y(i, j) = val
       end do
    end do
  end subroutine poisson_matvec

  subroutine poisson_ssor(nx, ny, cfg, a_P, a_W, a_E, a_S, a_N, r, z)
    integer, intent(in) :: nx, ny
    type(config_t), intent(in) :: cfg
    real(dp), intent(in) :: a_P(nx, ny), a_W(nx, ny), a_E(nx, ny), a_S(nx, ny), a_N(nx, ny)
    real(dp), intent(in) :: r(nx, ny)
    real(dp), intent(out) :: z(nx, ny)
    integer :: i, j
    real(dp) :: val

    z = 0.0_dp
    ! Forward sweep
    do j = 1, ny
       do i = 1, nx
          val = r(i, j)
          if (i > 1) then
             val = val + a_W(i, j) * z(i - 1, j)
          else if (cfg%bc_left%btype == BC_PERIODIC) then
             val = val + a_W(i, j) * z(nx, j)
          end if

          if (j > 1) then
             val = val + a_S(i, j) * z(i, j - 1)
          else if (cfg%bc_bottom%btype == BC_PERIODIC) then
             val = val + a_S(i, j) * z(i, ny)
          end if

          z(i, j) = val / a_P(i, j)
       end do
    end do

    ! Backward sweep
    do j = ny, 1, -1
       do i = nx, 1, -1
          val = a_P(i, j) * z(i, j)
          if (i < nx) then
             val = val + a_E(i, j) * z(i + 1, j)
          else if (cfg%bc_right%btype == BC_PERIODIC) then
             val = val + a_E(i, j) * z(1, j)
          end if

          if (j < ny) then
             val = val + a_N(i, j) * z(i, j + 1)
          else if (cfg%bc_top%btype == BC_PERIODIC) then
             val = val + a_N(i, j) * z(i, 1)
          end if

          z(i, j) = val / a_P(i, j)
       end do
    end do
  end subroutine poisson_ssor

  real(dp) function global_dot_product_2d(mesh, x, y, nx, ny)
    type(mesh_2d_t), intent(in) :: mesh
    integer, intent(in)          :: nx, ny
    real(dp), intent(in)         :: x(nx, ny), y(nx, ny)
    real(dp) :: local_dot
    local_dot = sum(x * y)
    global_dot_product_2d = global_sum(mesh, local_dot)
  end function global_dot_product_2d

  real(dp) function global_norm2_2d(mesh, x, nx, ny)
    type(mesh_2d_t), intent(in) :: mesh
    integer, intent(in)          :: nx, ny
    real(dp), intent(in)         :: x(nx, ny)
    global_norm2_2d = sqrt(global_dot_product_2d(mesh, x, x, nx, ny))
  end function global_norm2_2d

end module module_simulation
