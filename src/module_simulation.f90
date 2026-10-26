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
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: local_max_u, local_max_div, local_du_max, speed, du, div_val
    ! Déclaration de variable locale / paramètre de calcul
    integer  :: step, total_iters, uzawa_max, m, iters, i, j
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), allocatable :: state_vector(:), rhs(:), u_prev(:, :)

    ! Évaluation de la condition logique
    if (trim(this%config%mode) == "LAPLACIAN") then
       ! Appel de sous-programme modulaire
       call this%run_laplacian()
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

       ! 1. Transport de l'interface si écoulement diphasique
       if (this%config%enable_two_phase) then
          ! Appel de sous-programme modulaire
          call this%tracker%advect(this%mesh, this%fields, dt)
          ! Actualisation des conditions aux limites
          call apply_velocity_bc(this%fields, this%mesh, this%config, .false.)
          call apply_vof_bc(this%fields, this%mesh, this%config)
       end if

       ! 2. Assemblage des blocs du préconditionneur
       call this%ns_solver%assemble_preconditioner_blocks(this%mesh, this%config, this%fields, dt)

       ! Affectation et calcul numérique
       total_iters = 0
       ! Évaluation de la méthode de résolution Navier-Stokes
       if (this%config%solver_type == SOLVER_STANDARD_PROJECTION) then
          call this%solve_standard_projection(dt, total_iters)
       else
          if (this%config%solver_type == SOLVER_AUGMENTED_LAGRANGIAN) then
             ! Affectation et calcul numérique
             uzawa_max = this%config%uzawa_iter
          ! Branche par défaut si condition non satisfaite
          else
             ! Affectation et calcul numérique
             uzawa_max = 1
          end if

          ! Boucle d'itération numérique
          do m = 1, uzawa_max
             ! 3. Construction du second membre rhs
             call this%ns_solver%assemble_rhs(this%mesh, this%config, this%fields, dt, rhs)

             ! 4. Initialisation du vecteur d'état à partir des champs courants
             call this%ns_solver%pack_vector(this%mesh, this%config, this%fields, state_vector)

             ! 5. Résolution du système linéaire par Krylov (BiCGSTAB2 / PCG / GMRES)
             call solve_linear_system(this%ns_solver, this%mesh, this%config, this%fields, dt, &
                                      rhs, state_vector, iters)
             ! Affectation et calcul numérique
             total_iters = total_iters + iters

             ! 6. Déballage de la vitesse résolue
             call this%ns_solver%unpack_vector(this%mesh, this%config, state_vector, this%fields, .false.)

             ! 7. Mise à jour de la pression par Lagrangien augmenté (étape d'Uzawa)
             if (this%config%solver_type == SOLVER_AUGMENTED_LAGRANGIAN) then
                ! Appel de sous-programme modulaire
                call this%ns_solver%update_pressure_al(this%mesh, this%config, this%fields)
             end if
          end do
       end if

       ! 7.5 Résolution de l'énergie thermique si activée
       if (this%config%enable_energy) then
          ! Appel de sous-programme modulaire
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
       if (mod(step, this%config%vtk_interval) == 0 .or. t >= this%config%t_final) then
          ! Appel de sous-programme modulaire
          call this%vtk_writer%write_step(this%mesh, this%fields, step, t)
       end if
    end do

    ! Appel de sous-programme modulaire
    call this%vtk_writer%finalize()

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
  ! Cas de vérification 3.1 : Laplacien -Delta(Psi) = -1 avec solution analytique
  !-----------------------------------------------------------------------------
  subroutine sim_run_laplacian(this)
    ! Déclaration de variable locale / paramètre de calcul
    class(simulation_t), intent(inout) :: this

    double precision :: a, b, dx, dy, psi_ex, xp, yp, num, den, n1, n2, diff
    double precision :: err_l1, err_l2, err_linf, tol, rz_old, rz_new, alpha, beta, pAp, r_norm
    ! Déclaration de variable locale / paramètre de calcul
    integer :: nx, ny, n_int_x, n_int_y, n_dofs, i, j, k1, k2, row, iters, max_it, err_unit
    ! Déclaration de variable locale / paramètre de calcul
    character(len=256) :: err_file
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
    rhs = -1.0_dp
    ! Affectation et calcul numérique
    x   = 0.0_dp

    ! Résolution par Gradient Conjugué préconditionné par Jacobi
    r = rhs
    ! Affectation et calcul numérique
    z = r / A_diag
    ! Affectation et calcul numérique
    p = z
    ! Affectation et calcul numérique
    rz_old = sum(r * z)

    ! Affectation et calcul numérique
    max_it = 20000
    ! Affectation et calcul numérique
    tol = this%config%tol
    ! Affectation et calcul numérique
    if (tol <= 0.0_dp) tol = 1.0e-12_dp

    ! Boucle d'itération numérique
    do iters = 1, max_it
       ! Ap = A * p
       do j = 0, n_int_y - 1
          ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
          do i = 0, n_int_x - 1
             ! Affectation et calcul numérique
             row = j * n_int_x + i + 1
             ! Actualisation de la pression statique au centre de la maille
             Ap(row) = A_diag(row) * p(row)
             ! Actualisation de la pression statique au centre de la maille
             if (i > 0)           Ap(row) = Ap(row) + A_west(row)  * p(row - 1)
             ! Actualisation de la pression statique au centre de la maille
             if (i < n_int_x - 1) Ap(row) = Ap(row) + A_east(row)  * p(row + 1)
             ! Actualisation de la pression statique au centre de la maille
             if (j > 0)           Ap(row) = Ap(row) + A_south(row) * p(row - n_int_x)
             ! Actualisation de la pression statique au centre de la maille
             if (j < n_int_y - 1) Ap(row) = Ap(row) + A_north(row) * p(row + n_int_x)
          end do
       end do

       ! Affectation et calcul numérique
       pAp = sum(p * Ap)
       if (abs(pAp) < 1.0e-30_dp) exit

       ! Affectation et calcul numérique
       alpha = rz_old / pAp
       ! Affectation et calcul numérique
       x = x + alpha * p
       ! Affectation et calcul numérique
       r = r - alpha * Ap

       ! Évaluation de la norme du résidu ou de l'erreur numérique
       r_norm = sqrt(sum(r * r))
       if (r_norm < tol) exit

       ! Affectation et calcul numérique
       z = r / A_diag
       ! Affectation et calcul numérique
       rz_new = sum(r * z)
       ! Affectation et calcul numérique
       beta = rz_new / rz_old
       ! Affectation et calcul numérique
       rz_old = rz_new

       ! Affectation et calcul numérique
       p = z + beta * p
    end do

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
                num = 16.0_dp * sin(n1 * PI * xp / a) * sin(n2 * PI * yp / b)
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
       ! Évaluation de la norme du résidu ou de l'erreur numérique
       err_file = trim(this%config%output_dir) // "/laplacian_error.csv"
       ! Affectation et calcul numérique
       open(newunit=err_unit, file=trim(err_file), status="replace", action="write")
       write(err_unit, '(A)') "nx,ny,dx,dy,iters,err_l1,err_l2,err_linf"
       write(err_unit, '(I0,",",I0,",",ES14.6,",",ES14.6,",",I0,",",ES14.6,",",ES14.6,",",ES14.6)') &
         nx, ny, dx, dy, iters, err_l1, err_l2, err_linf
       close(err_unit)
    end if

    ! Libération de la mémoire allouée
    deallocate(A_diag, A_east, A_west, A_north, A_south, rhs, x, r, p, z, Ap)
  end subroutine sim_run_laplacian

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
                      is_u_face_in_obstacle(this%mesh, this%config, i, j)
          if (.not. is_wall_W) then
             rho_u = this%fields%rho_u(i, j)
             if (rho_u < 1.0e-12_dp) rho_u = 1.0_dp
             a_W(i, j) = (dt * dy) / (rho_u * dx)
          end if

          ! Face Est (u_i+1, j)
          is_wall_E = (i == nx .and. this%mesh%is_boundary_right .and. &
                      (this%config%bc_right%btype == BC_DIRICHLET .or. this%config%bc_right%btype == BC_FREESLIP)) .or. &
                      is_u_face_in_obstacle(this%mesh, this%config, i + 1, j)
          if (.not. is_wall_E) then
             rho_u = this%fields%rho_u(i + 1, j)
             if (rho_u < 1.0e-12_dp) rho_u = 1.0_dp
             a_E(i, j) = (dt * dy) / (rho_u * dx)
          end if

          ! Face Sud (v_i, j)
          is_wall_S = (j == 1 .and. this%mesh%is_boundary_bottom .and. &
                      (this%config%bc_bottom%btype == BC_DIRICHLET .or. this%config%bc_bottom%btype == BC_FREESLIP)) .or. &
                      is_v_face_in_obstacle(this%mesh, this%config, i, j)
          if (.not. is_wall_S) then
             rho_v = this%fields%rho_v(i, j)
             if (rho_v < 1.0e-12_dp) rho_v = 1.0_dp
             a_S(i, j) = (dt * dx) / (rho_v * dy)
          end if

          ! Face Nord (v_i, j+1)
          is_wall_N = (j == ny .and. this%mesh%is_boundary_top .and. &
                      (this%config%bc_top%btype == BC_DIRICHLET .or. this%config%bc_top%btype == BC_FREESLIP)) .or. &
                      is_v_face_in_obstacle(this%mesh, this%config, i, j + 1)
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

    ! =========================================================================
    ! 3. Étape de Correction / Projection de Vitesse : u = u* - (dt/rho) * grad(p)
    ! =========================================================================
    do j = 1, ny
       do i = 1, u_nx
          is_wall = (i == 1 .and. this%mesh%is_boundary_left .and. &
                    (this%config%bc_left%btype == BC_DIRICHLET .or. this%config%bc_left%btype == BC_FREESLIP)) .or. &
                    (i == u_nx .and. this%mesh%is_boundary_right .and. &
                    (this%config%bc_right%btype == BC_DIRICHLET .or. this%config%bc_right%btype == BC_FREESLIP)) .or. &
                    is_u_face_in_obstacle(this%mesh, this%config, i, j)
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
                    is_v_face_in_obstacle(this%mesh, this%config, i, j)
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
