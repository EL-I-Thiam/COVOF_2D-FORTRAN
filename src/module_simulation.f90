!===============================================================================
! Module: module_simulation
! Description: Moteur principal de simulation FUGU 2D en Fortran 2018 :
!              - Gestion du cycle de vie du calcul (init, run, finalize)
!              - Pas de temps adaptatif CFL + onde capillaire
!              - Boucle temporelle avec couplage Lagrangien Augmenté / Uzawa
!              - Mode de vérification analytique du Laplacien (Cas 3.1 Thétis)
!              - Mode de vérification du vortex de cisaillement réversible (Cas 3.19 Thétis)
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
     type(config_t)         :: config
     type(mesh_2d_t)        :: mesh
     type(flow_fields_t)    :: fields
     type(interface_vof_t)  :: tracker
     type(navier_stokes_t)  :: ns_solver
     type(vtk_writer_t)     :: vtk_writer
  contains
     procedure :: init => sim_init
     procedure :: run => sim_run
     procedure :: compute_time_step => sim_compute_time_step
     procedure :: run_laplacian => sim_run_laplacian
     procedure :: run_vortex_advection => sim_run_vortex_advection
  end type simulation_t

contains

  !-----------------------------------------------------------------------------
  ! Initialisation complète de la simulation à partir du fichier JSON
  !-----------------------------------------------------------------------------
  function sim_init(this, comm, config_file) result(success)
    class(simulation_t), intent(inout) :: this
    integer, intent(in)                :: comm
    character(len=*), intent(in)       :: config_file
    logical                            :: success

    logical  :: periodic_x, periodic_y
    integer  :: i, j, u_nx
    real(dp) :: frac, y_step, h, u_mean, y_val, u_val, x_val

    success = .false.

    ! 1. Chargement de la configuration JSON
    call load_config_from_json(config_file, this%config)

    ! 2. Initialisation du maillage cartésien partitionné MPI
    periodic_x = (this%config%bc_left%btype == BC_PERIODIC .or. this%config%bc_right%btype == BC_PERIODIC)
    periodic_y = (this%config%bc_bottom%btype == BC_PERIODIC .or. this%config%bc_top%btype == BC_PERIODIC)

    call init_mesh(this%mesh, comm, this%config%nx, this%config%ny, &
                   this%config%lx, this%config%ly, periodic_x, periodic_y)

    ! 3. Allocation des champs décalés MAC
    call this%fields%allocate(this%mesh)

    ! 4. Initialisation du module d'interface VOF et des conditions initiales
    call this%tracker%init(this%config)
    call this%tracker%set_initial_condition(this%mesh, this%fields)

    ! Initialisation des conditions aux limites initiales
    call apply_boundary_conditions(this%mesh, this%config, this%fields, .false.)

    ! Initialisation du canal amont pour profil POISEUILLE_STEP afin de garantir la continuité
    if (this%config%bc_left%btype == BC_DIRICHLET .and. this%config%obs_enabled .and. &
        this%config%obs_type == OBS_BOX) then
       y_step = this%config%obs_y1
       h = this%mesh%ly - y_step
       u_mean = this%config%bc_left%u
       u_nx = this%mesh%nx_local
       if (this%mesh%is_boundary_right) u_nx = u_nx + 1

       do j = 1, this%mesh%ny_local
          y_val = y_center(this%mesh, j)
          u_val = 0.0_dp
          if (y_val >= y_step) then
             u_val = 6.0_dp * u_mean * (y_val - y_step) * (this%mesh%ly - y_val) / (h * h)
          end if
          do i = 1, u_nx
             x_val = x_face(this%mesh, i)
             if (x_val <= 4.0_dp + 1.0e-9_dp) then
                this%fields%u(i, j) = u_val
             end if
          end do
       end do
       call sync_u(this%mesh, this%fields%u)
    end if

    ! Initialisation du champ thermique si activé
    if (this%config%enable_energy) then
       do j = 1, this%mesh%ny_local
          do i = 1, this%mesh%nx_local
             frac = x_center(this%mesh, i) / this%mesh%lx
             this%fields%T(i, j) = this%config%bc_left%T + frac * (this%config%bc_right%T - this%config%bc_left%T)
          end do
       end do
       call sync_cell_scalar(this%mesh, this%fields%T)
       call apply_thermal_boundary_conditions(this%mesh, this%config, this%fields)
    end if

    ! 5. Initialisation du solveur Navier-Stokes et préconditionneurs
    call this%ns_solver%init(this%mesh, this%config)
    call this%ns_solver%assemble_preconditioner_blocks(this%mesh, this%config, this%fields, this%config%dt)

    ! 6. Initialisation de l'exportateur VTK et diagnostics CSV
    call this%vtk_writer%init(this%mesh, this%config%output_dir, this%config%case_name)

    success = .true.
  end function sim_init

  !-----------------------------------------------------------------------------
  ! Calcul du pas de temps adaptatif (CFL advection + onde capillaire)
  !-----------------------------------------------------------------------------
  function sim_compute_time_step(this) result(dt)
    class(simulation_t), intent(in) :: this
    real(dp) :: dt

    real(dp) :: local_max_speed, global_max_speed, speed, u_val, v_val
    real(dp) :: dt_cfl, dt_sigma, rho_sum, h_min, dt_target
    integer  :: i, j

    if (.not. this%config%adaptive_dt) then
       dt = this%config%dt
       return
    end if

    local_max_speed = 1.0e-6_dp
    do j = 1, this%mesh%ny_local
       do i = 1, this%mesh%nx_local
          u_val = abs(this%fields%u(i, j))
          v_val = abs(this%fields%v(i, j))
          speed = sqrt(u_val * u_val + v_val * v_val)
          if (speed > local_max_speed) local_max_speed = speed
       end do
    end do
    global_max_speed = global_max(this%mesh, local_max_speed)

    ! 1. Critère CFL de convection
    dt_cfl = this%config%cfl_max * min(this%mesh%dx, this%mesh%dy) / global_max_speed

    ! 2. Critère capillaire si tension superficielle active
    dt_sigma = 1.0e10_dp
    if (this%config%enable_surface_tension .and. this%config%sigma > 0.0_dp) then
       rho_sum = this%config%rho1 + this%config%rho2
       h_min = min(this%mesh%dx, this%mesh%dy)
       dt_sigma = sqrt((rho_sum * h_min**3) / (4.0_dp * PI * this%config%sigma))
    end if

    dt_target = min(dt_cfl, dt_sigma)
    dt = min(this%config%dt, dt_target)
  end function sim_compute_time_step

  !-----------------------------------------------------------------------------
  ! Boucle principale d'intégration temporelle Navier-Stokes
  !-----------------------------------------------------------------------------
  subroutine sim_run(this)
    class(simulation_t), intent(inout) :: this

    real(dp) :: t, dt, global_max_u, global_max_div, global_du_max, nu_mean
    real(dp) :: local_max_u, local_max_div, local_du_max, speed, du, div_val
    integer  :: step, total_iters, uzawa_max, m, iters, i, j
    real(dp), allocatable :: state_vector(:), rhs(:), u_prev(:, :)

    if (trim(this%config%mode) == "LAPLACIAN") then
       call this%run_laplacian()
       return
    end if
    if (trim(this%config%mode) == "VORTEX_ADVECTION") then
       call this%run_vortex_advection()
       return
    end if

    t = 0.0_dp
    step = 0

    allocate(state_vector(this%ns_solver%n_owned_total))
    allocate(rhs(this%ns_solver%n_owned_total))
    allocate(u_prev(0:this%mesh%nx_local+1, 0:this%mesh%ny_local+1))

    ! Sauvegarde initiale à t = 0.0
    call this%vtk_writer%write_step(this%mesh, this%fields, step, t)
    call this%vtk_writer%write_diagnostics(this%mesh, this%fields, step, t, this%config%dt)

    if (this%mesh%rank == 0) then
       write(*, '(A)') ""
       write(*, '(A)') ">>> DEMARRAGE DU CALCUL FUGU 2D FORTRAN 2018 PARALLELE <<<"
       write(*, '(A8,A14,A12,A14,A16,A16,A14)') &
         "Pas", "Temps [s]", "dt [s]", "Iter Krylov", "Vitesse Max", "Div Max", "Nu Moyen"
       write(*, '(A)') repeat("-", 94)
    end if

    do while (t < this%config%t_final - 1.0e-12_dp * (abs(this%config%t_final) + 1.0_dp))
       step = step + 1
       u_prev = this%fields%u
       dt = this%compute_time_step()

       if (t + dt > this%config%t_final) then
          dt = this%config%t_final - t
       end if
       if (dt < 1.0e-12_dp) exit

       ! 1. Transport de l'interface si écoulement diphasique
       if (this%config%enable_two_phase) then
          call this%tracker%advect(this%mesh, this%fields, dt)
       end if

       ! 2. Assemblage des blocs du préconditionneur
       call this%ns_solver%assemble_preconditioner_blocks(this%mesh, this%config, this%fields, dt)

       total_iters = 0
       if (this%config%solver_type == SOLVER_AUGMENTED_LAGRANGIAN) then
          uzawa_max = this%config%uzawa_iter
       else
          uzawa_max = 1
       end if

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

          ! 7. Mise à jour de la pression par Lagrangien augmenté (étape d'Uzawa)
          if (this%config%solver_type == SOLVER_AUGMENTED_LAGRANGIAN) then
             call this%ns_solver%update_pressure_al(this%mesh, this%config, this%fields)
          end if
       end do

       ! 7.5 Résolution de l'énergie thermique si activée
       if (this%config%enable_energy) then
          call this%ns_solver%solve_energy(this%mesh, this%config, this%fields, dt)
       end if

       t = t + dt

       ! Calcul des métriques de contrôle
       local_max_u   = 0.0_dp
       local_max_div = 0.0_dp
       local_du_max  = 0.0_dp

       do j = 1, this%mesh%ny_local
          do i = 1, this%mesh%nx_local
             speed = sqrt(this%fields%u(i, j)**2 + this%fields%v(i, j)**2)
             if (speed > local_max_u) local_max_u = speed

             du = abs(this%fields%u(i, j) - u_prev(i, j))
             if (du > local_du_max) local_du_max = du

             div_val = abs((this%fields%u(i + 1, j) - this%fields%u(i, j)) / this%mesh%dx + &
                           (this%fields%v(i, j + 1) - this%fields%v(i, j)) / this%mesh%dy)
             if (div_val > local_max_div) local_max_div = div_val
          end do
       end do

       global_max_u   = global_max(this%mesh, local_max_u)
       global_max_div = global_max(this%mesh, local_max_div)
       global_du_max  = global_max(this%mesh, local_du_max)

       nu_mean = 0.0_dp
       if (this%config%enable_energy) then
          nu_mean = this%ns_solver%compute_mean_nusselt_hot_wall(this%mesh, this%config, this%fields)
       end if

       ! Détection de convergence stationnaire
       if (.not. this%config%enable_two_phase .and. &
           this%config%steady_state_tol > 0.0_dp .and. step > 10 .and. &
           global_du_max < this%config%steady_state_tol) then
          if (this%mesh%rank == 0) then
             write(*, '(A)') ""
             write(*, '(A,ES12.4,A,ES12.4,A,F10.4,A)') &
               ">>> CONVERGENCE STATIONNAIRE ATTEINTE : du_max = ", global_du_max, &
               " < ", this%config%steady_state_tol, " (t = ", t, " s) <<<"
             if (this%config%enable_energy) then
                write(*, '(A,F10.5,A)') ">>> NOMBRE DE NUSSELT FINAL : Nu = ", nu_mean, " <<<"
             end if
          end if
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
          call this%vtk_writer%write_diagnostics(this%mesh, this%fields, step, t, dt)
       end if

       ! Export des fichiers VTK
       if (mod(step, this%config%vtk_interval) == 0 .or. t >= this%config%t_final) then
          call this%vtk_writer%write_step(this%mesh, this%fields, step, t)
       end if
    end do

    call this%vtk_writer%finalize()

    if (this%mesh%rank == 0) then
       write(*, '(A)') repeat("-", 94)
       write(*, '(A,F10.4,A)') ">>> SIMULATION FUGU 2D TERMINEE AVEC SUCCES A t = ", t, " s <<<"
       write(*, '(A,A)') "Fichiers de visualisation VTK et CSV disponibles dans : ", trim(this%config%output_dir)
    end if

    deallocate(state_vector, rhs, u_prev)
  end subroutine sim_run

  !-----------------------------------------------------------------------------
  ! Cas de vérification 3.1 : Laplacien -Delta(Psi) = -1 avec solution analytique
  !-----------------------------------------------------------------------------
  subroutine sim_run_laplacian(this)
    class(simulation_t), intent(inout) :: this

    double precision :: a, b, dx, dy, psi_ex, xp, yp, num, den, n1, n2, diff
    double precision :: err_l1, err_l2, err_linf, tol, rz_old, rz_new, alpha, beta, pAp, r_norm
    integer :: nx, ny, n_int_x, n_int_y, n_dofs, i, j, k1, k2, row, iters, max_it, err_unit
    character(len=256) :: err_file
    double precision, allocatable :: A_diag(:), A_east(:), A_west(:), A_north(:), A_south(:)
    double precision, allocatable :: rhs(:), x(:), r(:), p(:), z(:), Ap(:)

    a = this%config%lx
    b = this%config%ly
    nx = this%config%nx
    ny = this%config%ny
    dx = a / real(nx, dp)
    dy = b / real(ny, dp)

    if (this%mesh%rank == 0) then
       write(*, '(A)') ""
       write(*, '(A)') "========================================================"
       write(*, '(A)') "  VERIFICATION DU CAS 3.1 : LAPLACIEN DELTA(PSI) = 1"
       write(*, '(A,F8.4,A,F8.4,A,I0,A,I0,A)') &
         "  Domaine [", a, " x ", b, "], Grille ", nx, " x ", ny
       write(*, '(A)') "========================================================"
    end if

    n_int_x = nx - 1
    n_int_y = ny - 1
    n_dofs = n_int_x * n_int_y

    allocate(A_diag(n_dofs), A_east(n_dofs), A_west(n_dofs), A_north(n_dofs), A_south(n_dofs))
    allocate(rhs(n_dofs), x(n_dofs), r(n_dofs), p(n_dofs), z(n_dofs), Ap(n_dofs))

    A_diag  = 2.0_dp / (dx * dx) + 2.0_dp / (dy * dy)
    A_east  = -1.0_dp / (dx * dx)
    A_west  = -1.0_dp / (dx * dx)
    A_north = -1.0_dp / (dy * dy)
    A_south = -1.0_dp / (dy * dy)

    rhs = -1.0_dp
    x   = 0.0_dp

    ! Résolution par Gradient Conjugué préconditionné par Jacobi
    r = rhs
    z = r / A_diag
    p = z
    rz_old = sum(r * z)

    max_it = 20000
    tol = this%config%tol
    if (tol <= 0.0_dp) tol = 1.0e-12_dp

    do iters = 1, max_it
       ! Ap = A * p
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

    ! Évaluation de la solution analytique exacte (formule de Thétis)
    err_l1   = 0.0_dp
    err_l2   = 0.0_dp
    err_linf = 0.0_dp

    do j = 0, n_int_y - 1
       do i = 0, n_int_x - 1
          xp = real(i + 1, dp) * dx
          yp = real(j + 1, dp) * dy
          psi_ex = 0.0_dp

          do k2 = 0, 59
             do k1 = 0, 59
                n1 = real(2 * k1 + 1, dp)
                n2 = real(2 * k2 + 1, dp)
                num = 16.0_dp * sin(n1 * PI * xp / a) * sin(n2 * PI * yp / b)
                den = (PI**4) * n1 * n2 * ((n1 / a)**2 + (n2 / b)**2)
                psi_ex = psi_ex - num / den
             end do
          end do

          row = j * n_int_x + i + 1
          diff = abs(x(row) - psi_ex)
          err_l1 = err_l1 + diff
          err_l2 = err_l2 + diff * diff
          if (diff > err_linf) err_linf = diff
       end do
    end do

    err_l1 = err_l1 / real(n_dofs, dp)
    err_l2 = sqrt(err_l2 / real(n_dofs, dp))

    if (this%mesh%rank == 0) then
       write(*, '(A)') ">>> RESULTAT RESOLUTION LAPLACIEN FUGU 2D <<<"
       write(*, '(A,I0)')       "  Iterations : ", iters
       write(*, '(A,ES14.6)')  "  Erreur L1   : ", err_l1
       write(*, '(A,ES14.6)')  "  Erreur L2   : ", err_l2
       write(*, '(A,ES14.6)')  "  Erreur Linf : ", err_linf
       write(*, '(A)') "========================================================"

       call execute_command_line("mkdir -p " // trim(this%config%output_dir))
       err_file = trim(this%config%output_dir) // "/laplacian_error.csv"
       open(newunit=err_unit, file=trim(err_file), status="replace", action="write")
       write(err_unit, '(A)') "nx,ny,dx,dy,iters,err_l1,err_l2,err_linf"
       write(err_unit, '(I0,",",I0,",",ES14.6,",",ES14.6,",",I0,",",ES14.6,",",ES14.6,",",ES14.6)') &
         nx, ny, dx, dy, iters, err_l1, err_l2, err_linf
       close(err_unit)
    end if

    deallocate(A_diag, A_east, A_west, A_north, A_south, rhs, x, r, p, z, Ap)
  end subroutine sim_run_laplacian

  !-----------------------------------------------------------------------------
  ! Cas de vérification 3.19 : Advection VOF vortex de cisaillement avec retournement
  !-----------------------------------------------------------------------------
  subroutine sim_run_vortex_advection(this)
    class(simulation_t), intent(inout) :: this

    integer  :: nx, ny, u_nx, v_ny, step, i, j, err_unit
    real(dp) :: T_period, dt, t, cos_factor, x_u, y_c, x_c, y_v, sx, sy
    real(dp) :: local_err_l1, global_err_l1, local_cells, global_cells, err_l1
    character(len=256) :: err_file
    real(dp), allocatable :: C_initial(:, :)

    nx = this%config%nx
    ny = this%config%ny
    T_period = this%config%t_final
    if (T_period <= 0.0_dp) T_period = 6.0_dp
    dt = this%config%dt

    if (this%mesh%rank == 0) then
       write(*, '(A)') ""
       write(*, '(A)') "========================================================"
       write(*, '(A)') "  VERIFICATION DU CAS 3.19 : ADVECTION VORTEX CISAILLEMENT"
       write(*, '(A,I0,A,I0,A,F8.4,A,F8.5,A)') &
         "  Grille : ", nx, " x ", ny, " | Periode T = ", T_period, " s | dt = ", dt, " s"
       write(*, '(A)') "========================================================"
    end if

    allocate(C_initial(0:this%mesh%nx_local+1, 0:this%mesh%ny_local+1))
    C_initial = this%fields%C

    call this%vtk_writer%write_step(this%mesh, this%fields, 0, 0.0_dp)

    t = 0.0_dp
    step = 0
    u_nx = this%mesh%nx_local
    if (this%mesh%is_boundary_right) u_nx = u_nx + 1
    v_ny = this%mesh%ny_local
    if (this%mesh%is_boundary_top) v_ny = v_ny + 1

    do while (t < T_period - 1.0e-12_dp)
       step = step + 1
       if (t + dt > T_period) dt = T_period - t

       cos_factor = cos(PI * t / T_period)

       ! Imposer le champ de vitesse analytique sur les faces décalées
       do j = 1, this%mesh%ny_local
          y_c = y_center(this%mesh, j)
          do i = 1, u_nx
             x_u = x_face(this%mesh, i)
             sx = sin(PI * x_u)
             sy = sin(2.0_dp * PI * y_c)
             this%fields%u(i, j) = - (sx * sx) * sy * cos_factor
          end do
       end do

       do j = 1, v_ny
          y_v = y_face(this%mesh, j)
          do i = 1, this%mesh%nx_local
             x_c = x_center(this%mesh, i)
             sx = sin(2.0_dp * PI * x_c)
             sy = sin(PI * y_v)
             this%fields%v(i, j) = sx * (sy * sy) * cos_factor
          end do
       end do

       call sync_u(this%mesh, this%fields%u)
       call sync_v(this%mesh, this%fields%v)

       ! Advection TVD Superbee de la fraction volumique C
       call this%tracker%advect(this%mesh, this%fields, dt)

       t = t + dt

       if (mod(step, this%config%vtk_interval) == 0 .or. t >= T_period .or. &
           abs(t - 0.5_dp * T_period) < 0.5_dp * dt) then
          call this%vtk_writer%write_step(this%mesh, this%fields, step, t)
       end if

       if (this%mesh%rank == 0 .and. (mod(step, this%config%log_interval) == 0 .or. t >= T_period)) then
          write(*, '(A,I0,A,F8.4,A,F8.4,A)') &
            "  [Vortex] Pas ", step, " | t = ", t, " s / ", T_period, " s"
       end if
    end do

    ! Calcul de l'erreur L1 entre C_initial et C(T)
    local_err_l1 = 0.0_dp
    local_cells = real(this%mesh%nx_local * this%mesh%ny_local, dp)
    do j = 1, this%mesh%ny_local
       do i = 1, this%mesh%nx_local
          local_err_l1 = local_err_l1 + abs(this%fields%C(i, j) - C_initial(i, j))
       end do
    end do

    global_err_l1 = global_sum(this%mesh, local_err_l1)
    global_cells  = global_sum(this%mesh, local_cells)
    err_l1 = global_err_l1 / global_cells

    if (this%mesh%rank == 0) then
       write(*, '(A)') ""
       write(*, '(A)') ">>> RESULTAT DU CAS 3.19 : ADVECTION VORTEX FUGU 2D <<<"
       write(*, '(A,I0,A,I0)') "  Maillage : ", nx, " x ", ny
       write(*, '(A,ES14.6)')  "  Erreur L1 : ", err_l1
       write(*, '(A)') "========================================================"

       call execute_command_line("mkdir -p " // trim(this%config%output_dir))
       err_file = trim(this%config%output_dir) // "/vortex_error.csv"
       open(newunit=err_unit, file=trim(err_file), status="replace", action="write")
       write(err_unit, '(A)') "nx,ny,dt,T,err_l1"
       write(err_unit, '(I0,",",I0,",",ES14.6,",",ES14.6,",",ES14.6)') &
         nx, ny, dt, T_period, err_l1
       close(err_unit)
    end if

    call this%vtk_writer%finalize()
    deallocate(C_initial)
  end subroutine sim_run_vortex_advection

end module module_simulation
