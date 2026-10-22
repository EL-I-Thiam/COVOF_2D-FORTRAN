!===============================================================================
! Module: module_vtk_writer
! Description: Exportateur de maillage rectiligne XML VTK parallèle (.pvtr / .vtr)
!              et enregistreur de séries temporelles CSV pour COVOF 2D.
!              - Écriture locale non-bloquante du fichier .vtr sur chaque processus
!              - Assemblage du fichier maître .pvtr sur le rang racine
!              - Diagnostics diphasiques : position du centre de masse, vitesse
!                d'ascension, volume de bulle, énergie cinétique, divergence max
!===============================================================================
module module_vtk_writer
  use module_precision, only: dp
  use module_mesh, only: mesh_2d_t, x_center, y_center
  use module_fields, only: flow_fields_t
  use module_mpi_comm, only: global_sum, global_max
  implicit none
  private
  include 'mpif.h'

  public :: vtk_writer_t

  type :: vtk_writer_t
     ! Déclaration de variable locale / paramètre de calcul
     character(len=256) :: output_dir = "output"
     ! Déclaration de variable locale / paramètre de calcul
     character(len=128) :: case_name  = "covof2d_case"
     ! Déclaration de variable locale / paramètre de calcul
     integer            :: csv_unit   = 0
     ! Déclaration de variable locale / paramètre de calcul
     logical            :: is_open    = .false.
  contains
     ! Affectation et calcul numérique
     procedure :: init => vtk_init
     ! Affectation et calcul numérique
     procedure :: write_step => vtk_write_step
     ! Affectation et calcul numérique
     procedure :: write_diagnostics => vtk_write_diagnostics
     ! Affectation et calcul numérique
     procedure :: finalize => vtk_finalize
  end type vtk_writer_t

contains

  !-----------------------------------------------------------------------------
  ! Initialisation du répertoire de sortie et du fichier de métriques CSV
  !-----------------------------------------------------------------------------
  subroutine vtk_init(this, mesh, out_dir, name)
    ! Déclaration de variable locale / paramètre de calcul
    class(vtk_writer_t), intent(inout) :: this
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)        :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    character(len=*), intent(in)       :: out_dir, name

    ! Déclaration de variable locale / paramètre de calcul
    character(len=512) :: csv_path, mkdir_cmd
    ! Déclaration de variable locale / paramètre de calcul
    integer :: io_stat

    ! Affectation et calcul numérique
    this%output_dir = trim(out_dir)
    ! Affectation et calcul numérique
    this%case_name  = trim(name)

    ! Évaluation de la condition logique
    if (mesh%rank == 0) then
       ! Affectation et calcul numérique
       mkdir_cmd = "mkdir -p " // trim(this%output_dir)
       ! Appel de sous-programme modulaire
       call execute_command_line(trim(mkdir_cmd))

       ! Affectation et calcul numérique
       csv_path = trim(this%output_dir) // "/" // trim(this%case_name) // "_history.csv"
       ! Affectation et calcul numérique
       open(newunit=this%csv_unit, file=trim(csv_path), status="replace", action="write", iostat=io_stat)
       ! Évaluation de la condition logique
       if (io_stat == 0) then
          ! Affectation et calcul numérique
          this%is_open = .true.
          write(this%csv_unit, '(A)') &
            "step,time,dt,center_of_mass_y,rise_velocity,bubble_volume,kinetic_energy,max_u,max_div"
          flush(this%csv_unit)
       end if
    end if
  end subroutine vtk_init

  !-----------------------------------------------------------------------------
  ! Écriture d'un pas de temps complet en format XML RectilinearGrid (.pvtr / .vtr)
  !-----------------------------------------------------------------------------
  subroutine vtk_write_step(this, mesh, fields, step, time)
    ! Déclaration de variable locale / paramètre de calcul
    class(vtk_writer_t), intent(in) :: this
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)     :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(in) :: fields
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(in)             :: step
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)            :: time

    ! Déclaration de variable locale / paramètre de calcul
    character(len=6)   :: step_str
    ! Déclaration de variable locale / paramètre de calcul
    character(len=512) :: vtr_file, pvtr_file, piece_str
    ! Déclaration de variable locale / paramètre de calcul
    integer :: x_start, x_end, y_start, y_end
    ! Déclaration de variable locale / paramètre de calcul
    integer :: vtr_unit, pvtr_unit, io_stat, i, j, p, ierr, dummy(1)
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: uc, vc, div_val
    ! Déclaration de variable locale / paramètre de calcul
    integer, allocatable :: all_x_start(:), all_x_end(:), all_y_start(:), all_y_end(:)

    write(step_str, '(I6.6)') step

    ! Affectation et calcul numérique
    x_start = mesh%i_start_global - 1
    ! Affectation et calcul numérique
    x_end   = x_start + mesh%nx_local
    ! Affectation et calcul numérique
    y_start = mesh%j_start_global - 1
    ! Affectation et calcul numérique
    y_end   = y_start + mesh%ny_local

    ! 1. Écriture du fichier .vtr propre à chaque processus MPI
    write(vtr_file, '(A,"/",A,"_p",I0,"_",A,".vtr")') &
      trim(this%output_dir), trim(this%case_name), mesh%rank, step_str

    ! Affectation et calcul numérique
    open(newunit=vtr_unit, file=trim(vtr_file), status="replace", action="write", iostat=io_stat)
    ! Évaluation de la condition logique
    if (io_stat == 0) then
       ! Affectation et calcul numérique
       write(vtr_unit, '(A)') '<?xml version="1.0"?>'
       ! Affectation et calcul numérique
       write(vtr_unit, '(A)') '<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">'
       write(vtr_unit, '(A,I0,A,I0,A,I0,A,I0,A)') &
         ! Affectation et calcul numérique
         '  <RectilinearGrid WholeExtent="', x_start, ' ', x_end, ' ', y_start, ' ', y_end, ' 0 0">'
       write(vtr_unit, '(A,I0,A,I0,A,I0,A,I0,A)') &
         ! Affectation et calcul numérique
         '    <Piece Extent="', x_start, ' ', x_end, ' ', y_start, ' ', y_end, ' 0 0">'
       ! Affectation et calcul numérique
       write(vtr_unit, '(A)') '      <CellData Scalars="PhaseIndicator,Pressure,Density,Viscosity,Divergence" Vectors="Velocity">'

       ! PhaseIndicator C
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="PhaseIndicator" format="ascii">'
       ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
       do j = 1, mesh%ny_local
          ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
          do i = 1, mesh%nx_local
             ! Affectation et calcul numérique
             write(vtr_unit, '(ES16.8," ")', advance="no") fields%C(i, j)
          end do
          write(vtr_unit, *)
       end do
       write(vtr_unit, '(A)') '        </DataArray>'

       ! Pressure p
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="Pressure" format="ascii">'
       ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
       do j = 1, mesh%ny_local
          ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
          do i = 1, mesh%nx_local
             ! Affectation et calcul numérique
             write(vtr_unit, '(ES16.8," ")', advance="no") fields%p(i, j)
          end do
          write(vtr_unit, *)
       end do
       write(vtr_unit, '(A)') '        </DataArray>'

       ! Density rho
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="Density" format="ascii">'
       ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
       do j = 1, mesh%ny_local
          ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
          do i = 1, mesh%nx_local
             ! Affectation et calcul numérique
             write(vtr_unit, '(ES16.8," ")', advance="no") fields%rho(i, j)
          end do
          write(vtr_unit, *)
       end do
       write(vtr_unit, '(A)') '        </DataArray>'

       ! Viscosity mu
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="Viscosity" format="ascii">'
       ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
       do j = 1, mesh%ny_local
          ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
          do i = 1, mesh%nx_local
             ! Affectation et calcul numérique
             write(vtr_unit, '(ES16.8," ")', advance="no") fields%mu(i, j)
          end do
          write(vtr_unit, *)
       end do
       write(vtr_unit, '(A)') '        </DataArray>'

       ! Divergence div(u)
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="Divergence" format="ascii">'
       ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
       do j = 1, mesh%ny_local
          ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
          do i = 1, mesh%nx_local
             ! Affectation et calcul numérique
             div_val = (fields%u(i + 1, j) - fields%u(i, j)) / mesh%dx + &
                       (fields%v(i, j + 1) - fields%v(i, j)) / mesh%dy
             ! Affectation et calcul numérique
             write(vtr_unit, '(ES16.8," ")', advance="no") div_val
          end do
          write(vtr_unit, *)
       end do
       write(vtr_unit, '(A)') '        </DataArray>'

       ! Temperature T
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="Temperature" format="ascii">'
       ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
       do j = 1, mesh%ny_local
          ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
          do i = 1, mesh%nx_local
             ! Affectation et calcul numérique
             write(vtr_unit, '(ES16.8," ")', advance="no") fields%T(i, j)
          end do
          write(vtr_unit, *)
       end do
       write(vtr_unit, '(A)') '        </DataArray>'

       ! Velocity Vector (au centre des mailles)
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="Velocity" NumberOfComponents="3" format="ascii">'
       ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
       do j = 1, mesh%ny_local
          ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
          do i = 1, mesh%nx_local
             ! Affectation et calcul numérique
             uc = 0.5_dp * (fields%u(i, j) + fields%u(i + 1, j))
             ! Affectation et calcul numérique
             vc = 0.5_dp * (fields%v(i, j) + fields%v(i, j + 1))
             ! Affectation et calcul numérique
             write(vtr_unit, '(ES16.8," ",ES16.8," 0.0 ")', advance="no") uc, vc
          end do
          write(vtr_unit, *)
       end do
       write(vtr_unit, '(A)') '        </DataArray>'

       write(vtr_unit, '(A)') '      </CellData>'
       write(vtr_unit, '(A)') '      <Coordinates>'

       ! Coordonnées X
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="XCoordinates" format="ascii">'
       ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
       do i = 0, mesh%nx_local
          ! Affectation et calcul numérique
          write(vtr_unit, '(ES16.8," ")', advance="no") real(mesh%i_start_global - 1 + i, dp) * mesh%dx
       end do
       write(vtr_unit, *)
       write(vtr_unit, '(A)') '        </DataArray>'

       ! Coordonnées Y
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="YCoordinates" format="ascii">'
       ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
       do j = 0, mesh%ny_local
          ! Affectation et calcul numérique
          write(vtr_unit, '(ES16.8," ")', advance="no") real(mesh%j_start_global - 1 + j, dp) * mesh%dy
       end do
       write(vtr_unit, *)
       write(vtr_unit, '(A)') '        </DataArray>'

       ! Coordonnées Z
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="ZCoordinates" format="ascii">'
       write(vtr_unit, '(A)') '          0.0'
       write(vtr_unit, '(A)') '        </DataArray>'

       write(vtr_unit, '(A)') '      </Coordinates>'
       write(vtr_unit, '(A)') '    </Piece>'
       write(vtr_unit, '(A)') '  </RectilinearGrid>'
       write(vtr_unit, '(A)') '</VTKFile>'
       close(vtr_unit)
    end if

    ! 2. Écriture du fichier maître parallèle .pvtr sur le rang 0
    if (mesh%rank == 0) then
       ! Allocation dynamique de tableau contigu en mémoire
       allocate(all_x_start(mesh%size), all_x_end(mesh%size))
       ! Allocation dynamique de tableau contigu en mémoire
       allocate(all_y_start(mesh%size), all_y_end(mesh%size))

       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Gather(x_start, 1, MPI_INTEGER, all_x_start, 1, MPI_INTEGER, 0, mesh%comm, ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Gather(x_end,   1, MPI_INTEGER, all_x_end,   1, MPI_INTEGER, 0, mesh%comm, ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Gather(y_start, 1, MPI_INTEGER, all_y_start, 1, MPI_INTEGER, 0, mesh%comm, ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Gather(y_end,   1, MPI_INTEGER, all_y_end,   1, MPI_INTEGER, 0, mesh%comm, ierr)

       write(pvtr_file, '(A,"/",A,"_",A,".pvtr")') &
         trim(this%output_dir), trim(this%case_name), step_str

       ! Affectation et calcul numérique
       open(newunit=pvtr_unit, file=trim(pvtr_file), status="replace", action="write", iostat=io_stat)
       ! Évaluation de la condition logique
       if (io_stat == 0) then
          ! Affectation et calcul numérique
          write(pvtr_unit, '(A)') '<?xml version="1.0"?>'
          ! Affectation et calcul numérique
          write(pvtr_unit, '(A)') '<VTKFile type="PRectilinearGrid" version="0.1" byte_order="LittleEndian">'
          write(pvtr_unit, '(A,I0,A,I0,A)') &
            ! Affectation et calcul numérique
            '  <PRectilinearGrid WholeExtent="0 ', mesh%nx_global, ' 0 ', mesh%ny_global, ' 0 0" GhostLevel="0">'
          ! Affectation et calcul numérique
          write(pvtr_unit, '(A)') '    <PCellData Scalars="PhaseIndicator,Pressure,Density,Viscosity,Divergence" Vectors="Velocity">'
          ! Affectation et calcul numérique
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="PhaseIndicator"/>'
          ! Affectation et calcul numérique
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="Pressure"/>'
          ! Affectation et calcul numérique
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="Density"/>'
          ! Affectation et calcul numérique
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="Viscosity"/>'
          ! Affectation et calcul numérique
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="Divergence"/>'
          ! Affectation et calcul numérique
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="Temperature"/>'
          ! Affectation et calcul numérique
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="Velocity" NumberOfComponents="3"/>'
          write(pvtr_unit, '(A)') '    </PCellData>'
          write(pvtr_unit, '(A)') '    <PCoordinates>'
          ! Affectation et calcul numérique
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="XCoordinates"/>'
          ! Affectation et calcul numérique
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="YCoordinates"/>'
          ! Affectation et calcul numérique
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="ZCoordinates"/>'
          write(pvtr_unit, '(A)') '    </PCoordinates>'

          ! Boucle d'itération numérique
          do p = 0, mesh%size - 1
             write(piece_str, '(A,"_p",I0,"_",A,".vtr")') trim(this%case_name), p, step_str
             write(pvtr_unit, '(A,I0,A,I0,A,I0,A,I0,A,A,A)') &
               ! Affectation et calcul numérique
               '    <Piece Extent="', all_x_start(p+1), ' ', all_x_end(p+1), ' ', &
               ! Affectation et calcul numérique
               all_y_start(p+1), ' ', all_y_end(p+1), ' 0 0" Source="', trim(piece_str), '"/>'
          end do

          write(pvtr_unit, '(A)') '  </PRectilinearGrid>'
          write(pvtr_unit, '(A)') '</VTKFile>'
          close(pvtr_unit)
       end if

       ! Libération de la mémoire allouée
       deallocate(all_x_start, all_x_end, all_y_start, all_y_end)
    ! Branche par défaut si condition non satisfaite
    else
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Gather(x_start, 1, MPI_INTEGER, dummy, 0, MPI_INTEGER, 0, mesh%comm, ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Gather(x_end,   1, MPI_INTEGER, dummy, 0, MPI_INTEGER, 0, mesh%comm, ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Gather(y_start, 1, MPI_INTEGER, dummy, 0, MPI_INTEGER, 0, mesh%comm, ierr)
       ! Appel à la bibliothèque de communication parallèle MPI
       call MPI_Gather(y_end,   1, MPI_INTEGER, dummy, 0, MPI_INTEGER, 0, mesh%comm, ierr)
    end if
  end subroutine vtk_write_step

  !-----------------------------------------------------------------------------
  ! Calcul et enregistrement des métriques physiques temporelles dans l'historique CSV
  !-----------------------------------------------------------------------------
  subroutine vtk_write_diagnostics(this, mesh, fields, step, time, dt)
    ! Déclaration de variable locale / paramètre de calcul
    class(vtk_writer_t), intent(in) :: this
    ! Déclaration de variable locale / paramètre de calcul
    type(mesh_2d_t), intent(in)     :: mesh
    ! Déclaration de variable locale / paramètre de calcul
    type(flow_fields_t), intent(in) :: fields
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(in)             :: step
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in)            :: time, dt

    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: local_vol, local_y_integral, local_v_integral, local_ke
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: local_max_u, local_max_div, dvol, y_pos, c_val, ind_bubble
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: uc, vc, speed2, speed, div_val
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: global_vol, global_y_integral, global_v_integral, global_ke
    ! Déclaration de variable locale / paramètre de calcul
    real(dp) :: global_max_u, global_max_div, yc, vc_mean
    ! Déclaration de variable locale / paramètre de calcul
    integer  :: i, j

    ! Affectation et calcul numérique
    local_vol        = 0.0_dp
    ! Affectation et calcul numérique
    local_y_integral = 0.0_dp
    ! Affectation et calcul numérique
    local_v_integral = 0.0_dp
    ! Affectation et calcul numérique
    local_ke         = 0.0_dp
    ! Affectation et calcul numérique
    local_max_u      = 0.0_dp
    ! Affectation et calcul numérique
    local_max_div    = 0.0_dp

    ! Affectation et calcul numérique
    dvol = mesh%dx * mesh%dy

    ! Boucle spatiale sur les cellules locales selon l'axe vertical Y (ordonnées)
    do j = 1, mesh%ny_local
       ! Affectation et calcul numérique
       y_pos = y_center(mesh, j)
       ! Boucle spatiale selon l'axe horizontal X (accès contigu colonne-major)
       do i = 1, mesh%nx_local
          ! Affectation et calcul numérique
          c_val = fields%C(i, j)
          ! Affectation et calcul numérique
          ind_bubble = 1.0_dp - c_val ! 1 si phase fluide 2 (bulle / goutte)

          ! Affectation et calcul numérique
          local_vol        = local_vol + ind_bubble * dvol
          ! Affectation et calcul numérique
          local_y_integral = local_y_integral + y_pos * ind_bubble * dvol

          ! Affectation et calcul numérique
          uc = 0.5_dp * (fields%u(i, j) + fields%u(i + 1, j))
          ! Affectation et calcul numérique
          vc = 0.5_dp * (fields%v(i, j) + fields%v(i, j + 1))

          ! Affectation et calcul numérique
          local_v_integral = local_v_integral + vc * ind_bubble * dvol

          ! Affectation et calcul numérique
          speed2 = uc * uc + vc * vc
          ! Affectation et calcul numérique
          local_ke = local_ke + 0.5_dp * fields%rho(i, j) * speed2 * dvol

          ! Affectation et calcul numérique
          speed = sqrt(speed2)
          ! Affectation et calcul numérique
          if (speed > local_max_u) local_max_u = speed

          ! Affectation et calcul numérique
          div_val = abs((fields%u(i + 1, j) - fields%u(i, j)) / mesh%dx + &
                        (fields%v(i, j + 1) - fields%v(i, j)) / mesh%dy)
          ! Affectation et calcul numérique
          if (div_val > local_max_div) local_max_div = div_val
       end do
    end do

    ! Affectation et calcul numérique
    global_vol        = global_sum(mesh, local_vol)
    ! Affectation et calcul numérique
    global_y_integral = global_sum(mesh, local_y_integral)
    ! Affectation et calcul numérique
    global_v_integral = global_sum(mesh, local_v_integral)
    ! Affectation et calcul numérique
    global_ke         = global_sum(mesh, local_ke)
    ! Affectation et calcul numérique
    global_max_u      = global_max(mesh, local_max_u)
    ! Affectation et calcul numérique
    global_max_div    = global_max(mesh, local_max_div)

    ! Évaluation de la condition logique
    if (global_vol > 1.0e-12_dp) then
       ! Affectation et calcul numérique
       yc = global_y_integral / global_vol
       ! Affectation et calcul numérique
       vc_mean = global_v_integral / global_vol
    ! Branche par défaut si condition non satisfaite
    else
       ! Affectation et calcul numérique
       yc = 0.0_dp
       ! Affectation et calcul numérique
       vc_mean = 0.0_dp
    end if

    ! Évaluation de la condition logique
    if (mesh%rank == 0 .and. this%is_open) then
       write(this%csv_unit, '(I0,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6,",",ES14.6)') &
         step, time, dt, yc, vc_mean, global_vol, global_ke, global_max_u, global_max_div
       flush(this%csv_unit)
    end if
  end subroutine vtk_write_diagnostics

  !-----------------------------------------------------------------------------
  ! Clôture du fichier de diagnostic CSV
  !-----------------------------------------------------------------------------
  subroutine vtk_finalize(this)
    ! Déclaration de variable locale / paramètre de calcul
    class(vtk_writer_t), intent(inout) :: this
    ! Évaluation de la condition logique
    if (this%is_open) then
       close(this%csv_unit)
       ! Affectation et calcul numérique
       this%is_open = .false.
    end if
  end subroutine vtk_finalize

end module module_vtk_writer
