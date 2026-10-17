!===============================================================================
! Module: module_vtk_writer
! Description: Exportateur de maillage rectiligne XML VTK parallèle (.pvtr / .vtr)
!              et enregistreur de séries temporelles CSV pour FUGU 2D.
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
     character(len=256) :: output_dir = "output"
     character(len=128) :: case_name  = "fugu2d_case"
     integer            :: csv_unit   = 0
     logical            :: is_open    = .false.
  contains
     procedure :: init => vtk_init
     procedure :: write_step => vtk_write_step
     procedure :: write_diagnostics => vtk_write_diagnostics
     procedure :: finalize => vtk_finalize
  end type vtk_writer_t

contains

  !-----------------------------------------------------------------------------
  ! Initialisation du répertoire de sortie et du fichier de métriques CSV
  !-----------------------------------------------------------------------------
  subroutine vtk_init(this, mesh, out_dir, name)
    class(vtk_writer_t), intent(inout) :: this
    type(mesh_2d_t), intent(in)        :: mesh
    character(len=*), intent(in)       :: out_dir, name

    character(len=512) :: csv_path, mkdir_cmd
    integer :: io_stat

    this%output_dir = trim(out_dir)
    this%case_name  = trim(name)

    if (mesh%rank == 0) then
       mkdir_cmd = "mkdir -p " // trim(this%output_dir)
       call execute_command_line(trim(mkdir_cmd))

       csv_path = trim(this%output_dir) // "/" // trim(this%case_name) // "_history.csv"
       open(newunit=this%csv_unit, file=trim(csv_path), status="replace", action="write", iostat=io_stat)
       if (io_stat == 0) then
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
    class(vtk_writer_t), intent(in) :: this
    type(mesh_2d_t), intent(in)     :: mesh
    type(flow_fields_t), intent(in) :: fields
    integer, intent(in)             :: step
    real(dp), intent(in)            :: time

    character(len=6)   :: step_str
    character(len=512) :: vtr_file, pvtr_file, piece_str
    integer :: x_start, x_end, y_start, y_end
    integer :: vtr_unit, pvtr_unit, io_stat, i, j, p, ierr, dummy(1)
    real(dp) :: uc, vc, div_val
    integer, allocatable :: all_x_start(:), all_x_end(:), all_y_start(:), all_y_end(:)

    write(step_str, '(I6.6)') step

    x_start = mesh%i_start_global - 1
    x_end   = x_start + mesh%nx_local
    y_start = mesh%j_start_global - 1
    y_end   = y_start + mesh%ny_local

    ! 1. Écriture du fichier .vtr propre à chaque processus MPI
    write(vtr_file, '(A,"/",A,"_p",I0,"_",A,".vtr")') &
      trim(this%output_dir), trim(this%case_name), mesh%rank, step_str

    open(newunit=vtr_unit, file=trim(vtr_file), status="replace", action="write", iostat=io_stat)
    if (io_stat == 0) then
       write(vtr_unit, '(A)') '<?xml version="1.0"?>'
       write(vtr_unit, '(A)') '<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">'
       write(vtr_unit, '(A,I0,A,I0,A,I0,A,I0,A)') &
         '  <RectilinearGrid WholeExtent="', x_start, ' ', x_end, ' ', y_start, ' ', y_end, ' 0 0">'
       write(vtr_unit, '(A,I0,A,I0,A,I0,A,I0,A)') &
         '    <Piece Extent="', x_start, ' ', x_end, ' ', y_start, ' ', y_end, ' 0 0">'
       write(vtr_unit, '(A)') '      <CellData Scalars="PhaseIndicator,Pressure,Density,Viscosity,Divergence" Vectors="Velocity">'

       ! PhaseIndicator C
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="PhaseIndicator" format="ascii">'
       do j = 1, mesh%ny_local
          do i = 1, mesh%nx_local
             write(vtr_unit, '(ES16.8," ")', advance="no") fields%C(i, j)
          end do
          write(vtr_unit, *)
       end do
       write(vtr_unit, '(A)') '        </DataArray>'

       ! Pressure p
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="Pressure" format="ascii">'
       do j = 1, mesh%ny_local
          do i = 1, mesh%nx_local
             write(vtr_unit, '(ES16.8," ")', advance="no") fields%p(i, j)
          end do
          write(vtr_unit, *)
       end do
       write(vtr_unit, '(A)') '        </DataArray>'

       ! Density rho
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="Density" format="ascii">'
       do j = 1, mesh%ny_local
          do i = 1, mesh%nx_local
             write(vtr_unit, '(ES16.8," ")', advance="no") fields%rho(i, j)
          end do
          write(vtr_unit, *)
       end do
       write(vtr_unit, '(A)') '        </DataArray>'

       ! Viscosity mu
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="Viscosity" format="ascii">'
       do j = 1, mesh%ny_local
          do i = 1, mesh%nx_local
             write(vtr_unit, '(ES16.8," ")', advance="no") fields%mu(i, j)
          end do
          write(vtr_unit, *)
       end do
       write(vtr_unit, '(A)') '        </DataArray>'

       ! Divergence div(u)
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="Divergence" format="ascii">'
       do j = 1, mesh%ny_local
          do i = 1, mesh%nx_local
             div_val = (fields%u(i + 1, j) - fields%u(i, j)) / mesh%dx + &
                       (fields%v(i, j + 1) - fields%v(i, j)) / mesh%dy
             write(vtr_unit, '(ES16.8," ")', advance="no") div_val
          end do
          write(vtr_unit, *)
       end do
       write(vtr_unit, '(A)') '        </DataArray>'

       ! Temperature T
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="Temperature" format="ascii">'
       do j = 1, mesh%ny_local
          do i = 1, mesh%nx_local
             write(vtr_unit, '(ES16.8," ")', advance="no") fields%T(i, j)
          end do
          write(vtr_unit, *)
       end do
       write(vtr_unit, '(A)') '        </DataArray>'

       ! Velocity Vector (au centre des mailles)
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="Velocity" NumberOfComponents="3" format="ascii">'
       do j = 1, mesh%ny_local
          do i = 1, mesh%nx_local
             uc = 0.5_dp * (fields%u(i, j) + fields%u(i + 1, j))
             vc = 0.5_dp * (fields%v(i, j) + fields%v(i, j + 1))
             write(vtr_unit, '(ES16.8," ",ES16.8," 0.0 ")', advance="no") uc, vc
          end do
          write(vtr_unit, *)
       end do
       write(vtr_unit, '(A)') '        </DataArray>'

       write(vtr_unit, '(A)') '      </CellData>'
       write(vtr_unit, '(A)') '      <Coordinates>'

       ! Coordonnées X
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="XCoordinates" format="ascii">'
       do i = 0, mesh%nx_local
          write(vtr_unit, '(ES16.8," ")', advance="no") real(mesh%i_start_global - 1 + i, dp) * mesh%dx
       end do
       write(vtr_unit, *)
       write(vtr_unit, '(A)') '        </DataArray>'

       ! Coordonnées Y
       write(vtr_unit, '(A)') '        <DataArray type="Float64" Name="YCoordinates" format="ascii">'
       do j = 0, mesh%ny_local
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
       allocate(all_x_start(mesh%size), all_x_end(mesh%size))
       allocate(all_y_start(mesh%size), all_y_end(mesh%size))

       call MPI_Gather(x_start, 1, MPI_INTEGER, all_x_start, 1, MPI_INTEGER, 0, mesh%comm, ierr)
       call MPI_Gather(x_end,   1, MPI_INTEGER, all_x_end,   1, MPI_INTEGER, 0, mesh%comm, ierr)
       call MPI_Gather(y_start, 1, MPI_INTEGER, all_y_start, 1, MPI_INTEGER, 0, mesh%comm, ierr)
       call MPI_Gather(y_end,   1, MPI_INTEGER, all_y_end,   1, MPI_INTEGER, 0, mesh%comm, ierr)

       write(pvtr_file, '(A,"/",A,"_",A,".pvtr")') &
         trim(this%output_dir), trim(this%case_name), step_str

       open(newunit=pvtr_unit, file=trim(pvtr_file), status="replace", action="write", iostat=io_stat)
       if (io_stat == 0) then
          write(pvtr_unit, '(A)') '<?xml version="1.0"?>'
          write(pvtr_unit, '(A)') '<VTKFile type="PRectilinearGrid" version="0.1" byte_order="LittleEndian">'
          write(pvtr_unit, '(A,I0,A,I0,A)') &
            '  <PRectilinearGrid WholeExtent="0 ', mesh%nx_global, ' 0 ', mesh%ny_global, ' 0 0" GhostLevel="0">'
          write(pvtr_unit, '(A)') '    <PCellData Scalars="PhaseIndicator,Pressure,Density,Viscosity,Divergence" Vectors="Velocity">'
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="PhaseIndicator"/>'
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="Pressure"/>'
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="Density"/>'
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="Viscosity"/>'
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="Divergence"/>'
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="Temperature"/>'
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="Velocity" NumberOfComponents="3"/>'
          write(pvtr_unit, '(A)') '    </PCellData>'
          write(pvtr_unit, '(A)') '    <PCoordinates>'
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="XCoordinates"/>'
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="YCoordinates"/>'
          write(pvtr_unit, '(A)') '      <PDataArray type="Float64" Name="ZCoordinates"/>'
          write(pvtr_unit, '(A)') '    </PCoordinates>'

          do p = 0, mesh%size - 1
             write(piece_str, '(A,"_p",I0,"_",A,".vtr")') trim(this%case_name), p, step_str
             write(pvtr_unit, '(A,I0,A,I0,A,I0,A,I0,A,A,A)') &
               '    <Piece Extent="', all_x_start(p+1), ' ', all_x_end(p+1), ' ', &
               all_y_start(p+1), ' ', all_y_end(p+1), ' 0 0" Source="', trim(piece_str), '"/>'
          end do

          write(pvtr_unit, '(A)') '  </PRectilinearGrid>'
          write(pvtr_unit, '(A)') '</VTKFile>'
          close(pvtr_unit)
       end if

       deallocate(all_x_start, all_x_end, all_y_start, all_y_end)
    else
       call MPI_Gather(x_start, 1, MPI_INTEGER, dummy, 0, MPI_INTEGER, 0, mesh%comm, ierr)
       call MPI_Gather(x_end,   1, MPI_INTEGER, dummy, 0, MPI_INTEGER, 0, mesh%comm, ierr)
       call MPI_Gather(y_start, 1, MPI_INTEGER, dummy, 0, MPI_INTEGER, 0, mesh%comm, ierr)
       call MPI_Gather(y_end,   1, MPI_INTEGER, dummy, 0, MPI_INTEGER, 0, mesh%comm, ierr)
    end if
  end subroutine vtk_write_step

  !-----------------------------------------------------------------------------
  ! Calcul et enregistrement des métriques physiques temporelles dans l'historique CSV
  !-----------------------------------------------------------------------------
  subroutine vtk_write_diagnostics(this, mesh, fields, step, time, dt)
    class(vtk_writer_t), intent(in) :: this
    type(mesh_2d_t), intent(in)     :: mesh
    type(flow_fields_t), intent(in) :: fields
    integer, intent(in)             :: step
    real(dp), intent(in)            :: time, dt

    real(dp) :: local_vol, local_y_integral, local_v_integral, local_ke
    real(dp) :: local_max_u, local_max_div, dvol, y_pos, c_val, ind_bubble
    real(dp) :: uc, vc, speed2, speed, div_val
    real(dp) :: global_vol, global_y_integral, global_v_integral, global_ke
    real(dp) :: global_max_u, global_max_div, yc, vc_mean
    integer  :: i, j

    local_vol        = 0.0_dp
    local_y_integral = 0.0_dp
    local_v_integral = 0.0_dp
    local_ke         = 0.0_dp
    local_max_u      = 0.0_dp
    local_max_div    = 0.0_dp

    dvol = mesh%dx * mesh%dy

    do j = 1, mesh%ny_local
       y_pos = y_center(mesh, j)
       do i = 1, mesh%nx_local
          c_val = fields%C(i, j)
          ind_bubble = 1.0_dp - c_val ! 1 si phase fluide 2 (bulle / goutte)

          local_vol        = local_vol + ind_bubble * dvol
          local_y_integral = local_y_integral + y_pos * ind_bubble * dvol

          uc = 0.5_dp * (fields%u(i, j) + fields%u(i + 1, j))
          vc = 0.5_dp * (fields%v(i, j) + fields%v(i, j + 1))

          local_v_integral = local_v_integral + vc * ind_bubble * dvol

          speed2 = uc * uc + vc * vc
          local_ke = local_ke + 0.5_dp * fields%rho(i, j) * speed2 * dvol

          speed = sqrt(speed2)
          if (speed > local_max_u) local_max_u = speed

          div_val = abs((fields%u(i + 1, j) - fields%u(i, j)) / mesh%dx + &
                        (fields%v(i, j + 1) - fields%v(i, j)) / mesh%dy)
          if (div_val > local_max_div) local_max_div = div_val
       end do
    end do

    global_vol        = global_sum(mesh, local_vol)
    global_y_integral = global_sum(mesh, local_y_integral)
    global_v_integral = global_sum(mesh, local_v_integral)
    global_ke         = global_sum(mesh, local_ke)
    global_max_u      = global_max(mesh, local_max_u)
    global_max_div    = global_max(mesh, local_max_div)

    if (global_vol > 1.0e-12_dp) then
       yc = global_y_integral / global_vol
       vc_mean = global_v_integral / global_vol
    else
       yc = 0.0_dp
       vc_mean = 0.0_dp
    end if

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
    class(vtk_writer_t), intent(inout) :: this
    if (this%is_open) then
       close(this%csv_unit)
       this%is_open = .false.
    end if
  end subroutine vtk_finalize

end module module_vtk_writer
