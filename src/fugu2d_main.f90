!===============================================================================
! Programme Principal: fugu2d_main
! Description: Point d'entrée exécutable du code CFD FUGU 2D en Fortran 2018 + OpenMPI.
!              Exécute les simulations 2D monophasiques, diphasiques, thermiques
!              et les cas de validation selon les spécifications Thétis.
! Usage:
!   mpirun -np <N> ./bin/fugu2d_fortran cases/<cas>.json
!===============================================================================
program fugu2d_main
  use module_simulation
  implicit none
  include 'mpif.h'

  integer :: rank, size, ierr, num_args, stat
  character(len=256) :: config_file, arg
  type(simulation_t) :: sim
  logical :: success

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

  config_file = "cases/cavity.json"
  num_args = command_argument_count()

  if (num_args >= 1) then
     call get_command_argument(1, arg, status=stat)
     if (trim(arg) == "-h" .or. trim(arg) == "--help") then
        if (rank == 0) then
           write(*, '(A)') "Usage: mpirun -np <N> ./bin/fugu2d_fortran [fichier_configuration.json]"
           write(*, '(A)') "Exemples :"
           write(*, '(A)') "  mpirun -np 2 ./bin/fugu2d_fortran cases/cas3_1_laplacian.json"
           write(*, '(A)') "  mpirun -np 2 ./bin/fugu2d_fortran cases/cas4_5_cavity.json"
           write(*, '(A)') "  mpirun -np 4 ./bin/fugu2d_fortran cases/cas7_5_rising_bubble.json"
           write(*, '(A)') "  mpirun -np 2 ./bin/fugu2d_fortran cases/cas6_1_natural_convection.json"
        end if
        call MPI_Finalize(ierr)
        stop
     end if
     config_file = trim(arg)
  end if

  if (rank == 0) then
     write(*, '(A)') "================================================================================"
     write(*, '(A)') "         FUGU 2D - FORTRAN 2018 + OPENMPI NAVIER-STOKES CFD SOLVER              "
     write(*, '(A)') "================================================================================"
     write(*, '(A,I0)')       "  Nombre de processus MPI : ", size
     write(*, '(A,A)')        "  Fichier de cas JSON     : ", trim(config_file)
     write(*, '(A)') "================================================================================"
  end if

  success = sim%init(MPI_COMM_WORLD, trim(config_file))
  if (.not. success) then
     if (rank == 0) then
        write(*, '(A)') "ERREUR CRITIQUE : Echec de l'initialisation de FUGU 2D Fortran."
     end if
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
     stop 1
  end if

  call sim%run()

  call MPI_Finalize(ierr)

end program fugu2d_main
