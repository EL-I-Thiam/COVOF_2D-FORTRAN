!===============================================================================
! Projet        : COVOF 2D (Convection & Volume-Of-Fluid 2D CFD Solver)
! Fichier       : covof2d_main.f90 (Point d'entrée principal Fortran 2018 / MPI)
! Description   : Programme principal pilotant l'initialisation de l'environnement
!                 OpenMPI, le traitement des arguments CLI, la lecture du cas
!                 JSON, l'orchestration du calcul Navier-Stokes / VOF / Énergie,
!                 et la libération finale des ressources système.
! Auteur        : Elhadji Ibrahima Thiam
! Standard      : Fortran 2018 / MPI-3
!===============================================================================
program covof2d_main
  ! Importation du module orchestrateur de la simulation CFD
  use module_simulation, only: simulation_t
  ! Typage explicite obligatoire pour toutes les variables locales
  implicit none
  ! Inclusion des définitions de constantes et prototypes MPI Fortran
  include 'mpif.h'

  ! --- Déclarations des Variables Locales ---
  integer            :: rank       ! Rang unique du processus courant (0 .. size-1)
  integer            :: size       ! Nombre total de processus alloués au calcul
  integer            :: ierr       ! Code d'erreur retourné par les appels MPI
  integer            :: num_args   ! Nombre d'arguments passés en ligne de commande
  integer            :: stat       ! Statut de lecture de l'argument CLI
  character(len=256) :: config_file! Chemin vers le fichier de configuration JSON
  character(len=256) :: arg        ! Tampon de lecture d'un argument textuel
  type(simulation_t) :: sim        ! Instance du type dérivé orchestrateur
  logical            :: success    ! Indicateur booléen de réussite de l'initialisation

  ! 1. Initialisation du runtime parallèle OpenMPI
  call MPI_Init(ierr)

  ! Récupération du rang du processus dans le communicateur global MPI_COMM_WORLD
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  ! Récupération de la taille du groupe de processus MPI
  call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

  ! Fichier de configuration par défaut si aucun argument n'est fourni
  config_file = "cases/cavity.json"

  ! Lecture du nombre d'arguments passés au binaire lors de l'appel mpirun
  num_args = command_argument_count()

  ! 2. Analyse des arguments passés en ligne de commande
  if (num_args >= 1) then
     ! Récupération du premier argument textuel
     call get_command_argument(1, arg, status=stat)

     ! Détection de la demande d'affichage d'aide (-h ou --help)
     if (trim(arg) == "-h" .or. trim(arg) == "--help") then
        ! Seul le processus maître (rang 0) écrit sur la sortie standard
        if (rank == 0) then
           write(*, '(A)') "================================================================================"
           write(*, '(A)') "   COVOF 2D — Solveur CFD Navier-Stokes Incompressible Parallèle (Fortran 2018)"
           write(*, '(A)') "================================================================================"
           write(*, '(A)') "Usage :"
           write(*, '(A)') "  mpirun -np <N> ./bin/covof2d_fortran [fichier_configuration.json]"
           write(*, '(A)') ""
           write(*, '(A)') "Exemples d'exécution :"
           write(*, '(A)') "  mpirun -np 2 ./bin/covof2d_fortran cases/cas3_1_laplacian.json"
           write(*, '(A)') "  mpirun -np 2 ./bin/covof2d_fortran cases/cas4_5_cavity.json"
           write(*, '(A)') "  mpirun -np 4 ./bin/covof2d_fortran cases/cas7_5_rising_bubble.json"
           write(*, '(A)') "  mpirun -np 2 ./bin/covof2d_fortran cases/cas6_1_natural_convection.json"
           write(*, '(A)') "================================================================================"
        end if
        ! Clôture propre du communicateur MPI avant terminaison
        call MPI_Finalize(ierr)
        stop
     end if

     ! Affectation du chemin de fichier spécifié
     config_file = trim(arg)
  end if

  ! Affichage du récapitulatif de lancement par le rang racine
  if (rank == 0) then
     write(*, '(A)') "================================================================================"
     write(*, '(A)') "   COVOF 2D — Solveur CFD Navier-Stokes Incompressible Parallèle (Fortran 2018)"
     write(*, '(A)') "================================================================================"
     write(*, '(A,I0)')       "  Nombre de processus MPI  : ", size
     write(*, '(A,A)')        "  Fichier de cas JSON      : ", trim(config_file)
     write(*, '(A)') "================================================================================"
  end if

  ! 3. Initialisation de la simulation CFD (lecture JSON, partitionnement MAC, allocations)
  success = sim%init(MPI_COMM_WORLD, trim(config_file))
  if (.not. success) then
     ! En cas d'erreur de parsing ou d'incohérence des données physiques
     if (rank == 0) then
        write(*, '(A)') "[ERREUR CRITIQUE] Échec de l'initialisation de COVOF 2D Fortran."
        write(*, '(A)') "                  Vérifiez le fichier de configuration : " // trim(config_file)
     end if
     ! Interruption immédiate et coordonnée de tous les rangs MPI
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
     stop 1
  end if

  ! 4. Exécution de la boucle temporelle principale de simulation
  call sim%run()

  ! 5. Libération des ressources MPI et terminaison normale
  call MPI_Finalize(ierr)

end program covof2d_main
