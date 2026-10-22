!===============================================================================
! Module : module_precision
! Description : Définition standardisée de la précision arithmétique flottante
!               et entière pour l'ensemble du solveur CFD COVOF 2D (Fortran 2018).
! Conforme aux recommandations IEEE 754 et ISO_FORTRAN_ENV.
!===============================================================================
module module_precision
  ! Obligation stricte de déclarer tous les types de variables (évite le typage implicite)
  implicit none
  ! Encapsulation par défaut : aucune entité n'est visible de l'extérieur sans export
  private

  ! Précision double (Double Precision 64-bit IEEE 754) :
  ! - Au moins 15 chiffres significatifs en base 10 (mantisse de 53 bits)
  ! - Plage dynamique d'exposant de 10^-307 à 10^+307 (exposant de 11 bits)
  ! - Utilisée pour tous les calculs hydrodynamiques, de pression et de transport
  integer, parameter, public :: dp = selected_real_kind(15, 307)
  
  ! Précision simple (Single Precision 32-bit IEEE 754) :
  ! - Au moins 6 chiffres significatifs en base 10 (mantisse de 24 bits)
  ! - Plage dynamique d'exposant de 10^-37 à 10^+37 (exposant de 8 bits)
  ! - Utilisée optionnellement pour les exports VTK légers ou diagnostics rapides
  integer, parameter, public :: sp = selected_real_kind(6, 37)

  ! Entier standard 32 bits signé (-2*10^9 à +2*10^9) :
  ! - Utilisé pour l'indexation spatiale des mailles (i, j) et les compteurs de boucle
  integer, parameter, public :: i4 = selected_int_kind(9)
  ! Alias d'entier par défaut pour le maillage et les solveurs
  integer, parameter, public :: ip = i4
  ! Entier 64 bits signé (-9*10^18 à +9*10^18) pour les grands volumes de données
  integer, parameter, public :: i8 = selected_int_kind(18)

end module module_precision
