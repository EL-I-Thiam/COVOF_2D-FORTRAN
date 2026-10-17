module module_precision
  implicit none
  private

  ! Double precision (64-bit IEEE 754): 15-17 significant digits, exponent range 10^-307 to 10^307
  integer, parameter, public :: dp = selected_real_kind(15, 307)
  
  ! Single precision (32-bit):
  integer, parameter, public :: sp = selected_real_kind(6, 37)

  ! Standard integer kinds
  integer, parameter, public :: i4 = selected_int_kind(9)
  integer, parameter, public :: ip = i4
  integer, parameter, public :: i8 = selected_int_kind(18)

end module module_precision
