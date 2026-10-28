!===============================================================================
! Module        : module_eos
! Description   : Équations d'état (EOS) et propriétés thermodynamiques
!                 pour la formulation All-Mach / Compressible diphasique.
!                 Référence : Fluids 2024, 9, 90 (El Ouafa et al.)
! Auteur        : Elhadji Ibrahima Thiam
! Standard      : Fortran 2018
!===============================================================================
module module_eos
  use module_precision, only: dp
  use module_constants, only: EOS_INCOMPRESSIBLE, EOS_LINEAR, EOS_IDEAL_GAS, EOS_NASG
  implicit none
  private

  public :: eos_parse_type, eos_compute_chi_T, eos_compute_density, &
            eos_compute_temperature, eos_compute_sound_speed

contains

  !-----------------------------------------------------------------------------
  ! Convertit une chaîne de caractères en entier d'énumération EOS
  !-----------------------------------------------------------------------------
  pure function eos_parse_type(str) result(t)
    character(len=*), intent(in) :: str
    integer :: t

    if (trim(str) == "LINEAR" .or. trim(str) == "ACOUSTIC") then
       t = EOS_LINEAR
    else if (trim(str) == "IDEAL_GAS" .or. trim(str) == "PERFECT_GAS") then
       t = EOS_IDEAL_GAS
    else if (trim(str) == "NASG" .or. trim(str) == "STIFFENED_GAS") then
       t = EOS_NASG
    else
       t = EOS_INCOMPRESSIBLE
    end if
  end function eos_parse_type

  !-----------------------------------------------------------------------------
  ! Calcule la compressibilité isotherme chi_T = (1/rho)*(drho/dp)_T [1/Pa]
  !-----------------------------------------------------------------------------
  pure function eos_compute_chi_T(eos_type, chi_T_param, p, rho, p_inf, b_nasg) result(chi_T)
    integer,  intent(in) :: eos_type
    real(dp), intent(in) :: chi_T_param, p, rho
    real(dp), intent(in), optional :: p_inf, b_nasg
    real(dp) :: chi_T
    real(dp) :: p_pos, rho_pos, p_tot, b_val

    select case (eos_type)
    case (EOS_INCOMPRESSIBLE)
       chi_T = 0.0_dp
    case (EOS_LINEAR)
       chi_T = chi_T_param
    case (EOS_IDEAL_GAS)
       p_pos = max(1.0_dp, p)
       chi_T = 1.0_dp / p_pos
    case (EOS_NASG)
       p_tot = p
       if (present(p_inf)) p_tot = p_tot + p_inf
       p_tot = max(1.0_dp, p_tot)
       rho_pos = max(1.0e-3_dp, rho)
       b_val = 0.0_dp
       if (present(b_nasg)) b_val = b_nasg
       chi_T = (1.0_dp - rho_pos * b_val) / (rho_pos * p_tot)
    case default
       chi_T = 0.0_dp
    end select
  end function eos_compute_chi_T

  !-----------------------------------------------------------------------------
  ! Calcule la masse volumique rho(p, T) [kg/m^3]
  !-----------------------------------------------------------------------------
  pure function eos_compute_density(eos_type, rho0, chi_T_param, p, T, r_gas, p_ref) result(rho)
    integer,  intent(in) :: eos_type
    real(dp), intent(in) :: rho0, chi_T_param, p, T
    real(dp), intent(in), optional :: r_gas, p_ref
    real(dp) :: rho
    real(dp) :: pref, rgas, T_pos

    pref = 1.01325e5_dp
    if (present(p_ref)) pref = p_ref
    rgas = 287.05_dp
    if (present(r_gas)) then
       if (r_gas > 0.0_dp) rgas = r_gas
    end if

    select case (eos_type)
    case (EOS_INCOMPRESSIBLE)
       rho = rho0
    case (EOS_LINEAR)
       rho = rho0 * (1.0_dp + chi_T_param * (p - pref))
    case (EOS_IDEAL_GAS)
       T_pos = max(1.0_dp, T)
       rho = p / (rgas * T_pos)
    case default
       rho = rho0
    end select
  end function eos_compute_density

  !-----------------------------------------------------------------------------
  ! Calcule la température à partir de l'énergie spécifique et vitesse [K]
  !-----------------------------------------------------------------------------
  pure function eos_compute_temperature(eos_type, e, v_sq, gamma, r_gas) result(T)
    integer,  intent(in) :: eos_type
    real(dp), intent(in) :: e, v_sq, gamma, r_gas
    real(dp) :: T
    real(dp) :: u, gam, rgas

    if (eos_type == EOS_INCOMPRESSIBLE .or. eos_type == EOS_LINEAR) then
       T = 300.0_dp
       return
    end if

    gam = 1.4_dp
    if (gamma > 1.0_dp) gam = gamma
    rgas = 287.05_dp
    if (r_gas > 0.0_dp) rgas = r_gas

    u = max(1.0_dp, e - 0.5_dp * v_sq)
    T = (gam - 1.0_dp) / rgas * u
  end function eos_compute_temperature

  !-----------------------------------------------------------------------------
  ! Calcule la vitesse du son c_s [m/s]
  !-----------------------------------------------------------------------------
  pure function eos_compute_sound_speed(eos_type, rho, chi_T_param, p, T, gamma, r_gas) result(cs)
    integer,  intent(in) :: eos_type
    real(dp), intent(in) :: rho, chi_T_param, p, T, gamma, r_gas
    real(dp) :: cs
    real(dp) :: gam, rgas

    select case (eos_type)
    case (EOS_INCOMPRESSIBLE)
       cs = 1.0e5_dp
    case (EOS_LINEAR)
       if (chi_T_param > 0.0_dp .and. rho > 0.0_dp) then
          cs = 1.0_dp / sqrt(rho * chi_T_param)
       else
          cs = 1500.0_dp
       end if
    case (EOS_IDEAL_GAS)
       gam = 1.4_dp
       if (gamma > 1.0_dp) gam = gamma
       rgas = 287.05_dp
       if (r_gas > 0.0_dp) rgas = r_gas
       cs = sqrt(gam * rgas * max(1.0_dp, T))
    case default
       cs = 1500.0_dp
    end select
  end function eos_compute_sound_speed

end module module_eos
