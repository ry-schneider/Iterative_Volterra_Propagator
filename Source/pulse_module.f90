!! This module defines various types of pulses at time T using
!! the parameters loaded globally by the parameters module.
!!
!!    procedure(pulse_at_t_func), pointer :: my_pulse
!!    my_pulse => select_pulse_type('smooth')
!!
!! This module also makes the assumption that there is exactly one pulse,
!! with the same parameters being pulled from the parameters module.
!!
module pulse_module
  use parameters, only: E_0, omega, phase, pi, t_intv, t_on, E_d, omega_d, phase_d, peak_d, E_p, omega_p, phase_p, peak_p
  implicit none

  real(8) :: envelope_time

  interface 
    function pulse_at_t_func(t) result(E_t)
      real(8) :: t
      real(8) :: E_t
    end function pulse_at_t_func
  end interface

contains

  !> Function to enable selecting a pulse type appropite for the problem
  function select_pulse_type(pulse_name) result(f_ptr)
    character(len=*), intent(in) :: pulse_name
    procedure(pulse_at_t_func), pointer :: f_ptr 

    envelope_time = t_intv !t_on

    f_ptr => null()

    select case (trim(pulse_name))
      case("square_pulse")
        f_ptr => square_length
      case("smooth_pulse")
         f_ptr => smooth_length
      case("linearly_polarized")
         f_ptr => linearly_polarized
    end select

    ! if we didn't find a function we liked in the switch, then undefined behavior
    if (.NOT. ASSOCIATED(f_ptr)) then
      print*, pulse_name
      stop "invalid pulse type"
   end if
   
  end function select_pulse_type

  function zero_pulse(t) result(E_t)
    real(8) :: E_t, t
    E_t = 0d0
  end function zero_pulse

  function square_length(t) result(E_t)
    real(8) :: E_t, t

    E_t = E_0 * sin(omega*t + phase)
  end function

  function smooth_length(t) result(E_t)
    real(8) :: E_t, t

    E_t = E_0 * cos(omega*t + phase) * &
     sin(pi * t / envelope_time) * &
     sin(pi * t / envelope_time)

  end function smooth_length

  function linearly_polarized(t) result(E_t)
    real(8) :: E_t, t

    ! linearly polarized pulse of Douguet et al. (with sin^2 envelope)
    E_t = E_d*sin(omega_d*(t-peak_d)+phase_d)*sin(pi*t/(2d0*peak_d))**2

    ! E_t = E_d*exp(-(t-peak_d)**2/(t_intv**2))*sin(omega_d*t+phase_d) + &
    !     E_p*exp(-(t-peak_p)**2/(t_intv**2))*sin(omega_p*t+phase_p)
    
  end function linearly_polarized

end module pulse_module
