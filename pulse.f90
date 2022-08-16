!> This module defines various types of pulses at time T using
!! the parameters loaded globally by the parameters module.
!! 
!! The functions are named like so <pulse_type>_<gauge>
!!
!! Where the available pulse types are `square` or `smooth`
!! and the gauge is either length or velocity.
!!
!! It is recommended to use the convience function `select_pulse_type`
!! to get your desired pulse, as it will get the required pulse function
!! automagically. See example below
!!
!!    procedure(pulse_at_t_func), pointer :: my_pulse
!!    my_pulse => select_pulse_type('smooth', 'length')
!!
!! This module also makes the assumption that there is exactly one pulse,
!! with the same parameters being pulled from the parameters module.
!!
module pulse
  use parameters, only: E_0, omega, phase, pi, t_on
  implicit none

  ! public square_length, square_velocity, smooth_length, smooth_velocity

  real(8) :: envelope_time

  interface 
    function pulse_at_t_func(t) result(E_t)
      real(8) :: t
      real(8) :: E_t
    end function pulse_at_t_func
  end interface

contains

  !> Function to enable selecting a pulse type appropite for the problem
  function select_pulse_type(pulse_name, gauge) result(f_ptr)
    character(len=*), intent(in) :: pulse_name, gauge
    procedure(pulse_at_t_func), pointer :: f_ptr 

    envelope_time = t_on

    f_ptr => null()

    if (trim(gauge) == 'length') then
      select case (trim(pulse_name))
      case("square_pulse")
        f_ptr => square_length
      case("smooth_pulse")
        f_ptr => smooth_length
      end select
    else if (trim(gauge) == 'velocity') then
      select case (trim(pulse_name)) 
      case("square_pulse")
        f_ptr => square_velocity
      case("smooth_pulse")
        f_ptr => smooth_velocity
      end select
    end if

    ! if we didn't find a function we liked in the switch, then undefined behavior
    if (.NOT. ASSOCIATED(f_ptr)) then
      print*, pulse_name, gauge
      stop "function2::select_pulse_type - Invalid pulse type or gauge requested"
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

  function square_velocity(t) result(E_t)
    real(8) :: E_t, t

    E_t = E_0 * &
      (cos(omega*t + phase) - cos(phase)) &
      / omega
  end function square_velocity


  function smooth_velocity(t) result(E_t)
    real(8) :: E_t, t

    E_t = -E_0/4.0 *(envelope_time *&
      cos(2.0* pi * t / envelope_time + omega*t)/(envelope_time*omega + 2.0*pi) + &
      envelope_time*                                                         &
      cos(2.0* pi * t / envelope_time - omega*t)/(envelope_time*omega - 2.0*pi) - &
      2.0*cos(omega*t)/omega-                                  &
      (8*pi*pi/(omega*(envelope_time**2*omega**2-4*pi**2)))) 
  end function smooth_velocity


end module pulse
