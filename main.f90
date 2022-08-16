program main
  use parameters
  implicit none

  !reading inputs
  call parameter_reader

  !selet problem to run
  select case(trim(example_problem))
     
  case default
     print *, example_problem, ' is not programmed -try harmonic_oscillator.'
     stop
 
  end select

end program main
