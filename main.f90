program main
  use parameters
  use examples
  implicit none

  !reading inputs
  call parameter_reader

  !selet problem to run
  select case(trim(example_problem))

  case('harmonic_oscillator')
     call harmonic_oscillator
     
  case default
     print *, example_problem, ' is not programmed -try harmonic_oscillator.'
     stop
 
  end select

end program main
