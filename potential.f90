module potential
  use parameters, only: soft_core_eps,gauss_pot_amplitude, &
                        gauss_pot_exponent
implicit none


  interface 
    function potential_func(x) result(pot_x)
      real(8) :: x(:)
      real(8), allocatable :: pot_x(:)
    end function potential_func
  end interface

  contains 

  !> function to enable selecting a potential type appropite for the problem
  function select_potential_type(potential_type) result(f_ptr)
    character(len=*), intent(in) :: potential_type
    procedure(potential_func), pointer :: f_ptr 


    f_ptr => null()

      select case (trim(potential_type))
      case("soft_coulomb")
        f_ptr => soft_coulomb_potential
      case("gaussian")
        f_ptr => gaussian_potential
      case default
        print*,trim(potential_type)
        stop 'type of potential is not defined!'        
      end select

  end function select_potential_type

  function soft_coulomb_potential(x) result(potential)
     real(8) :: x(:)
     real(8),allocatable :: potential(:)
     
     allocate(potential, source=x)
     
     potential(:) = - 1.d0 / sqrt(x(:)*x(:) + soft_core_eps * soft_core_eps)

     ! potential(:) = 1.d0 - 1.d0/sqrt(x(:)*x(:) + soft_core_eps*soft_core_eps)
  end function

  function gaussian_potential(x) result(potential)
     real(8) :: x(:)
     real(8), allocatable :: potential(:)

     allocate(potential, source=x)

     potential(:) = - gauss_pot_amplitude * exp(- gauss_pot_exponent * x(:)*x(:))

  end function

end module  
