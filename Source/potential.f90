module potential
  use parameters
  implicit none

  interface 
    function potential_func(x) result(pot_x)
      real(8) :: x(:)
      real(8), allocatable :: pot_x(:)
    end function potential_func
  end interface

contains 

  ! function to enable selecting a potential type appropite for the problem
  function select_potential_type(potential_type) result(f_ptr)
    character(len=*), intent(in) :: potential_type
    procedure(potential_func), pointer :: f_ptr 

    f_ptr => null()

    select case (trim(potential_type))
      case("zero_potential")
        f_ptr => zero_potential
      case("soft_coulomb")
        f_ptr => soft_coulomb_potential
      case("gaussian")
         f_ptr => gaussian_potential
      case("gobbler")
         f_ptr => gobbler
      case("quad_absorb")
         f_ptr => quad_absorb
      case("log_absorb")
         f_ptr => log_absorb
      case("quart_absorb")
         f_ptr => quart_absorb
      case default
        print*,trim(potential_type)
        stop 'type of potential is not defined!'        
      end select

   end function select_potential_type

  ! zero potential (necessary since banded_matrices assumes a potential is used)
  function zero_potential(x) result(potential)
     real(8)              :: x(:)
     real(8), allocatable :: potential(:)

     allocate(potential, source=x)
     potential(:) = 0
  end function zero_potential
    
  ! soft coulomb potential
  function soft_coulomb_potential(x) result(potential)
     real(8) :: x(:)
     real(8),allocatable :: potential(:)
     
     allocate(potential, source=x)
     
     potential(:) = - 1.d0 / sqrt(x(:)*x(:) + soft_core_eps * soft_core_eps)
   ! potential(:) = 1.d0 - 1.d0/sqrt(x(:)*x(:) + soft_core_eps*soft_core_eps)
     
   end function soft_coulomb_potential
   
  ! gaussian potential
  function gaussian_potential(x) result(potential)
     real(8) :: x(:)
     real(8), allocatable :: potential(:)

     allocate(potential, source=x)

     potential(:) = - gauss_pot_amplitude * exp(- gauss_pot_exponent * x(:)*x(:))

  end function gaussian_potential

  ! absorbing potential (gobbler) from Kazansky J. Phys. B: At. Mol. Opt. Phys. 31 (1998)
  function gobbler(x) result(potential)
    real(8)                :: x(:)
    real(8), allocatable   :: potential(:)
    integer                :: i

    allocate(potential, source=x)

    do i = 1,size(x)
       if (x(i) < sqrt(r_0)) then
          potential(i) = 0d0
       else
          potential(i) = G_0*(x(i)*x(i)/r_0 - 1)**3
       end if
    end do
    
  end function gobbler

  ! quadratic absorber
  function quad_absorb(x) result(potential)
    real(8)                :: x(:)
    real(8), allocatable   :: potential(:)
    integer                :: i

    allocate(potential, source=x)

    do i = 1,size(x)
       if (x(i) < r_0) then
          potential(i) = 0d0
       else
          potential(i) = G_0*(x(i)-r_0)**2
       end if
    end do
    
  end function quad_absorb

  ! logarithmic absorber
  function log_absorb(x) result(potential)
    real(8)               :: x(:)
    real(8), allocatable  :: potential(:)
    integer               :: i

    allocate(potential, source=x)

    do i = 1,size(x)
       if (x(i) < r_0) then
          potential(i) = 0d0
       else
          potential(i) = G_0*log(cos((x(i)-r_0)/(r_max-r_0)))
       end if
    end do
    
  end function log_absorb

  ! quartic absorber
  function quart_absorb(x) result(potential)
    real(8)               :: x(:)
    real(8), allocatable  :: potential(:)
    integer               :: i

    allocate(potential, source=x)

    do i = 1,size(x)
       if (x(i) < r_0) then
          potential(i) = 0d0
       else
          potential(i) = G_0*5d0/2d0*((x(i)-r_0)/(r_max-r_0))**4
       end if
    end do

  end function quart_absorb
  
end module  
