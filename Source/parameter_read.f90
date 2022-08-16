module parameter_read
  use parameters
  use fconfig
  implicit none
  private
  public harmonic_oscillator_read

contains

  subroutine harmonic_oscillator_read
    type(config)        :: conf

    call conf%read_file('harmonic_oscillator_input')

    ! read spatial parameters
    call conf%value_from_key("box_size", l)
    call conf%value_from_key("spatial_step_size", h)
    call conf%value_from_key('grid_type', grid_type)

    ! read potential parameters
    call conf%value_from_key('potential_type', potential_type)
    select case(potential_type)
        case('soft_coulomb')
           call conf%value_from_key('soft_core_eps',soft_core_eps)
        
        case('gaussian')
           call  conf%value_from_key('gauss_pot_amplitude',gauss_pot_amplitude)   
           call  conf%value_from_key('gauss_pot_exponent',gauss_pot_exponent)
           
        case default
           stop 'potential_type is not a valid key'
    end select

    ! read propagation parameters
    call conf%value_from_key('total_time', t_intv)
    call conf%value_from_key('time_step_size', dt)
    t_tot = int(t_intv/dt)
    call conf%value_from_key('time_field_on_percent', prcntg_on)
    t_on  = int(t_tot * dt * prcntg_on)

    ! read pulse parameters
    call conf%value_from_key('pulse_amp', E_0)
    call conf%value_from_key('pulse_freq', omega)
    call conf%value_from_key('pulse_phase', phase)
    call conf%value_from_key('pulse_type', pulse_name)

    ! read Lanczos parameters
    call conf%value_from_key('lancz_iterations', lancz_itnum)
    call conf%value_from_key('lancz_threshold', lanc_threshold)
    
    ! read banded matrix parameters
    call conf%value_from_key('band_num_sym_mat', band_num_sym_mat)

    ! read Chebyshev parameters
    call conf%value_from_key('chebyshev_terms', chebyshev_terms)
    call conf%value_from_key('chebyshev_threshold', chebyshev_threshold)

    ! read iterative parameters
    call conf%value_from_key('prop_method', prop_method)
    call conf%value_from_key('it_type', it_type)
    call conf%value_from_key('quad_pt', quad_pt)
    call conf%value_from_key('it_tolerance', it_tolerance)
    call conf%value_from_key('states', states)

    call make_datafile_name

    print *, '************************************'
    print *, 'Beginning computations!'
    print *, '************************************'
    
  end subroutine harmonic_oscillator_read


  !--------------------------------------------------------------------
  ! DESCRIPTION:
  !> makes the datafile name from the parameters passed into the
  !> program
  !>
  !> Data file name is formatted such that:
  !>
  !> method_name pulse_name box_size spatial_step dt total_time
  !>
  !> Values that are small are expected to be small are multiplied by
  !> a 1000. All values are formatted as ints
  !--------------------------------------------------------------------
  subroutine make_datafile_name
    character(len=64)     :: info
    character(len=32)     :: lanczos
    character (len=32)    :: sb, sbcnl

    write(info, 1001) h, l, dt, t_intv


    datafilename = trim(example_problem) // '_' // trim(pulse_name) // &
      '_' // trim(info)

    datafilename = trim(datafilename) // '.dat'

1001 format (G0.2, '_', G0.4, '_', G0.4, '_', G0.4)
1002 format ('_', I0.3, '_', G0.4)
1003 format ('_', I1,'_', I0.3, '_', G0.4)
1004 format ('_b',I1)
    
  end subroutine

  
end module parameter_read

