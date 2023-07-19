module parameter_read
  use parameters
  use fconfig
  implicit none
  private
  public harmonic_oscillator_read, two_level_atom_read, two_channel_read, model_ode_read, hydrogen_read

contains


  subroutine hydrogen_read
    type(config)        :: conf
    character(len=255)  :: conf_file_name
    integer             :: ierr

    call GET_COMMAND_ARGUMENT(1, conf_file_name, status=ierr)

    if (ierr .gt. 0) then
           conf_file_name = '../../Input/hydrogen_input.in'
    end if

    call conf%read_file(conf_file_name)

    ! read expansion parameters
    call conf%value_from_key('box_size', r_max)
    call conf%value_from_key('spatial_step_size', dr)
    r_size = ceiling(r_max/dr)
    call conf%value_from_key('grid_type', grid_type)
    call conf%value_from_key('grid_tol', grid_tol)
    call conf%value_from_key('l_max', l_max)
    call conf%value_from_key('coul_num', coul_num)
    call conf%value_from_key('e_min', e_min)
    call conf%value_from_key('dE', dE)

    ! read propagation parameters
    call conf%value_from_key('total_time', t_intv)
    call conf%value_from_key('time_step_size', dt)
    t_tot = int(t_intv/dt)

    ! read pulse parameters
    call conf%value_from_key('pulse_type', pulse_name)
    call conf%value_from_key('dressing_amplitude', E_d)
    call conf%value_from_key('probing_amplitude', E_p)
    call conf%value_from_key('dressing_frequency', omega_d)
    call conf%value_from_key('probing_frequency', omega_p)
    call conf%value_from_key('dressing_phase', phase_d)
    call conf%value_from_key('probing_phase', phase_p)
    call conf%value_from_key('dressing_peak', peak_d)
    call conf%value_from_key('probing_peak', peak_p)

    ! read absorber parameters
    call conf%value_from_key('absorber_type', potential_type)
    call conf%value_from_key('radial_cutoff', r_0)
    call conf%value_from_key('absorber_amp', G_0)

    ! read exponential parameters
    call conf%value_from_key('arnoldi_iterations', arnoldi_itnum)
    call conf%value_from_key('arnoldi_threshold', arnoldi_threshold)
    call conf%value_from_key('arnoldi_reortho', arnoldi_reortho)
    call conf%value_from_key('cn_gmres_max', cn_gmres_max)
    call conf%value_from_key('cn_gmres_tol', cn_gmres_tol)
    call conf%value_from_key('lancz_iterations', lancz_itnum)
    call conf%value_from_key('lancz_threshold', lanc_threshold)
    call conf%value_from_key('lancz_reortho', lancz_reortho)
    call conf%value_from_key('chebyshev_terms', chebyshev_terms)
    call conf%value_from_key('chebyshev_threshold',chebyshev_threshold)
    
    ! read banded matrix parameters
    call conf%value_from_key('band_num_sym_mat', band_num_sym_mat)

    ! read GMRES parameters
    call conf%value_from_key('gmres_tol', gmres_tol)
    call conf%value_from_key('gmres_max', gmres_max)

    ! read quadrature/iterative parameters
    call conf%value_from_key('prop_method', prop_method)
    call conf%value_from_key('it_type', it_type)
    call conf%value_from_key('quad_pt', quad_pt)
    call conf%value_from_key('quad_type', quad_type)
    call conf%value_from_key('it_tolerance', it_tolerance)
    call conf%value_from_key('it_cap', it_cap)

    ! read harmonic oscillator parameters
    call conf%value_from_key('example_problem', example_problem)
    call conf%value_from_key('soln_method', soln_method)
    call conf%value_from_key('itvolt_version', itvolt_version)

    call make_datafile_name

  end subroutine hydrogen_read

  
  subroutine harmonic_oscillator_read
    type(config)        :: conf
    character(len=255)  :: conf_file_name
    integer             :: ierr

    call GET_COMMAND_ARGUMENT(1, conf_file_name, status=ierr)

    if (ierr .gt. 0) then
           conf_file_name = '../../Input/harmonic_oscillator_input.in'
    end if

    call conf%read_file(conf_file_name)

    ! read spatial parameters
    call conf%value_from_key("box_size", r_max)
    call conf%value_from_key("spatial_step_size", dr)
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
    call conf%value_from_key('lancz_reortho', lancz_reortho)
    
    ! read banded matrix parameters
    call conf%value_from_key('band_num_sym_mat', band_num_sym_mat)

    ! read Chebyshev parameters
    call conf%value_from_key('chebyshev_terms', chebyshev_terms)
    call conf%value_from_key('chebyshev_threshold', chebyshev_threshold)

    ! read GMRES parameters
    call conf%value_from_key('gmres_tol', gmres_tol)
    call conf%value_from_key('gmres_max', gmres_max)

    ! read quadrature/iterative parameters
    call conf%value_from_key('prop_method', prop_method)
    call conf%value_from_key('exp_it_type', exp_it_type)
    call conf%value_from_key('it_type', it_type)
    call conf%value_from_key('quad_pt', quad_pt)
    call conf%value_from_key('quad_type', quad_type)
    call conf%value_from_key('it_tolerance', it_tolerance)
    call conf%value_from_key('it_cap', it_cap)

    ! read harmonic oscillator parameters
    call conf%value_from_key('states', states)
    call conf%value_from_key('example_problem', example_problem)
    call conf%value_from_key('soln_method', soln_method)

    call make_datafile_name
    
  end subroutine harmonic_oscillator_read


  subroutine two_level_atom_read
    type(config)        :: conf
    character(len=255)  :: conf_file_name
    integer             :: ierr

    call GET_COMMAND_ARGUMENT(1, conf_file_name, status=ierr)

    if (ierr .gt. 0) then
       conf_file_name = '../../Input/two_level_atom_input.in'
    end if

    call conf%read_file(conf_file_name)

    ! read propagation parameters
    call conf%value_from_key('prop_method', prop_method)
    call conf%value_from_key('total_time', t_intv)
    call conf%value_from_key('time_step_size', dt)
    t_tot = int(t_intv/dt)

    ! read pulse parameters
    call conf%value_from_key('pulse_type', pulse_name)
    call conf%value_from_key('pulse_freq', omega)
    call conf%value_from_key('pulse_phase', phase)
    call conf%value_from_key('pulse_amp', E_0)
    t_on = t_intv

    ! read quadrature/iterative parameters
    call conf%value_from_key('it_type', it_type)
    call conf%value_from_key('quad_pt', quad_pt)
    call conf%value_from_key('quad_type', quad_type)
    call conf%value_from_key('it_tolerance', it_tolerance)
    call conf%value_from_key('it_cap', it_cap)
    call conf%value_from_key('gmres_max', gmres_max)
    call conf%value_from_key('gmres_tol', gmres_tol)

    ! read example problem parameters
    call conf%value_from_key('example_problem', example_problem)

    call make_datafile_name

  end subroutine two_level_atom_read


  subroutine two_channel_read
    type(config)        :: conf
    character(len=255)  :: conf_file_name
    integer             :: ierr

    call GET_COMMAND_ARGUMENT(1, conf_file_name, status=ierr)

    if (ierr .gt. 0) then
       conf_file_name = '../../Input/two_channel_input.in'
    end if
    
    call conf%read_file(conf_file_name)

    ! read propagation parameters
    call conf%value_from_key('total_time', t_intv)
    call conf%value_from_key('time_step_size', dt)
    t_tot = int(t_intv/dt)

    ! read quadrature/iterative parameters
    call conf%value_from_key('it_type', it_type)
    call conf%value_from_key('quad_pt', quad_pt)
    call conf%value_from_key('quad_type', quad_type)
    call conf%value_from_key('it_tolerance', it_tolerance)
    call conf%value_from_key('it_cap', it_cap)

    ! read example problem parameters
    call conf%value_from_key('example_problem', example_problem)

    call make_datafile_name

  end subroutine two_channel_read


  subroutine model_ode_read
    type(config)        :: conf
    character(len=255)  :: conf_file_name
    integer             :: ierr

    call GET_COMMAND_ARGUMENT(1, conf_file_name, status=ierr)

    if (ierr .gt. 0) then
       conf_file_name = '../../Input/model_ode_input.in'
    end if

    call conf%read_file(conf_file_name)

    ! read propagation parameters
    call conf%value_from_key('total_time', t_intv)
    call conf%value_from_key('time_step_size', dt)
    t_tot = int(t_intv/dt)

    ! read quadrature/iterative parameters
    call conf%value_from_key('it_type', it_type)
    call conf%value_from_key('quad_pt', quad_pt)
    call conf%value_from_key('quad_type', quad_type)
    call conf%value_from_key('it_tolerance', it_tolerance)
    call conf%value_from_key('it_cap', it_cap)

    ! read example problem parameters
    call conf%value_from_key('example_problem', example_problem)
    call conf%value_from_key('soln_method', soln_method)
    call conf%value_from_key('add', add)

    call make_datafile_name
    
  end subroutine model_ode_read

  
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

    write(info, 1001) dt, t_intv


    datafilename = trim(example_problem) // '_' // trim(info)

    datafilename = trim(datafilename) // '.dat'

1001 format (G0.2, '_', G0.4, '_', G0.4, '_', G0.4)
1002 format ('_', I0.3, '_', G0.4)
1003 format ('_', I1,'_', I0.3, '_', G0.4)
1004 format ('_b',I1)
    
  end subroutine

  
end module parameter_read

