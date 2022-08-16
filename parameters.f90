!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!
!> MODULE: parameters
!
!> This module is used to load the parameters for the simulation, and
!> provide them to other modules through globals.
!---------------------------------------------------------------------

module parameters
  use fconfig

  implicit none

  !>***********************************************************************
  !> time variables
  !>***********************************************************************
  !> time 'box' size
  real(8)              :: t_intv
  !> size of time step 
  real(8)              :: dt
  !> total time to propagate the simulation
  integer              :: t_tot
  !> Number of steps pulse is on
  real(8)              :: t_on
  !> Percent of time the pulse is on
  real(8)              :: prcntg_on
  !> sqrt(-1) * dt
  complex(8)           :: idt ! sqrt(-1)*dt

  !>***********************************************************************
  !> spatial variables
  !>***********************************************************************
  !> size of spatial box
  real(8)              :: l
  !> size of steps in spatial box
  real(8)              :: h

  !>**********************************************************************
  !> matrix variables
  !>**********************************************************************
  !> soft core potential
  real(8)              :: soft_core_eps
  real(8)              :: gauss_pot_amplitude
  real(8)              :: gauss_pot_exponent

  !***********************************************************************
  !lanczos variables
  !***********************************************************************
  !> number of iterations of lanzos
  integer              :: lancz_itnum
  !> lanczos convergence threshold
  real(8)              :: lanc_threshold

  !***********************************************************************
  !chebyshev variables
  !***********************************************************************
  !> number of terms in chebyshev expansion
  integer              :: chebyshev_terms
  !> threshold for truncating chebyshev expansion
  real(8)              :: chebyshev_threshold
  
  !***********************************************************************
  !symetric banded matrices variable
  !***********************************************************************
  !> Number of bands of a symmetric banded matrix
  integer              :: band_num_sym_mat
  
  !***********************************************************************
  ! pulse variables
  !***********************************************************************
  !> pulse amplitude
  real(8)              :: E_0   
  !> frequency of pulse
  real(8)              :: omega 
  !> phase of pulse
  real(8)              :: phase
  !> type of pulse to simulate
  character(:), allocatable  :: pulse_name
  character(:), allocatable  :: potential_type
  character(:), allocatable  :: grid_type

  !**********************************************************************
  ! propagation variables
  !**********************************************************************
  !> propagation method
  character(:), allocatable  :: prop_method
  !> iteration type
  character(:), allocatable  :: it_type
  !> number of quadrature points
  integer                    :: quad_pt
  !> convergence criteria for iteration
  real(8)                    :: it_tolerance

  !**********************************************************************
  ! example problem parameters
  !**********************************************************************
  !> example problem to run
  character(:), allocatable  :: example_problem
  !> number of states in harmonic oscillator basis expansion
  integer                    :: states

  !**********************************************************************
  ! other complex variables and parameters
  !**********************************************************************
  !> 0 + sqrt(-1)
  complex*16, parameter      :: ii     = cmplx(0.d0,1.d0)
  !> zero
  complex*16, parameter      :: z_zero = cmplx(0.d0,0.d0)
  !> 1 + 0i
  complex*16, parameter      :: z_one  = cmplx(1.d0,0.d0)
  !> pi
  real*8, parameter          :: pi = 4*atan(1.0d0)
  !> name of data file to write results out to
  character(len=512)         :: datafilename
contains

  !-----------------------------------------------------------------------
  ! DESCRIPTION:
  !> Reads in the parameters from file. The file is either the default or
  !> passed in on the command line. If the file passed on the command line
  !> does not exist then the program terminates.
  !
  !> The loaded parameters are placed into the global variables of this
  !> module.
  !-----------------------------------------------------------------------
  subroutine parameter_reader
    type(config)        :: conf
    character(len=255)  :: conf_file_name, default_conf_file = 'default.conf'
    integer             :: ierr

    call GET_COMMAND_ARGUMENT(1, conf_file_name, status=ierr)

    if (ierr .gt. 0) then
      conf_file_name = default_conf_file
      print*, "Using default configuration file"
    end if

    call conf%read_file(conf_file_name)
    print*, "Using config file = ", trim(conf_file_name)


    !     read size of box on each side
    call conf%value_from_key("box_size", l)

    !     read spatial step size
    call conf%value_from_key("spatial_step_size", h)
    !     number of steps needed to fill the box from 0 to l (look in grid.f90)
    !      n = int(l/h)

    !>    total size of the matrix (number of steps to fill the box from -l to l
    !>     and NOT include point 0 in the middle) look in grid module
    !>      m_size = 2*n
    !>      ...| * | * | * | ... where * is the grid point and | corresponds to
    !>      the normal division points

    !     read grid type
    call conf%value_from_key('grid_type', grid_type)
    
    call conf%value_from_key('potential_type', potential_type)
    select case(potential_type)
        case('soft_coulomb')
    
    !>     soft core potential factor
           call conf%value_from_key('soft_core_eps',soft_core_eps)
        
        case('gaussian')
    !>     gaussian potential
           call  conf%value_from_key('gauss_pot_amplitude',gauss_pot_amplitude)   
           call  conf%value_from_key('gauss_pot_exponent',gauss_pot_exponent)
    case default
         stop 'potential_type is not a valid key'
    end select

    !>     total time in au that it takes for the propagation of the wave
    call conf%value_from_key('total_time', t_intv)

    !     time step size
    call conf%value_from_key('time_step_size', dt)

    !     total time steps needed to be taken
    t_tot = int(t_intv/dt)
    !     percentage of the time that the laser pulse is on

    call conf%value_from_key('time_field_on_percent', prcntg_on)

    !     time (in a.u.) for which the pulse is on (T in Eberly et. al.)
    t_on  = int(t_tot * dt * prcntg_on)
    !     amplitude of the laser puls
    call conf%value_from_key('pulse_amp', E_0)

    !     frequency of the laser pulse
    call conf%value_from_key('pulse_freq', omega)

    !     phase of the laser pulse
    call conf%value_from_key('pulse_phase', phase)

    !     shape of pulse, currently keywords smooth_pulse and square_pulse work
    call conf%value_from_key('pulse_type', pulse_name)

    !     number of lanczos iterations to be performed
    call conf%value_from_key('lancz_iterations', lancz_itnum)
    call conf%value_from_key('lancz_threshold', lanc_threshold)
    
    !     numbers of bands of a symetric banded matrix
    call conf%value_from_key('band_num_sym_mat', band_num_sym_mat)

    !    read in chebyshev parameters
    call conf%value_from_key('chebyshev_terms', chebyshev_terms)
    call conf%value_from_key('chebyshev_threshold', chebyshev_threshold)

    !    read in remaining propagation and example problem parameters
    call conf%value_from_key('prop_method', prop_method)
    call conf%value_from_key('it_type', it_type)
    call conf%value_from_key('quad_pt', quad_pt)
    call conf%value_from_key('it_tolerance', it_tolerance)
    call conf%value_from_key('example_problem', example_problem)
    call conf%value_from_key('states', states)

    call make_datafile_name

    print *, '************************************'
    print *, 'Beginning computations!'
    print *, '************************************'

    !Print relavant data
    ! print *, ''
    ! print* , potential_type
    ! print *, 'Spatial Box ranges from (a.u.)  ', -l, ' to ', l
    ! print *, 'Step size 		= ', h
    ! print *, 'Matrix size               = ', int(2 * real(l) / h) + 1
    ! print *, 'Propagation time in a.u. ', t_intv
    ! print *, 'Time step size 		= ', dt
    ! print *, 'Time pulse is on (a.u.) = ', t_on
    ! print *, 'Amplitude of pulse	= ', E_0
    ! print *, 'Frequency 		= ', omega
    ! print *, 'Phase                 = ', phase

    ! print *, 'type of pulse:', pulse_name
    ! print *, 'Saving results to: ', trim(datafilename)
    ! print *, ''
    !     initialize some parameters
    idt = 0.5d0*ii*dt


  end subroutine

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
  subroutine make_datafile_name()
    character(len=64)     :: info
    character(len=32)     :: lanczos
    character (len=32)    :: sb, sbcnl

    write(info, 1001) h, l, dt, t_intv


    datafilename = trim(example_problem) // '_' // trim(pulse_name) // &
      '_' // trim(info)

    ! if (index(trim(method_name), "lanczos") > 0) then
    !  write(lanczos, 1002) lancz_itnum, lanc_threshold

    !  datafilename = trim(datafilename) // trim(lanczos)
    ! end if

    ! if (index(trim(method_name), "spectrum_glncz") > 0 .or. &
    !    index(trim(method_name), "sb_lanczos") > 0 )    then
    !  write(sb, 1003)  band_num_sym_mat, lancz_itnum, lanc_threshold

    !  datafilename = trim(datafilename) // trim(sb)
    ! end if

    ! if (index(trim(method_name), "spectrum_sb_cnl") > 0 )  then
    !   write(sbcnl, 1004)  band_num_sym_mat

    !  datafilename = trim(datafilename) // trim(sbcnl)
    ! end if


    datafilename = trim(datafilename) // '.dat'

1001 format (G0.2, '_', G0.4, '_', G0.4, '_', G0.4)
1002 format ('_', I0.3, '_', G0.4)
1003 format ('_', I1,'_', I0.3, '_', G0.4)
1004 format ('_b',I1)
  end subroutine


end module
