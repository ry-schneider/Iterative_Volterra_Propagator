module parameters

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
  real(8)              :: r_max
  !> size of steps in spatial box
  real(8)              :: dr

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
  !> number of vectors to reorthogonalize against
  integer              :: lancz_reortho

  !***********************************************************************
  !chebyshev variables
  !***********************************************************************
  !> number of terms in chebyshev expansion
  integer              :: chebyshev_terms
  !> threshold for truncating chebyshev expansion
  real(8)              :: chebyshev_threshold

  !***********************************************************************
  !arnoldi variables
  !***********************************************************************
  !> number of arnoldi iterations
  integer              :: arnoldi_itnum
  !> arnoldi convergence threshold
  real(8)              :: arnoldi_threshold
  !> number of vectors to reorthogonalize against
  integer              :: arnoldi_reortho
  
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
  !> linearly polarized pulse variables (assume both have durations t_intv)
  real(8)       :: E_d ! dressing field amplitude
  real(8)       :: E_p ! probe field amplitude
  real(8)       :: omega_d ! dressing field frequency
  real(8)       :: omega_p ! probe field frequency
  real(8)       :: phase_d ! dressing field phase
  real(8)       :: phase_p ! probe field phase
  real(8)       :: peak_d ! time of dressing field peak
  real(8)       :: peak_p ! time of probing field peak

  !**********************************************************************
  ! complex absorber variables
  !**********************************************************************
  !> radial cut off value
  real(8)                     :: r_0
  !> absorber amplitude
  real(8)                     :: G_0

  !**********************************************************************
  ! propagation variables
  !**********************************************************************
  !> propagation method
  character(:), allocatable  :: prop_method
  !> iteration type for itvolt exponential
  character(:), allocatable  :: exp_it_type
  !> iteration type
  character(:), allocatable  :: it_type
  !> number of quadrature points
  integer                    :: quad_pt
  !> type of quadrature points
  character(:), allocatable  :: quad_type
  !> convergence criteria for iteration
  real(8)                    :: it_tolerance
  !> limit on the number of iterations allowed
  integer                    :: it_cap

  !**********************************************************************
  ! GMRES parameters
  !**********************************************************************
  !> convergence criteria
  real(8)                    :: gmres_tol
  !> maximum number of iterations
  integer                    :: gmres_max

  !**********************************************************************
  ! example problem parameters
  !**********************************************************************
  !> example problem type
  character(:), allocatable  :: example_problem
  !> number of states in harmonic oscillator basis expansion
  integer                    :: states
  !> number of samples for exponential comparison
  integer                    :: samples
  !> solution method for model ode
  character(:), allocatable  :: soln_method
  !> add/subtract midpoint in model ode?
  integer                    :: add
  !> number of angular momentum terms in hydrogen expansion
  integer                    :: l_max
  !> number of Coulomb functions to project on
  integer                    :: coul_num

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

end module
