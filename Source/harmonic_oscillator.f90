program harmonic_oscillator
  use parameters
  use parameter_read
  use general_utility
  use banded_matrices
  use timing
  use integral_method
  implicit none

    type(banded_sym_mat)        :: h_zero, v, test
    real(8)                     :: I1, I2, x_zero, alpha, p_zero
    complex(8)                  :: g_t, h_t
    complex(8), allocatable     :: psi(:), b(:)
    real(8), allocatable        :: eigenvectors(:,:)
    real(8), allocatable        :: eigenvalues(:), h_diagonal(:)
    real(8), allocatable        :: h_off_diagonal(:), off_diagonal_copy(:)
    real(8), allocatable        :: prop_off_diagonal(:,:), factorial(:)
    logical                     :: converged
    real(8)                     :: t, E, E_diff, norm_error, max_norm_error
    real(8), allocatable        :: pop_prob(:), prob_error(:), max_prob_error(:)
    real(8), allocatable        :: work(:)
    complex(8), allocatable     :: gs_diag(:)
    complex(8), allocatable     :: gs_offd1(:), gs_offd2(:)
    integer                     :: i, j, k, m, l, info
    integer                     :: iteration_count, max_iter, step_count

    call harmonic_oscillator_read
    m = states
    
    allocate(psi(1:m), b(1:m), eigenvectors(1:m,1:m), eigenvalues(1:m), h_diagonal(1:m), h_off_diagonal(1:m-1), &
         off_diagonal_copy(1:m-1), prop_off_diagonal(1:m-1,1:2), factorial(1:m-1), pop_prob(1:m), prob_error(1:m),&
         max_prob_error(1:m), work(1:2*m-2), gs_diag(1:m), gs_offd1(1:m-1), gs_offd2(1:m-1))

    alpha = 2*pi/t_intv
    max_iter = 0
    step_count = 0
    max_prob_error(:) = 0
    max_norm_error = 0


    factorial(:) = factorial_array(m-1)
    
    ! set up Hamiltonian, V, and initial wavefunction
    do i=1,m
       h_diagonal(i) = 0.5d0*( 2*(i-1) + 1)
    end do

    do j = 1,m-1
       h_off_diagonal(j) = sqrt(0.5d0*j)
    end do

    off_diagonal_copy(:) = h_off_diagonal(:)

    psi = 0
    psi(1) = 1

    prop_off_diagonal(:,1) = h_off_diagonal(:)
    call h_zero%initialize(m, 1, h_diagonal, prop_off_diagonal)
    
    v = h_zero
    v%diagonal(:) = 0
    h_zero%offdiagonal(:,:) = 0

    open(unit=53, file=datafilename)

    ! start propagation
    call begin_timing
    t = 0

    do while (t < t_intv)

       call iterative_loop(h_zero, t, psi, v, max_iter)
     
       step_count = step_count + 1

       ! compute analytic population probabilities at time t+dt
       if (omega == 0d0) then
          I1 = -E_0/4d0*(2*sin(t+dt) - sin((alpha-1)*(t+dt))/(alpha-1) - sin((alpha+1)*(t+dt))/(alpha+1))
         
          I2 = -E_0/4d0*(-2*cos(t+dt) + 2 + cos((alpha+1)*(t+dt))/(alpha+1) - 1/(alpha+1)&
               - cos((alpha-1)*(t+dt))/(alpha-1) + 1/(alpha-1))
         
       else if (omega == 1d0) then
          I1 = -E_0/8d0*(2*(t+dt) - 2*sin(alpha*(t+dt))/alpha + sin(2*(t+dt))&
               - sin((alpha-2)*(t+dt))/(alpha-2) - sin((alpha+2)*(t+dt))/(alpha+2))
         
          I2 = -E_0/8d0*(-cos(2*(t+dt)) + 1 + cos((alpha+2)*(t+dt))/(alpha+2) - 1/(alpha+2)&
               - cos((alpha-2)*(t+dt))/(alpha-2) + 1/(alpha-2))

       else
          print *, 'Analytic solution for choice of omega_zero not programmed.'
          stop
       end if

       x_zero = sin(t+dt)*I1 - cos(t+dt)*I2
       p_zero = cos(t+dt)*I1 + sin(t+dt)*I2

       g_t = 0.25d0*(x_zero + ii*p_zero)*(x_zero + ii*p_zero) - 0.5d0*x_zero*x_zero
       h_t = 0.5d0*(x_zero + ii*p_zero)
       
       pop_prob(1) = abs(zexp(g_t))*abs(zexp(g_t))
       do l = 2,m
          pop_prob(l) = (2d0**(l-1))/factorial(l-1)*pop_prob(1)*(abs(h_t)**(2*(l-1)))
       end do

       ! run error comparison
       do k = 1,m
          prob_error(k) = abs(pop_prob(k) - abs(psi(k))*abs(psi(k)))
   
          if (max_prob_error(k) < prob_error(k)) then
             max_prob_error(k) = prob_error(k)
          end if
       end do

       norm_error = abs(1 - dot_product(psi, psi))
       if (max_norm_error < norm_error) then
          max_norm_error = norm_error
       end if

       if (step_count == 10) then
          write(53,*) t+dt, prob_error(1), norm_error
          step_count = 0
       end if

       t = t+dt
 
    end do
    
    close(53)

    print *, 'Problem: driven harmonic oscillator'
    print *, 'Method: ', it_type, ' with ', prop_method
    print *, '************************************'
    print *, 'Worst case population probability errors:'
    print *, 'Ground state:', max_prob_error(1)
    print *, 'First excited state:', max_prob_error(2)
    print *, 'Second excited state:', max_prob_error(3)
    print *, 'Third excited state:', max_prob_error(4)
    print *, '************************************'

    print *, 'Maximum norm error:', max_norm_error
    print *, '************************************'

    if (it_type == 'jacobi' .or. it_type == 'gauss_seidel') then
       print *, 'Maximum number of iterations:', max_iter
       print *, '************************************'
    end if

    call stop_timing
    call print_timing

    
end program harmonic_oscillator

