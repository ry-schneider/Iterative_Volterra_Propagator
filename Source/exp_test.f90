program exp_test
  use parameters
  use parameter_read
  use general_utility
  use banded_matrices
  use propagator
  use timing
  implicit none

    type(banded_sym_mat)        :: A
    real(8), allocatable        :: diag(:), offd(:,:), re(:), im(:)
    complex(8), allocatable     :: v(:), soln(:), cheby_soln(:), lancz_soln(:), iv_soln(:)
    real(8)                     :: cheby_error, lancz_error, iv_error, cpu_ti, cpu_tf
    real(8)                     :: max_cheby_error, max_lancz_error, max_iv_error
    real(8)                     :: tot_cheby_error, tot_lancz_error, tot_iv_error
    real(8)                     :: diag_time, cheby_time, lancz_time, iv_time, dx
    real(8), allocatable        :: x(:)
    integer                     :: j, m, k

    call exp_test_read
    m = states

    print *, 'Problem: Matrix Exponential Test'
    print *, 'Size:', m
    ! print *, 'Number of samples:', samples 

    ! max_cheby_error = 0
    ! max_lancz_error = 0
    ! max_iv_error = 0
    
    diag_time = 0
    cheby_time = 0
    lancz_time = 0
    iv_time = 0
    
    ! tot_cheby_error = 0
    ! tot_lancz_error = 0
    ! tot_iv_error = 0

    allocate(diag(1:m), offd(1:m-1,1:2), re(1:m), im(1:m), v(1:m), soln(1:m), cheby_soln(1:m), &
         lancz_soln(1:m), iv_soln(1:m), x(1:m))

    ! do j = 1,samples
       ! draw random (banded symmetric) matrix and vector
       ! call random_number(diag)
       ! call random_number(offd(:,1))
       
       ! if (A%init) then
       !    A%diagonal(:) = diag(:)
       !    A%offdiagonal(:,:) = offd(:,:)
       ! else   
       !    call A%initialize(m, 1, 10d0*diag, 10d0*offd)
       ! end if

       ! call random_number(re)
       ! call random_number(im)
       ! v(:) = re(:) + ii*im(:)

    ! construct grid on [-10, 10]
    dx = 20d0/(m-1)
    do k = 1,m
       x(k) = -10+(k-1)*dx
    end do

    ! construct time independent hamiltonian and starting vector
    ! offd(:,1) = -0.5d0/dx**2
    ! diag(:) = 1d0/dx**2 + 0.5d0*x(:)*x(:)
    
    offd(1:m-1,1) = -2d0/3d0/dx**2
    offd(1:m-2,2) = 1d0/24d0/dx**2
    diag(:) = 5d0/4d0/dx**2 + 0.5d0*x(:)*x(:)

    call A%initialize(m, 2, diag, offd)

    v(:) = 1d0/sqrt(sqrt(pi))*exp(-0.5d0*x(:)*x(:))

    ! compute exact answer via diagonalization
    call cpu_time(cpu_ti)
    
    call sb_eigensolve(A)
    soln = diag_prop(A, dt, v)

    print *, '************************************'
    print *, 'First six eigenvalues:'
    print *, A%eigenvalues(1:6)

    call cpu_time(cpu_tf)
    diag_time = diag_time + (cpu_tf - cpu_ti)

    ! compute exponential via chebyshev
    call cpu_time(cpu_ti)

    call sb_eigenvalues(A)
    cheby_soln = cheby_prop(A, dt, v)

    call cpu_time(cpu_tf)
    cheby_time = cheby_time + (cpu_tf - cpu_ti)

    ! compute exponential via lanczos
    call cpu_time(cpu_ti)
       
    lancz_soln = lanczos_prop(A, dt, v)

    call cpu_time(cpu_tf)
    lancz_time = lancz_time + (cpu_tf - cpu_ti)

    ! compute exponential via itvolt
    call cpu_time(cpu_ti)
       
    iv_soln = itvolt_exp(A, dt, v)

    call cpu_time(cpu_tf)
    iv_time = iv_time + (cpu_tf - cpu_ti)

    ! error comparison
    cheby_error = sqrt(dot_product(soln-cheby_soln, soln-cheby_soln))
    lancz_error = sqrt(dot_product(soln-lancz_soln, soln-lancz_soln))
    iv_error = sqrt(dot_product(soln-iv_soln, soln-iv_soln))

       ! tot_cheby_error = tot_cheby_error + cheby_error
       ! tot_lancz_error = tot_lancz_error + lancz_error
       ! tot_iv_error = tot_iv_error + iv_error

       ! if (max_cheby_error < cheby_error) then
       !    max_cheby_error = cheby_error
       ! end if

       ! if (max_lancz_error < lancz_error) then
       !    max_lancz_error = lancz_error
       ! end if

       ! if (max_iv_error < iv_error) then
       !    max_iv_error = iv_error
       ! end if
       
    ! end do

    print *, '************************************'
    print *, 'Error:'
    print *, 'Chebyshev:', cheby_error
    print *, 'Lanczos:', lancz_error
    print *, 'ITVOLT:', iv_error
    print *, '************************************'
    
    ! print *, 'Worst-case error:'
    ! print *, 'Chebyshev:', max_cheby_error
    ! print *, 'Lanczos:', max_lancz_error
    ! print *, 'ITVOLT:', max_iv_error
    ! print *, '************************************'
    ! print *, 'Average error:'
    ! print *, 'Chebyshev:', tot_cheby_error/samples
    ! print *, 'Lanczos:', tot_lancz_error/samples
    ! print *, 'ITVOLT:', tot_iv_error/samples
    ! print *, '************************************'
    
    print *, 'CPU time (in milliseconds):'
    print *, 'Diagonalization:', diag_time*1000
    print *, 'Chebyshev:', cheby_time*1000
    print *, 'Lanczos:', lancz_time*1000
    print *, 'ITVOLT:', iv_time*1000
    
end program exp_test

