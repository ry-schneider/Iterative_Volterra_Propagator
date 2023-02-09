program two_level_atom
  use parameters
  use parameter_read
  use general_utility
  use banded_matrices
  use timing
!  use pulse_module
!  use propagator
!  use Lagrange_weights
  use integral_method
  implicit none

    type(banded_sym_mat)       :: h_zero, v
    complex(8), dimension(2)   :: tl_psi, analytic_soln, soln_error
    real(8), dimension(2)      :: h_diagonal
    real(8), dimension(1,1)    :: h_off_diagonal
    real(8)                    :: time
    real(8)                    :: max_soln_error_g, solution_error_g
    real(8)                    :: max_soln_error_e, solution_error_e
    integer                    :: n, max_iter

!    real(8)                    :: E_mid, theta
!    real(8), allocatable       :: pt(:), RWORK(:)
!    real(8), allocatable       :: wt(:,:), comp_wt(:,:)
!    complex(8), dimension(2)   :: v_1, v_2
!    complex(8), allocatable    :: A(:,:), eigs(:), VL(:,:), VR(:,:), WORK(:)
!    real(8)                    :: rho, rho_max
!    integer                    :: i, j, k, info

!    procedure(pulse_at_t_func), pointer :: pulse
!    procedure(propagator_func), pointer :: propagator
!    procedure(initialize_variables), pointer :: initialize

    ! read in parameters
    call two_level_atom_read
    n = quad_pt

    print *, 'Problem: Driven Two-Level Atom'
    print *, 'Iteration type: ', it_type
    print *, 'Pulse amplitude:', E_0
    print *, 'Total propagation time:', t_intv
    print *, 'Propagation step size:', dt
    print *, 'Quadrature type: ', quad_type
    print *, 'Number of quadrature points:', n
    print *, '************************************'

!    pulse => select_pulse_type(pulse_name)
!    call select_propagator_type(prop_method, propagator, initialize)
!    allocate(pt(1:quad_pt), wt(1:quad_pt,1:quad_pt-2), comp_wt(1:quad_pt,1:quad_pt-1), &
!         A(2*quad_pt-2,2*quad_pt-2), eigs(1:2*quad_pt-2), VL(1:2*quad_pt-2,1:2*quad_pt-2), &
!         VR(1:2*quad_pt-2,1:2*quad_pt-2), WORK(1:4*quad_pt-4), RWORK(1:4*quad_pt-4))
!    rho_max = 0

    ! initialize the problem
    tl_psi(1) = 1
    tl_psi(2) = 0
    max_soln_error_g = 0
    max_soln_error_e = 0
    max_iter = 0

    h_diagonal = 0
    h_off_diagonal = 0

    call h_zero%initialize(2, 1, h_diagonal, h_off_diagonal)
    v = h_zero
    v%offdiagonal(1,1) = 1

    open(unit=73, file=datafilename)

    ! start propagation
    call begin_timing
    time = 0

    do while (time < t_intv)
       if (it_type == 'gmres') then
          call linear_solve(h_zero, time, tl_psi, v, max_iter)
       else
          call iterative_loop(h_zero, time, tl_psi, v, max_iter)
       end if

       ! compare converged result with analytic solution
       analytic_soln(1) = cos(0.5d0*E_0*((time+dt) - (t_intv/(2*pi))*sin(2*pi*(time+dt)/t_intv)))
       analytic_soln(2) = -ii*sin(0.5d0*E_0*((time+dt) - (t_intv/(2*pi))*sin(2*pi*(time+dt)/t_intv)))

       solution_error_g = abs(abs(analytic_soln(1))*abs(analytic_soln(1)) - abs(tl_psi(1))*abs(tl_psi(1)))
       solution_error_e = abs(abs(analytic_soln(2))*abs(analytic_soln(2)) - abs(tl_psi(2))*abs(tl_psi(2)))

       if (max_soln_error_g < solution_error_g) then
          max_soln_error_g = solution_error_g
       end if

       if (max_soln_error_e < solution_error_e) then
          max_soln_error_e = solution_error_e
       end if

       write(73,*) time+dt, solution_error_g, solution_error_e
       
       time = time + dt

!       ! compute spectral radius
!       E_mid = pulse(time+0.5d0*dt)
!       h_zero%diagonal(:) = h_zero%diagonal(:) + E_mid*v%diagonal(:)
!       h_zero%offdiagonal(:,:) = h_zero%offdiagonal(:,:) + E_mid*v%offdiagonal(:,:)

!       call initialize(h_zero)

!       call lgr_weights(time, time+dt, pt, wt, n-2, quad_type)
       
!       comp_wt(:,1) = wt(:,1)
!       do i = 2,n-1
!          comp_wt(:,i) = comp_wt(:,i-1) + wt(:,i)
!       end do

!       do j = 1,n-1
!          theta = pulse(pt(j+1))-E_mid
!          v_1(1) = 0
!          v_1(2) = theta
!          v_2(1) = theta
!          v_2(2) = 0

!          do k = 1,j-1
!             A(2*k-1:2*k,2*j-1) = -ii*comp_wt(j+1,k)*propagator(h_zero, pt(k+1)-pt(j+1), v_1)
!             A(2*k-1:2*k,2*j) = -ii*comp_wt(j+1,k)*propagator(h_zero, pt(k+1)-pt(j+1), v_2) 
!          end do

!          A(2*j-1:2*j,2*j-1) = -ii*comp_wt(j+1,j)*v_1
!          A(2*j-1:2*j,2*j) = -ii*comp_wt(j+1,j)*v_2
          
!          do k = j+1,n-1
!             A(2*k-1:2*k,2*j-1) = -ii*comp_wt(j+1,k)*propagator(h_zero, pt(k+1)-pt(j+1), v_1)
!             A(2*k-1:2*k,2*j) = -ii*comp_wt(j+1,k)*propagator(h_zero, pt(k+1)-pt(j+1), v_2)
!          end do
!       end do

!       call zgeev('N', 'N', 2*n-2, A, 2*n-2, eigs, VL, 2*n-2, VR, 2*n-2, WORK, 4*n-4, RWORK, info)

!       if (info /= 0) then
!          print *, 'Eiegenvalue solve for spectral radius failed'
!          stop
!       end if

!       rho = maxval(abs(eigs))

!       if (rho_max < rho) then
!          rho_max = rho
!       end if
       
    end do

    print *, 'Results:'
    print *, 'Maximum ground state error:', max_soln_error_g
    print *, 'Maximum excited state error:', max_soln_error_e
    print *, 'Maximum number of iterations:', max_iter
    print *, '**************************************************'
!    print *, 'Spectral radius bound:', rho_max

    close(73)

    call stop_timing
    call print_timing

end program two_level_atom

