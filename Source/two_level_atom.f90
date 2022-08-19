program two_level_atom
  use parameters, only: E_0, t_intv, dt, it_type, quad_pt, it_tolerance, it_cap, ii, pi, datafilename, quad_type
  use parameter_read
  use timing
  use integral_method
  use Lagrange_weights
  implicit none

    complex(8), dimension(2,2) :: H, A
    complex(8), dimension(2)   :: tl_psi, tl_integral, cn_solution, b, analytic_soln, soln_error, dummy
    complex(8), allocatable    :: tl_inhomogeneity(:,:), tl_iterate(:,:), tl_phi(:,:), model_soln(:,:), model_diff(:,:)
    real(8), allocatable       :: nc_points(:)
    real(8), allocatable       :: nc_weights(:,:), weights(:,:)
    real(8)                    :: time, cn_solution_error, max_soln_error_g, max_cn_error, theta, E_t, solution_error_g
    real(8)                    :: max_soln_error_e, solution_error_e, tolerance
    complex(8)                 :: tl_pulse, tl_pulse_midpoint
    integer                    :: n, i, j, k, l, m, r, s, p, q, u, info, count, max_count, it_max
    logical                    :: tl_converged
    integer, dimension(2,2)    :: IPIV

    ! read in parameters
    call two_level_atom_read
    n = quad_pt

    allocate(tl_inhomogeneity(1:2,1:n), tl_iterate(1:2,1:n), tl_phi(1:2,1:n), model_soln(1:2,1:n), model_diff(1:2,1:n), &
         nc_points(1:n), nc_weights(1:n,1:n-1), weights(1:n,1:n-1))
    
    ! initialize the problem
    H = 0
    tl_psi(1) = 1
    tl_psi(2) = 0
    cn_solution(1) = 1
    cn_solution(2) = 0
    max_soln_error_g = 0
    max_soln_error_e = 0
    max_cn_error = 0
    count = 0
    max_count = 0

    open(unit=73, file=datafilename)

    print *, '***Begin Time Loop***'
    call begin_timing
    time = 0

    do while (time < t_intv)
       ! compute pulse at midpoint and add to hamiltonian
       tl_pulse_midpoint = 0.5d0*E_0*sin(pi*(time + 0.5d0*dt)/t_intv)*sin(pi*(time + 0.5d0*dt)/t_intv)
       H(1,2) = tl_pulse_midpoint
       H(2,1) = tl_pulse_midpoint

       ! short-time approximation
       if (it_type == 'short_time') then
          theta = tl_pulse_midpoint*dt
          dummy(1) = cos(theta)*tl_psi(1) - ii*sin(theta)*tl_psi(2)
          dummy(2) = -ii*sin(theta)*tl_psi(1) + cos(theta)*tl_psi(2)
          
          tl_psi(:) = dummy(:)

       else
          ! compute quadrature points and weights
          call lgr_weights(time, time+dt, nc_points, nc_weights, n-2, quad_type)

          ! compute composite weights
          weights(:,1) = nc_weights(:,1)
          do q=2,n-1
             do s = 1,n
                weights(s,q) = weights(s,q-1) + nc_weights(s,q)
             end do
          end do

          ! compute inhomogeneity at the quadrature points
          tl_inhomogeneity(:,1) = tl_psi(:)
          do i=2,n
             theta = tl_pulse_midpoint*(nc_points(i) - nc_points(1))
             tl_inhomogeneity(1,i) = cos(theta)*tl_psi(1) - ii*sin(theta)*tl_psi(2)
             tl_inhomogeneity(2,i) = -ii*sin(theta)*tl_psi(1) + cos(theta)*tl_psi(2)
          end do

          ! iterate starting with the inhomogeneity
          tl_converged = .FALSE.
          tl_iterate(:,1:n) = tl_inhomogeneity(:,1:n)

          ! Jacobi iteration
          if (it_type == 'jacobi') then
             do while(.not. tl_converged)
                tl_phi(:,1:n) = tl_iterate(:,1:n)

                do u=2,n
                   b(:) = tl_inhomogeneity(:,u)
                
                   do p=1,n
                      theta = tl_pulse_midpoint*(nc_points(u) - nc_points(p))
                      E_t = 0.5d0*E_0*sin(pi*nc_points(p)/t_intv)*sin(pi*nc_points(p)/t_intv) - tl_pulse_midpoint
                      b(1) = b(1) - ii*weights(p,u-1)*(cos(theta)*E_t*tl_phi(2,p) - ii*sin(theta)*E_t*tl_phi(1,p))
                      b(2) = b(2) - ii*weights(p,u-1)*(cos(theta)*E_t*tl_phi(1,p) - ii*sin(theta)*E_t*tl_phi(2,p))
                   end do

                   tl_iterate(:,u) = b(:)
                end do

                tl_converged = convergence_test(tl_phi, tl_iterate, it_tolerance)
                count = count + 1

                if (count >= it_cap) then
                   tl_converged = .TRUE.
                end if
          
             end do

             tl_psi(:) = tl_iterate(:,n)
          

          ! Gauss-Seidel iteration
          else if (it_type == 'gauss_seidel') then
             do while(.not. tl_converged)
                tl_phi(:,1:n) = tl_iterate(:,1:n)
         
                do j=2,n
                   A(1,1) = 1
                   A(1,2) = ii*weights(j,j-1)*(0.5d0*E_0*sin(pi*nc_points(j)/t_intv)*sin(pi*nc_points(j)/t_intv)&
                        - tl_pulse_midpoint)
                   A(2,1) = A(1,2)
                   A(2,2) = A(1,1)
             
                   b(:) = tl_inhomogeneity(:,j)
             
                   do k=1,j-1
                      theta = tl_pulse_midpoint*(nc_points(j) - nc_points(k))
                      E_t = 0.5d0*E_0*sin(pi*nc_points(k)/t_intv)*sin(pi*nc_points(k)/t_intv) - tl_pulse_midpoint
                      b(1) = b(1) - ii*weights(k,j-1)*(cos(theta)*E_t*tl_iterate(2,k) - ii*sin(theta)*E_t*tl_iterate(1,k))
                      b(2) = b(2) - ii*weights(k,j-1)*(cos(theta)*E_t*tl_iterate(1,k) - ii*sin(theta)*E_t*tl_iterate(2,k))
                   end do

                   do l=j+1,n
                      theta = tl_pulse_midpoint*(nc_points(j) - nc_points(l))
                      E_t = 0.5d0*E_0*sin(pi*nc_points(l)/t_intv)*sin(pi*nc_points(l)/t_intv) - tl_pulse_midpoint
                      b(1) = b(1) - ii*weights(l,j-1)*(cos(theta)*E_t*tl_phi(2,l) - ii*sin(theta)*E_t*tl_phi(1,l))
                      b(2) = b(2) - ii*weights(l,j-1)*(cos(theta)*E_t*tl_phi(1,l) - ii*sin(theta)*E_t*tl_phi(2,l))
                   end do

                   call zgesv(2,1,A,2,IPIV,b,2,info)
                   tl_iterate(:,j) = b(:)
                
                end do

                tl_converged = convergence_test(tl_phi, tl_iterate, it_tolerance)
                count = count + 1

                if (count >= it_cap) then
                   tl_converged = .TRUE.
                end if
            
             end do
   
             tl_psi(:) = tl_iterate(:,n)

          end if

       end if

       ! compare converged result with analytic solution
       if (it_type == 'short_time') then
          analytic_soln(1) = cos(0.25d0*E_0*((time+dt) - (t_intv/(2*pi))*sin(2*pi*(time+dt)/t_intv)))
          analytic_soln(2) = -ii*sin(0.25d0*E_0*((time+dt) - (t_intv/(2*pi))*sin(2*pi*(time+dt)/t_intv)))

          solution_error_g = abs(abs(analytic_soln(1))*abs(analytic_soln(1)) - abs(tl_psi(1))*abs(tl_psi(1)))
          solution_error_e = abs(abs(analytic_soln(2))*abs(analytic_soln(2)) - abs(tl_psi(2))*abs(tl_psi(2)))

          if (max_soln_error_g < solution_error_g) then
             max_soln_error_g = solution_error_g
          end if

          if (max_soln_error_e < solution_error_e) then
             max_soln_error_e = solution_error_e
          end if

       else
          do r = 1,n
             analytic_soln(1) = cos(0.25d0*E_0*(nc_points(r) - (t_intv/(2*pi))*sin(2*pi*nc_points(r)/t_intv)))
             analytic_soln(2) = -ii*sin(0.25d0*E_0*(nc_points(r) - (t_intv/(2*pi))*sin(2*pi*nc_points(r)/t_intv)))

             solution_error_g = abs(abs(analytic_soln(1))*abs(analytic_soln(1)) - abs(tl_iterate(1,r))*abs(tl_iterate(1,r)))
             solution_error_e = abs(abs(analytic_soln(2))*abs(analytic_soln(2)) - abs(tl_iterate(2,r))*abs(tl_iterate(2,r)))
     
             if (max_soln_error_g < solution_error_g) then
                max_soln_error_g = solution_error_g
             end if

             if (max_soln_error_e < solution_error_e) then
                max_soln_error_e = solution_error_e
             end if
          end do
          
       end if

       write(73,*) time+dt, solution_error_g, solution_error_e

       if (max_count < count) then
          max_count = count
       end if
       
       ! reset for next step
       H(1,2) = 0
       H(2,1) = 0
       time = time + dt
       count = 0
       
    end do

    print *, 'Maximum ground state error:', max_soln_error_g
    print *, 'Maximum excited state error:', max_soln_error_e
    print *, 'Maximum number of iterations:', max_count
    print *, '**************************************************'

    close(73)

    call stop_timing
    call print_timing

end program two_level_atom

