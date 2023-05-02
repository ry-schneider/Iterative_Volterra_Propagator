program model_ode
  use parameters, only: t_intv, dt, it_type, quad_pt, it_tolerance, it_cap, ii, quad_type, datafilename, soln_method, add
  use parameter_read
  use timing
  use Lagrange_weights
  implicit none

    complex(8)                     :: mp_psi, exact_solution, b, integral
    complex(8), allocatable        :: mp_inhomogeneity(:), mp_iterate(:), mp_phi(:), mp_dummy(:), mp_vpsi(:)
    complex(8), allocatable        :: A(:,:)
    complex(8), allocatable        :: RHS(:)
    complex(8), allocatable        :: X(:,:)
    complex(8), allocatable        :: iterative_difference(:), y(:)
    real(8), allocatable           :: dvr_weights(:), mp_points(:)
    real(8)                        :: mp_time, theta, mp_max_diff, diff_dummy, mp_error, max_mp_error, alpha, theta1
    real(8), allocatable           :: mp_weights(:,:), comp_wt(:,:)
    integer                        :: n, i, j, k, l, m, p, q, t, r, step_count, iteration_number, max_iter, info
    integer, allocatable           :: IPIV(:)
    integer, allocatable           :: IPIV2(:)
    logical                        :: mp_converged

    ! read in parameters
    call model_ode_read
    n = quad_pt

    print *, '*************************************'
    print *, 'Problem: Model ODE'
    print *, 'Solution method: ', soln_method
    if (soln_method == 'volterra_iteration') then
       print *, 'Iteration type: ', it_type
    end if
    print *, 'Propagation step size:', dt
    print *, 'Quadrature type: ', quad_type
    print *, 'Number of quadrature points:', n
    print *, '************************************'

    allocate(mp_inhomogeneity(1:n), mp_iterate(1:n), mp_phi(1:n), mp_dummy(1:n), mp_vpsi(1:n), iterative_difference(1:n), &
         A(1:n-1,1:n-1), RHS(1:n-1), X(1:n,1:n), y(1:n), dvr_weights(1:n), mp_points(1:n), mp_weights(1:n,1:n-1), &
         comp_wt(1:n,1:n-1), IPIV(1:n-1), IPIV2(1:n))

    ! initialize the problem
    mp_error = 0
    max_mp_error = 0
    mp_psi = cmplx(1,0,8)
    exact_solution = cmplx(1,0,8)
    step_count = 1
    max_iter = 0

    alpha = 0d0

    open(unit=81, file=datafilename)
    
    call begin_timing
    mp_time = 0

    do while (mp_time < t_intv)

       if (add == 1) then
          alpha = mp_time + 0.5d0*dt
       end if

       if (soln_method == 'dvr') then
          ! compute points and weights
          call lgr_weights(mp_time, mp_time+dt, mp_points, mp_weights, n-2, quad_type)

          dvr_weights(:) = 0
          do j = 2,n
             dvr_weights(:) = dvr_weights(:) + mp_weights(:,j-1)
          end do

          ! set up system to solve
          theta1 = alpha*mp_time
          do l = 2,n
             RHS(l-1) = sqrt(dvr_weights(l))*(mp_points(l)-alpha)*cmplx(cos(theta1), sin(theta1), 8)*mp_psi
          end do

          do k = 2,n
             do i = 2,k-1
                A(k-1,i-1) = ii*sqrt(dvr_weights(k)/dvr_weights(i))*(1d0/(mp_points(i) - mp_points(k)))
                do m = 1,n
                   if ((m /= k) .and. (m /= i)) then
                      A(k-1,i-1) = A(k-1,i-1)*(mp_points(k) - mp_points(m))/(mp_points(i) - mp_points(m))
                   end if
                end do
             end do

             A(k-1,k-1) = -(mp_points(k)-alpha)
             do p = 1,k-1
                A(k-1,k-1) = A(k-1,k-1) + ii*1d0/(mp_points(k) - mp_points(p))
             end do

             do p = k+1,n
                A(k-1,k-1) = A(k-1,k-1) + ii*1d0/(mp_points(k) - mp_points(p))
             end do

             do i = k+1,n
                A(k-1,i-1) = ii*sqrt(dvr_weights(k)/dvr_weights(i))*(1d0/(mp_points(i) - mp_points(k)))
                do m = 1,n
                   if ((m /= k) .and. (m /= i)) then
                      A(k-1,i-1) = A(k-1,i-1)*(mp_points(k) - mp_points(m))/(mp_points(i) - mp_points(m))
                   end if
                end do
             end do
             
          end do

          info = 1
          call zgesv(n-1, 1, A, n-1, IPIV, RHS, n-1, info)

          if (info /= 0) then
             print *, 'DVR system solve failed'
             stop
          end if

          mp_iterate(1) = mp_psi
          do q = 2,n
             theta = alpha*mp_points(q)
             mp_iterate(q) = cmplx(cos(-theta), sin(-theta), 8)*(RHS(q-1)/sqrt(dvr_weights(q)) &
                  + cmplx(cos(theta1), sin(theta1), 8)*mp_psi)
          end do

          mp_psi = mp_iterate(n)

       else if (soln_method == 'power_expansion') then
          ! find interpolation points
          call lgr_weights(mp_time, mp_time+dt, mp_points, mp_weights, n-2, quad_type)

          ! evaluate inhomogeneity at quadrature points
          mp_inhomogeneity(1) = mp_psi
          do k = 2,n
             theta = alpha*(mp_points(k) - mp_time)
             mp_inhomogeneity(k) = cmplx(cos(-theta), sin(-theta), 8)*mp_psi
          end do

          mp_converged = .FALSE.
          mp_iterate(:) = mp_inhomogeneity(:)
          iteration_number = 0

          ! iterate until converged
          do while ((.not. mp_converged) .and. (iteration_number < it_cap))
             mp_phi(:) = mp_iterate(:)

             ! set up a system to solve for power expansion
             do l = 2,n
                do i = 1,n
                   do j = 1,n
                      X(i,j) = mp_points(i)**(j-1)
                   end do

                   theta = alpha*(mp_points(l) - mp_points(i))
                   y(i) = cmplx(cos(-theta), sin(-theta), 8)*(mp_points(i) - alpha)*mp_phi(i)
  
                end do

                info = 1
                call zgesv(n, 1, X, n, IPIV2, y, n, info)

                if (info /= 0) then
                   print *, 'Power expansion system solve failed'
                   stop
                end if

                ! analytically evaluate integral of power series
                mp_iterate(l) = mp_inhomogeneity(l)
                do m = 1,n
                   mp_iterate(l) = mp_iterate(l) - ii*y(m)*(1d0/m)*(mp_points(l)**m - mp_time**m)
                end do
                
             end do

             ! compare new iterate with previous
             iterative_difference(:) = mp_iterate(:) - mp_phi(:)
             mp_max_diff = 0
             
             do  p=1,n
                diff_dummy = sqrt(iterative_difference(p) * conjg(iterative_difference(p)))
                if (mp_max_diff < diff_dummy) then
                   mp_max_diff = diff_dummy
                end if
             end do
        
             if(mp_max_diff  <= it_tolerance) then
                mp_converged = .TRUE.
             end if

             iteration_number = iteration_number + 1

          end do

          if (max_iter < iteration_number) then
             max_iter = iteration_number
          end if

          mp_psi = mp_iterate(n)

       else
          ! compute quadrature points and weights
          call lgr_weights(mp_time, mp_time+dt, mp_points, mp_weights, n-2, quad_type)

          comp_wt(:,1) = mp_weights(:,1)
          do i = 2,n-1
             comp_wt(:,i) = comp_wt(:,i-1) + mp_weights(:,i)
          end do
 
          ! compute inhomogeneity at the quadrature points
          mp_inhomogeneity(1) = mp_psi
          do j = 2,n
             theta = alpha*(mp_points(j) - mp_time)
             mp_inhomogeneity(j) = cmplx(cos(-theta), sin(-theta), 8)*mp_psi
          end do

          mp_converged = .FALSE.
          mp_iterate(:) = mp_inhomogeneity(:)
          iteration_number = 0

          ! iterate starting with the inhomogeneity until convergence is reached
          do while ((.not. mp_converged) .and. (iteration_number < it_cap))
             mp_phi(:) = mp_iterate(:)
             if (it_type == 'gauss_seidel') then
                do k = 2,n
                
                   do l = 1,n
                      mp_vpsi(l) = (mp_points(l) - alpha)*mp_iterate(l)
                   end do

                   b = mp_inhomogeneity(k)
                
                   do m = 1,k-1
                      theta = alpha*(mp_points(k) - mp_points(m))
                      b = b - ii*comp_wt(m,k-1)*cmplx(cos(-theta), sin(-theta), 8)*mp_vpsi(m)
                   end do

                   do m = k+1,n
                      theta = alpha*(mp_points(k) - mp_points(m))
                      b = b - ii*comp_wt(m,k-1)*cmplx(cos(-theta), sin(-theta), 8)*mp_vpsi(m)
                   end do

                   mp_iterate(k) = b/(1 + ii*comp_wt(k,k-1)*(mp_points(k) - alpha))
               
                end do

             else if (it_type == 'jacobi') then
                do l = 1,n
                   mp_vpsi(l) = (mp_points(l) - alpha)*mp_iterate(l)
                end do

                do k = 2,n
                   mp_iterate(k) = mp_inhomogeneity(k)
                   
                   do m = 1,n
                      theta = alpha*(mp_points(k) - mp_points(m))
                      mp_iterate(k) = mp_iterate(k) - ii*comp_wt(m,k-1)*cmplx(cos(-theta), sin(-theta), 8)*mp_vpsi(m)
                   end do
                   
                end do

             end if
                
             ! compare new iterate with previous
             iterative_difference(:) = mp_iterate(:) - mp_phi(:)
             mp_max_diff = 0
             
             do  p=1,n
                diff_dummy = sqrt(iterative_difference(p) * conjg(iterative_difference(p)))
                if (mp_max_diff < diff_dummy) then
                   mp_max_diff = diff_dummy
                end if
             end do

             if(mp_max_diff  <= it_tolerance) then
                mp_converged = .TRUE.
             end if

             iteration_number = iteration_number + 1
          
          end do

          if (max_iter < iteration_number) then
             max_iter = iteration_number
          end if

          mp_psi = mp_iterate(n)
          
       end if

       do r = 1,n
          exact_solution = cmplx(cos(-0.50d0*mp_points(r)*mp_points(r)), sin(-0.50d0*mp_points(r)*mp_points(r)), 8)
          mp_error = abs(mp_iterate(r) - exact_solution)
          
          if (max_mp_error < mp_error) then
             max_mp_error = mp_error
          end if
       end do

       write(81,*) mp_time+dt, mp_psi, exact_solution
       
       mp_time = mp_time + dt
       step_count = step_count + 1
        
    end do

    close(81)

    print *, 'Maximum solution error (at any quadrature point):', max_mp_error
    if ((soln_method == 'volterra_iteration') .or. (soln_method == 'power_expansion')) then
       print *, 'Maximum number of iterations:', max_iter
    end if
    print *, '**************************************'
    
    call stop_timing
    call print_timing

end program model_ode

