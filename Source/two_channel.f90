program two_channel
  use parameters, only: t_intv, dt, it_type, quad_pt, it_tolerance, it_cap, ii, quad_type
  use parameter_read
  use timing
  use Lagrange_weights
  implicit none

    real(8), allocatable              :: trig_ih(:,:), trig_iterate(:,:), trig_phi(:,:), soln_error(:,:)
    real(8), dimension(2)             :: difference
    real(8), allocatable              :: points(:), weight_sum(:)
    real(8), allocatable              :: weights(:,:)
    real(8), dimension(2)             :: b
    real(8), dimension(2,2)           :: A
    real(8)                           :: time, norm, max_diff
    integer                           :: n,i,j,k,l,m,p,r,info,iterations,max_iter
    logical                           :: converged
    integer, dimension(2,2)           :: IPIV

    ! read in parameters
    call two_channel_read
    n = quad_pt

    print *, 'Problem: Two-channel'
    print *, 'Iteration type: ', it_type
    print *, 'Quadrature type: ', quad_type
    print *, 'Number of quadrature points:', n
    print *, '************************************'
    
    allocate(trig_ih(1:2,1:n), trig_iterate(1:2,1:n), trig_phi(1:2,1:n), soln_error(1:2,1:n), &
         points(1:n), weight_sum(1:n), weights(1:n,1:n-1))

    ! initialize the problem
    iterations = 0
    max_iter = 0

    call begin_timing
    time = 0

    do while (time < t_intv)
       ! compute quadrature points and weights
       call lgr_weights(time, time+dt, points, weights, n-2, quad_type)

       ! compute inhomogeneity at the quadrature points
       do i=1,n
          trig_ih(1,i) = (-cos(points(i))*cos(points(i)) - 0.5d0*sin(points(i) - 1))*points(i) + 2*cos(points(i))&
               + sin(points(i))*cos(points(i)) - 0.25d0*cos(1 + points(i)) + 0.25d0*cos(1 - points(i)) - 1
          trig_ih(2,i) = sin(points(i)) - points(i)
       end do
  
       converged = .FALSE.
       trig_iterate(:,1:n) = trig_ih(:,1:n)

       do while ((.not. converged) .and. (iterations < it_cap))
          iterations = iterations + 1
          
          ! make a copy of iterate to compate
          trig_phi(:,1:n) = trig_iterate(:,1:n)
          weight_sum(:) = 0

          if (it_type == 'gauss_seidel') then

             ! set up and solve system of equations at each point
             do j = 2,n

                weight_sum(:) = weight_sum(:) + weights(:,j-1)

                A(1,1) = 1 - weight_sum(j)*sin(-1.d0)
                A(1,2) = weight_sum(j)*(points(j)*cos(points(j)) - 1)
                A(2,1) = -weight_sum(j)
                A(2,2) = 1

                b(:) = trig_ih(:,j)

                do l = 1,j-1
                   b(1) = b(1) + weight_sum(l)*(sin((points(j)-points(l))-1)*trig_iterate(1,l)&
                        + (1 - points(l)*cos(points(j)))*trig_iterate(2,l))
                   b(2) = b(2) + weight_sum(l)*(trig_iterate(1,l) + (points(j)-points(l))*trig_iterate(2,l))
                end do

                do m = j+1,n
                   b(1) = b(1) + weight_sum(m)*(sin((points(j)-points(m))-1)*trig_phi(1,m)&
                        + (1 - points(m)*cos(points(j)))*trig_phi(2,m))
                   b(2) = b(2) + weight_sum(m)*(trig_phi(1,m) + (points(j)-points(m))*trig_phi(2,m))
                end do

                call dgesv(2, 1, A, 2, IPIV, b, 2, info)

                if (info /= 0) then
                   print *, 'Gauss-Seidel system solve failed'
                   stop
                end if

                trig_iterate(:,j) = b(:)
             
             end do

          else if (it_type == 'jacobi') then
             do j = 2,n

                weight_sum(:) = weight_sum(:) + weights(:,j-1)

                trig_iterate(:,j) = trig_ih(:,j)

                do l = 1,n
                   trig_iterate(1,j) = trig_iterate(1,j) + weight_sum(l)*(sin((points(j)-points(l))-1)*trig_phi(1,l)&
                        + (1 - points(l)*cos(points(j)))*trig_phi(2,l))
                   trig_iterate(2,j) = trig_iterate(2,j) + weight_sum(l)*(trig_phi(1,l) + (points(j)-points(l))*trig_phi(2,l))
                end do

             end do
          end if

          ! test convergence
          max_diff = 0
          do p=1,n
             difference = trig_iterate(:,p) - trig_phi(:,p)
             norm = sqrt(dot_product(difference, difference))
             
             if (max_diff < norm) then
                max_diff = norm
             end if
          end do

          print *, 'Iteration:', iterations, 'Max Error:', max_diff

          if (max_diff <= it_tolerance) then
             converged = .TRUE.
          end if
          
       end do

       print *, '************************************'
       print *, 'Error at each quadrature point: (point, first state error, second state error)'

       ! compute error at each point
       do r=1,n
          soln_error(1,r) = abs(cos(points(r)) - trig_iterate(1,r))
          soln_error(2,r) = abs(sin(points(r)) - trig_iterate(2,r))

          print *, points(r), soln_error(:,r)
       end do

       time = time + dt

       if (max_iter < iterations) then
          max_iter = iterations
       end if

       iterations = 0
       
    end do

    print *, '************************************'

    print *, 'Maximum number of iterations:', max_iter

    call stop_timing
    call print_timing

  
end program two_channel

