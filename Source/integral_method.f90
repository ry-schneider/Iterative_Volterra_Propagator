 module integral_method
   use parameters
   use general_utility
   use banded_matrices
   use pulse_module
   use timing
   use propagator
   use grid, only: x
   use Lagrange_weights
   implicit none
   private
   public iterative_loop, convergence_test, linear_solve

   
contains


   ! runs Jacobi/Gauss Seidel iterations for the time step [t, t+dt]
   ! assumes V(t) = E(t)V for some fixed (time independent) symmetric banded matrix V
   subroutine iterative_loop(mat, t, psi, v, max_iter)
     type(banded_sym_mat)                     :: mat, v
     real(8)                                  :: t
     complex(8)                               :: psi(:)
     integer                                  :: max_iter 
     complex(8)                               :: inhomogeneity(size(psi),quad_pt)
     complex(8)                               :: iterative_ans(size(psi),quad_pt), phi(size(psi),quad_pt)
     complex(8)                               :: v_psi(size(psi),quad_pt), b(size(psi))
     complex(8)                               :: gs_diag(size(psi)), gs_off1(size(psi)-1), gs_off2(size(psi)-1)
     integer                                  :: IPIV(size(psi))
     complex(8)                               :: AB(3*mat%bsz+1, size(psi))
     complex(8)                               :: alpha
     integer                                  :: it_num
     real(8)                                  :: pt(quad_pt)
     real(8)                                  :: wt(quad_pt,quad_pt-1), comp_wt(quad_pt,quad_pt-1)
     logical                                  :: converged
     real(8)                                  :: E_diff
     integer                                  :: i, j, k, l, r, s, u, p, q, a, c, n, info

     procedure(pulse_at_t_func), pointer      :: pulse
     procedure(propagator_func), pointer      :: propagator
     procedure(initialize_variables), pointer :: initialize

     pulse => select_pulse_type(pulse_name)
     call select_propagator_type(prop_method, propagator, initialize)

     it_num = 0
     n = quad_pt

     mat%diagonal(:) = mat%diagonal(:) + pulse(t + 0.5d0*dt)*v%diagonal(:)
     mat%offdiagonal(:,:) = mat%offdiagonal(:,:) + pulse(t + 0.5d0*dt)*v%offdiagonal(:,:)

     call initialize(mat)
     
     if (it_type == 'short_time') then
        psi(:) = propagator(mat, dt, psi)
        
     else
        ! compute quadrature points and weights
        call lgr_weights(t, t+dt, pt, wt, n-2, quad_type)

        comp_wt(:,1) = wt(:,1)
        do i = 2,n-1
           comp_wt(:,i) = comp_wt(:,i-1) + wt(:,i)
        end do

        ! evaluate inhomogeneity at the quadrature points
        inhomogeneity(:,1) = psi(:)
        do j = 2,n
           inhomogeneity(:,j) = propagator(mat, pt(j)-pt(1), psi)
        end do

        converged = .FALSE.
        iterative_ans(:,:) = inhomogeneity(:,:)

        ! iterate until converged
        do while (.not. (converged .or. it_num >= it_cap))
           phi(:,:) = iterative_ans(:,:)

           ! jacobi iteration
           if (it_type == 'jacobi') then
              do k = 1,n
                 E_diff = pulse(pt(k)) - pulse(t+0.5d0*dt)
                 v_psi(:,k) = E_diff*sb_matvec_mul(v, iterative_ans(:,k))
              end do

              do l = 2,n
                 b(:) = inhomogeneity(:,l)

                 do r = 1,n
                    b(:) = b(:) - ii*comp_wt(r,l-1)*propagator(mat, pt(l)-pt(r), v_psi(:,r))
                 end do

                 iterative_ans(:,l) = b(:)
              end do

           ! gauss-seidel iteration
           else if (it_type == 'gauss_seidel') then
              do u = 2,n
                 do s = 1,n
                    E_diff = pulse(pt(s)) - pulse(t+0.5d0*dt)
                    v_psi(:,s) = E_diff*sb_matvec_mul(v, iterative_ans(:,s))
                 end do

                 b(:) = inhomogeneity(:,u)

                 do p = 1,u-1
                    b(:) = b(:) - ii*comp_wt(p,u-1)*propagator(mat, pt(u)-pt(p), v_psi(:,p))
                 end do
                 
                 do p = u+1,n
                    b(:) = b(:) - ii*comp_wt(p,u-1)*propagator(mat, pt(u)-pt(p), v_psi(:,p))
                 end do

                 ! tridiagonal system solve
                 if (mat%bsz == 1) then
                    gs_diag(:) = 1
                    alpha = ii*comp_wt(u,u-1)*(pulse(pt(u)) - pulse(t+0.5d0*dt))
                    
                    gs_diag(:) = gs_diag(:) + alpha*v%diagonal(:)
                    gs_off1(:) = alpha*v%offdiagonal(:,1)
                    gs_off2(:) = gs_off1(:)
                    
                    call zgtsv(mat%mat_size, 1, gs_off1, gs_diag, gs_off2, b, mat%mat_size, info)

                    iterative_ans(:,u) = b(:)

                 ! general system solve   
                 else
                    alpha = ii*comp_wt(u,u-1)*(pulse(pt(u))-pulse(t+0.5d0*dt))
                    AB = 0

                    ! construct the matrix AB for Lapack call
                    do a = 1,mat%bsz
                       do c = 1,a-1
                          AB(2*mat%bsz+1+c-a, a) = alpha*v%offdiagonal(c,a-c)
                       end do

                       AB(2*mat%bsz+1+c-a,a) = 1 + alpha*v%diagonal(a)

                       do c = a+1,a+mat%bsz
                          AB(2*mat%bsz+1+c-a, a) = alpha*v%offdiagonal(c,c-a)
                       end do
                    end do

                    do a = mat%bsz+1,mat%mat_size-mat%bsz-1
                       do c = a-mat%bsz,a-1
                          AB(2*mat%bsz+1+c-a,a) = alpha*v%offdiagonal(c,a-c)
                       end do

                       AB(2*mat%bsz+1+c-a,a) = 1 + alpha*v%diagonal(a)
                       
                       do c = a+1,a+mat%bsz
                          AB(2*mat%bsz+1+c-a,a) = alpha*v%offdiagonal(c,c-a)
                       end do
                    end do

                    do a = mat%mat_size-mat%bsz,mat%mat_size
                       do c = a-mat%bsz,a-1
                          AB(2*mat%bsz+1+c-a,a) = alpha*v%offdiagonal(c,a-c)
                       end do

                       AB(2*mat%bsz+1+c-a,a) = 1 + alpha*v%diagonal(a)

                       do c = a+1,mat%mat_size
                          AB(2*mat%bsz+1+c-a,a) = alpha*v%offdiagonal(c,c-a)
                       end do
                    end do
                       
                    call zgbsv(mat%mat_size, mat%bsz, mat%bsz, 1, AB, 3*mat%bsz+1, IPIV, b, mat%mat_size, info)
                    
                    iterative_ans(:,u) = b(:)
                    
                 end if

                 if (info /= 0) then
                    print *, 'system solve failed in iteration'
                    stop
                 end if
                 
              end do

           else
              print *, 'iteration type not programmed'
              stop
           end if

           it_num = it_num + 1
           converged = convergence_test(phi, iterative_ans, it_tolerance)
           
        end do

        psi(:) = iterative_ans(:,n)

        if (max_iter < it_num) then
           max_iter = it_num
        end if
        
     end if

     mat%diagonal(:) = mat%diagonal(:) - pulse(t + 0.5d0*dt)*v%diagonal(:)
     mat%offdiagonal(:,:) = mat%offdiagonal(:,:) - pulse(t + 0.5d0*dt)*v%offdiagonal(:,:)

   end subroutine iterative_loop
 

   ! compares two vector valued functions (i.e. wavefunctions) at a set of points 
   function convergence_test(wave_one, wave_two, tolerance) result(ans)
     complex(8), intent(in)           :: wave_one(:,:), wave_two(:,:)
     real(8), intent(in)              :: tolerance
     logical                          :: ans
     real(8)                          :: max_diff, dummy
     integer                          :: d, i, n
     complex(8), allocatable          :: difference(:)

     max_diff = 0.d0

     ans = .FALSE.

     d = size(wave_one(:,1))
     n = size(wave_one(1,:))

     allocate(difference(1:d))

     ! compare the two waves at each point, saving the maximum difference as you go
     do i = 1,n
        difference = wave_one(:,i) - wave_two(:,i)
        dummy = sqrt(dot_product(difference, difference))
        if (max_diff < dummy) then
           max_diff = dummy
        end if
     end do

     if (max_diff <= tolerance) then
        ans = .TRUE.
     end if

   end function convergence_test

   
   ! solves linear system on [t, t+dt] via GMRES
   ! assumes V(t) = E(t)V for some fixed (time independent) symmetric banded matrix V
   subroutine linear_solve(mat, t, psi, v, max_iter)
     type(banded_sym_mat)                          :: mat, v
     real(8)                                       :: t
     complex(8)                                    :: psi(:)
     complex(8)                                    :: RHS(size(psi)*(quad_pt-1)), v_psi(size(psi))
     complex(8)                                    :: soln(size(psi)*(quad_pt-1)), A_v(size(psi)*(quad_pt-1))
     complex(8)                                    :: v_soln(size(psi)*(quad_pt-1))
     real(8)                                       :: pt(quad_pt)
     real(8)                                       :: wt(quad_pt, quad_pt-1), comp_wt(quad_pt, quad_pt-1)
     real(8)                                       :: cntl(1:5), rinfo(1:2)
     complex(8), allocatable                       :: work(:)
     integer                                       :: irc(1:5), icntl(1:8), info(1:3), lwork, max_iter
     integer                                       :: n, d, i, j, k, p

     procedure(pulse_at_t_func), pointer           :: pulse
     procedure(propagator_func), pointer           :: propagator
     procedure(initialize_variables), pointer      :: initialize


     pulse => select_pulse_type(pulse_name)
     call select_propagator_type(prop_method, propagator, initialize)

     ! dimensions of the problem
     n = quad_pt
     d = size(psi)

     ! replace mat with the value of the time-dependent hamiltonian at the midpoint of [t, t+dt]
     mat%diagonal(:) = mat%diagonal(:) + pulse(t+0.5d0*dt)*v%diagonal(:)
     mat%offdiagonal(:,:) = mat%offdiagonal(:,:) + pulse(t+0.5d0*dt)*v%offdiagonal(:,:)

     ! initialize method for handling exponentials
     call initialize(mat)

     ! find quadrature points and weights
     call lgr_weights(t, t+dt, pt, wt, n-2, quad_type)

     comp_wt(:,1) = wt(:,1)
     do i = 2,n-1
        comp_wt(:,i) = comp_wt(:,i-1) + wt(:,i)
     end do

     ! construct right hand side vector
     v_psi = (pulse(pt(1))-pulse(t+0.5d0*dt))*sb_matvec_mul(v,psi)
     do j = 1,n-1
        RHS(1+(j-1)*d:j*d) = propagator(mat, pt(j+1)-pt(1), psi-ii*comp_wt(1,j)*v_psi)
     end do

     ! initialize GMRES parameters
     call init_zgmres(icntl, cntl)
     icntl(4) = 0
     icntl(6) = 1
     icntl(7) = d*(n-1)
     cntl(1) = gmres_tol

     ! construct input (i.e. RHS vector and initial guess)
     lwork = gmres_max*gmres_max + gmres_max*(d*(n-1)+5) + 5*d*(n-1) + 2
     allocate(work(1:lwork))

     do j = 1,n-1
        work(1+(j-1)*d:j*d) = propagator(mat, pt(j+1)-pt(1), psi)
     end do
     work(d*(n-1)+1:2*d*(n-1)) = RHS

     ! call GMRES driver 
10   call drive_zgmres(d*(n-1), d*(n-1), gmres_max, lwork, work, irc, icntl, cntl, info, rinfo)

     ! perform matrix/vector multiplication if necessary
     if (irc(1) == 1) then
        soln = work(irc(2):irc(2)+d*(n-1)-1)
        do k = 1,n-1
           v_soln(1+(k-1)*d:k*d) = (pulse(pt(k+1))-pulse(t+0.5d0*dt))*sb_matvec_mul(v,soln(1+(k-1)*d:k*d))
        end do
       
        do k = 1,n-1
           A_v(1+(k-1)*d:k*d) = soln(1+(k-1)*d:k*d) + ii*comp_wt(k+1,k)*v_soln(1+(k-1)*d:k*d)
          
           do j = 1,k-1
              A_v(1+(k-1)*d:k*d) = A_v(1+(k-1)*d:k*d) + ii*comp_wt(j+1,k)*propagator(mat, pt(k+1)-pt(j+1), v_soln(1+(j-1)*d:j*d))
           end do

           do j = k+1,n-1
              A_v(1+(k-1)*d:k*d) = A_v(1+(k-1)*d:k*d) + ii*comp_wt(j+1,k)*propagator(mat, pt(k+1)-pt(j+1), v_soln(1+(j-1)*d:j*d))
           end do
        end do
       
        work(irc(4):irc(4)+d*(n-1)-1) = A_v
       
        go to 10

     ! perform dot products if necessary
     else if (irc(1) == 4) then
        do j = 1,irc(5)
           work(irc(4)+(j-1)) = dot_product(work(irc(2)+(j-1)*d*(n-1):irc(2)+j*d*(n-1)-1), work(irc(3):irc(3)+d*(n-1)-1))
        end do
       
        go to 10
       
     else
        if (info(1) /= 0) then
           print *, 'GMRES failed: info(1) = ', info(1)
           stop
        end if
       
     end if

     psi = work(d*(n-2)+1:d*(n-1))
    
     ! reset for next step
     mat%diagonal(:) = mat%diagonal(:) - pulse(t+0.5d0*dt)*v%diagonal(:)
     mat%offdiagonal(:,:) = mat%offdiagonal(:,:) - pulse(t+0.5d0*dt)*v%offdiagonal(:,:)
     deallocate(work)

     ! track maximum number of iterations
     if (max_iter < info(2)) then
        max_iter = info(2)
     end if

   end subroutine linear_solve

   
end module integral_method

