module system_solve
  use parameters
  use general_utility
  use banded_matrices
  use pulse_module
  use propagator
  use Lagrange_weights
  implicit none
  private
  public linear_solve


contains

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
    integer                                       :: irc(1:5), icntl(1:8), info(1:3), lwork, max_iter, it_count
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

    it_count = 0

    ! call GMRES driver 
10  call drive_zgmres(d*(n-1), d*(n-1), gmres_max, lwork, work, irc, icntl, cntl, info, rinfo)

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

       it_count = it_count + 1
       
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
    if (max_iter < it_count) then
       max_iter = it_count
    end if

  end subroutine linear_solve

end module system_solve
