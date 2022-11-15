module propagator
use parameters !, only : idt, ii, dt, t_intv, datafilename, z_one, z_zero, lancz_itnum, lanc_threshold, chebyshev_terms, chebyshev_threshold
use banded_matrices
use pulse_module
use timing
use general_utility
use Lagrange_weights
use grid, only: x
implicit none
private
public  propagator_func, initialize_variables, select_propagator_type

  interface 
     function propagator_func(mat, local_dt, psi) result(ans)
      use banded_matrices, only: banded_sym_mat
      type(banded_sym_mat), intent(in) :: mat
      real(8), intent(in)              :: local_dt
      complex(8),intent(in)            :: psi(:)
      complex(8)                       :: ans(size(psi))
    end function propagator_func
    
    subroutine initialize_variables(mat)
      use banded_matrices, only: banded_sym_mat
      type(banded_sym_mat)            :: mat
    end subroutine initialize_variables  
  end interface

 contains

  !Function to enable selecting a propagator type appropite for the problem
  subroutine select_propagator_type(propagator_name, f_ptr, init_ptr)
    character(len=*), intent(in) :: propagator_name
    procedure(propagator_func), pointer :: f_ptr 
    procedure(initialize_variables), pointer :: init_ptr

    f_ptr => null()
    init_ptr => null()

    select case (trim(propagator_name))
      
      case('diagonalization')
         init_ptr => sb_eigensolve
         f_ptr => diag_prop

      case('chebyshev')
         init_ptr => sb_eigenvalues
         f_ptr => cheby_prop

      case('lanczos')
         init_ptr => init_lanczos
         f_ptr => lanczos_prop

      case('itvolt')
         init_ptr => init_itvolt_exp
         f_ptr => itvolt_exp
     
      case default 
        print *, propagator_name
        stop 'above method is not programmed'
        
    end select

    ! if we didn't find a function we liked in the switch, then undefined behavior
    if (.NOT. ASSOCIATED(f_ptr)) then
      print *, propagator_name
      stop "invalid pulse type"
   end if
   
  end subroutine select_propagator_type


  !applies a matrix exponential to a vector via ITVOLT
  function itvolt_exp(mat, local_dt, v) result(ans)
    type(banded_sym_mat), intent(in)                :: mat
    real(8), intent(in)                             :: local_dt
    complex(8), intent(in)                          :: v(:)
    type(banded_sym_mat)                            :: W_j
    complex(8)                                      :: expH_0(size(v), quad_pt)
    complex(8)                                      :: inv_expH_0(size(v), quad_pt) 
    complex(8)                                      :: ans(size(v))
    real(8)                                         :: sign
    real(8)                                         :: pt(quad_pt)
    real(8)                                         :: wt(quad_pt, quad_pt-1), comp_wt(quad_pt, quad_pt-1)
    complex(8)                                      :: iterate(size(v), quad_pt), comp(size(v), quad_pt)
    complex(8)                                      :: inhomogeneity(size(v), quad_pt)
    complex(8)                                      :: integral(size(v))
    complex(8)                                      :: b(size(v)), gs_diag(size(v)), gs_off1(size(v)-1), gs_off2(size(v)-1)
    logical                                         :: converged
    real(8)                                         :: diff, max_diff
    integer                                         :: n, i, j, k, r, p, u, s, a, it_count, info

    ! check sign of local_dt
    if (local_dt < 0d0) then
       sign = -1d0
    else
       sign = 1d0
    end if

    ! select number of quadrature points
    n = quad_pt

    ! construct matrix minus its diagonal (i.e. H_j - H_0)
    call W_j%initialize(mat%mat_size, mat%bsz, sign*mat%diagonal, sign*mat%offdiagonal)
    W_j%diagonal(:) = 0d0

    ! find quadrature points and weights
    call lgr_weights(0d0, sign*local_dt, pt, wt, n-2, quad_type)
    
    ! compute the exponentials of the diagonal piece H_0
    do a = 1,n
       do s = 1,size(v)
          expH_0(s,a) = cmplx(cos(sign*mat%diagonal(s)*pt(a)), sin(sign*mat%diagonal(s)*pt(a)), 8)
       end do
    end do
    inv_expH_0(:,:) = conjg(expH_0(:,:))

    ! add up composite weights
    comp_wt(:,1) = wt(:,1)
    do i=2,n-1
       comp_wt(:,i) = comp_wt(:,i-1) + wt(:,i)
    end do

    ! iterate beginning with e^{-i*H_0*t}v 
    inhomogeneity(:,1) = v(:)
    do j=2,n       
       inhomogeneity(:,j) = inv_expH_0(:,j)*v(:)
    end do

    it_count = 0
    converged = .FALSE.
    iterate(:,:) = inhomogeneity(:,:)
    
    do while (.not. (converged .or. it_count > n))
       comp(:,:) = iterate(:,:)
       
       ! Jacobi iteration
       if (it_type == 'jacobi') then
          do r = 2,n
             integral(:) = 0

             do p = 1,n
                integral(:) = integral(:) + comp_wt(p,r-1)*expH_0(:,p)*sb_matvec_mul(W_j, comp(:,p)) 
             end do

             iterate(:,r) = inhomogeneity(:,r) - ii*inv_expH_0(:,r)*integral(:)
          end do

       ! Gauss-Seidel iteration
       else if (it_type == 'gauss_seidel') then
          do r = 2,n
             integral(:) = 0
             
             do p = 1,r-1
                integral(:) = integral(:) + comp_wt(p,r-1)*expH_0(:,p)*sb_matvec_mul(W_j, iterate(:,p))
             end do

             do p = r+1,n
                integral(:) = integral(:) + comp_wt(p,r-1)*expH_0(:,p)*sb_matvec_mul(W_j, iterate(:,p))
             end do

             b(:) = inhomogeneity(:,r) - ii*inv_expH_0(:,r)*integral(:)

             ! tridiagonal system solve
             if (mat%bsz == 1) then
                gs_diag(:) = 1
                gs_off1(:) = ii*comp_wt(r,r-1)*W_j%offdiagonal(:,1)
                gs_off2(:) = gs_off1(:)
                
                call zgtsv(mat%mat_size, 1, gs_off1, gs_diag, gs_off2, b, mat%mat_size, info)

                if (info /= 0) then
                   print *, 'exponential system solve failed'
                   stop
                end if

                iterate(:,r) = b(:)

             else
                print *, 'iterative system solve for band size not programmed'
             end if
             
          end do

       else
          print *, 'iteration type for exponential not programmed'
          stop
          
       end if

       it_count = it_count + 1

       ! check for convergence
       max_diff = 0
       do u = 2,n
          diff = sqrt(dot_product(iterate(:,u)-comp(:,u), iterate(:,u)-comp(:,u)))
          if (max_diff < diff) then
             max_diff = diff
          end if
       end do

       if (max_diff <= 1.d-10) then
          converged = .TRUE.
       end if
       
    end do

    ans(:) = iterate(:,n)
    
  end function itvolt_exp

  
  ! applies a matrix exponential to a vector via full diagonalization
  ! must initialize with a full eigensolve before using
  function diag_prop(mat, local_dt, v) result(ans)
    type(banded_sym_mat), intent(in)                :: mat
    real(8), intent(in)                             :: local_dt
    complex(8), intent(in)                          :: v(:)
    complex(8)                                      :: ans(size(v))
    integer                                         :: i

    ans = matmul(transpose(mat%eigenvectors), v)
    
    do i = 1,size(v)
       ans(i) = cmplx(cos(-local_dt*mat%eigenvalues(i)), sin(-local_dt*mat%eigenvalues(i)), 8)*ans(i)
    end do

    ans = matmul(mat%eigenvectors, ans)

  end function diag_prop


  ! applies a matrix exponential to a vector via the chebyshev propagator
  ! must initialize with an eigenvalue solve before using
  function cheby_prop(mat, local_dt, psi) result(ans)
    type(banded_sym_mat), intent(in)        :: mat
    complex(8), intent(in)                  :: psi(:)
    real(8), intent(in)                     :: local_dt
    type(banded_sym_mat)                    :: mat_norm
    complex(8)                              :: ans(size(psi)), phi_old(size(psi)), phi_new(size(psi)), phi_copy(size(psi))
    real(8)                                 :: eigenvalues(size(psi)), ones(size(psi))
    real(8)                                 :: off_diagonal(size(psi)-1)
    real(8)                                 :: eig(size(psi), size(psi))
    real(8)                                 :: work(2*size(psi)-2)
    integer                                 :: m, info, num_cheby, i, j
    real(8)                                 :: delta, e_min
    complex(8)                              :: coeff(chebyshev_terms)
    logical                                 :: truncated

    e_min = minval(mat%eigenvalues)
    delta = maxval(mat%eigenvalues) - e_min

    ! find expansion coefficients (truncate once below tolerance)
    coeff(:) = 0
    coeff(1) = zexp(-ii*(0.50d0*delta + e_min)*local_dt)*bessel_jn(0, 0.50d0*delta*local_dt)
    
    num_cheby = 1
    i = 2
    truncated = .FALSE.
    do while (.not. truncated)
       coeff(i) = 2.0d0*zexp(-ii*(0.50d0*delta + e_min)*local_dt)*bessel_jn(i-1, 0.50d0*delta*local_dt)
       if ((abs(coeff(i)) <= chebyshev_threshold) .or. (i == chebyshev_terms)) then
          num_cheby = i
          truncated = .TRUE.
       end if
       
       i = i+1
    end do
    
    ! construct normalized version of mat
    mat_norm = mat
    ones(:) = 1
    mat_norm%diagonal(:) = (2.0d0/delta)*(mat_norm%diagonal(:) - e_min*ones(:)) - ones(:)
    mat_norm%offdiagonal(:,1) = (2.0d0/delta)*mat_norm%offdiagonal(:,1)

    ! recursively apply expansion to input vector
    phi_old(:) = psi(:)
    phi_new(:) = -ii*sb_matvec_mul(mat_norm, phi_old)
    ans(:) = coeff(1)*phi_old(:) + coeff(2)*phi_new(:)
    
    if (num_cheby > 2) then
       do j = 3, num_cheby
          phi_copy(:) = phi_new(:)
          phi_new(:) = -2.0d0*ii*sb_matvec_mul(mat_norm, phi_new) + phi_old(:)
          phi_old(:) = phi_copy(:)

          ans(:) = ans(:) + coeff(j)*phi_new(:)
       end do
    end if
    
  end function cheby_prop


  ! applies a matrix exponential to a vector via Lanczos iteration
  function lanczos_prop(mat, local_dt, psi) result(ans)
    type(banded_sym_mat), intent(in)             :: mat
    real(8), intent(in)                          :: local_dt
    complex(8), intent(in)                       :: psi(:)
    complex(8)                                   :: ans(size(psi))
    complex(8)                                   :: Q(size(psi), lancz_itnum)
    real(8)                                      :: beta(lancz_itnum), alpha(lancz_itnum)
    complex(8)                                   :: w(size(psi)), phi(size(psi))
    real(8)                                      :: eigenvalues(lancz_itnum), eigenvectors(lancz_itnum, lancz_itnum)
    real(8)                                      :: off_diagonal(lancz_itnum), work(2*lancz_itnum-2)
    complex(8)                                   :: dummy(lancz_itnum)
    integer                                      :: k, info, i, m
    real(8)                                      :: error
    logical                                      :: converged

    ! if starting vector is zero, skip iteration and output zero  
    if (sqrt(dot_product(psi,psi)) == 0.d0) then
       ans(:) = 0
       
    else
       ! first vector is normalized psi
       Q(:,1) = psi(:)/sqrt(dot_product(psi, psi))

       ! compute second lanczos vector
       w(:) = sb_matvec_mul(mat, Q(:,1))
       alpha(1) = dot_product(Q(:,1), w)
       w(:) = w(:) - alpha(1)*Q(:,1)
       beta(1) = sqrt(dot_product(w,w))
       Q(:,2) = w(:)/beta(1)

       ! re-orthogonalize against first to ensure convergence
       call zschmab(Q(:,1), Q(:,2), 1.d-10, size(psi), 1, 1, m)
       
       ! project on first lanczos vector
       phi(:) = cmplx(cos(-local_dt*alpha(1)), sin(-local_dt*alpha(1)), 8)*dot_product(Q(:,1), psi)*Q(:,1)

       ! start iteration, adding an extra vector until convergence is reached 
       k = 2
       converged = .FALSE.
       do while ( .not. (converged .or. k > lancz_itnum))
          w(:) = sb_matvec_mul(mat, Q(:,k))
          alpha(k) = dot_product(Q(:,k), w)

          ! find eigenvalues/eigenvectors of tridiagonal matrix 
          eigenvalues(:) = alpha(:)
          off_diagonal(:) = beta(:)
          call dstev('V', k, eigenvalues(1:k), off_diagonal(1:k-1), eigenvectors(1:k,1:k), k, work(1:2*k-2), info)

          if (info /= 0) then
             print *, 'eigensolve in lanczos iteration failed at step', k-1
             stop
          end if

          ! use diagonalization to apply the exponential to psi
          dummy(1:k) = matmul(transpose(eigenvectors(1:k,1:k)), matmul(transpose(conjg(Q(:,1:k))), psi))
          do i = 1,k
             dummy(i) = cmplx(cos(-local_dt*eigenvalues(i)), sin(-local_dt*eigenvalues(i)), 8)*dummy(i)
          end do
          ans(:) = matmul(Q(:,1:k), matmul(eigenvectors(1:k,1:k), dummy(1:k)))

          ! check for convergence
          error = sqrt(dot_product(phi-ans, phi-ans))
          if (error < lanc_threshold) then
             converged = .TRUE.

          else if (k /= lancz_itnum) then
             w(:) = w(:) - beta(k-1)*Q(:,k-1) - alpha(k)*Q(:,k)
             beta(k) = sqrt(dot_product(w,w))
             Q(:,k+1) = w(:)/beta(k)

             call zschmab(Q(:,k-1:k), Q(:,k+1), 1.d-10, size(psi), 2, 1, m)
          
          end if

          k = k+1
          phi(:) = ans(:)
       end do

       ans(:) = phi(:)

    end if

  end function lanczos_prop


  subroutine init_lanczos(mat)
    type(banded_sym_mat)                    :: mat
  end subroutine init_lanczos

  subroutine init_itvolt_exp(mat)
    type(banded_sym_mat)                    :: mat
  end subroutine init_itvolt_exp
  
  
end module propagator
