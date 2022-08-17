module propagator
use parameters !, only : idt, ii, dt, t_intv, datafilename, z_one, z_zero, lancz_itnum, lanc_threshold, chebyshev_terms, chebyshev_threshold
use banded_matrices
use pulse
use timing
use general_utility
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
     
      case default 
        print*, propagator_name
        stop 'above method is not programmed'
        
    end select

    ! if we didn't find a function we liked in the switch, then undefined behavior
    if (.NOT. ASSOCIATED(f_ptr)) then
      print *, propagator_name
      stop "invalid pulse type"
    end if
   
  end subroutine select_propagator_type

  
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
       ans(i) = cmplx(cos(-local_dt*mat%eigenvalues(i)), sin(-local_dt*mat%eigenvalues(i)))*ans(i)
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
       phi(:) = cmplx(cos(-local_dt*alpha(1)), sin(-local_dt*alpha(1)))*dot_product(Q(:,1), psi)*Q(:,1)

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
             dummy(i) = cmplx(cos(-local_dt*eigenvalues(i)), sin(-local_dt*eigenvalues(i)))*dummy(i)
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
  
  
end module propagator
