module propagator
use parameters 
use banded_matrices
use pulse_module
use timing
use general_utility
use Lagrange_weights
use grid, only: x
use omp_lib
implicit none
private
public  propagator_func, initialize_variables, select_propagator_type, itvolt_exp, diag_prop, &
     cheby_prop, lanczos_prop, arnoldi_prop, split_operator_cn

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
         init_ptr => d_sb_eigensolve
         f_ptr => diag_prop

      case('chebyshev')
         init_ptr => d_sb_eigenvalues
         f_ptr => cheby_prop

      case('lanczos')
         init_ptr => init_lanczos
         f_ptr => lanczos_prop

      case('itvolt')
         init_ptr => init_itvolt_exp
         f_ptr => itvolt_exp

      case('two_level')
         init_ptr => init_two_level
         f_ptr => two_level_prop

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
    complex(8), allocatable                         :: expH_0(:,:)
    complex(8), allocatable                         :: inv_expH_0(:,:) 
    complex(8)                                      :: ans(size(v))
    real(8)                                         :: sign
    real(8), allocatable                            :: pt(:)
    real(8), allocatable                            :: wt(:,:), comp_wt(:,:)
    complex(8), allocatable                         :: iterate(:,:), comp(:,:)
    complex(8), allocatable                         :: inhomogeneity(:,:)
    complex(8)                                      :: integral(size(v))
    complex(8)                                      :: b(size(v)), gs_diag(size(v)), gs_off1(size(v)-1), gs_off2(size(v)-1)
    logical                                         :: converged
    real(8)                                         :: diff, max_diff
    complex(8)                                      :: AB(3*mat%bsz+1,mat%mat_size)
    integer                                         :: IPIV(size(v))
    integer                                         :: n, i, j, k, r, p, u, s, a, c, d, it_count, info

    ! check sign of local_dt
    if (local_dt < 0d0) then
       sign = -1d0
    else
       sign = 1d0
    end if

    ! select number of quadrature points
    n = min(quad_pt, max(4, ceiling(quad_pt*local_dt/dt)))

    ! allocate arrays based on choice of n
    allocate(expH_0(1:size(v),1:n), inv_expH_0(1:size(v),1:n), pt(1:n), wt(1:n,1:n-1), comp_wt(1:n,1:n-1),&
         iterate(1:size(v),1:n), comp(1:size(v),1:n), inhomogeneity(1:size(v),1:n))

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
    
    do while (.not. (converged .or. it_count > it_cap-1))
       comp(:,:) = iterate(:,:)
       
       ! Jacobi iteration
       if (exp_it_type == 'jacobi') then
          do r = 2,n
             integral(:) = 0

             do p = 1,n
                integral(:) = integral(:) + comp_wt(p,r-1)*expH_0(:,p)*sb_matvec_mul(W_j, comp(:,p)) 
             end do

             iterate(:,r) = inhomogeneity(:,r) - ii*inv_expH_0(:,r)*integral(:)
          end do

       ! Gauss-Seidel iteration
       else if (exp_it_type == 'gauss_seidel') then
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
                ! comstruct AB for general lapack system solve
                AB = 0
                do c = 1,mat%mat_size
                   do d = max(1, c-mat%bsz),  min(mat%mat_size,c+mat%bsz)
                      if (c == d) then
                         AB(2*mat%bsz+1+d-c,c) = 1d0
                      else if (d < c) then 
                         AB(2*mat%bsz+1+d-c,c) = ii*comp_wt(r,r-1)*W_j%offdiagonal(d,c-d)
                      else
                         AB(2*mat%bsz+1+d-c,c) = ii*comp_wt(r,r-1)*W_j%offdiagonal(d,d-c)
                      end if
                   end do
                end do
                
                call zgbsv(mat%mat_size, mat%bsz, mat%bsz, 1, AB, 3*mat%bsz+1, IPIV, b, mat%mat_size, info)

                if (info /= 0) then
                   print *, 'exponential system solve failed'
                   stop
                end if

                iterate(:,r) = b(:)
             end if
             
          end do

       else
          print *, 'iteration type for exponential not programmed'
          stop
          
       end if

       it_count = it_count + 1
       
       ! check for convergence
       diff = sqrt(dot_product(iterate(:,n)-comp(:,n), iterate(:,n)-comp(:,n)))
       if (diff <= it_tolerance) then
          converged = .TRUE.
       end if

       ! max_diff = 0
       ! do u = 2,n
       !    diff = sqrt(dot_product(iterate(:,u)-comp(:,u), iterate(:,u)-comp(:,u)))
       !    if (max_diff < diff) then
       !       max_diff = diff
       !    end if
       ! end do

       ! if (max_diff <= 1.d-15) then
       !    converged = .TRUE.
       ! end if
       
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
    mat_norm%offdiagonal(:,:) = (2.0d0/delta)*mat_norm%offdiagonal(:,:)
    ! mat_norm%offdiagonal(:,1) = (2.0d0/delta)*mat_norm%offdiagonal(:,1)

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
    integer                                      :: k, info, i, m, j
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
       do while ( .not. (converged .or. k >= lancz_itnum))
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

             j = min(k, lancz_reortho)
             call zschmab(Q(:,k-j+1:k), Q(:,k+1), 1.d-10, size(psi), j, 1, m)
          
          end if

          k = k+1
          phi(:) = ans(:)
       end do

    end if

  end function lanczos_prop


  ! applies a (complex) matrix exponential to a vector via Arnoldi
  function arnoldi_prop(mat, local_dt, psi) result(ans)
    type(complex_sb_mat), intent(in)        :: mat
    real(8), intent(in)                     :: local_dt
    complex(8), intent(in)                  :: psi(:)
    complex(8)                              :: ans(size(psi))
    complex(8)                              :: Q(size(psi), arnoldi_itnum)
    complex(8)                              :: w(size(psi)), phi(size(psi))
    complex(8)                              :: H(arnoldi_itnum, arnoldi_itnum), A(arnoldi_itnum, arnoldi_itnum)
    complex(8)                              :: eigenvalues(arnoldi_itnum), eigenvectors(arnoldi_itnum, arnoldi_itnum)
    complex(8)                              :: eig_copy(arnoldi_itnum)
    complex(8)                              :: b(arnoldi_itnum)
    complex(8)                              :: VL(1,arnoldi_itnum), work(11*arnoldi_itnum), Z(1,arnoldi_itnum)
    complex(8)                              :: work2(arnoldi_itnum*arnoldi_itnum)
    real(8)                                 :: rwork(2*arnoldi_itnum)
    real(8)                                 :: error
    integer                                 :: k, m, i, j, l, info, ifaill, ifailr
    integer                                 :: IPIV(arnoldi_itnum)
    logical                                 :: converged
    logical                                 :: select_array(arnoldi_itnum)

    ! if starting vector is zero, skip iteration and output zero
    if (sqrt(dot_product(psi,psi)) == 0.d0) then
       ans(:) = 0

    else
       ! initialize hessenberg matrix as zero
       H(:,:) = 0

       ! first vector is normalized psi
       Q(:,1) = psi(:)/sqrt(dot_product(psi,psi))

       ! compute second arnoldi vector
       w(:) = sb_matvec_mul(mat,Q(:,1))
       H(1,1) = dot_product(Q(:,1),w)
       w(:) = w(:) - H(1,1)*Q(:,1)

       H(2,1) = sqrt(dot_product(w,w))
       Q(:,2) = w(:)/H(2,1)

       ! re-orthogonalize against first to ensure convergence
       call zschmab(Q(:,1), Q(:,2), 1.d-10, size(psi), 1, 1, m)

       ! project on first arnoldi vector (i.e., evaluate the exponential with only Q(:,1))
       phi(:) = exp(-ii*local_dt*H(1,1))*dot_product(Q(:,1),psi)*Q(:,1)

       ! start iteration, adding an extra vector until convergence is reached
       select_array(1) = .TRUE.
       select_array(2) = .TRUE.
       k = 2
       converged = .FALSE.
       do while ( .not. (converged .or. k >= arnoldi_itnum))
          ! print *, k
          w(:) = sb_matvec_mul(mat, Q(:,k))
          do i=1,k
             H(i,k) = dot_product(Q(:,i),w)
             w(:) = w(:) - H(i,k)*Q(:,i)
          end do

          ! find eigenvalues/eigenvectors of hessenberg matrix H(1:k,1:k)
          A(1:k,1:k) = H(1:k,1:k)
          call zgeev('N', 'V', k, A(1:k,1:k), k, eigenvalues(1:k), VL, 1, eigenvectors(1:k,1:k), k, work, &
               10*arnoldi_itnum,  rwork(1:2*k), info)
          if (info /= 0) then
             print *, 'eigensolve in arnoldi iteration failed at step ', k-1
          end if

          ! call Hessenberg routine for eigenvalues
          ! call zhseqr('E', 'N', k, 1, k, A(1:k,1:k), k, eigenvalues(1:k), Z, 1, work, 11*arnoldi_itnum, info)
          ! if (info /= 0) then
          !    print *, 'eigenvalue solve in arnoldi failed at step', k-1
          !    stop
          ! end if

          ! use inverse iteration to find corresponding eigenvectors
          ! A(1:k,1:k) = H(1:k,1:k)
          ! eig_copy(1:k) = eigenvalues(1:k)
          ! call zhsein('R', 'Q', 'N', select_array(1:k), k, A(1:k,1:k), k, eig_copy(1:k), VL(1,1:k), 1, &
          !      eigenvectors(1:k,1:k), k, k, m, work2(1:k*k), rwork(1:k), ifaill, ifailr, info)
          ! if (info /= 0) then
          !    print *, 'eiegenvector solve in arnoldi iteration failed at step', k-1
          !    stop
          ! end if

          ! use diagonalization to apply the exponential to psi (note eigenvectors are not orthonormal)
          b(1:k) = matmul(transpose(conjg(Q(:,1:k))),psi)
          A(1:k,1:k) = eigenvectors(1:k,1:k)
          call zgesv(k, 1, A(1:k,1:k), k, IPIV(1:k), b(1:k), k, info)
          do j = 1,k
             b(j) = exp(-ii*local_dt*eigenvalues(j))*b(j)
          end do
          ans(:) = matmul(Q(:,1:k), matmul(eigenvectors(1:k,1:k), b(1:k)))

          ! check for convergence
          error = sqrt(dot_product(phi-ans, phi-ans))
          ! print *, error
          if (error < arnoldi_threshold) then
             converged = .TRUE.

          else if (k /= arnoldi_itnum) then
             ! add next vector and reorthogonalize
             H(k+1,k) = sqrt(dot_product(w,w))
             Q(:,k+1) = w(:)/H(k+1,k)
             l = min(k, arnoldi_reortho)
             call zschmab(Q(:,k-l+1:k), Q(:,k+1), 1.d-10, size(psi), l, 1, m)
             
          end if

          k = k+1
          select_array(k) = .TRUE.
          phi(:) = ans(:)
          
       end do
       
    end if
    
  end function arnoldi_prop


  ! applies linearly polarized hydrogen hamiltonian exponential to psi via split operator and Crank-Nicolson (CN)
  ! h_grid is the r_reach x r_reach (tridiagonal) radial piece, to be updated for each l
  ! v is the l_max x l_max angular piece (tridiagonal and real), to be updated at each grid point
  ! psi is the wavefuction stored in a r_reach x l_max array
  function split_operator_cn(time_step, h_grid, v, r_grid, local_dt, psi, split_type) result(ans)
    type(complex_sb_mat), intent(in)        :: h_grid
    type(banded_sym_mat), intent(in)        :: v
    real(8), intent(in)                     :: r_grid(:)
    real(8), intent(in)                     :: local_dt
    complex(8), intent(in)                  :: psi(:,:)
    integer, intent(in)                     :: time_step, split_type
    complex(8)                              :: ans(size(r_grid),l_max)
    complex(8), allocatable                 :: l_diagonal(:), l_rhs(:)
    complex(8), allocatable                 :: l_offdiag(:), l_offdiag2(:)
    complex(8), allocatable                 :: h_diagonal(:)
    complex(8), allocatable                 :: r_diagonal(:), r_rhs(:)
    complex(8), allocatable                 :: r_offdiag(:), r_offdiag2(:)
    integer                                 :: i, j, k, d, info

    d = size(r_grid)

    ! Nico's split
    if (split_type == 1) then
       ! first exponential: e^{-i(dt)v}*psi
       !$OMP parallel do private(i,j,l_rhs,l_diagonal,l_offdiag,l_offdiag2,info)
       do i = 1,d
          allocate(l_rhs(1:l_max), l_diagonal(1:l_max), l_offdiag(1:l_max-1), l_offdiag2(1:l_max-1))
          ! construct RHS vector and coefficient diagonal/off diagonal
          l_rhs(:) = psi(i,:) - ii*r_grid(i)*(local_dt/2d0)*sb_matvec_mul(v,psi(i,:))
          do j = 1,l_max
             l_diagonal(j) = 1d0
          end do
          l_offdiag(:) = ii*r_grid(i)*(local_dt/2d0)*v%offdiagonal(:,1)
          l_offdiag2(:) = l_offdiag(:)

          ! call tridiagonal system solve
          call zgtsv(l_max, 1, l_offdiag, l_diagonal, l_offdiag2, l_rhs, l_max, info)
          if (info /= 0) then
             print *, 'system solve in split operator failed'
             stop
          end if

          ans(i,:) = l_rhs(:)
          deallocate(l_rhs, l_diagonal, l_offdiag, l_offdiag2)
       end do
       !$OMP end parallel do

       ! second exponential: e^{-i(dt)h_grid}*ans
       if (time_step <= int(t_intv/dt)-1) then
          !$OMP parallel do private(i,j,h_diagonal,r_rhs,r_diagonal,r_offdiag,r_offdiag2,info)
          do i = 1,l_max
             allocate(h_diagonal(1:d), r_rhs(1:d), r_diagonal(1:d), r_offdiag(1:d-1), r_offdiag2(1:d-1))
             ! update h_grid with l
             do j = 1,d
                h_diagonal(j) = h_grid%diagonal(j) + dble((i-1)*i)/(2d0*r_grid(j)**2)
             end do

             ! construct RHS vector and coefficient diagonal/off diagonal
             r_rhs(:) = ans(:,i) - ii*(local_dt/2d0)*h_diagonal(:)*ans(:,i)
             r_rhs(1:d-1) = r_rhs(1:d-1) - ii*(local_dt/2d0)*h_grid%offdiagonal(:,1)*ans(2:d,i)
             r_rhs(2:d) = r_rhs(2:d) - ii*(local_dt/2d0)*h_grid%offdiagonal(:,1)*ans(1:d-1,i)
             do j = 1,d
                r_diagonal(j) = 1d0 + ii*(local_dt/2d0)*h_diagonal(j)
             end do
             r_offdiag(:) = ii*(local_dt/2d0)*h_grid%offdiagonal(:,1)
             r_offdiag2(:) = r_offdiag(:)

             ! call tridiagonal system solve
             call zgtsv(d, 1, r_offdiag, r_diagonal, r_offdiag2, r_rhs, d, info)
             if (info /= 0) then
                print *, 'system solve in split operator failed'
                stop
             end if

             ans(:,i) = r_rhs(:)
             deallocate(h_diagonal, r_rhs, r_diagonal, r_offdiag, r_offdiag2)
          end do
          !$OMP end parallel do

       end if

    ! h_grid x v x h_grid split
    else if (split_type == 2) then
       ! first exponential: e^{-i(dt/2)h_grid}*psi
       allocate(h_diagonal(1:d), r_rhs(1:d), r_diagonal(1:d), r_offdiag(1:d-1), r_offdiag2(1:d-1))
       do i = 1,l_max
          ! update h_grid with l
          do j = 1,d
             h_diagonal(j) = h_grid%diagonal(j) + dble((i-1)*i)/(2d0*r_grid(j)**2)
          end do

         ! construct RHS vector and coefficient diagonal/off diagonal
         r_rhs(:) = psi(:,i) - ii*(local_dt/4d0)*h_diagonal(:)*psi(:,i)
         r_rhs(1:d-1) = r_rhs(1:d-1) - ii*(local_dt/4d0)*h_grid%offdiagonal(:,1)*psi(2:d,i)
         r_rhs(2:d) = r_rhs(2:d) - ii*(local_dt/4d0)*h_grid%offdiagonal(:,1)*psi(1:d-1,i)
         do j = 1,d
            r_diagonal(j) = 1d0 + ii*(local_dt/4d0)*h_diagonal(j)
         end do
         r_offdiag(:) = ii*(local_dt/4d0)*h_grid%offdiagonal(:,1)
         r_offdiag2(:) = r_offdiag(:)

         ! call tridiagonal system solve
         call zgtsv(d, 1, r_offdiag, r_diagonal, r_offdiag2, r_rhs, d, info)
         if (info /= 0) then
            print *, 'system solve in split operator failed'
            stop
         end if

         ans(:,i) = r_rhs(:)
      end do

      ! second exponential: e^{-i(dt)v}*ans
      allocate(l_rhs(1:l_max), l_diagonal(1:l_max), l_offdiag(1:l_max-1), l_offdiag2(1:l_max-1))
      do i = 1,d
         ! construct RHS vector and coefficient diagonal/off diagonal
         l_rhs(:) = ans(i,:) - ii*r_grid(i)*(local_dt/2d0)*sb_matvec_mul(v,ans(i,:))
         do j = 1,l_max
            l_diagonal(j) = 1d0
         end do
         l_offdiag(:) = ii*r_grid(i)*(local_dt/2d0)*v%offdiagonal(:,1)
         l_offdiag2(:) = l_offdiag(:)

         ! call tridiagonal system solve
         call zgtsv(l_max, 1, l_offdiag, l_diagonal, l_offdiag2, l_rhs, l_max, info)
         if (info /= 0) then
            print *, 'system solve in split operator failed'
            stop
         end if

         ans(i,:) = l_rhs(:)
      end do

      ! third exponential: e^{-i(dt/2)h_grid}*ans
      do i = 1,l_max
         ! update h_grid with l
         do j = 1,d
            h_diagonal(j) = h_grid%diagonal(j) + dble((i-1)*i)/(2d0*r_grid(j)**2)
         end do

         ! construct RHS vector and coefficient diagonal/off diagonal
         r_rhs(:) = ans(:,i) - ii*(local_dt/4d0)*h_diagonal(:)*ans(:,i)
         r_rhs(1:d-1) = r_rhs(1:d-1) - ii*(local_dt/4d0)*h_grid%offdiagonal(:,1)*ans(2:d,i)
         r_rhs(2:d) = r_rhs(2:d) - ii*(local_dt/4d0)*h_grid%offdiagonal(:,1)*ans(1:d-1,i)
         do j = 1,d
            r_diagonal(j) = 1d0 + ii*(local_dt/4d0)*h_diagonal(j)
         end do
         r_offdiag(:) = ii*(local_dt/4d0)*h_grid%offdiagonal(:,1)
         r_offdiag2(:) = r_offdiag(:)

         ! call tridiagonal system solve
         call zgtsv(d, 1, r_offdiag, r_diagonal, r_offdiag2, r_rhs, d, info)
         if (info /= 0) then
            print *, 'system solve in split operator failed'
            stop
         end if

         ans(:,i) = r_rhs(:)
      end do 

    ! v x h_grid x v split   
    else
      ! first exponential: e^{-i(dt/2)v}*psi
      allocate(l_rhs(1:l_max), l_diagonal(1:l_max), l_offdiag(1:l_max-1), l_offdiag2(1:l_max-1))
      do i = 1,d
         ! construct RHS vector and coefficient diagonal/off diagonal
         l_rhs(:) = psi(i,:) - ii*r_grid(i)*(local_dt/4d0)*sb_matvec_mul(v,psi(i,:))
         do j = 1,l_max
            l_diagonal(j) = 1d0
         end do
         l_offdiag(:) = ii*r_grid(i)*(local_dt/4d0)*v%offdiagonal(:,1)
         l_offdiag2(:) = l_offdiag(:)

         ! call tridiagonal system solve
         call zgtsv(l_max, 1, l_offdiag, l_diagonal, l_offdiag2, l_rhs, l_max, info)
         if (info /= 0) then
            print *, 'system solve in split operator failed'
            stop
         end if

         ans(i,:) = l_rhs(:)
      end do

      ! second exponential: e^{-i(dt)h_grid}*ans
      allocate(h_diagonal(1:d), r_rhs(1:d), r_diagonal(1:d), r_offdiag(1:d-1), r_offdiag2(1:d-1))
      do i = 1,l_max
         ! update h_grid with l
         do j = 1,d
            h_diagonal(j) = h_grid%diagonal(j) + dble((i-1)*i)/(2d0*r_grid(j)**2)
         end do

         ! construct RHS vector and coefficient diagonal/off diagonal
         r_rhs(:) = ans(:,i) - ii*(local_dt/2d0)*h_diagonal(:)*ans(:,i)
         r_rhs(1:d-1) = r_rhs(1:d-1) - ii*(local_dt/2d0)*h_grid%offdiagonal(:,1)*ans(2:d,i)
         r_rhs(2:d) = r_rhs(2:d) - ii*(local_dt/2d0)*h_grid%offdiagonal(:,1)*ans(1:d-1,i)
         do j = 1,d
            r_diagonal(j) = 1d0 + ii*(local_dt/2d0)*h_diagonal(j)
         end do
         r_offdiag(:) = ii*(local_dt/2d0)*h_grid%offdiagonal(:,1)
         r_offdiag2(:) = r_offdiag(:)

         ! call tridiagonal system solve
         call zgtsv(d, 1, r_offdiag, r_diagonal, r_offdiag2, r_rhs, d, info)
         if (info /= 0) then
            print *, 'system solve in split operator failed'
            stop
         end if

         ans(:,i) = r_rhs(:)
      end do

      ! final exponential: e^{-i(dt/2)v}*ans
      do i = 1,d
         ! construct RHS vector and coefficient diagonal/off diagonal
         l_rhs(:) = ans(i,:) - ii*r_grid(i)*(local_dt/4d0)*sb_matvec_mul(v,ans(i,:))
         do j = 1,l_max
            l_diagonal(j) = 1d0
         end do
         l_offdiag(:) = ii*r_grid(i)*(local_dt/4d0)*v%offdiagonal(:,1)
         l_offdiag2(:) = l_offdiag(:)

         ! call tridiagonal system solve
         call zgtsv(l_max, 1, l_offdiag, l_diagonal, l_offdiag2, l_rhs, l_max, info)
         if (info /= 0) then
            print *, 'system solve in split operator failed'
            stop
         end if

         ans(i,:) = l_rhs(:)
      end do

    end if
    
  end function split_operator_cn


  ! propagator for the two level atom
  ! exploits the fact that the Hamiltonian is always anti-diagonal
  function two_level_prop(mat, local_dt, psi) result(ans)
    type(banded_sym_mat), intent(in)        :: mat
    real(8), intent(in)                     :: local_dt
    complex(8), intent(in)                  :: psi(:)
    complex(8)                              :: ans(size(psi))
    real(8)                                 :: theta

    theta = mat%offdiagonal(1,1)*local_dt
    ans(1) = cos(theta)*psi(1) - ii*sin(theta)*psi(2)
    ans(2) = -ii*sin(theta)*psi(1) + cos(theta)*psi(2) 
    
  end function two_level_prop


  subroutine init_two_level(mat)
    type(banded_sym_mat)                    :: mat
  end subroutine init_two_level

  subroutine init_lanczos(mat)
    type(banded_sym_mat)                    :: mat
  end subroutine init_lanczos

  subroutine init_itvolt_exp(mat)
    type(banded_sym_mat)                    :: mat
  end subroutine init_itvolt_exp
  
  
end module propagator
