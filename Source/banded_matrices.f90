module banded_matrices
  use potential
  use parameters
  implicit none
  private
  public banded_sym_mat, complex_sb_mat, sb_matvec_mul, d_sb_eigensolve, d_sb_eigenvalues

  ! real symmetric banded matrix
  type banded_sym_mat
    integer                 :: mat_size
    integer                 :: bsz
    real(8), allocatable    :: diagonal(:)
    real(8), allocatable    :: offdiagonal(:,:)
    real(8), allocatable    :: eigenvalues(:)
    real(8), allocatable    :: eigenvectors(:,:)
    logical                 :: init = .false.

  contains

    procedure :: d_initialize !(mat_size, band_size, diagonal(mat_size), offdiagonal(mat_size-1))
    generic   :: initialize => d_initialize
    procedure :: d_make_banded_matrix_on_grid !(x(mat_size),dx,band_size)
    generic   :: make_banded_matrix_on_grid => d_make_banded_matrix_on_grid

  end type banded_sym_mat

  ! complex banded symmetric matrix
  type complex_sb_mat
    integer                   :: mat_size ! matrix size
    integer                   :: bsz ! number of nonzero bands
    complex(8), allocatable   :: diagonal(:) ! main diagonal
    complex(8), allocatable   :: offdiagonal(:,:) ! nonzero off-diagonal bands
    complex(8), allocatable   :: eigenvalues(:)
    complex(8), allocatable   :: eigenvectors(:,:)
    logical                   :: init = .false.

  contains

    procedure :: z_initialize !(mat_size, band_size, diagonal(mat_size), offdiagonal(mat_size-1,band_size))
    generic   :: initialize => z_initialize
    procedure :: z_make_banded_matrix_on_grid !(x(mat_size), dx, band_size)
    generic   :: make_banded_matrix_on_grid => z_make_banded_matrix_on_grid
     
  end type complex_sb_mat

  ! generic call to real/complex matrix vector multiplication routines
  interface sb_matvec_mul
     procedure d_sb_matvec_mul
     procedure z_sb_matvec_mul
  end interface sb_matvec_mul

  
contains

  
  subroutine d_initialize(self, m_size, b_size, diag, offd)
    !******************************************************************
    ! Initializes and allocates the components of the
    ! real banded symmetric matrix (self).  
    !  m_size : symmetric matrix size 
    !  b_size : number of bands in matrix
    !  diag   : the diagonal vector of the matrix dimension(m_size)
    !  offd   : the off diagonal matrix has dimensions (m_size, b_size)
    !*******************************************************************           
    class(banded_sym_mat) :: self
    integer, intent(in)   :: m_size, b_size
    real(8), intent(in)   :: diag(:), offd(:,:)
    integer               :: i, j

    self%mat_size = m_size
    self%bsz      = b_size
    allocate(self%diagonal(m_size))
    allocate(self%offdiagonal(m_size-1, b_size))
    allocate(self%eigenvalues(m_size))
    allocate(self%eigenvectors(m_size, m_size))

    do j = 1, b_size
      do i = 1, m_size - 1
        self%offdiagonal(i, j) = offd(i, j)
      end do
    end do

    do i = 1, m_size
      self%diagonal(i) = diag(i)
    end do

    self%init = .true.

  end subroutine d_initialize


  subroutine z_initialize(self, m_size, b_size, diag, offd)
    !******************************************************************
    ! Initializes and allocates the components of the complex banded matrix (self).  
    !  m_size  : symmetric matrix size 
    !  b_size  : number of nonzero bands in matrix
    !  diag    : the diagonal vector of the matrix (m_size)
    !  offd    : the nonzero off diagonal bands (m_size-1, b_size) 
    !*******************************************************************           
    class(complex_sb_mat)    :: self
    integer, intent(in)      :: m_size, b_size
    complex(8), intent(in)   :: diag(:), offd(:,:)
    integer                  :: i, j
    
    self%mat_size = m_size
    self%bsz      = b_size
    allocate(self%diagonal(m_size))
    allocate(self%offdiagonal(m_size-1, b_size))
    allocate(self%eigenvalues(m_size))
    ! allocate(self%eigenvectors(m_size, m_size))

    self%diagonal(:) = diag(:)
    self%offdiagonal(:,:) = 0
    do j = 1,b_size
       do i = 1,m_size-j
         self%offdiagonal(i,j) = offd(i,j)
       end do
    end do

    self%init = .true.

  end subroutine z_initialize  

  
  function d_sb_matvec_mul(mat, vec) result(ans)
    !***************************************************************
    !  Multiplies symmetric_banded matrix into a complex vector 
    !  mat : real symmetric banded matrix (NxN)
    !  vec : complex vector the matrix is multiplied into (N)
    !  ans : the outcome complex vector of multiplication (N)
    !**************************************************************
    type(banded_sym_mat), intent(in)  :: mat
    complex(8), intent(in)            :: vec(:)
    complex(8)                        :: ans(size(vec))
    integer :: i, j
    integer :: b, m_size
    complex(8) :: z_zero = (0.d0, 0.d0)
    
    if (mat%init) then

      b = mat%bsz
      m_size = mat%mat_size

      ans(1:m_size) = mat%diagonal(1:m_size) * vec(1:m_size)

      do i = 1,b
        ans(1:m_size) = ans(1:m_size) + [mat%offdiagonal(1:m_size-i,i) * vec(i+1:m_size),(z_zero, j=1,i)]
        ans(1:m_size) = ans(1:m_size) + [(z_zero, j=1,i), mat%offdiagonal(1:m_size-i,i) * vec(1:m_size-i)]
      end do

     else
      stop 'type banded_sym_mat is not initilized!'
     end if

  end function d_sb_matvec_mul


  function z_sb_matvec_mul(mat, vec) result(ans)
    !***************************************************************
    !  Multiplies symmetric_banded matrix into a complex vector 
    !  mat : complex symmetric banded matrix (NxN)
    !  vec : complex vector the matrix is multiplied into (N)
    !  ans : the outcome complex vector of multiplication (N)
    !  CAUTION: ONLY APPLICABLE TO LINEARLY POLARIZED HYDROGEN
    !**************************************************************
    type(complex_sb_mat), intent(in)  :: mat
    complex(8), intent(in)            :: vec(:)
    complex(8)                        :: ans(size(vec))
    integer                           :: i, j, m
    complex(8)                        :: z_zero = (0.d0, 0.d0)
    
    if (mat%init) then
       
      m = mat%mat_size

      ans(1:m) = mat%diagonal(1:m)*vec(1:m)

      do i = 1,mat%bsz-1
         ans(1:m) = ans(1:m) + [mat%offdiagonal(1:m-i,i)*vec(i+1:m),(z_zero, j=1,i)]
         ans(1:m) = ans(1:m) + [(z_zero, j=1,i), mat%offdiagonal(1:m-i,i)*vec(1:m-i)]
      end do

      ans(1:m) = ans(1:m) + [mat%offdiagonal(1:m-r_size,mat%bsz)*vec(r_size+1:m),(z_zero, j=1,r_size)]
      ans(1:m) = ans(1:m) + [(z_zero, j=1,r_size), mat%offdiagonal(1:m-r_size,mat%bsz)*vec(1:m-r_size)]

     else
        stop 'type complex_sb_mat is not initilized!'
     end if

  end function z_sb_matvec_mul


  !function sb_matvec_mul_real(mat, vec) result(ans)
  !  !***************************************************************
  !  !  Multiplies symmetric_banded matrix into a real vector 
  !  !  mat : matrix of symmetrix_banded_type (NxN)
  !  !  vec : real vector the matrix is multiplied into (N)
  !  !  ans : the outcome real vector of multiplication (N)
  !  !**************************************************************
  !  type(banded_sym_mat), intent(in)  :: mat
  !  real(8), intent(in)               :: vec(:)
  !  real(8)                           :: ans(size(vec))
  !  integer :: i, j
  !  integer :: b, m_size
    
  !  if (mat%init) then

  !    b = mat%bsz
  !    m_size = mat%mat_size

  !    ans(1:m_size) = mat%diagonal(1:m_size) * vec(1:m_size)

  !    do i = 1,b
  !      ans(1:m_size) = ans(1:m_size) + [mat%offdiagonal(1:m_size-i,i) * vec(i+1:m_size),(0.d0, j=1,i)]
  !      ans(1:m_size) = ans(1:m_size) + [(0.d0, j=1,i), mat%offdiagonal(1:m_size-i,i) * vec(1:m_size-i)]
  !    end do

  !   else
  !    stop 'type banded_sym_mat is not initilized!'
  !   end if

  !end function sb_matvec_mul_real 


  ! multiplies two banded symmetric matrices (mat1*mat2)
  ! requires that mat1 and mat2 have the same number of bands
  !function sb_matmul(mat1, mat2) result(ans)
  !  type(banded_sym_mat), intent(in)         :: mat1, mat2
  !  real(8)                                  :: ans(mat1%mat_size, mat1%mat_size)
  !  real(8)                                  :: diagonal(mat1%mat_size), off_diagonal(mat1%mat_size,mat1%bsz)
  !  real(8)                                  :: dummy(mat1%mat_size)
  !  integer                                  :: i, j, k

  !  if (mat1%mat_size /= mat2%mat_size) then
  !     print *, 'incorrect dimensions for matrix multiplication'
  !     stop
  !  end if

  !  if (mat1%bsz /= mat2%bsz) then
  !     print *, 'matrices must have the same number of bands to multiply'
  !     stop
  !  end if

  !  off_diagonal(:,:) = 0
  !  diagonal(:) = 0

  !  do i = 1, mat1%bsz
  !     dummy(:) = 0
  !     do j = 1,i-1
  !        dummy(j) = mat2%offdiagonal(j,i-j)
  !     end do

  !     dummy(i) = mat2%diagonal(i)

  !     do j = i+1, i+mat1%bsz
  !        dummy(j) = mat2%offdiagonal(i,j-i)
  !     end do

  !     ans(:,i) = sb_matvec_mul_real(mat1, dummy)
       
  !  end do

  !  do i = mat1%bsz+1, mat1%mat_size-mat1%bsz
  !     dummy(:) = 0
  !     do j = i-mat1%bsz, i-1
  !        dummy(j) = mat2%offdiagonal(j,i-j)
  !     end do

  !     dummy(i) = mat2%diagonal(i)

  !     do j = i+1, i+mat1%bsz
  !        dummy(j) = mat2%offdiagonal(i,j-i)
  !     end do

  !     ans(:,i) = sb_matvec_mul_real(mat1, dummy)
             
  !  end do

  !  do i = mat1%mat_size-mat1%bsz+1, mat1%mat_size
  !     dummy(:) = 0
  !     do j = i-mat1%bsz, i-1
  !        dummy(j) = mat2%offdiagonal(j,i-j)
  !     end do

  !     dummy(i) = mat2%diagonal(i)

  !     do j = i+1, mat1%mat_size
  !        dummy(j) = mat2%offdiagonal(i,j-i)
  !     end do

  !     ans(:,i) = sb_matvec_mul_real(mat1, dummy)
       
  !  end do
    
  !end function sb_matmul


  subroutine d_make_banded_matrix_on_grid(self, x_grid, dx, bandnum)
    !***********************************************************************************************   
    ! Constructs the matrix of the Hamiltonian d^2/dx^2 + potential on
    ! as a real banded symmetric matrix (self)
    !   xgrid      : on input, the already filled grid points vector
    !   dx         : the size of the uniform grid increments
    !   bandnum    : determines the size of finite-difference scheme to be used
    !                (i.e. the size of the bands in matrix)
    !***********************************************************************************************
    class(banded_sym_mat)  :: self
    real(8),intent(in)     :: x_grid(:)
    real(8)                :: dx
    integer                :: bandnum
    integer                :: ndim
    real(8),allocatable    :: diag(:), offdiag(:,:)
    procedure(potential_func), pointer :: pot


    pot => select_potential_type(potential_type)

    ndim = size(x_grid)

    allocate(diag(ndim), offdiag(ndim-1,bandnum))

    select case (bandnum)

    case(1)

      diag(:) = (1.d0 / dx**2) +pot(x_grid(:))
      offdiag(1:ndim-1,1) = -0.5d0/dx**2


    case(2)

      diag(:) = (2.5d0 /2.d0 /dx**2) + pot(x_grid(:))
      offdiag(1:ndim-1,1) = -2d0/3.d0/dx**2
      offdiag(1:ndim-2,2) = 1d0/24.d0/dx**2

    case(3)

      diag(:) = (49.d0 / 36.d0/dx**2) &
        + pot(x_grid(:))
      offdiag(1:ndim-1,1) = -1.5d0/2.d0/dx**2
      offdiag(1:ndim-2,2) = 3.d0/40.d0/dx**2
      offdiag(1:ndim-3,3) = -1.d0/180.d0/dx**2

    case(4)

      diag(:) = (205.d0 / 144.d0/dx**2) &
        + pot(x_grid(:))
      offdiag(1:ndim-1,1) = -4.d0/5.d0/dx**2
      offdiag(1:ndim-2,2) = 1.d0/10.d0/dx**2
      offdiag(1:ndim-3,3) = -4.d0/315.d0/dx**2
      offdiag(1:ndim-4,4) = 1.d0/1120.d0/dx**2


    case default

      stop 'only matcies with 1-4 bands are coded'

    end select

    call self%initialize(ndim, bandnum, diag, offdiag)

  end subroutine d_make_banded_matrix_on_grid


  subroutine z_make_banded_matrix_on_grid(self, x_grid, dx, bandnum)
    !***********************************************************************************************   
    ! Constructs the matrix of the Hamiltonian d^2/dx^2 + potential on
    ! as a complex banded symmetric matrix (self)
    !   xgrid      : on input, the already filled grid points vector
    !   dx         : the size of the uniform grid increments
    !   bandnum    : determines the size of finite-difference scheme to be used
    !                (i.e. the size of the bands in matrix)
    !***********************************************************************************************
    class(complex_sb_mat)     :: self
    real(8), intent(in)       :: x_grid(:)
    real(8)                   :: dx
    integer                   :: bandnum
    integer                   :: ndim
    complex(8), allocatable   :: diag(:), offdiag(:,:)
    procedure(potential_func), pointer :: pot


    pot => select_potential_type(potential_type)

    ndim = size(x_grid)

    allocate(diag(ndim), offdiag(ndim-1,bandnum))

    select case (bandnum)

    case(1)

      diag(:) = (1.d0 / dx**2) - ii*pot(x_grid(:))
      offdiag(1:ndim-1,1) = -0.5d0/dx**2

    case(2)

      diag(:) = (2.5d0 /2.d0 /dx**2) - ii*pot(x_grid(:))
      offdiag(1:ndim-1,1) = -2d0/3.d0/dx**2
      offdiag(1:ndim-2,2) = 1d0/24.d0/dx**2

    case(3)

      diag(:) = (49.d0 / 36.d0/dx**2) - ii*pot(x_grid(:))
      offdiag(1:ndim-1,1) = -1.5d0/2.d0/dx**2
      offdiag(1:ndim-2,2) = 3.d0/40.d0/dx**2
      offdiag(1:ndim-3,3) = -1.d0/180.d0/dx**2

    case(4)

      diag(:) = (205.d0 / 144.d0/dx**2) - ii*pot(x_grid(:))
      offdiag(1:ndim-1,1) = -4.d0/5.d0/dx**2
      offdiag(1:ndim-2,2) = 1.d0/10.d0/dx**2
      offdiag(1:ndim-3,3) = -4.d0/315.d0/dx**2
      offdiag(1:ndim-4,4) = 1.d0/1120.d0/dx**2

    case default

      stop 'only matrices with 1-4 bands are coded'

    end select

    call self%initialize(ndim, bandnum, diag, offdiag)

  end subroutine z_make_banded_matrix_on_grid


  ! computes a  full set of eigenvalues and eigenvectors for real symmetric banded matrices
  subroutine d_sb_eigensolve(mat)
    type(banded_sym_mat)            :: mat 
    real(8)                         :: diagonal(mat%mat_size)
    real(8)                         :: off_diagonal(mat%mat_size - 1)
    real(8)                         :: work(2*mat%mat_size - 2)
    real(8)                         :: AB(mat%bsz + 1, mat%mat_size)
    real(8)                         :: Q(mat%mat_size, mat%mat_size)
    real(8)                         :: work2(7*mat%mat_size)
    real(8)                         :: VL, VU
    integer                         :: info, IL, IU, i, j, m
    integer                         :: iwork(5*mat%mat_size), ifail(mat%mat_size)

    ! tridiagonal case
    if (mat%bsz == 1) then
       mat%eigenvalues(:) = mat%diagonal(:)
       off_diagonal(:) = mat%offdiagonal(:,1)

       call dstev('V', mat%mat_size, mat%eigenvalues, off_diagonal, mat%eigenvectors, mat%mat_size, work, info)

    ! general case   
    else
       AB = 0
       do j = 1,mat%mat_size
          AB(1,j) = mat%diagonal(j)
          do i = j+1,min(mat%mat_size, j+mat%bsz)
             AB(1+i-j,j) = mat%offdiagonal(i, i-j)
          end do
       end do

       call dsbevx('V', 'A', 'L', mat%mat_size, mat%bsz, AB, mat%bsz+1, Q, mat%mat_size, VL, VU, IL, &
            IU, 0.d0, m, mat%eigenvalues, mat%eigenvectors, mat%mat_size, work2, iwork, ifail, info)
            
    end if

    ! check for success
    if (info /= 0) then
       print *, 'Eigensolver failed.'
       stop
    end if
    
  end subroutine d_sb_eigensolve

  
  ! computes only a full set of eigenvalues for real symmetric banded matrices
  subroutine d_sb_eigenvalues(mat)
    type(banded_sym_mat)               :: mat
    real(8)                            :: diagonal(mat%mat_size)
    real(8)                            :: off_diagonal(mat%mat_size - 1)
    real(8)                            :: work(2*mat%mat_size - 2)
    real(8)                            :: AB(mat%bsz + 1, mat%mat_size)
    real(8)                            :: Q(mat%mat_size, mat%mat_size)
    real(8)                            :: VL, VU
    real(8)                            :: work2(7*mat%mat_size)
    integer                            :: iwork(5*mat%mat_size)
    integer                            :: ifail(mat%mat_size)
    integer                            :: info, IL, IU, i, j, m

    ! tridiagonal case
    if (mat%bsz == 1) then
       mat%eigenvalues(:) = mat%diagonal(:)
       off_diagonal(:) = mat%offdiagonal(:,1)

       call dstev('N', mat%mat_size, mat%eigenvalues, off_diagonal, mat%eigenvectors, mat%mat_size, work, info)

    ! general case
    else
       AB = 0
       do j = 1,mat%mat_size
          AB(1,j) = mat%diagonal(j)
          do i = j+1,min(mat%mat_size, j+mat%bsz)
             AB(1+i-j, j) = mat%offdiagonal(i, i-j)
          end do
       end do
       
       call dsbevx('N', 'A', 'L', mat%mat_size, mat%bsz, AB, mat%bsz + 1, Q, mat%mat_size, VL, VU, &
            IL, IU, 0d0, m, mat%eigenvalues, mat%eigenvectors, mat%mat_size, work2, iwork, ifail, info) 
       
    end if

    ! check for success
    if (info /= 0) then
       print *, 'Eigensolver failed.'
       stop
    end if
    
  end subroutine d_sb_eigenvalues


end module banded_matrices
