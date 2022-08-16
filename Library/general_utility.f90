module general_utility
  implicit none
  public sym_tri_matmul, schmab, zschmab, factorial_array


  contains


  ! mutliplies a real symmetric tridiagonal matrix and a complex vector
  function sym_tri_matmul(diagonal, off_diagonal, v) result(ans)
    real(8), intent(in)               :: diagonal(:)
    real(8), intent(in)               :: off_diagonal(:)
    complex(8), intent(in)            :: v(:)
    complex(8), allocatable           :: ans(:)
    integer                           :: m, i

    m = size(diagonal)
    allocate(ans(1:m))

    ans(1) = diagonal(1)*v(1) + off_diagonal(1)*v(2)
    do i = 2,m-1
       ans(i) = off_diagonal(i-1)*v(i-1) + diagonal(i)*v(i) + off_diagonal(i)*v(i+1)
    end do
    ans(m) = off_diagonal(m-1)*v(m-1) + diagonal(m)*v(m)
    
  end function sym_tri_matmul


  
  subroutine schmab(va,vb,thresh,n,na,nb,nout,prnt)
    ! c***begin prologue     schmab
    ! c***date written       960801  (yymmdd)
    ! c***revision date      170820  (yymmdd)
    ! c
    ! c***keywords           gram-schmidt, orthogonalization
    ! c***author             barry schneider(nsf), revised heman gharibnejad (nist)
    ! c***source
    ! c***purpose            gram-schmidt orthogonalization.
    ! c***description        a set of non-orthonormal vectors, vb, are orthogonalized 
    ! c                      to another set of vectors, va, using a gram-schmidt process 
    ! c                      that checks for linear dependencies. the set of vectors,
    ! c                      va, are already assumed to be internally orthonormal. 
    ! c
    ! c                          va(n,*) = input vectors
    ! c                          vb(n,*) = input as non-orthogonal vectors and output
    ! c                                    as orthonormal set.
    ! c                          thresh  = acceptance tolerance on overlap
    ! c                          n       = dimension of vector space
    ! c                          na      = number of va vectors
    ! c                          nb      = number of initial vb vectors
    ! c                          nout    = number of final vb vectors
    ! c***references
    ! c***routines called    saxpy(clams), sdot(clams), sscal(clams)
    ! c
    ! c***end prologue       gschmt
    integer :: i, j, k
    integer :: n, na, nb
    integer :: nout

    real(8), dimension(n,na) :: va
    real(8), dimension(n,nb) :: vb
    real(8) :: ddot, norm, thresh, ovrlap
    logical, optional :: prnt
    ! print*,"inside"
    ! print*,"va"
    !      do i=1,na       
    !       print*, i, va(:,i)
    ! print*,'******************'
    !       end do
    !        print*, 'vb:'
    !        print*,vb
    nout=0
    do  i=1,nb
      do  j=1,na
        ovrlap=sum(va(:,j)*vb(:,i))!ddot(n,va(:,j),1,vb(:,i),1)
        !             if(na==15 .and. j==4)then
        !                do k=1,n
        !                  write(*,*)'va(i,4):',k, va(k,j)
        !                  write(*,*)'vb(i,15):',k,i, vb(k,i)
        !                  write(*,*) 'mult:',va(k,j)*vb(k,i)
        !                end do  
        !                stop
        !             endif   
        !            print*,'overlap:',j,ovrlap
        call daxpy(n,-ovrlap,va(:,j),1,vb(:,i),1)
        !             print*,'vb after saxpy:'
        !             print*, vb(:,i)
        !            print*,'****************************************************'
      end do

      norm=sqrt(sum(vb(:,i)*vb(:,i)))
      !sqrt(ddot(n,vb(:,i),1,vb(:,i),1))
      !          print*, 'norm',norm
      if(norm.gt.thresh) then
        !call dscal(n,1.0d+00/norm,vb(:,i),1)
        vb(:,i)=vb(:,i)/norm
        nout=nout+1
        !call copy(vb(1,i),vb(1,nout),n)
        vb(:,nout) = vb(:,i)
      endif
    end do
    if(present(prnt).and.prnt) then
      write(*,1) na, nb
      write(*,2) nout
    endif             
    if(nout.eq.0) then
      stop 'no vectors from schmab'
    endif
    !
    !
    !       print*,'vb being passed on'
    !       print*,vb
    !       print*,"*******************"
    return
    1    format(/,1x,'schmidt orthogonalization of two sets of vectors', &
      /,1x,                                                     &
      'set a already assumed to be orthonormal',/,1x,     &
      'number of vectors in set a = ',i4,/,1x,            &
      'number of vectors in set b = ',i4)                 
    2    format(/,1x,'number of set b vectors after orthogonalization '  &
      'to set a vectors = ',i4)
  end subroutine schmab

  
  
  !*************************************************************************
  !******************************************************************************************
  subroutine zschmab(va,vb,thresh,n,na,nb,nout)!,prnt)
    ! c***begin prologue     schmab
    ! c***date written       960801  (yymmdd)
    ! c***revision date      170820  (yymmdd)
    ! c
    ! c***keywords           gram-schmidt, orthogonalization
    ! c***author             barry schneider(nsf), revised heman gharibnejad (nist)
    ! c***source
    ! c***purpose            gram-schmidt orthogonalization.
    ! c***description        a set of non-orthonormal vectors, vb, are orthogonalized 
    ! c                      to another set of vectors, va, using a gram-schmidt process 
    ! c                      that checks for linear dependencies. the set of vectors,
    ! c                      va, are already assumed to be internally orthonormal. 
    ! c
    ! c                          va(n,*) = input vectors
    ! c                          vb(n,*) = input as non-orthogonal vectors and output
    ! c                                    as orthonormal set.
    ! c                          thresh  = acceptance tolerance on overlap
    ! c                          n       = dimension of vector space
    ! c                          na      = number of va vectors
    ! c                          nb      = number of initial vb vectors
    ! c                          nout    = number of final vb vectors
    ! c***references
    ! c***routines called    saxpy(clams), sdot(clams), sscal(clams)
    ! c
    ! c***end prologue       gschmt
    integer :: i, j, k
    integer :: n, na, nb
    integer :: nout

    complex(8), dimension(n,na) :: va
    complex(8), dimension(n,nb) :: vb
    complex(8) :: zdotc, ovrlap
    real(8)    :: norm, thresh
    !logical, optional :: prnt
    ! print*,"inside"
    ! print*,"va"
    !      do i=1,na       
    !       print*, i, va(:,i)
    ! print*,'******************'
    !       end do
    !        print*, 'vb:'
    !        print*,vb
    nout=0
    do  i=1,nb
      do  j=1,na
        ovrlap= dot_product(va(:,j),vb(:,i))!!zdotc(n,va(:,j),1,vb(:,i),1)!sum(va(:,j)*vb(:,i))!ddot(n,va(:,j),1,vb(:,i),1)
        !             if(na==15 .and. j==4)then
        !                do k=1,n
        !                  write(*,*)'va(i,4):',k, va(k,j)
        !                  write(*,*)'vb(i,15):',k,i, vb(k,i)
        !                  write(*,*) 'mult:',va(k,j)*vb(k,i)
        !                end do  
        !                stop
        !             endif   
        !            print*,'overlap:',j,ovrlap
        call zaxpy(n,-ovrlap,va(:,j),1,vb(:,i),1)
        !             print*,'vb after saxpy:'
        !             print*, vb(:,i)
        !            print*,'****************************************************'
      end do

      norm=real(sqrt(dot_product(vb(:,i),vb(:,i))))!zdotc(n,vb(:,i),1,vb(:,i),1)))
      !sqrt(ddot(n,vb(:,i),1,vb(:,i),1))
      !          print*, 'norm',norm
      if(norm.gt.thresh) then
        !call dscal(n,1.0d+00/norm,vb(:,i),1)
        vb(:,i)=vb(:,i)/norm
        nout=nout+1
        !call copy(vb(1,i),vb(1,nout),n)
        vb(:,nout) = vb(:,i)
      end if
    end do
   
    !print*, present(prnt).and. prnt
  ! if(present(prnt).and.prnt) then
  !    write(*,1) na, nb
  !    write(*,2) nout
  ! end if             
  ! if(nout.eq.0) then
  !     stop 'no vectors from schmab'
  ! end if
    
  !       print*,'vb being passed on'
  !       print*,vb
  !       print*,"*******************"
  !  return
  !  1    format(/,1x,'schmidt orthogonalization of two sets of vectors', &
  !    /,1x,                                                     &
  !    'set a already assumed to be orthonormal',/,1x,     &
  !    'number of vectors in set a = ',i4,/,1x,            &
  !    'number of vectors in set b = ',i4)                 
  !  2    format(/,1x,'number of set b vectors after orthogonalization '  &
  !    'to set a vectors = ',i4)
    
  end subroutine zschmab

  

  !*************************************************************************
  ! Returns factorials up to the integer argument as an array of 
  ! double precision numbers.
  !*************************************************************************
  function factorial_array(max_number) result(array)
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, intent(in)            :: max_number
      integer                        :: i
      real(8), dimension(max_number) :: array
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(max_number == 0) array(1)=1.d0
      if(max_number < 0 ) stop 'argument of factorial must be positive'

      array(1) = 1.d0
      if(max_number > 1) then
        do i = 2,max_number
          array(i)=array(i-1)*real(i)
        end do
      end if

  end function factorial_array
  !**************************************************************************
  
 end module general_utility
