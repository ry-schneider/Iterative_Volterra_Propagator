module grid

  implicit none
  private
  !box variables
  public  x, n, m_size, l_indx, h_indx, &
    spatial_grid, spatial_grid_fill, &
    spatial_grid_fill_sym_z

  !***********************************************************************
  ! Grid variables
  !***********************************************************************    
  integer             :: n  !number of spatial steps in each direction
  !> lowest index of allocated arrays and matrices
  integer             :: l_indx
  !> highest index of allocated arrays and matrices
  integer             :: h_indx
  integer             :: m_size ! 2n+1 size of matrices
  real(8),allocatable :: x(:) !spatial grid
  !************************************************************************
contains

  subroutine spatial_grid(l,h)
    ! Generating the grid fo box of 2n+1 gird points centered at zero
    !**************************************************************************************
    integer             :: i
    integer             :: ierr
    real(8),intent(in)  :: l  !size of box
    real(8),intent(in)  :: h  !size of spatial steps
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	   
    n = int(l/h)
    m_size = n
    l_indx = 1
    h_indx = m_size

    ! x is a grid from -n to n
    if(.not. allocated(x))then
      allocate(x(1:n),stat = ierr)
      if ( ierr /= 0 ) then
        write(*,*) 'Oops, x allocation in grid failed!'
      endif
    endif
    do i = 1,n 
      x(i) = i * h - (l / 2d0)
    end do


    call print_grid
  end subroutine spatial_grid
  !***************************************************************************************

  subroutine print_grid
    integer    :: num
    
    num = min(10, size(x))
    write(*,*) "Lowest values in grid:"
    write(*,*) x(1:num)
    write(*,*) "Highest Values in grid:"
    write(*,*) x(size(x) - num + 1: size(x))
    
  end subroutine print_grid

  !> Creates a spatial grid starting at start, and goes for length l of stepping
  !! forward by h
  !! \deprecated
  subroutine spatial_grid_fill(l, h, start)
    real(8), intent(in)   :: l, h, start
    integer               :: i

    n = l / h
    write(*,*) n
    m_size = n
    l_indx = 1
    h_indx = m_size
    allocate(x(1:n))
    do i = 1, n
      x(i) = (i-1)*h + start
    end do
    call print_grid
  end subroutine spatial_grid_fill


  !> Fills the grid up symetrically.
  !! 
  subroutine spatial_grid_fill_sym_z(xmax,h)
    real(8), intent(in) :: h
    real(8), intent(in) :: xmax
    integer             :: ngrid
    integer             :: i, center

    n = int(xmax/h) * 2
    if (mod(n, 2) .eq. 0) then
      !stop "Not making a grid with zero"
      ngrid = n + 1
    else
      ngrid = n  
    end if


    ! We need to allocate X
    if (allocated(x)) then
      deallocate(x)
    end if

    allocate(x(ngrid))
    n = ngrid

    ! Set up the low high index thing
    m_size = n
    l_indx = 1
    h_indx = m_size

    !h = xmax / dble((ngrid - 1) / 2d0)
    write(*,*) h

    center = (ngrid / 2) + 1
    x(center) = 0d0
    do i = 1, center - 1
      x(center + i) = h * i
      x(center - i) = -x(center + i)
    end do



    call print_grid
  end subroutine spatial_grid_fill_sym_z

  
  
end module
