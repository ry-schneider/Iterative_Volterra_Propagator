!********************************************************************************
  Module  LaGrange_Module
!********************************************************************************
  USE precin
  USE Gauss_Quadrature
  ! USE Matrix_Print
!********************************************************************************
!********************************************************************************
  LOGICAL,           DIMENSION(5)                     :: prnt_lagrange
  CHARACTER(LEN=80), DIMENSION(5)                     :: print_title_lagrange                                    
  CHARACTER(LEN=80)                                   :: type_points!='newton-cotes'
  ! DATA print_title_lagrange /'print_interpolation=polynomials',                                        &             
  !                            'print_interpolation=derivatives_of_polynomials',                         &             
  !                            'print_interpolation=grid_points',                                        &             
  !                            'print_interpolation=interpolation_points',                               &             
  !                            'print_interpolation=lagrange_weights' /
!********************************************************************************
!********************************************************************************
                            TYPE  LaGrange
!********************************************************************************
!********************************************************************************              
                            Contains
!********************************************************************************
!********************************************************************************
                                Procedure  :: LaGrange_Polynomials
                                Procedure  :: Taylor_Polynomials
                                Procedure  :: LaGrange_Integration_Weights
                                Procedure  :: LaGrange_Integration_weights_2
!********************************************************************************
                            END TYPE  LaGrange
!********************************************************************************
!********************************************************************************    
!
                                  CONTAINS
!********************************************************************************
!********************************************************************************
!deck Taylor_Polynomials
!***begin prologue     Taylor_Polynomials
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           Interpolation_Points
!***author             schneider, b. i.(nsf)
!***source             Taylor_Polynomials
!***purpose            Taylor_Polynomials
!***                   
!***                   
!***                   
!***                   
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Taylor_Polynomials
  SUBROUTINE Taylor_Polynomials(lgr,p,dp,ddp,y,prnt_lagrange)
  IMPLICIT NONE
!
  CLASS(LaGrange)                                    :: lgr
  REAL(idp),            DIMENSION(:,:)               :: p
  REAL(idp),  OPTIONAL, DIMENSION(:,:)               :: dp
  REAL(idp),  OPTIONAL, DIMENSION(:,:)               :: ddp
  REAL(idp),            DIMENSION(:)                 :: y
  REAL(idp)                                          :: y_i
  INTEGER                                            :: n_x
  INTEGER                                            :: n_y
  INTEGER                                            :: i
  INTEGER                                            :: j
  LOGICAL,              DIMENSION(:)                 :: prnt_lagrange
  
!
  n_y=size(p,1)
  n_x=size(p,2)
  p(1:n_y,1) = 1.d0
  DO j = 1, n_y
     DO i = 2, n_x
        p(j,i) = y(j) * p(j,i-1)
     END DO
  END DO
  IF ( prnt_lagrange(1) ) THEN
      ! call Print_Matrix(type_real_matrix,p,n_y,n_x,title='Taylor Polynomials')
  END IF
  IF ( PRESENT(dp) ) THEN
       DO j = 1, n_y
          y_i = 1.d0 / y(j)
          DO i = 1, n_x
             dp(j,i)  = ( i - 1 ) * p(j,i) * y_i
             ddp(j,i) = ( i - 1 ) * ( i - 2 ) * p(j,i) * y_i * y_i
          END DO
       END DO
       IF ( prnt_lagrange(2) ) THEN
            ! call Print_Matrix(type_real_matrix,dp,n_y,n_x,title='First Derivative of Taylor Polynomials')
            ! call Print_Matrix(type_real_matrix,ddp,n_y,n_x,title='Second Derivative of Taylor Polynomials')
       END IF
  END IF
!
  END SUBROUTINE Taylor_Polynomials
!********************************************************************************
!********************************************************************************
!deck Lagrange_Polynomials
!***begin prologue     Lagrange_Polynomials
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           Interpolation_Points
!***author             schneider, b. i.(nsf)
!***source             Lagrange_Polynomials
!***purpose            Lagrange_Polynomials
!***                   
!***                   
!***                   
!***                   
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Lagrange_Polynomials
  SUBROUTINE Lagrange_Polynomials(lgr,p,dp,ddp,x,y,prnt_lagrange)
  IMPLICIT NONE
!
  CLASS(LaGrange)                            :: lgr
  REAL(idp),            DIMENSION(:,:)       :: p
  REAL(idp),  OPTIONAL, DIMENSION(:,:)       :: dp
  REAL(idp),  OPTIONAL, DIMENSION(:,:)       :: ddp
  REAL(idp),            DIMENSION(:)         :: x
  REAL(idp),            DIMENSION(:)         :: y
  REAL(idp),            DIMENSION(2)         :: sn
  INTEGER                                    :: n_x
  INTEGER                                    :: n_y
  INTEGER                                    :: i
  INTEGER                                    :: j
  INTEGER                                    :: k
  INTEGER                                    :: zerfac
  INTEGER                                    :: first
  INTEGER                                    :: second
  REAL(idp)                                  :: fac
  LOGICAL,              DIMENSION(:)         :: prnt_lagrange
!
  n_y=size(p,1)
  n_x=size(p,2)
  p(1:n_y,1:n_x) = 1.d0
  DO i = 1, n_y
     zerfac = 0
     DO j = 1, n_x
       fac =abs ( y(i) - x(j) )
        IF ( fac < 1.d-10 ) THEN
             zerfac = j
             exit
        END IF
     END DO
     DO j = 1, n_x
       DO k = 1, j - 1
          p(i,j) = p(i,j) * ( y(i) - x(k) ) / ( x(j) - x(k) )
       END DO
       DO k = j + 1, n_x
          p(i,j) = p(i,j) * ( y(i) - x(k) ) / ( x(j) - x(k) )
       END DO
       IF ( PRESENT(dp) ) THEN
            IF ( abs(p(i,j)) > 1.d-10) THEN
                 sn(1:2) = 0.0d0
                 DO k = 1, j - 1
                    fac = 1.d0 / ( y(i) - x(k) )
                    sn(1) = sn(1) + fac
                    sn(2) = sn(2) + fac * fac
                 END DO
                 DO k = j + 1, n_x
                   fac = 1.d0 / ( y(i) - x(k) )
                   sn(1) = sn(1) + fac
                   sn(2) = sn(2) + fac * fac
                 END DO
                 dp(i,j)  = sn(1) * p(i,j)
                 ddp(i,j)  = sn(1) * dp(i,j) - sn(2) * p(i,j)
            ELSE
                  first = j
                  second = zerfac
                  IF ( first > second ) THEN
                       first = zerfac
                       second = j
                  END IF
                  sn(1) = 1.d0
                  sn(2) = 0.d0
                  DO k = 1, first-1
                     fac = 1.d0 / ( x(j) - x(k) )
                     sn(1) = sn(1) * fac * ( y(i) - x(k) )
                     sn(2) = sn(2) + 1.d0 / ( y(i) - x(k) )
                  END DO
                  DO k = first + 1, second-1
                     fac = 1.d0 / ( x(j) - x(k) )
                     sn(1) = sn(1) * fac * ( y(i) - x(k) )
                     sn(2) = sn(2) + 1.d0 / ( y(i) - x(k) )
                  END DO
                  DO k = second + 1, n_x
                     fac = 1.d0 / ( x(j) - x(k) )
                     sn(1) = sn(1) * fac * ( y(i) - x(k) )
                     sn(2) = sn(2) + 1.d0 / ( y(i) - x(k) )
                  END DO
!               write(*,*) ' sums',sn(1), sn(2)
                  dp(i,j)  = sn(1) / ( x(j) - x(zerfac) )
                  ddp(i,j) =  2.d0 * sn(2) * dp(i,j)
            END IF
       END IF
     END DO
  END DO
  IF ( prnt_lagrange(1) ) THEN
      ! call Print_Matrix(type_real_matrix,p,n_y,n_x,title='LaGrange Polynomials')
  END IF
  IF ( prnt_lagrange(2) ) THEN
      IF ( PRESENT(dp) ) THEN
          !  call Print_Matrix(type_real_matrix,dp,n_y,n_x,title='First Derivative of LaGrange Polynomials')
          !  call Print_Matrix(type_real_matrix,ddp,n_y,n_x,title='Second Derivative of LaGrange Polynomials')
      END IF
  END IF
!
!
  END SUBROUTINE Lagrange_Polynomials
!***********************************************************************
!********************************************************************************
!deck Lagrange_Integration_Weights
!***begin prologue     Integration_Weights
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           Interpolation_Points
!***author             schneider, b. i.(nsf)
!***source             Lagrange_Weights
!***purpose            Compute Lagrange integration weights
!***                   
!***                   
!***                   
!***                   
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Lagrange_Weights
  SUBROUTINE Lagrange_Integration_Weights(lgr,x,nc_wt,edge,step_size,prnt_lagrange)
  IMPLICIT NONE

  CLASS(LaGrange)                           :: lgr
  REAL(idp),  DIMENSION(:)                  :: x
  REAL(idp),  DIMENSION(:,:)                :: nc_wt
  REAL(idp),  DIMENSION(2)                  :: edge
  REAL(idp),  DIMENSION(2)                  :: edge_std
  REAL(idp),  DIMENSION(:),   ALLOCATABLE   :: q_std
  REAL(idp),  DIMENSION(:),   ALLOCATABLE   :: wt_std
  REAL(idp),  DIMENSION(:),   ALLOCATABLE   :: wt
  REAL(idp),  DIMENSION(:),   ALLOCATABLE   :: y
  REAL(idp),  DIMENSION(:,:), ALLOCATABLE   :: p
  REAL(idp)                                 :: step_size
  INTEGER                                   :: n_x
  INTEGER                                   :: n_y
  INTEGER                                   :: degree
  INTEGER                                   :: n_1
  INTEGER                                   :: n_2    
  INTEGER                                   :: i  
  INTEGER                                   :: int  
  INTEGER                                   :: pt
  INTEGER                                   :: gpt
  LOGICAL,    DIMENSION(:)                  :: prnt_lagrange
  Data edge_std / -1.d0, 1.d0 /
!
!
!
!                  Get  the LG points.
!
  n_x=size(x,1)
  degree = ( n_x - 1 ) / 2
  n_y = degree  + 1
  IF ( type_points == 'newton-cotes' ) THEN
       x(1) = edge(1)
       DO pt = 2, n_x
          x(pt) = x(pt-1) + step_size
       END DO
  ELSE
       ALLOCATE( wt(1:n_x) )      
       Call Gauss(q=x, wt=wt, edge=edge_std, type_quadrature='gauss', n=n_x)
       call cnvtpt(x,wt,edge,n_x)
       DEALLOCATE( wt )      
  END IF
!
  IF ( prnt_lagrange(3)  ) THEN
       write(*,2) x(:)
  END IF
!                  To get the weights we need to integrate the LG polynomials exactly over sub-intervals.
!                  Use a Gauss quadrature of proper size on each sub-interval to compute it exactly.
  
  ALLOCATE( q_std(1:n_y), wt_std(1:n_y) )
  Call Gauss(q=q_std, wt=wt_std, edge=edge_std, type_quadrature='gauss', n=n_y,    &
                      print_points_weights=prnt_lagrange(1))
  IF (  prnt_lagrange(4) ) THEN
     write (*,3) n_y
     write (*,4) ( q_std(i), wt_std(i), i = 1, n_y )     
  END IF
  Allocate( y(1:n_y), wt(1:n_y), p(1:n_y,1:n_x) )
  edge(1) = x(1)
  DO int = 1, n_x - 1
     edge(2) = x(int+1)           
!
!    Convert the standard interval to the actual one.
!
     y(1:n_y)  = q_std(1:n_y)
     wt(1:n_y) = wt_std(1:n_y)
     call cnvtpt(y,wt,edge,n_y)
     IF (  prnt_lagrange(4) ) THEN
          write(*,5) int 
          write(*,4) ( y(i), wt(i), i = 1,n_y)
     END IF
!     Call LaGrange_Polynomials(p=p,x=x,y=y,n_x=n_x,n_y=n_y)
     Call lgr%LaGrange_Polynomials(p=p,x=x,y=y,prnt_lagrange=prnt_lagrange(1:))
     nc_wt(1:n_x,int) = 0.d0
     DO pt = 1, n_x
        DO gpt = 1, n_y
           nc_wt(pt,int) = nc_wt(pt,int) + p(gpt,pt) * wt(gpt)
        END DO
     END DO              
     edge(1) = edge(2)
  END DO
  Deallocate( q_std, wt_std, y, wt, p )
  IF (  prnt_lagrange(5) ) THEN 
      ! call Print_Matrix(type_real_matrix,nc_wt,n_x,n_x-1,title='LaGrange Weights')
  END IF
!
!
1 Format(/10x,'Number of Points = ',i3,2x,'Lagange Step Size = ', e15.8)
2 Format(/10x,'Points = ',5e15.8)
3 Format(/10x,'Size of Legendre Quadrature = ',i3,/,25x,'Standard (Point Weight) = ' )
4 Format(50x,'(',e15.8,')','(',e15.8,')' )
5 Format(/,25x,'Interval = ',i3,1x,'(Point Weight) = ' )  
END SUBROUTINE Lagrange_Integration_Weights
!***********************************************************************
!***********************************************************************
!deck Lagrange_Integration_Weights
!***begin prologue     Integration_Weights
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           Interpolation_Points
!***author             schneider, b. i.(nsf)
!***source             Lagrange_Weights
!***purpose            Compute Lagrange integration weights
!***                   
!***                   
!***                   
!***                   
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Lagrange_Weights
  SUBROUTINE Lagrange_Integration_Weights_2(lgr,x,z,nc_wt,prnt_lagrange)
  IMPLICIT NONE
!
  CLASS(LaGrange)                           :: lgr
  REAL(idp),  DIMENSION(:)                  :: x
  REAL(idp),  DIMENSION(:)                  :: z
  REAL(idp),  DIMENSION(:,:)                :: nc_wt
  REAL(idp),  DIMENSION(2)                  :: edge_std
  REAL(idp),  DIMENSION(:),   ALLOCATABLE   :: q_std
  REAL(idp),  DIMENSION(:),   ALLOCATABLE   :: wt_std
  REAL(idp),  DIMENSION(:),   ALLOCATABLE   :: y_i
  REAL(idp),  DIMENSION(:),   ALLOCATABLE   :: wt
  REAL(idp),  DIMENSION(:,:), ALLOCATABLE   :: p
  INTEGER                                   :: n_x
  INTEGER                                   :: n_z
  INTEGER                                   :: ny_i
  INTEGER                                   :: degree
  INTEGER                                   :: n_1
  INTEGER                                   :: n_2    
  INTEGER                                   :: i  
  INTEGER                                   :: int  
  INTEGER                                   :: pt
  INTEGER                                   :: gpt
  LOGICAL,    DIMENSION(:)                  :: prnt_lagrange
  Data edge_std / -1.d0, 1.d0 /
!
!
!
!                  Get  the LG points.
!
  n_x=size(x,1)   ! LaGrange pivot points
  n_z=size(z,1)   ! Integration grid
!
!  IF ( prnt_lagrange(3)  ) THEN
!       write(iout,1) n_x
!       write(iout,2) x(:)
!       write(iout,3) n_z
!       write(iout,4) z(:)
!  END IF
  degree = ( n_x - 1 ) / 2  ! Need a Gauss grid of ny_i points to integrate a polynomial of degree n_x
  ny_i = degree  + 1
!                             To get the weights we need to integrate the LG polynomials exactly over sub-intervals.
!                             Use a Gauss quadrature of proper size on each sub-interval to compute it exactly.
  ALLOCATE( q_std(1:ny_i), wt_std(1:ny_i) )
  Call Gauss(q=q_std, wt=wt_std, edge=edge_std, type_quadrature='gauss', n=ny_i,    &
                      print_points_weights=prnt_lagrange(1))
!  IF (  prnt_lagrange(4) ) THEN
!     write (iout,5) ny_i
!     write (iout,6) ( q_std(i), wt_std(i), i = 1, ny_i )     
!  END IF
  Allocate( y_i(1:ny_i), wt(1:ny_i), p(1:ny_i,1:n_x) )
  DO int = 1, n_z - 1
!
!    Convert the standard interval to the integration grid interval.
!
     y_i(1:ny_i)  = q_std(1:ny_i)
     wt(1:ny_i) = wt_std(1:ny_i)
     call cnvtpt(y_i,wt,z(int:int+1),ny_i)
!     IF (  prnt_lagrange(4) ) THEN
!          write(iout,5) int 
!          write(iout,6) ( y_i(i), wt(i), i = 1,ny_i)
!     END IF
     Call lgr%LaGrange_Polynomials(p=p,x=x,y=y_i,prnt_lagrange=prnt_lagrange(1:))  ! Compute polynomials at Gauss points
     nc_wt(1:n_x,int) = 0    ! Compute the integration weight for the integration grid.
     DO pt = 1, n_x
        DO gpt = 1, ny_i
           nc_wt(pt,int) = nc_wt(pt,int) + p(gpt,pt) * wt(gpt)
        END DO
     END DO              
  END DO

  Deallocate( q_std, wt_std, y_i, wt, p )
!  IF ( prnt_lagrange(5) ) THEN 
!      call Print_Matrix(type_real_matrix,nc_wt,n_x,n_x-1,title='LaGrange Weights for Each Subinterval')
!  END IF
!
!

1 Format(/10x,'Number of Pivot Points = ',i4)
2 Format(/10x,'Pivot Points = ',5e15.8)
3 Format(/10x,'Number of Integration Points = ',i4)
4 Format(/,10x,'Integration Points = ',5e15.8)
5 Format(/10x,'Number of Points to Integrate Lagrange Polynomial = ',i4,/,25x,'Standard (Point Weight) = ' )
6 Format(50x,'(',e15.8,')','(',e15.8,')' )
7 Format(/,25x,'Integration Interval = ',i4,1x,'(Point Weight) = ' )  
END SUBROUTINE Lagrange_Integration_Weights_2
!***********************************************************************
!***********************************************************************
           END MODULE LaGrange_Module
!***********************************************************************
!***********************************************************************

