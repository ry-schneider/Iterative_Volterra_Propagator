program test_LagrnLib
    use precin
    use Lagrange_weights
    implicit NONE
    
    real(idp), allocatable :: integrand(:), integral(:,:)
    real(idp), allocatable :: x(:)
    real(idp), allocatable :: nc_wt(:,:)
    real(idp)              :: edge(2), original_edges(2),exponent

    integer                :: n_x, int, pt

    print*, 'test Lagrange Module Call'

    !edges:
    original_edges = [0.d0,1.d0]
    edge = original_edges
    
    n_x = 50
    exponent = 1.d0

    call lgr_weights(edge,n_x,x,nc_wt)
    allocate(integrand,source=x)
    allocate(integral(n_x,2))
    integral = 0.d0

    Integrand(:) = exp(exponent * x(:) )
    edge(1) = x(1)
    DO int = 1, n_x
       edge(2) = x(int+1)
       integral(int,2) = integral(int,2)   + ( exp(exponent*edge(2))  - exp(exponent*edge(1)) ) /exponent
       DO pt = 1, size(x)
          integral(int,1) = integral(int,1) + nc_wt(pt,int) * Integrand(pt) 
       END DO
       write(*,2) integral(int,1),integral(int,2), edge(:)

            edge(1) = edge(2)

    END DO
write (*,*) "******************Final comparison****************************"

 write(*,2) sum(integral(:,1)), sum(integral(:,2)), original_edges

 2 Format(/10x,'Numerical Integral = ',e15.8,2x,'Exact Integral = ',e15.8,2x,' on interval (',e15.8,',',e15.8,')' )

print*, "diff. Error:" ,abs(sum(integral(:,1))- sum(integral(:,2)))
    
end    

! Program Test_Lagrange
!     USE LaGrange_Module
!     TYPE(LaGrange)                                :: lgr
!     INTEGER                                       :: intkey
!     CHARACTER(LEN=80)                             :: chrkey
!     LOGICAL                                       :: logkey
!     REAL(idp)                                     :: fpkey
!     REAL(idp),      DIMENSION(2)                  :: edge
!     REAL(idp)                                     :: step_size
!     CHARACTER(LEN=5)                              :: itoc
!     REAL(idp),      DIMENSION(:),   ALLOCATABLE   :: x
!     REAL(idp),      DIMENSION(:,:), ALLOCATABLE   :: p
!     REAL(idp),      DIMENSION(:,:), ALLOCATABLE   :: dp
!     REAL(idp),      DIMENSION(:,:), ALLOCATABLE   :: ddp
!     REAL(idp),      DIMENSION(:,:), ALLOCATABLE   :: nc_wt
!     ! Call Drum
!     ! Call Get_Environment
!     ! call posinp('$lagrange',card)
!     ! call cardin(card)
!     ! call fparr(card,'region_edges',edge,2,' ') 
!     edge(1) = 0.d0
!     edge(2) = 1.d0
!     ! prnt_lagrange=.true.
!     ! n_x = intkey(card,'number_of_points',3,' ')
!     n_x = 3
!     step_size= ( edge(2) - edge(1) ) / n_x
!     n_x = n_x + 1
!     ALLOCATE( x(1:n_x), p(1:n_x,1:n_x), dp(1:n_x,1:n_x), ddp(1:n_x,1:n_x), nc_wt(1:n_x,1:n_x-1) )
!   !  ALLOCATE( x(0:n_x), p(0:n_x,0:n_x), dp(0:n_x,0:n_x), ddp(0:n_x,0:n_x) )
!     x(1) = edge(1)
!   !  x(0) = edge(1)
!     DO i =2, n_x
!   !  DO i =1, n_x
!        x(i) = x(i-1) + step_size
!     End DO
!     write(*,*) n_x, edge, step_size
!   !  write(*,*) x(0:)
!   !  write(*,*) x(:)
!     write(*,*)  
!     write(*,*) 'Computing Polynomials'
!     write(*,*)  
!     call lgr%LaGrange_Polynomials(p=p,dp=dp,ddp=ddp,x=x,y=x,prnt_lagrange=prnt_lagrange)
!     write(*,*)  
!     write(*,*) 'Computing Weights'
!     write(*,*)  
!     call lgr%Lagrange_Integration_Weights(x=x,nc_wt=nc_wt,edge=edge,step_size=step_size,prnt_lagrange=prnt_lagrange)
!     DEALLOCATE( x, p, dp, ddp, nc_wt )
!     n_x = n_x - 1
!     ALLOCATE( x(0:n_x), p(0:n_x,0:n_x), dp(0:n_x,0:n_x), ddp(0:n_x,0:n_x), nc_wt(0:n_x,0:n_x-1) )
!   !  n_x=n_x+1
!      write(*,*) 'Computing Polynomials'
!     call lgr%LaGrange_Polynomials(p=p(0:,0:),dp=dp(0:,0:),ddp=ddp(0:,0:),x=x(0:),y=x(0:),prnt_lagrange=prnt_lagrange)
!   !  
!     write(*,*)  
!     write(*,*) 'Computing Weights'
!     write(*,*)  
!     call lgr%Lagrange_Integration_Weights(x=x(0:),nc_wt=nc_wt(0:,0:),edge=edge,step_size=step_size,prnt_lagrange=prnt_lagrange)
!     print*, nc_wt
!   !  allocate(x(0:n_x),nc_wt(0:n_x,0:n_x-1))
!   !  Call chainx(0)
!   !  call exit
! end program Test_LaGrange
  
