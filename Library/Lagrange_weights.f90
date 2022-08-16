module Lagrange_weights
use precin
use LaGrange_Module
use Gauss_Quadrature

implicit none

private
public lgr_weights
contains

subroutine lgr_weights(left_edge,right_edge,x,nc_wt,n,quad_type)
    TYPE(LaGrange)                                            :: lgr
    REAL(idp),        INTENT(IN)                              :: left_edge, right_edge
    CHARACTER(len=*), INTENT(IN)                              :: quad_type
    REAL(idp)                                                 :: step_size
    REAL(idp)                                                 :: edge(2)
    REAL(idp),   DIMENSION(2)                                 :: edge_std
    REAL(idp),   DIMENSION(:,:), ALLOCATABLE                  :: p
    REAL(idp),   DIMENSION(:,:), ALLOCATABLE                  :: dp
    REAL(idp),   DIMENSION(:,:), ALLOCATABLE                  :: ddp
    REAL(idp),   DIMENSION(n+2,n+1), INTENT(OUT)              :: nc_wt
    REAL(idp),   DIMENSION(n+2),   INTENT(OUT)                :: x
    INTEGER  , INTENT(IN)                                     :: n
    INTEGER  :: i, n_x
    Data edge_std / -1.d0, 1.d0 /

    edge = [left_edge, right_edge]
    prnt_lagrange = .false.
    n_x = n + 2
    step_size = (edge(2) - edge(1))/(n+1)

    select case(trim(quad_type))

    case('newton_cotes')
       x(1) = edge(1)

       DO i=2,n_x
          x(i) = x(i-1) + step_size
       End DO

    case('gauss')
       call gauss(q=x, wt=nc_wt(:,1), edge=edge_std, type_quadrature='lobatto', n=n_x)
       call cnvtpt(x, nc_wt(:,1), edge, n_x)

    end select
       
!    write(*,*) n_x, edge, step_size
  !  write(*,*) x(0:)
  !  write(*,*) x(:)
!    write(*,*)  
!    write(*,*) 'Computing Polynomials'
!    write(*,*)  
!    call lgr%LaGrange_Polynomials(p=p,dp=dp,ddp=ddp,x=x,y=x,prnt_lagrange= prnt_lagrange)
!    write(*,*)  
!    write(*,*) 'Computing Weights'
!    write(*,*)  
     call lgr%Lagrange_Integration_Weights_2(x, x, nc_wt, prnt_lagrange)
    
end subroutine

end module
