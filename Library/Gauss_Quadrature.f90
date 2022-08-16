      module Gauss_Quadrature
         use precin

         implicit none

      contains

      SUBROUTINE gauss(q,wt,edge,p,dp,ddp,type_quadrature,fixed_point,n,   &
                        print_points_weights,print_polynomials)
      IMPLICIT NONE
      REAL(idp), DIMENSION(:)                :: q
      REAL(idp), DIMENSION(:)                :: wt
      REAL(idp), DIMENSION(:)                :: edge
      REAL(idp), OPTIONAL, DIMENSION(:,:)    :: p
      REAL(idp), OPTIONAL, DIMENSION(:,:)    :: dp
      REAL(idp), OPTIONAL, DIMENSION(:,:)    :: ddp
      REAL(idp), DIMENSION(:), ALLOCATABLE   :: b
      REAL(idp)                              :: h
      CHARACTER(LEN=*)                       :: type_quadrature
      LOGICAL, DIMENSION(2), OPTIONAL        :: fixed_point
      INTEGER                                :: n
      INTEGER                                :: i
      REAL(idp), DIMENSION(2)                :: endpts
      REAL(idp), DIMENSION(2)                :: ptfix
      LOGICAL, OPTIONAL                      :: print_points_weights
      LOGICAL, OPTIONAL                      :: print_polynomials
      DATA ptfix / -1.d0, 1.d0 /
      !
      !
      ! If the weight function is a one of the classical weight functions the
      ! points and weights are known analytically and after computing them we
      ! go directly to getting the coordinate functions.
      !
      IF (type_quadrature == 'trapezoidal' ) THEN
            wt(1) = ( edge(2) - edge(1) ) / ( n - 1)
            wt(2:n) = wt(1)
            wt(1) = 0.5_idp * wt(1)
            wt(n) = 0.5_idp * wt(n)
            q(1) = edge(1)
            DO i = 2, n
               q(i) = q(i-1) + wt(2)
            END DO
      ELSE IF ( type_quadrature == 'simpson' ) THEN
            IF ( 2*(n/2) == n ) THEN
               Stop 'n must be odd for simpsons rule'
            END IF
            IF(n <= 1) THEN
               q(1) = edge(1)
               wt(1) = 2.d+00
               return
            END IF
            DO i = 2, n, 2
               wt(i) = 4_idp
            ENd DO
            DO i = 3, n, 2
               wt(i) = 2_idp
            END DO
            wt(1) = 1_idp
            wt(n) = 1_idp
            h = (edge(2) - edge(1))/(n-1)
            q(1) = edge(1)
            DO i = 2, n
               q(i) = q(i-1) + h
            END DO
            h = h/3_idp
            wt(1:n) = h * wt(1:n)
      ELSE           
            endpts(:)=edge(:)
            ALLOCATE( b(1:n) )
            IF (type_quadrature == "gauss") THEN
               CALL gaussq('one',n,0.d0,0.d0,0,ptfix,b,q,wt)
            ELSE IF ( type_quadrature == "radau") THEN
               IF (fixed_point(1) .eqv. .true.) THEN
                  ptfix(1) = -1.d0
               ELSE IF ( fixed_point(2) .eqv. .true.) THEN
                  ptfix(1) = 1.d0
               END IF
               CALL gaussq('one',n,0.d0,0.d0,1,ptfix,b,q,wt)
            ELSE IF ( type_quadrature == "lobatto" ) THEN
               ptfix(1) = -1.d0
               ptfix(2) = 1.d0
               CALL gaussq('one',n,0.d0,0.d0,2,ptfix,b,q,wt)
            END IF
            CALL cnvtpt(q,wt,edge,n)
            DEALLOCATE(b)
      END IF
      IF(PRESENT(print_points_weights) .eqv. .true.) THEN
          IF ( print_points_weights .eqv. .true. ) THEN
         !       call Print_Matrix(type_real_vector,q,title='Final Nodes from Gauss')
         !       call Print_Matrix(type_real_vector,wt,title='Final Weights from Gauss')
         END IF
      END IF
      !  
      !  
      IF ( PRESENT (p) .eqv. .true.) THEN
            Call cpoly(p,dp,ddp,q,n,print_polynomials)

      ! The DVR library assumes that the polynomials are $\delta$
      ! functions at the quadrature points.  Convert to this normalization
      !
      !
         DO i=1,n
            h = 1_idp/sqrt(wt(i))     
            p(:,i) = h * p(:,i)
            dp(:,i) = h * dp(:,i)
            ddp(:,i) = h * ddp(:,i)
         END DO
         IF(PRESENT(print_polynomials) .eqv. .true.) THEN
            IF(print_polynomials .eqv. .true.) THEN       
               ! call Print_Matrix(type_real_matrix,p,n,n,title='Normalized Coordinate Function')
               ! call Print_Matrix(type_real_matrix,dp,n,n,title=                                  &
               !                   'Normalized First Derivative of Coordinate Function')
               ! call Print_Matrix(type_real_matrix,ddp,n,n,title=                                 &
               !                   'Normalized Second Derivative of Coordinate Function')
            END IF
         END IF
      END IF
      END SUBROUTINE gauss

      !deck cpoly.f
!***begin prologue     cpoly
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            generate coordinate functions.
!***
!***description
!***references
!***routines called
!***end prologue       cpoly

  SUBROUTINE cpoly(cp,dcp,ddcp,pt,n,print)
   IMPLICIT NONE
   INTEGER                             :: n
   REAL(idp), DIMENSION(:,:)           :: cp
   REAL(idp), DIMENSION(:,:)           :: dcp
   REAL(idp), DIMENSION(:,:)           :: ddcp
   REAL(idp), DIMENSION(:)             :: pt
   LOGICAL, OPTIONAL                   :: print
   CHARACTER (LEN=80)                  :: title
   CALL lgngr(cp,dcp,ddcp,pt,pt,n,n,drctv='on')
   ! IF(PRESENT(print) == .true.) THEN
   !    ! call Print_Matrix(type_real_matrix,cp,n,n,title='Coordinate Function')
   !    ! call Print_Matrix(type_real_matrix,dcp,n,n,title='First Derivative of Coordinate Function')
   !    ! call Print_Matrix(type_real_matrix,ddcp,n,n,title='Second Derivative of Coordinate Function')
   !  END IF
 END SUBROUTINE cpoly

 !***********************************************************************
!*deck lgngr
 SUBROUTINE LGNGR(p,dp,ddp,x,y,nx,ny,type,drctv,print) 
   !***begin prologue     lgngr
   !***date written       940504   (yymmdd)
   !***revision date               (yymmdd)
   !***keywords
   !***author             schneider, barry (nsf)
   !***source             %W% %G% 
   !***purpose            lagrange polynomials at arbitrary points.
   !***description
   !***            
   !               
   !               
   !***references
   !
   !***routines called
   !
   !***end prologue       lgngr
   !
     IMPLICIT NONE
     REAL*8, DIMENSION (ny,nx)           :: p
     REAL*8, DIMENSION (ny,nx)           :: dp
     REAL*8, DIMENSION (ny,nx)           :: ddp
     REAL*8, DIMENSION(nx)               :: x
     REAL*8, DIMENSION(ny)               :: y
     REAL*8, DIMENSION(:), ALLOCATABLE   :: xt
     REAL*8, DIMENSION(:), ALLOCATABLE   :: yt
     REAL*8                              :: sn
     REAL*8                              :: ssn
     REAL*8                              :: fac
     LOGICAL, DIMENSION(:), OPTIONAL     :: print
     CHARACTER (LEN = 80)                :: title 
     CHARACTER (LEN = *), OPTIONAL       :: drctv
     CHARACTER (LEN = *), OPTIONAL       :: type
     INTEGER                             :: nx
     INTEGER                             :: ny
     INTEGER                             :: i
     INTEGER                             :: j
     INTEGER                             :: k
     INTEGER                             :: first
     INTEGER                             :: second
     INTEGER                             :: zerfac
     INTEGER                             :: inp

   !
   !     generate polynomials and derivatives with respect to x
   !
     p(:,:) = 1.d0
     IF (present(type) ) THEN 
         ALLOCATE(xt(nx),yt(ny))
         xt(:) = x(:)
         yt(:) = y(:)
         x(:) = x(:) * x(:)
         y(:) = y(:) * y(:)
     END IF
     DO i=1,ny
        zerfac = 0
        DO j=1,nx
           fac =  y(i) - x(j) 
           IF(abs(fac) <= 1.d-10) THEN
              zerfac = j
           ENDIF  
        END DO
        DO j=1,nx
           DO k = 1, j-1
              p(i,j) = p(i,j) * ( y(i) - x(k) )   &
                              / ( x(j) - x(k) )
           END DO
           DO k=j+1,nx
              p(i,j) = p(i,j) * ( y(i) - x(k) )   &
                              / ( x(j) - x(k) )
           END DO
           IF(present(drctv) ) THEN
               IF ( abs(p(i,j)) > 1.d-10) THEN
                    sn = 0.d0
                    ssn = 0.d0
                    DO k=1,j-1
                       fac = 1.d0/( y(i) - x(k) )
                       sn = sn + fac
                       ssn = ssn + fac*fac
                    END DO
                    DO k=j+1,nx
                       fac = 1.d0/( y(i) - x(k) )
                       sn = sn + fac
                       ssn = ssn + fac*fac
                    END DO                                 
                    dp(i,j) = sn*p(i,j)               
                    ddp(i,j) = sn*dp(i,j) - ssn*p(i,j)
               ELSE
                    first=j
                    second=zerfac
                    IF(first > second) THEN
                       first=zerfac
                       second=j
                    END IF
                    sn = 1.d0
                    ssn = 0.d0
                    DO k=1,first-1
                       fac = 1.d0/( x(j) - x(k) )
                       sn = sn*fac*( y(i) - x(k) )
                       ssn = ssn + 1.d0/(y(i) - x(k))
                    END DO
                    DO k=first+1,second-1
                       fac = 1.d0/( x(j) - x(k) )
                       sn = sn*fac*( y(i) - x(k) )
                       ssn = ssn + 1.d0/( y(i) - x(k) )             
                    END DO
                    DO k=second+1,nx
                       fac = 1.d0/( x(j) - x(k) )
                       sn = sn*fac*( y(i) - x(k) )
                       ssn = ssn + 1.d0/( y(i) - x(k) )             
                    END DO
                    dp(i,j) = sn/( x(j) - x(zerfac) )
                    ddp(i,j) = 2.d0*ssn*dp(i,j)
               END IF                    
           END IF
   !
        END DO
     END DO
   !
     IF (present(type)) THEN 
         DO i=1,ny
            ddp(i,:) = 2.d0*dp(i,:) + 4.d0 * yt(i) * yt(i) * ddp(i,:) 
            dp(i,:) = 2.d0 * yt(i) * dp(i,:)
         END DO
         x(:) = xt(:)
         y(:) = yt(:)
         DEALLOCATE(xt,yt)
   !
     END IF
     IF(present(print)) THEN
      !   call Print_Matrix(type_real_matrix,p,ny,nx,title='Polynomials')
        IF(present(drctv)) THEN
         !   call Print_Matrix(type_real_matrix,dp,ny,nx,title='First Derivative of Polynomials')
         !   call Print_Matrix(type_real_matrix,ddp,ny,nx,title='Second Derivative of Polynomials')
        END IF
     END IF
     END SUBROUTINE Lgngr

   end module
