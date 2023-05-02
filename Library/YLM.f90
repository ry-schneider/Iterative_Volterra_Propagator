!***********************************************************************                
!
                     MODULE Spherical_Harmonics
!***************************************************************************
!***************************************************************************
  USE iso_fortran_env, only:  real64
  IMPLICIT NONE
  INTEGER,     PARAMETER                                  :: dp= real64
  REAL(dp),    PARAMETER                                  :: pi      = atan(1.0d0)*4                                 
  REAL(dp),    PARAMETER                                  :: two_pi  = 2.d0*pi                                   
  REAL(dp),    PARAMETER                                  :: four_pi = 4.d0*pi       
  INTEGER                                                 :: iout=6
  REAL(dp),             PARAMETER                         :: inv_sqrt_2=sqrt(1.d0/2.d0)                             
  REAL(dp)                                                :: norm_0 = sqrt(1.d0/two_pi)                             
  REAL(dp)                                                :: norm_m = sqrt(1.d0/pi)  
  REAL(dp),    DIMENSION(:),     ALLOCATABLE              :: scratch
  REAL(dp),    DIMENSION(10)                              :: a_lm
                              Contains
!***************************************************************************                
!***********************************************************************
!deck Upward_Recursion
!***begin prologue     Upward_Recursion
!***date written       20140706  (yyyymmdd)                                                
!***revision date                (yyyymmdd)                                                
!***keywords           Initialize
!***author             schneider, b. i.(nist)                                            
!***source                                                                        
!***purpose      Compute the UN-NORMALIZED or normalized P_lm.             
!***             A direct calculation of the normalized functions is a bit more      
!***             complex but avoides the calculation of factorials which become
!***             unstable for large arguments.
!***description  This uses the standard approach but modified to incude the normalization
!***             factors.  It must be called separately for each m value so there is a minor      
!***             loss of efficiency in the computation of the starting values.      
!***                   
!***                   
!***references                                                                          
!***routines called                                                                     
!***end prologue       Upward_Recursion
  Subroutine Upward_Recursion ( p_lm, x, m, l_max, normalized )
  IMPLICIT NONE
!**********************************************************************
!**********************************************************************
  REAL(dp),   DIMENSION(:)                              :: x
  REAL(dp),   DIMENSION(m:,:)                           :: p_lm
  INTEGER                                               :: l_max
  INTEGER                                               :: m
  INTEGER                                               :: ipt
  INTEGER                                               :: l
  INTEGER                                               :: n
  INTEGER                                               :: two_m
  LOGICAL,    OPTIONAL                                  :: normalized
!
  write(iout,*) '                    P_LM for M = ',m  
  ALLOCATE(scratch(1:size(x,1)))
  scratch = sqrt( 1.d0 - x*x )
  two_m = m + m
  IF ( PRESENT(normalized) ) Then
!
!             The recursion for the normailzed P_lm was derived from the un-normalized version by
!             using the normalization values and then re-arranging all of the factors to relpace
!             all the factorials by recurrances  of previous values.  This approach is far more
!             stable numerically, avoiding the cancellation of large factors in the numerator and
!             denominator.  It was used in the CPC program published ~2010.
!
!          Initialize first value.                                                                   
!                                                                                                    
       p_lm(m,:) = inv_sqrt_2                                                                      
       a_lm(1) = 3                                                                               
       a_lm(2) = 2                                                                               
       DO l = 1, m                                                                               
          a_lm(3) = sqrt(a_lm(1)/a_lm(2))                                                        
          p_lm(m,:) = - a_lm(3) * scratch(:) * p_lm(m,:)                                              
          a_lm(1) = a_lm(1) + 2                                                                   
          a_lm(2) = a_lm(2) + 2                                                                   
       END DO                                                                                     
!                                                                                                     
!               Calculate second value.                                                             
!                                                                                                    
       IF (l_max /= m) THEN                                                                      
!                                                                                                    
!          Now calculate:  P_(M+1)M                                                       
!                                                                                                    
           a_lm(1) = two_m + 3                                                                    
           p_lm(m+1,:) = sqrt(a_lm(1)) * x(:) * p_lm(m,:)                                          
           a_lm(1) = a_lm(1) + 2                                                            
       END IF
       a_lm(1) = 2            ! Starting value for l-m+1 for l =  m + 2                       
       a_lm(2) = two_m + 2    ! Starting value for l+m+1 for l =  m + 2                       
       a_lm(3) = two_m + 1    ! Starting value for 2*l-1 for l =  m + 2                       
       a_lm(4) = two_m + 3    ! Starting value for 2*l+1 for l =  m + 2                       
       a_lm(5) = two_m + 5    ! Starting value for 2*l+3 for l =  m + 2                       
       a_lm(6) = two_m + 1    ! Starting value for l+m   for l =  m + 2                       
       a_lm(7) = 1            ! Starting value for l-m   for l =  m + 2                       
       DO l = m + 2, l_max                                                                    
          a_lm(8)  = a_lm(4) * a_lm(5) / ( a_lm(1) * a_lm(2) )                                    
          a_lm(8)  = sqrt ( a_lm(8) )                                                            
          a_lm(9)  =  a_lm(5) * a_lm(6) * a_lm(7) / ( a_lm(1) * a_lm (2) * a_lm(3) )             
          a_lm(9)  = sqrt ( a_lm(9) )                                                             
          p_lm(l,:) = a_lm(8) * x(:) * p_lm(l-1,:)                &               
                                       -                          &                
                      a_lm(9) * p_lm(l-2,:)                                    
          a_lm(1) = a_lm(1) + 1                                                               
          a_lm(2) = a_lm(2) + 1                                                               
          a_lm(3) = a_lm(3) + 2                                                               
          a_lm(4) = a_lm(4) + 2                                                               
          a_lm(5) = a_lm(5) + 2                                                               
          a_lm(6) = a_lm(6) + 1                                                               
          a_lm(7) = a_lm(7) + 1                                                               
       END DO
  ELSE
!                                                                                                   
!         The upward recursion.  Stable for all values of z                                          
!         This is just the standard reecursion in textbooks.                                          
!         Begin by computing th first two p_lm                                              
!       
          p_lm(m,:) = 1.d0
          a_lm(1) = 1
          DO l = 1, m
             p_lm(m,:) = - a_lm(1) * scratch(:) * p_lm(m,:)
             a_lm(1) = a_lm(1) + 2
          END DO
          IF ( l_max /= m ) Then
               p_lm(m+1,:) = a_lm(1) * x(:) * p_lm(m,:)
          END IF
!
!         Now recur
!
          a_lm(1) = two_m + 3                                                                
          a_lm(2) = two_m + 1                                                                  
          a_lm(3) = 2                                                                          
          DO l = m + 2, l_max                                                                   
             p_lm(l,:) = ( a_lm(1) * x(:) * p_lm(l - 1,:)         &             
                                   -                              &              
                           a_lm(2) * p_lm(l - 2,:) ) / a_lm(3)                      
             a_lm(1) = a_lm(1) + 2
             a_lm(2) = a_lm(2) + 1
             a_lm(3) = a_lm(3) + 1
          END DO                                                                                     
  END IF
  DEALLOCATE ( scratch )
!***********************************************************************
!***********************************************************************
  END Subroutine Upward_Recursion
!***********************************************************************
!***********************************************************************
!deck  Phi_Functions 
!***begin prologue     Phi_Functions 
!***date written       20140706  (yyyymmdd)                                                
!***revision date                (yyyymmdd)                                                
!***keywords           Initialize
!***author             schneider, b. i.(nist)                                            
!***source                                                                        
!***purpose      Compute the NORMALIZED spherical harmonics.             
!***             A direct calculation of the normalized functions is a bit more      
!***             complex but avoides the calculation of factorials which become
!***             unstable for large arguments.
!***description  This uses the standard approach but modified to incude the normalization
!***             factors.      
!***                   
!***                   
!***                   
!***references                                                                          
!***routines called                                                                     
!***end prologue       phi_functions
  Subroutine Phi_Functions ( phi_m, x, m )
  IMPLICIT NONE
!**********************************************************************
!**********************************************************************
  REAL(dp),   DIMENSION(:,:)                            :: phi_m
  REAL(dp),   DIMENSION(:)                              :: x
  INTEGER                                               :: m
  INTEGER                                               :: n
  INTEGER                                               :: ipt
  n= size(x,1)
  IF ( m == 0 ) Then
      phi_m(1,:) = norm_0
  ELSE
      phi_m(1,:) = norm_m * sin ( m * x(:) )  
      phi_m(2,:) = norm_m * cos ( m * x(:) )  
      DO ipt = 1, n
         write(iout,*) x(ipt), phi_m(1,ipt),  phi_m(2,ipt)
      END DO
  END IF
!***********************************************************************
!***********************************************************************
  END Subroutine Phi_Functions
!*****************************************************************************
!*****************************************************************************
  END MODULE Spherical_Harmonics
!***********************************************************************
!***********************************************************************

