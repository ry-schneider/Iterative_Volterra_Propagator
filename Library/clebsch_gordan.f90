module clebsch_gordan

  contains

!==================================================================== 
      Real*8 FUNCTION F_3J (J1,M1,J2,M2,J,M,clebsch)
!==================================================================== 
! 
!     determines the 3J or Clebsh-Gordon coefficients 
!     (the momenta are used in (2J+1)-representation) 
! 
!     Call:  F_3J 
! 
!     determines the value of the 3j-symbols without direct using of 
!     factorials. The following expression for the 3j-symbols is used: 
!         (A.P.JUCYS, A.A.BANDZAITIS, 1977) 
! 
!     3j{j1,m1,j2,m2,j3,m3} = delta(m1+m2,m3) * (2j3+1)^1/2 * {j1,j2,j3} * 
!       sqrt[ (j1+m1)!*(j1-m1)!*(j2+m2)!*(j2-m2)!*(j3+m3)!*(j3-m3)! ] 
!                         SUM(z) {   (-1)^z  / 
!          [ z! *  (j1+j2-j3-z)! * (j1-m1-z)! * (j2-m2-z)! * 
!                  (j3-j2+m1+z)! * (j3-j1-m2+z)! ] } 
! 
!     where {a,b,c}=sqrt[ (a+b-c)! * (a-b+c)! * (b-a+c)! / (a+b+c+1)! ] 
! 
!     If we introduce the auxiliary values a(i) and b(i) 
!     (see below the text of program) then 
! 
!     3j =         (-1) ^ Sum[a(i)] 
!          sqrt{ Pr[ (b(j)-a(i))! ] / [ Sum (b(j)-a(i))+1 ] } 
!                  Sum(z) { (-1)^z  / 
!          [  Pr[ (z-a(i))! ]  * Pr[ (b(j)-z)! ]   ] } 
! 
!     (below the moments are used in (2J+1)-representation) 
! 
!-------------------------------------------------------------------- 
!-------------------------------------------------------------------- 
      Implicit None 
      INTEGER, intent(in) :: J1, M1, J2, M2, J, M 
      INTEGER             :: J_a, M_a, J_b, M_b, J_c, M_c 
      LOGICAL             :: clebsch 
      ! Real*8, EXTERNAL    :: Z_3j 
      J_a = J1 + J1 + 1
      M_a = M1 + M1 + 1
      J_b = J2 + J2 + 1
      M_b = M2 + M2 + 1
      J_c = J + J + 1
      M_c = M + M + 1
      IF (clebsch) THEN
          F_3J = (-1)**((J_a-J_b+M_c-1)/2)*sqrt(DBLE(J_c))*   & 
                  Z_3J(J_a,M_a,J_b,M_b,J_c,-M_c+2) 
      ELSE
          F_3J = Z_3J(J_a,M_a,J_b,M_b,J_c,M_c)
      END IF
      END FUNCTION F_3J
!-------------------------------------------------------------------- 
!-------------------------------------------------------------------- 


!-------------------------------------------------------------------- 
!-------------------------------------------------------------------- 
      REAL*8 FUNCTION Z_3J(J1,M1,J2,M2,J3,M3)
      Implicit None 
      Integer(4), intent(in) :: j1,m1,j2,m2,j3,m3 
      Integer(4) :: i,i_max,k,kk,m,iz,iz_min,iz_max 
      Real(8) :: x,y,z 
      Integer(4) a(3),b(3),J(16) 
      Z_3j=0.0 
      IF(M1+M2+M3-3.ne.0) RETURN    ! check of conservation rules 
      J(1)= J1+J2-J3-1 
      J(2)= J1-J2+J3-1 
      J(3)= J2-J1+J3-1 
      J(4)= J1+M1-2 
      J(5)= J1-M1 
      J(6)= J2-M2 
      J(7)= J2+M2-2 
      J(8)= J3+M3-2 
      J(9)= J3-M3 
      Do I=1,9 
       IF(J(i).lt.0.or.mod(J(i),2).eq.1) RETURN 
      End do 
      a(1) = 0                         ! auxiliary values 
      a(2) = (j2-j3-m1+1)/2 
      a(3) = (j1-j3+m2-1)/2 
      b(1) = (j1+j2-j3-1)/2 
      b(2) = (j1-m1)/2 
      b(3) = (j2+m2-2)/2 
      IZ_min=MAX0(a(1),a(2),a(3))      ! limits of the sum 
      IZ_max=MIN0(b(1),b(2),b(3)) 
      IF(IZ_max.LT.IZ_min) Return 
      Do I=1,3                         ! constant factorial parameters 
      Do K=1,3 
       J(I+3*K-3)=b(i)-a(k) 
      End do 
      End do 
      J(10)=(j1+j2+j3-3)/2+1 
      Do I=1,3 
       J(I+10)=IZ_min-a(i)               ! initial factorial parameters 
       J(I+13)=b(i)-IZ_min               ! in the sum 
      End do 
      Z=0.0 
      DO IZ=IZ_min,IZ_max                 ! summation 
       I_max=0                            ! max. factorial 
       Do I=1,16 
        if(J(i).gt.I_max) I_max=J(i) 
       End do 
       Y=1.0 
       DO I=2,I_max         ! estimation of one term in sum 
        K=0                 ! K - the extent of the integer I in term 
        DO M=1,9 
         IF(J(M).GE.I) K=K+1 
        End do 
        IF(J(10).GE.I) K=K-1 
        DO M=11,16 
         IF(J(M).GE.I) K=K-2 
        End do 
        IF(K.EQ.0) Cycle 
        X=DBLE(I)                   ! Y = Y * I ** K/2 
        KK=IABS(K)/2 
        IF(KK.GT.0) THEN 
         DO M=1,KK 
          IF(K.GT.0) Y=Y*X 
          IF(K.LT.0) Y=Y/X 
         END DO 
        END IF 
        IF(mod(K,2).EQ.+1) Y=Y*SQRT(X) 
        IF(mod(K,2).EQ.-1) Y=Y/SQRT(X) 
       End do 
       IF(mod(IZ,2).eq.1) Y=-Y 
       Z=Z+Y 
       Do I=11,13                  ! new factorial parameters in sum 
        J(I)=J(I)+1 
       End do 
       DO I=14,16 
        J(I)=J(I)-1 
       End do 
      End do                       ! end of summation 
      K=a(1)+a(2)+a(3) 
      if(mod(k,2).ne.0) Z=-Z 
      Z_3j=Z 
      END FUNCTION Z_3j 
      
end module clebsch_gordan
