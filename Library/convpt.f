      !deck cnvtpt.f
      SUBROUTINE cnvtpt(pt,wt,endpts,n)
        !use precin
        IMPLICIT NONE
        INTEGER                                :: n
        REAL(8), DIMENSION(n)                :: pt
        REAL(8), DIMENSION(n)                :: wt
        REAL(8), DIMENSION(2)                :: endpts
        REAL(8)                              :: f1
        REAL(8)                              :: f2
        f1 = ( endpts(2)-endpts(1) )*.5D0
        f2 = ( endpts(1) + endpts(2) )*.5D0
        pt =  f1*pt + f2
        wt = wt*f1
      END SUBROUTINE cnvtpt

