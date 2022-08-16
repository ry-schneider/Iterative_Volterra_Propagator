*deck @(#)gaussq.f	1.1 9/9/91
      subroutine gaussq(ckind, n, alpha, beta, kpts, endpts, b, t, w)
c
c           this set of routines computes the nodes x(i) and weights
c        c(i) for gaussian-type quadrature rules with pre-assigned
c        nodes.  these are used when one wishes to approximate
c
c                 integral (from a to b)  f(x) w(x) dx
c
c                              n
c        by                   sum c  f(x )
c                             i=1  i    i
c
c        here w(x) is one of six possible non-negative weight
c        functions (listed below), and f(x) is the
c        function to be integrated.  gaussian quadrature is particularly
c        useful on infinite intervals (with appropriate weight
c        functions), since then other techniques often fail.
c
c           associated with each weight function w(x) is a set of
c        orthogonal polynomials.  the nodes x(i) are just the zeroes
c        of the proper n-th degree polynomial.
c
c     input parameters
c
c        kind     an integer between 0 and 6 giving the type of
c                 quadrature rule
c
c        kind = 0=  simpson's rule w(x) = 1 on (-1, 1) n must be odd.
c        kind = 1=  legendre quadrature, w(x) = 1 on (-1, 1)
c        kind = 2=  chebyshev quadrature of the first kind
c                   w(x) = 1/dsqrt(1 - x*x) on (-1, +1)
c        kind = 3=  chebyshev quadrature of the second kind
c                   w(x) = dsqrt(1 - x*x) on (-1, 1)
c        kind = 4=  hermite quadrature, w(x) = exp(-x*x) on
c                   (-infinity, +infinity)
c        kind = 5=  jacobi quadrature, w(x) = (1-x)**alpha * (1+x)**
c                   beta on (-1, 1), alpha, beta .gt. -1.
c                   note= kind=2 and 3 are a special case of this.
c        kind = 6=  generalized laguerre quadrature, w(x) = exp(-x)*
c                   x**alpha on (0, +infinity), alpha .gt. -1
c
c        n        the number of points used for the quadrature rule
c        alpha    real parameter used only for gauss-jacobi and gauss-
c                 laguerre quadrature (otherwise use 0.).
c        beta     real parameter used only for gauss-jacobi quadrature--
c                 (otherwise use 0.).
c        kpts     (integer) normally 0, unless the left or right end-
c                 point (or both) of the interval is required to be a
c                 node (this is called gauss-radau or gauss-lobatto
c                 quadrature).  then kpts is the number of fixed
c                 endpoints (1 or 2).
c        endpts   real array of length 2.  contains the values of
c                 any fixed endpoints, if kpts = 1 or 2.
c        b        real scratch array of length n
c
c     output parameters (both arrays of length n)
c
c        t        will contain the desired nodes x(1),,,x(n)
c        w        will contain the desired weights c(1),,,c(n)
c
c     subroutines required
c
c        gbslve, class, and gbtql2 are provided. underflow may sometimes
c        occur, but it is harmless if the underflow interrupts are
c        turned off as they are on this machine.
c
c     accuracy
c
c        the routine was tested up to n = 512 for legendre quadrature,
c        up to n = 136 for hermite, up to n = 68 for laguerre, and up
c        to n = 10 or 20 in other cases.  in all but two instances,
c        comparison with tables in ref. 3 showed 12 or more significant
c        digits of accuracy.  the two exceptions were the weights for
c        hermite and laguerre quadrature, where underflow caused some
c        very small weights to be set to zero.  this is, of course,
c        completely harmless.
c
c     method
c
c           the coefficients of the three-term recurrence relation
c        for the corresponding set of orthogonal polynomials are
c        used to form a symmetric tridiagonal matrix, whose
c        eigenvalues (determined by the implicit ql-method with
c        shifts) are just the desired nodes.  the first components of
c        the orthonormalized eigenvectors, when properly scaled,
c        yield the weights.  this technique is much faster than using a
c        root-finder to locate the zeroes of the orthogonal polynomial.
c        for further details, see ref. 1.  ref. 2 contains details of
c        gauss-radau and gauss-lobatto quadrature only.
c
c     references
c
c        1.  golub, g. h., and welsch, j. h.,  calculation of gaussian
c            quadrature rules,  mathematics of computation 23 (april,
c            1969), pp. 221-230.
c        2.  golub, g. h.,  some modified matrix eigenvalue problems,
c            siam review 15 (april, 1973), pp. 318-334 (section 7).
c        3.  stroud and secrest, gaussian quadrature formulas, prentice-
c            hall, englewood cliffs, n.j., 1966.
c
c     ..................................................................
c
      implicit real*8 (a-h,o-z)
      real*8  muzero
      character *(*) ckind
      dimension  b(n),t(n),w(n),endpts(2)
      
      if (ckind.eq.'simpson') then
          kind=0
      elseif(ckind.eq.'legendre'.or.
     1       ckind.eq.'one'.or.
     2       ckind.eq.'theta'.or.
     3       ckind.eq.'cylindrical'.or.
     4       ckind.eq.'spherical') then
          kind=1
       elseif(ckind.eq.'chebyshev-1') then
          kind=2
       elseif(ckind.eq.'chebyshev-2') then
          kind=3
       elseif(ckind.eq.'hermite') then
          kind=4
       elseif(ckind.eq.'jacobi') then
          kind=5
       elseif(ckind.eq.'laguerre') then
          kind=6
       else
          stop 'error in quadrature type'
       endif
       if(kind.eq.0) then
       if(2*(n/2).eq.n) then
        stop 'n must be odd for simpson rule'
      endif
        if(n.le.1) then
        t(1) = 0.d+00
        w(1) = 2.d+00
        return
      endif
       h = 2.d+00/(n-1)
       t(1) = -1.d+00
       t(n) = 1.d+00
       w(1) = h/3.d+00
       w(n) = h/3.d+00
       nm1 = n-1
       do 801 i=2,nm1
       t(i) = t(i-1) + h
       w(i) = 4.d+00 - 2.d+00*(i-2*(i/2))
       w(i) = w(i)*h/3.d+00
801      continue 
       return
      endif
c
      call class (kind, n, alpha, beta, b, t, muzero)
c
c           the matrix of coefficients is assumed to be symmetric.
c           the array t contains the diagonal elements, the array
c           b the off-diagonal elements.
c           make appropriate changes in the lower right 2 by 2
c           submatrix.
c
      if (kpts.eq.0)  go to 100
      if (kpts.eq.2)  go to  50
c
c           if kpts=1, only t(n) must be changed
c
      t(n) =gbslve(endpts(1), n, t, b)*b(n-1)**2 + endpts(1)
      go to 100
c
c           if kpts=2, t(n) and b(n-1) must be recomputed
c
50    gam =gbslve(endpts(1), n, t, b)
      t1 = ((endpts(1) - endpts(2))/(gbslve(endpts(2), n, t, b) - gam))
      b(n-1) =  sqrt(t1)
      t(n) = endpts(1) + gam*t1
c
c           note that the indices of the elements of b run from 1 to n-1
c           and thus the value of b(n) is arbitrary.
c           now compute the eigenvalues of the symmetric tridiagonal
c           matrix, which has been modified as necessary.
c           the method used is a ql-type method with origin shifting
c
100   w(1) = 1.0d0
      do 105 i = 2, n
         w(i) = 0.0d0
105   continue      
c
      call gbtql2 (n, t, b, w, ierr)
      do 110 i = 1, n
         w(i) = muzero * w(i) * w(i)
110   continue         
c
      return
      end







