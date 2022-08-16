*deck @(#)gbslve.f	1.1 9/9/91
      function gbslve(shift, n, a, b)
c
c       this procedure performs elimination to solve for the
c       n-th component of the solution delta to the equation
c
c             (jn - shift*identity) * delta  = en,
c
c       where en is the vector of all zeroes except for 1 in
c       the n-th position.
c
c       the matrix jn is symmetric tridiagonal, with diagonal
c       elements a(i), off-diagonal elements b(i).  this equation
c       must be solved to obtain the appropriate changes in the lower
c       2 by 2 submatrix of coefficients for orthogonal polynomials.
c
c
      implicit real*8 (a-h,o-z)
      dimension  a(n),b(n)
c
      alpha = a(1) - shift
      nm1 = n - 1
      do 10 i = 2, nm1
        alpha = a(i) - shift - b(i-1)**2/alpha
10    continue  
      gbslve = 1.0d0  /alpha
      return
      end

