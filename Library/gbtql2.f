*deck @(#)gbtql2.f	1.1 9/9/91
      subroutine gbtql2(n, d, e, z, ierr)
c
c     this subroutine is a translation of the algol procedure imtql2,
c     num. math. 12, 377-383(1968) by martin and wilkinson,
c     as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c
c     this subroutine finds the eigenvalues and first components of the
c     eigenvectors of a symmetric tridiagonal matrix by the implicit ql
c     method, and is adapted from the eispak routine imtql2
c
c     on input=
c
c        n is the order of the matrix;
c
c        d contains the diagonal elements of the input matrix;
c
c        e contains the subdiagonal elements of the input matrix
c          in its first n-1 positions.  e(n) is arbitrary;
c
c        z contains the first row of the identity matrix.
c
c      on output=
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1, 2, ..., ierr-1;
c
c        e has been destroyed;
c
c        z contains the first components of the orthonormal eigenvectors
c          of the symmetric tridiagonal matrix.  if an error exit is
c          made, z contains the eigenvectors associated with the stored
c          eigenvalues;
c
c        ierr is set to
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     ------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      integer i, j, k, l, m, n, ii, mml, ierr
      real*8  machep
      dimension  d(n),e(n),z(n)
c
c     ========== machep is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c                machep = 16.0d0**(-13) for long form arithmetic
c                on s360 ==========
       machep=1.0e-14
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      e(n) = 0.0d0
      do 240 l = 1, n
         j = 0
c     ========== look for small sub-diagonal element ==========
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            if ( abs(e(m)) .le. machep * ( abs(d(m)) +  abs(d(m+1))))
     x         go to 120
  110    continue
c
  120    p = d(l)
         if (m .eq. l) go to 240
         if (j .eq. 30) go to 1000
         j = j + 1
c     ========== form shift ==========
         g = (d(l+1) - p) / (2.0d0   * e(l))
         r =  sqrt(g*g+1.0d0  )
         g = d(m) - p + e(l) / (g +  sign(r, g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
c     ========== for i=m-1 step -1 until l do -- ==========
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            if ( abs(f) .lt.  abs(g)) go to 150
            c = g / f
            r =  sqrt(c*c+1.0d0  )
            e(i+1) = f * r
            s = 1.0d0   / r
            c = c * s
            go to 160
  150       s = f / g
            r =  sqrt(s*s+1.0d0  )
            e(i+1) = g * r
            c = 1.0d0   / r
            s = s * c
  160       g = d(i+1) - p
            r = (d(i) - g) * s + 2.0d0   * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
c     ========== form first component of vector ==========
            f = z(i+1)
            z(i+1) = s * z(i) + c * f
            z(i) = c * z(i) - s * f
c
  200    continue
c
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0d0
         go to 105
  240 continue
c     ========== order eigenvalues and eigenvectors ==========
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         p = z(i)
         z(i) = z(k)
         z(k) = p
c
  300 continue
c
      go to 1001
c     ========== set error -- no convergence to an
c                eigenvalue after 30 iterations ==========
 1000 ierr = l
 1001 return
c     ========== last card of gbtql2 ==========
      end

