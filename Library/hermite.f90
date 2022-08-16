module hermite

  implicit none

  contains

subroutine h_polynomial_value ( m, n, x, p )

!*****************************************************************************80
!
!! H_POLYNOMIAL_VALUE evaluates H(i,x).
!
!  Discussion:
!
!    H(i,x) is the physicist's Hermite polynomial of degree I.
!
!  Differential equation:
!
!    Y'' - 2 X Y' + 2 N Y = 0
!
!  First terms:
!
!      1
!      2 X
!      4 X^2     -  2
!      8 X^3     - 12 X
!     16 X^4     - 48 X^2     + 12
!     32 X^5    - 160 X^3    + 120 X
!     64 X^6    - 480 X^4    + 720 X^2    - 120
!    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
!    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
!    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
!   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
!
!  Recursion:
!
!    H(0,X) = 1,
!    H(1,X) = 2*X,
!    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
!
!  Norm:
!
!    Integral ( -oo < X < oo ) exp ( - X^2 ) * H(N,X)^2 dX
!    = sqrt ( PI ) * 2^N * N!
!
!    H(N,X) = (-1)^N * exp ( X^2 ) * dn/dXn ( exp(-X^2 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Larry Andrews,
!    Special Functions of Mathematics for Engineers,
!    Second Edition, 
!    Oxford University Press, 1998.
!
!  Parameters:
!
!    Input, integer M, the number of evaluation points.
!
!    Input, integer N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = rk ) X(M), the evaluation points.
!
!    Output, real ( kind = rk ) P(M,0:N), the values of the first N+1 Hermite
!    polynomials at the point X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer j
  real ( kind = rk ) p(m,0:n)
  real ( kind = rk ) x(m)

  if ( n < 0 ) then
    return
  end if

  p(1:m,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  p(1:m,1) = 2.0D+00 * x(1:m)
 
  do j = 2, n
    p(1:m,j) = 2.0D+00 * x(1:m) * p(1:m,j-1) &
             - 2.0D+00 * real ( j - 1, kind = rk ) * p(1:m,j-2)
  end do
 
  return
end
subroutine h_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! H_POLYNOMIAL_VALUES: tabulated values of H(i,x).
!
!  Discussion:
!
!    H(i,x) is the physicist's Hermite polynomial of degree I.
!
!    In Mathematica, the function can be evaluated by:
!
!      HermiteH[n,x]
!
!  Differential equation:
!
!    Y'' - 2 X Y' + 2 N Y = 0
!
!  First terms:
!
!      1
!      2 X
!      4 X^2     -  2
!      8 X^3     - 12 X
!     16 X^4     - 48 X^2     + 12
!     32 X^5    - 160 X^3    + 120 X
!     64 X^6    - 480 X^4    + 720 X^2    - 120
!    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
!    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
!    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
!   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
!
!  Recursion:
!
!    H(0,X) = 1,
!    H(1,X) = 2*X,
!    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
!
!  Norm:
!
!    Integral ( -oo < X < +oo ) exp ( - X^2 ) * H(N,X)^2 dX
!    = sqrt ( PI ) * 2^N * N!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer N, the order of the polynomial.
!
!    Output, real ( kind = rk ) X, the point where the polynomial is evaluated.
!
!    Output, real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 18

  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
      0.1000000000000000D+01, &
      0.1000000000000000D+02, &
      0.9800000000000000D+02, &
      0.9400000000000000D+03, &
      0.8812000000000000D+04, &
      0.8060000000000000D+05, &
      0.7178800000000000D+06, &
      0.6211600000000000D+07, &
      0.5206568000000000D+08, &
      0.4212712000000000D+09, &
      0.3275529760000000D+10, &
      0.2432987360000000D+11, &
      0.1712370812800000D+12, &
      0.0000000000000000D+00, &
      0.4100000000000000D+02, &
     -0.8000000000000000D+01, &
      0.3816000000000000D+04, &
      0.3041200000000000D+07 /)
  integer n
  integer n_data
  integer, save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12,  5,  5, &
     5,  5, 5 /)
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    0.0D+00, &
    0.5D+00, &
    1.0D+00, &
    3.0D+00, &
    1.0D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine h_polynomial_zeros ( nt, z )

!*****************************************************************************80
!
!! H_POLYNOMIAL_ZEROS: zeros of H(i,x).
!
!  Discussion:
!
!    H(i,x) is the physicist's Hermite polynomial of degree I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NT, the degree of the polynomial.
!
!    Output, real ( kind = rk ) Z(NT), the zeros of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nt

  real ( kind = rk ) bj(nt)
  integer i
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) wts(nt)
  real ( kind = rk ) z(nt)

  z(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = sqrt ( real ( i, kind = rk ) / 2.0D+00 )
  end do

  wts(1:nt) = 0.0D+00
  wts(1) = sqrt ( sqrt ( r8_pi ) )

  call imtqlx ( nt, z, bj, wts )

  return
end
subroutine h_quadrature_rule ( nt, t, wts )

!*****************************************************************************80
!
!! H_QUADRATURE_RULE: quadrature for H(i,x).
!
!  Discussion:
!
!    H(i,x) is the physicist's Hermite polynomial of degree I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NT, the order of the rule.
!
!    Output, real ( kind = rk ) T(NT), WTS(NT), the points and weights 
!    of the rule.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nt

  real ( kind = rk ) bj(nt)
  integer i
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) t(nt)
  real ( kind = rk ) wts(nt)

  t(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = sqrt ( real ( i, kind = rk ) / 2.0D+00 )
  end do

  wts(1:nt) = 0.0D+00
  wts(1) = sqrt ( sqrt ( r8_pi ) )

  call imtqlx ( nt, t, bj, wts )

  wts(1:nt) = wts(1:nt) ** 2

  return
end subroutine h_quadrature_rule

subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to 
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine. 
!
!    It has been modified to produce the product Q' * Z, where Z is an input 
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
!    The changes consist (essentially) of applying the orthogonal 
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real ( kind = rk ) D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, real ( kind = rk ) E(N), the subdiagonal entries of the 
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, real ( kind = rk ) Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d(n)
  real ( kind = rk ) e(n)
  real ( kind = rk ) f
  real ( kind = rk ) g
  integer i
  integer ii
  integer, parameter :: itn = 30
  integer j
  integer k
  integer l
  integer m
  integer mml
  real ( kind = rk ) p
  real ( kind = rk ) prec
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) z(n)

  prec = epsilon ( prec )

  if ( n == 1 ) then
    return
  end if

  e(n) = 0.0D+00

  do l = 1, n

    j = 0

    do

      do m = l, n

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        exit
      end if

      if ( itn <= j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IMTQLX - Fatal error!'
        write ( *, '(a)' ) '  Iteration limit exceeded.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  L = ', l
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  N = ', n
        stop
      end if

      j = j + 1
      g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
      r =  sqrt ( g * g + 1.0D+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0D+00
      c = 1.0D+00
      p = 0.0D+00
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r =  sqrt ( c * c + 1.0D+00 )
          e(i+1) = f * r
          s = 1.0D+00 / r
          c = c * s
        else
          s = f / g
          r =  sqrt ( s * s + 1.0D+00 )
          e(i+1) = g * r
          c = 1.0D+00 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0D+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0D+00

    end do

  end do
!
!  Sorting.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end if

  end do

  return
end subroutine imtqlx


end module hermite

