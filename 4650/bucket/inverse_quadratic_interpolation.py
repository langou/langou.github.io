def inverse_quadratic_interpolation( f, a, b, c, numiter, rr ):

  fa = f(a)
  fb = f(b)

  for i in range(0,numiter):

    fc = f(c)

    q = fa / fb
    r = fc / fb
    s = fc / fa

    x = c - ( r * ( r - q ) * ( c - b ) + ( 1 - r ) * s * ( c - a ) ) \
          / ( ( q - 1. ) * ( r - 1. ) * ( s - 1. ) )

    true_fwd_rel_error = abs( rr - x ) / abs( rr )

    print( f"{i+1:2d}", f"{x:20.16f}", f"{true_fwd_rel_error:.2e}" )

    a = b
    fa = fb
    b = c
    fb = fc
    c = x

  return c
