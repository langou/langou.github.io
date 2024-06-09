def method_of_false_position( f, a, b, numiter, r ):

  fa = f(a)
  fb = f(b)

  for i in range(0,numiter):
    c = ( b*f(a) - a*fb ) / ( fa - fb )
    true_fwd_rel_error = abs( r - c ) / abs( r )
    print( f"{i+1:2d}", f"{c:.16f}", f"{true_fwd_rel_error:.2e}" )
    fc = f(c)
    if ( ( fb > 0 ) and ( fc > 0 ) ) or ( ( fb < 0 ) and ( fc < 0 ) ):
      b = c
      fb = fc
    else:
      a = c
      fa = fc

  return c
