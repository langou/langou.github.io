def secant_method( f, x0, x1, numiter, r ):

  true_fwd_rel_error = abs( r - x0 ) / abs( r )
  print( f"{0:2d}", f"{x0:.16f}", f"{true_fwd_rel_error:.2e}" )

  true_fwd_rel_error = abs( r - x1 ) / abs( r )
  print( f"{1:2d}", f"{x1:.16f}", f"{true_fwd_rel_error:.2e}" )

  fx0 = f( x0 )
  for i in range(1,numiter):
    fx1 = f( x1 )
    x = x1 - fx1 * ( x1 - x0 ) / ( fx1 - fx0 )
    true_fwd_rel_error = abs( r - x ) / abs( r )
    print( f"{i+1:2d}", f"{x:.16f}", f"{true_fwd_rel_error:.2e}" )
    x0 = x1
    fx0 = fx1
    x1 = x
  return x
