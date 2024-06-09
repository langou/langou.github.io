def newton_method( f, dfdx, x, numiter, r ):
  for i in range(0,numiter):
    x = x - f(x) / dfdx(x)
    true_fwd_rel_error = abs( r - x ) / abs( r )
    print( f"{i+1:2d}", f"{x:.16f}", f"{true_fwd_rel_error:.2e}" )
  return x
