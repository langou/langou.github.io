# Bisection Method
# Computes an approximate solution of f(x)=0 using the bisection method
#
# Input:  f:        function: ( f: R -> R )
#         a and b:     float: numbers a and b such that f(a) * f(b) < 0
#         tol:         float: tolerance
#         verbose:   boolean: True=print, False = no print
#         
# Output: c:           float: an approximate solution c of f(x)=0 
#         err_bound:   float: bound on the error
#         numeval:       int: number of function evaluations
#
# Assumption: f is continuous on [a,b]

def bisect( f, a, b, tol, verbose ):

  numeval = 0

  err_bound = []

  k = 0

  fa = f(a)
  numeval += 1

  if ( fa == 0 ): # `a` is a solution, done
    if ( verbose ): print( 'f(a) = 0, returning `a` as a solution ' )
    return a, 0., numeval

  fb = f(b)
  numeval += 1

  if ( fb == 0 ): # `b` is a solution, done
    if ( verbose ): print( 'f(b) = 0, returning `b` as a solution ' )
    return b, 0., numeval

  if ( ( fa > 0 ) and ( fb > 0 ) ) or ( ( fa < 0 ) and ( fb < 0 ) ):
    print( '[ error ] f(a) and f(b) have the same sign' )
    return float("NaN"), numeval, float("NaN")

  if ( ( b - a ) / 2 <=  tol ):
    c = ( a + b ) / 2
    err_bound.append( ( b - a ) / 2 )
    if ( verbose ): print( '( b - a ) / 2 <=  tol, f(a) and f(b) opposite signs' )
    if ( verbose ): print( 'returning c = ( b + a ) / 2 as an approximate solution' )
    return c, numeval, err_bound
  
  while ( ( b - a ) / 2 > tol ):

    c = ( a + b ) / 2

    err_bound.append( ( b - a ) / 2 )

    if ( verbose ): print( f"{k:3d}", f"{c:20.16f}", '  ', f"{err_bound[k]:3.1e}" )
    
    k = k + 1

    fc = f(c)
    numeval += 1

    if ( fc == 0 ): # `c` is a solution, done
      return c, 0., numeval

    if ( ( fb > 0 ) and ( fc > 0 ) ) or ( ( fb < 0 ) and ( fc < 0 ) ):
#     `a` and `c` make the new interval
#     so replace `b` with `c`  
      b = c
      fb = fc
      
    else:
#     `c` and `b` make the new interval
#     so replace `a` with `c` 
      a = c
      fa = fc

    c = ( a + b ) / 2

  err_bound.append( (b - a ) / 2 )

  if ( verbose ): print( f"{k:3d}", f"{c:20.16f}", '  ', f"{err_bound[k]:3.1e}" )

  return c, err_bound, numeval
