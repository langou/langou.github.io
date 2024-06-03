### This work is licensed under a Creative Commons Attribution 4.0 International License
### https://creativecommons.org/licenses/by/4.0/
### Copyright (c) 2021, Julien Langou. All rights reserved.

### def polyval_basic( c, x ):
###   return px

### purpose
# `polyval_basic` evaluates the polynomial `p` (as represented in `c`) at the
# data points `x`

### interface
# input:  `c`   a numpy array of size n+1 containing the coefficients 
#               of polynomial `p` ordered from degree 0 to degree n
# input:  `x`   a numpy array of m data points at which we want to evaluate `p`
# output: `px`  a numpy array of the m values of `p` at `x`

### useage
# to evaluate `p(x) = 3x^2 - 2x + 5` at x = -7 and x = 2, do
#     p = np.array( [ 5., -2., 3. ]);
#     x = np.array( [ -7., 2. ]);
#     px = polyval( p, x);
# you should get px = [ 166., 13. ]
# indeed p(-7) = 166 and p(2) = 13

### notes
# (1) "same" interface as numpy `np.polyval( np.flip(p), x )`
#     see: https://numpy.org/doc/stable/reference/generated/numpy.polyval.html
#     see: https://numpy.org/doc/stable/reference/generated/numpy.flip.html
# (2) using a Horner's method (nested multiplication)

def polyval_nested( c, x ):
  d = np.size(c)
  px = c[ d-1 ]
  for i in range( d-2, -1, -1 ):
    px = px * x + c[i]
  return px
