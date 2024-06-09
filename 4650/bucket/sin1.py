# creating our sine function: sin1

# input: x: (x can be an array)
# output: y: y is an array of the same length of x such that, for each entry i of the vector x,
#      y[i] = sin( x[i] )
#
# the function sin1 is accurate about (absolute error) 2.5e-3
# ( so, for all x, | sin1(x) - sin(x) | < 2.5e-3 )

# method: we approximate sin(x) with the interpolating polynomial of degree 3 that 
# passes by (x,f(x)) for x = [0., pi/6., 2.*pi/6., pi/2.]
# The interpolation uses Newton's nested form to store the interpolating polynomial.

# the coefficients `c` are pre-computed and hard coded. We only store
# 5 digits after the point (as opposed to 16 for example) since ``5 digits``
# seems to be enough to give an error below the interpolation error (about 2.5e-3)

# there is an in-place and an out-of-place version

# for out-of-place
#      y = sin1( x )

# for in-place
#      sin1_inplace( x )

# the bulk of the work is done in the function sin1_inplace

# for an example of use, see: 
# https://colab.research.google.com/drive/1re2cmh8443yRIuqkU_uUkWyzBWt1innZ

import numpy as np
from math import pi
import copy

def sin1( x ):

    y = copy.deepcopy( x )
    sin1_inplace( y )
    return y

def sin1_inplace( y ):
    b = np.array([0., pi/6., 2.*pi/6., pi/2.])
    c = np.array([ 0.00000,
                   0.95493,
                  -0.24434,
                  -0.11387 ] )
    for i in range(0,len(y)):
        y[i] = np.mod( y[i], 2*pi )
        s = 1
        if y[i] > pi:
            y[i] = 2 * pi - y[i]
            s = -1
        if y[i] > pi / 2.:
            y[i] = pi - y[i]
        y[i] = s * polyval_nested_w_base_points(c,b,y[i])

def polyval_nested_w_base_points( c, b, x ):
  d = np.size(c)
  px = c[ d-1 ]
  for i in range( d-2, -1, -1 ):
    px = px * ( x - b[i] ) + c[i]
  return px
