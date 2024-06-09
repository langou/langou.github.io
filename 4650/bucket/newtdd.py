# Newton Divided Difference algorithm

# input:  x and y: the points to be interpolated ( x_i, y_i) by the interpolating polynomial
# output: c:       the coefficients of the interpolating polynomial in nested form (with base nodes x_i)
 
# dimensions:
#    x is of size n
#    y must be of at least size n and the entries 0 to n-1 of y are considered
#    c is of size n-1
#    the polynomial p has degree n-1

# there are an in-place and an out-of-place versions

# for out-of-place
#      c = newtdd( x, y )

# for in-place
#      newtdd_inplace( x, y )
# then `y` in input are the y-values to be interpolated, 
# and `y` in output is `c`

# the bulk of the work is done in the function newtdd_inplace

# once `c` is constructed, to evaluate the polynomial p, at the points xx, 
# to get yy such that yy = p(xx), one can use the accompanying polyval_nested.py
# function, in the following manner
#     yy = polyval_nested( c, x, xx )
# with the `x`, the input x of newtdd, and the `y`, the output c of newtdd.

# a check that things might work correctly is to check that, for some x and y,
#    np.isclose( polyval_nested( *newtdd(x,y), x, x ), y ) 
# is true.

# example of use at: 
# https://colab.research.google.com/drive/1qi5fDnFfQDfR_Iy4NSfsUHHpcUnRPVD8

import copy
import numpy as np

def newtdd_inplace( x, y ):
    n = len( x )
    for i in range(1,n):
        for j in range(n-i-1,-1,-1):
            y[j+i] = ( y[j+i] - y[j+i-1] ) / ( x[j+i] - x[j] )
    return y

def newtdd( x, y ):
    c = copy.deepcopy( y )
    newtdd_inplace( x, c )
    return c
