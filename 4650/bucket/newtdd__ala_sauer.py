# Newton Divided Difference algorithm, following Sauer's algorithm

# input:  x and y: the points to be interpolated ( x_i, y_i) by the interpolating polynomial
# output: c:       the coefficients of the interpolating polynomial in nested form (with base nodes x_i)
 
# dimensions:
#    x is of size n
#    y must be of at least size n and the entries 0 to n-1 of y are considered
#    c is of size n-1
#    the polynomial p has degree n-1

# there are an in-place and an out-of-place versions

# for out-of-place
#      c = newtdd__ala_sauer( x, y )

# once `c` is constructed, to evaluate the polynomial p, at the points xx, 
# to get yy such that yy = p(xx), one can use the accompanying 
# polyval_nested_w_base_points.py
# function, in the following manner
#     yy = polyval_nested_w_base_points( c, x, xx )
# with the `x`, the input x of newtdd__ala_sauer, and the `y`, 
# the output c of newtdd__ala_sauer.

# a check that things might work correctly is to check that, for some x and y,
#    np.isclose( polyval_nested_w_base_points( *newtdd__ala_sauer(x,y), x, x ), y ) 
# is true.

# example of use at: 
# https://colab.research.google.com/drive/174ssa5vC6scoSSb1Dni2Z6lsDfc6MNuq

import copy
import numpy as np

def newtdd__ala_sauer( x, y ):

    n = len( x )
    v = np.zeros([n,n])
    c = np.zeros([n])

    for j in range(0,n):
        v[j,0] = y[j]

    for i in range(1,n):
        for j in range(0,n-i):
            v[j,i] = ( v[j+1,i-1] - v[j,i-1] ) / ( x[j+i] - x[j] )

    for i in range(0,n):
        c[i] = v[0,i]

    return c
