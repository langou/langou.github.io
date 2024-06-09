#
##############################################################################################
#
# backward_substitution
#
#
# perform backward substitution Ax = b where A is upper triangular
#
# Note that we do not check whether A is upper triangular, it is assumed to be,
# Only the upper part of A is referenced. The zeros in the lower part are assumed.
# You can store whatever you want in the lower part of A.
#
# useage: 
# => for the ``standard interface``, do
#           x = backward_substitution( A, b )
#
# => for the ``in-place`` interface, do
#           backward_substitution( A, b, inplace = True )
#
# input:  A and b: represent the upper triangular linear system of equations 
#                  Ax=b we want to solve
# output:          solution of x of Ax=b, 
#    
def backward_substitution( A, b_, inplace = False ):
#
  import numpy as np
  import copy
#    
  if ( inplace ):
    b = b_
  else: 
    b = copy.deepcopy(b_)
#    
  n = A.shape[0]
  for i in range(n-1,-1,-1):
    for j in range(i+1,n):
      b[i] = b[i] - A[i,j] * b[j]
    b[i] = b[i] / A[i,i]
#    
  if not inplace:
    return b
#
