#
##############################################################################################
#
# forward__substitution
#
# This work is licensed under a Creative Commons Attribution 4.0 International License
# https://creativecommons.org/licenses/by/4.0/
# Copyright (c) 2021, Julien Langou. All rights reserved.
#
# useage: 
# => for the ``standard interface``, do
#           x = forward__substitution___x( A, b )
#
# => for the ``in-place`` interface, do
#           forward__substitution___x( A, b, inplace = True )
#
# the flag ``unit`` is to say whether there are one on the diagonal of A.
# If unit is ``True`` then 1 are assumed on the diagonal of A. We do not reference the diagonal
# and you can store whatever in it.
#
# Note that we do not check whether A is lower triangular, it is assumed to be.
# Only the lower part of A is referenced. The zeros in the upper part are assumed.
# You can store whatever you want in the upper part of A.
#
def forward__substitution( A, b_, inplace = False, unit = False ):
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
  for i in range(0,n,):
    for j in range(0,i):
      b[i] = b[i] - A[i,j] * b[j]
    if not unit:
      b[i] = b[i] / A[i,i]
#    
  if not inplace:
    return b
