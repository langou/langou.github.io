#
# useage:   
# => for the ``standard interface``, do
#           L, U = lu_no_pivoting___x( A )
#
# => for the ``in-place`` interface, do
#           lu_no_pivoting___x( A, inplace = True )
#
# performs LU_no_pivoting on `A`
#
# input:  A:       n-by-n matrix to be factored
# output: L an U:  the L and U factors defined by:
#                  (1) L is n-by-n lower unit triangular
#                  (2) U is n-by-n upper triangular
#                  (3) A = L * U
#
# in the ``in-place`` interface L and U are stored one on top of the other, 
# in the array A. So, in input, you give A; in ouput, you get L and U stacked 
# one on top of the other.
#
# warning: no pivoting is used, this algorithm will exit with error message if 
# an exact-zero pivot is encountered, this code might be quite unstable when 
# small pivots are encountered during the factorization.
#
def lu_no_pivoting___x( A_, inplace = False ):
#    
  if ( inplace ):
    A = A_
  else: 
    A = copy.deepcopy(A_)
# begin main code
  n = A.shape[0]
  for k in range(0,n-1):
    if ( A[k,k] == 0. ): print("oops, zero pivot encountered\n"); break
    for i in range(k+1,n):
      A[i,k] = A[i,k] / A[k,k]
      for j in range(k+1,n):
        A[i,j] = A[i,j] - A[k,j] * A[i,k]
# end main code   
  if not inplace:
    L = np.tril(A,-1)+np.eye(n)
    U = np.triu(A)
    return L, U
