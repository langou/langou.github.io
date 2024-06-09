# perform row reduction on `A` using the Gaussian elimination
# apply same transformations to `b`
# algorithm as per Section 2.1
# no pivoting is used
#
# input:  A and b: represent linear system of equations Ax=b we want to solve
# output: A and b: equivalent linear system of equations Ax=b 
#                  such that A is upper triangular.
#
# there is a deep copy of A and b in input to mimic a matlab-like interface
# this might be more user-friendly for python beginers
def gaussian_elimination__section_2_1( A_, b_ ):
  A = copy.deepcopy(A_)
  b = copy.deepcopy(b_)
  n = A.shape[0]
  for k in range(0,n-1):
    if ( A[k,k] == 0. ): print("oops, zero pivot encountered\n"); break
    for i in range(k+1,n):
      mult = A[i,k] / A[k,k]
      for j in range(k+1,n):
        A[i,j] = A[i,j] - A[k,j] * mult
      b[i] = b[i] - b[k] * mult
      A[i,k] = 0.
  return A, b
