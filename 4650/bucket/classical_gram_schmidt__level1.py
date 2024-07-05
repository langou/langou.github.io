# perform a QR factorization of `A` using the Classical Gram-Schmidt algorithm
#
# input:  Q and R, (R is "memory space" where you want R, only upper part is referenced)
# R is in/out. Q is input only
# output: Q and R
def classical_gram_schmidt__level1( Q, R ):
  n = Q.shape[1]
  for j in range(0,n):
    for i in range(0,j):
      R[i,j] = Q[:,i].T @ Q[:,j]
    for i in range(0,j):
      Q[:,j] = Q[:,j] - Q[:,i] * R[i,j]
    R[j,j] = np.linalg.norm(Q[:,j])
    Q[:,j] = Q[:,j]/R[j,j]
  return
