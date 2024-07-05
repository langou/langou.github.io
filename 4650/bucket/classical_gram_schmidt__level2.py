# perform a QR factorization of `A` using the Classical Gram-Schmidt algorithm
#
# input:  Q and R, (R is "memory space" where you want R, only upper part is referenced)
# R is in/out. Q is input only
# output: Q and R
def classical_gram_schmidt__level2( Q, R ):
  n = Q.shape[1]
  for j in range(0,n):
    R[0:j,j] = Q[:,0:j].T @ Q[:,j]
    Q[:,j] = Q[:,j] - Q[:,0:j] @ R[0:j,j]
    R[j,j] = np.linalg.norm(Q[:,j])
    Q[:,j] = Q[:,j]/R[j,j]
  return

