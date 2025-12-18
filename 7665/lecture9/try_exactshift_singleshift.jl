
   using LinearAlgebra
   using Printf
   using Random
   rng = MersenneTwister()

   n = 10
   Λ = diagm(0 => Float32[2.0^(-i) for i=0:n-1])
#  Λ[1,1] = 10.
   X = rand(rng, n, n)
   A = X * Λ / X
#  A = rand(rng, n, n)

#####################################################
## this is an exact shift when we know one eigenvalue

   A0 = copy(A)
   A0, tau = LAPACK.gehrd!(A0)
   A0 = triu(A0,-1)
   H2 = copy(A0)
#  e = eigvals(A0)
#  mu = Real(e[n])
#  mu = Real(e[1])
#  mu = Λ[1,1]
   mu = Λ[n,n]
   #u = 
   Q, R = qr(H2 - mu*I)
   H2 = R * Q + mu * I
   display( H2 )
#  Q, R = qr(H2 - mu*I)
#  H2 = R * Q + mu * I
;

#####################################################
## this is an exact shift when we know one eigenvector

#  A0 = copy(A)
#  A0, tau = LAPACK.gehrd!(A0)
#  Q = tril(A0,-2) 
#  Q = Q[2:n,1:n-1]
#  LAPACK.orgqr!(Q,tau)
#  Q = [ 1 zeros(1,n-1); zeros(n-1,1) Q]
#  A0 = triu(A0,-1)
## check that
#  Q*A0*Q'-A
 
## we have
#  A0 = (Q'*X) * Λ / (Q*X) 
#  V = (Q'*X) 
## we have
#  A0 - V * Λ / V

#  v = Q'*X[:,1]
#  F = qr(v)
#  W = collect(F.Q)
#  W[:,[1 n]] = W[:,[n 1]]
#  display( W'*A0*W )

;
