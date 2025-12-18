
   using LinearAlgebra
   using Printf
   using Random
   rng = MersenneTwister()

#  the goal is to show that implicit single shift work
#  we are not claiming that this is a good shift!

   n = 8 
#  Λ = diagm(0 => Float32[2.0^(-i) for i=0:n-1])
   Λ = diagm(0 => Float32[i+1 for i=1:n])
#  Λ[n-1:n,n-1:n] = 0.25*[1. 2;-2 1]
   X = rand(rng, n, n)
   A = X * Λ / X
   A0 = copy(A)
   A0, tau = LAPACK.gehrd!(A0)
   A0 = triu(A0,-1)

#  display(A0)

## step of Schur iteration step (no implicit Q)
   B = A0 + 2.0 * Matrix{Float64}(I,n,n)
   Q,R = qr(B)
#  A1 = Q' * A0 * Q
   A1 = triu(Q' * A0 * Q,-1)
   display(A1)

## step of QR iteration without implicit Q 
   A2 = copy(A0)
   Q,R = qr(A2 + 2.0 * Matrix{Float64}(I,n,n))
   A2 = R*Q - 2.0 * Matrix{Float64}(I,n,n)
   display(A2)

## step of QR iteration with implicit Q 
   A3 = copy(A0)
   B = A3 + 2.0 * Matrix{Float64}(I,n,n)
   b = B[1:n,1]
   F = qr(b)
   Q = collect(F.Q)
   A3 = Q*A3*Q' 
   A3, tau = LAPACK.gehrd!(A3)
   A3 = triu(A3,-1)
   display(A3)

;



