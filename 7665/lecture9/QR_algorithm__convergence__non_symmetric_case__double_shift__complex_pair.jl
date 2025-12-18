
   using LinearAlgebra
   using Printf
   using Random
   rng = MersenneTwister()

   n = 8 
#  Λ = diagm(0 => Float32[2.0^(-i) for i=0:n-1])
   Λ = diagm(0 => Float32[i for i=1:n])
   Λ[n-1:n,n-1:n] = 0.25*[1. 2;-2 1]
   X = rand(rng, n, n)
   A = X * Λ / X
   A0 = copy(A)
   A0, tau = LAPACK.gehrd!(A0)
   A0 = triu(A0,-1)

#  here we do 10 step of QR algorithm to "prep" the lower 2x2 block to be our complex conjugate pair
#  we could have done the experiment with a :
   A1 = copy(A0)
   for i=1:10
   global A1, Q, R
   Q,R = qr(A1)
   A1 = R*Q
   end

function double_shift_st(A)
    a = A[1,1]
    b = A[1,2]
    c = A[2,1]
    d = A[2,2]
    s = a+d       # sum
    t = a*d - b*c # product
    return s,t
end

  display(A1)
  display(eigvals(A1[n-1:n,n-1:n]))

   A2 = copy(A1)
   for i=1:3
      global A2, Q, R, B, s, t, b, F, tau
      s,t = double_shift_st(A2[n-1:n,n-1:n])

##    step of Schur iteration step (no implicit Q)
#     B = A2^2 - s * A2 + t * Matrix{Float64}(I,n,n)
#     Q,R = qr(B)
#     A2 = triu(Q' * A2 * Q,-1)

##    step of QR iteration with implicit Q 
      B = A2^2 - s * A2 + t * Matrix{Float64}(I,n,n)
      b = B[1:n,1]
      F = qr(b)
      Q = collect(F.Q)
      A2 = Q*A2*Q' 
      A2, tau = LAPACK.gehrd!(A2)
      A2 = triu(A2,-1)

      display(A2)
      display(eigvals(A2[n-1:n,n-1:n]))
   end

;



