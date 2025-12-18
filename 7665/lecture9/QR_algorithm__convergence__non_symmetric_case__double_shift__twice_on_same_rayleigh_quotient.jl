
   using LinearAlgebra
   using Printf
   using Random
   rng = MersenneTwister()

   n = 8 
   Λ = diagm(0 => Float32[i for i=1:n])
   X = rand(rng, n, n)
   A = X * Λ / X
   A0 = copy(A)
   A0, tau = LAPACK.gehrd!(A0)
   A0 = triu(A0,-1)

#  here we can do a few steps of QR algorithm to "prep" the Hessenberg to have "ordered" eigenvalues
   A1 = copy(A0)
#  for i=1:10
#  global A1, Q, R
#  Q,R = qr(A1)
#  A1 = R*Q
#  end

   display(A1)

   residual = [] 
   push!( residual, abs(A1[n,n-1]) )

   A2 = copy(A1)
   for i=1:6

      global A2, Q, R, B, s, t, b, F, tau

##    either case 1: step of Schur iteration step

##    either case 1-a: double shift on the real eigenvalue
      s = 2*A2[n,n]
      t = A2[n,n]^2
      B = A2^2 - s * A2 + t * Matrix{Float64}(I,n,n)

##    or     case 1-b: single shift on the real eigenvalue
#     μ = A2[n,n]
#     B = A2 - μ * Matrix{Float64}(I,n,n)

      Q,R = qr(B)
      A2 = triu(Q' * A2 * Q,-1)

##    step of QR iteration with implicit Q 
#     B = A2^2 - s * A2 + t * Matrix{Float64}(I,n,n)
#     b = B[1:n,1]
#     F = qr(b)
#     Q = collect(F.Q)
#     A2 = Q*A2*Q' 
#     A2, tau = LAPACK.gehrd!(A2)
#     A2 = triu(A2,-1)

      push!( residual, abs(A2[n,n-1]) )
      display(A2)

   end

   @printf("\n")
   for i in 1:length(residual)
      @printf("%6.2e ", residual[i] )
   end
   @printf("\n")
;



