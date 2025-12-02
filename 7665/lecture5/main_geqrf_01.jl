   using LinearAlgebra
   using Printf
   using Random
   rng = MersenneTwister()

   m = 10; n = 7; A = randn(rng,Float64,m,n)


# “reduced” QR factorization

   Q, R = qr(A)
   Q = Matrix(Q)

   @printf("|| I - Qᴴ Q ||₁             = %6.2e\n",
      opnorm( I(n) - Q' * Q, 1 ) )
   @printf("|| A - Q * R ||₁ / || A ||₁ = %6.2e\n",
      opnorm( A - Q * UpperTriangular(R), 1) / opnorm(A,1))


# To retrieve the “full” Q factor, an m×m orthogonal matrix, use F.Q*I or collect(F.Q)

   Q, R = qr(A)
   Q = collect(Q)
   R = [R;zeros(m-n,n)]

   @printf("|| I - Qᴴ Q ||₁             = %6.2e\n",
      opnorm( I(m) - Q' * Q, 1 ) )
   @printf("|| A - Q * R ||₁ / || A ||₁ = %6.2e\n",
      opnorm( A - Q * R, 1) / opnorm(A,1))

