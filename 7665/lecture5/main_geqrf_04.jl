   using LinearAlgebra
   using Printf
   using Random
   rng = MersenneTwister()

   m = 10; n = 7; A = randn(rng,Float64,m,n)

# “reduced” QR factorization 

   Q = Matrix{Float64}(undef,(m,n))
   R = Matrix{Float64}(undef,(n,n))
   tau = Vector{Float64}(undef,n)
 
   Q .= A                                                      # do not forget this dot in the .=
   LAPACK.geqrf!( Q, tau )
   R .= triu(Q[1:n,1:n])
   LAPACK.orgqr!( Q, tau, n )
 
   @printf("|| I - Qᴴ Q ||₁             = %6.2e\n",
      opnorm( I(n) - Q' * Q, 1 ) )                             # you would hope that SYRK is being called
   @printf("|| A - Q * R ||₁ / || A ||₁ = %6.2e\n",
      opnorm( A - Q * UpperTriangular(R), 1) / opnorm(A,1))    # you would hope that TRSM is being called

# “full” QR factorization 

#  Q = Matrix{Float64}(undef,(m,m))
#  R = Matrix{Float64}(undef,(m,n))
#  tau = Vector{Float64}(undef,n)
#
#  Q[1:m,1:n] .= A
#  LAPACK.geqrf!( Q[1:m,1:n], tau )
#  LAPACK.orgqr!( Q, tau, m )
#  R[1:n,1:n] .= triu(Q[1:n,1:n])
#  R = [R;zeros(m-n,n)]
#
#  @printf("|| I - Qᴴ Q ||₁             = %6.2e\n",
#     opnorm( I(m) - Q' * Q, 1 ) )
#  @printf("|| A - Q * R ||₁ / || A ||₁ = %6.2e\n",
#     opnorm( A - Q * UpperTriangular(R), 1) / opnorm(A,1))

