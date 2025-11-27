   using LinearAlgebra
   using Printf
   using Random
   rng = MersenneTwister()

   function house(x)
      """Computes the Householder transformation for input vector x"""
      sigma = dot(x[2:end],x[2:end])
      v = copy(x)

      if sigma == 0
          beta = 0
          return beta, v
      end

      sq = sqrt(x[1]^2 + sigma)
      if x[1] > 0
          v[1] += sq
      else
          v[1] -= sq
      end

      beta = 2.0 / (v[1]^2 + sigma)

      return beta, v
   end
   
   function geqrf!(A)
      m = size(A,1)
      n = size(A,2)
      vA = zeros(n)
      kend = (m > n ? n : m-1)
      for k=1:kend
          beta, v = house(A[k:end,k])
          for j=k:n
              # vA = beta * v^T * A
              vA[j] = 0.0
              for i=k:m
                  vA[j] += v[i-k+1] * A[i,j]
              end
              vA[j] *= beta
          end
          # A - beta v (v^T A)
          for j=k:n, i=k:m
              A[i,j] -= v[i-k+1] * vA[j]
          end
          A[k+1:end,k] = v[2:end] 
          # Saving v in the lower triangular part of A.
          # This was not done here but one can always
          # divide v by v[1] such that v[1] = 1 is always true.
          # In that case, v[1] does not need to be stored.
      end
      # Lower triangular part of A: sequence of v vectors
      # Upper triangular part: factor R
   end

   # Testing QR factorization using Householder transformations
   n = 4

   # Building an orthogonal matrix Q
   Q = zeros(n,n)
   for j=0:n-1, i=0:n-1
      Q[i+1,j+1] = cos(π*(2i+1)*j/2n)
   end
   for j=1:n
      Q[:,j] /= norm(Q[:,j])
   end

   # Initializing an upper triangular matrix R
   R = triu(Float64[ i/j for i=1:n, j=1:n ])

   # Matrix A
   A = Q*R

   geqrf!(A)

   @show R,triu(A)

#  @printf("|| I - Qᴴ Q ||₁             = %6.2e\n",
#     opnorm( I(n) - Q' * Q, 1 ) )
#  @printf("|| A - Q * R ||₁ / || A ||₁ = %6.2e\n",
#     opnorm( A - Q * UpperTriangular(R), 1) / opnorm(A,1))


