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

   m = 10; x = randn(rng,Float64,m,1)

   τ, v = house(x)

   n = length(x)
   y = zeros(n); y[1] = norm(x)
   Px = x - v * τ * v' * x
   if Px[1] < 0
      Px = -Px
   end

   @show Px
   @show y
   @show norm(Px - y)
