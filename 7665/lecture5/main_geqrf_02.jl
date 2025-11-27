   using LinearAlgebra
   using Printf
   using Random
   rng = MersenneTwister()

   m = 10; n = 1; x = randn(rng,Float64,m,n)

   beta, v = house(x)

n = length(x)
y = zeros(n); y[1] = norm(x)
Px = x - beta * dot(v,x) * v
if Px[1] < 0
    Px = -Px
end
@show Px
@show y
@show norm(Px - y)
