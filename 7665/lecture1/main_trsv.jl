using Random
using LinearAlgebra
using Printf

rng = MersenneTwister(2025);

# Size of the matrix
n = 9;

# Lower triangular matrix
L = zeros(Float64,n,n)
# Filling the matrix with random integer entries
for j=1:n # Column j
    L[j,j] = 1 # Should be non-zero
    L[j+1:n,j] = rand(rng, -2:2, n-j) 
    #L[j+1:n,j] = rand(n-j) 
    #L[j+1:n,j] = L[j+1:n,j] + randn(n-j)*1e-15; 
end
display(L[1:6,1:6])

# Initializing the right-hand side
xe = rand(rng, 0:9, n) # This will be our solution
#xe = xe + 1e-15 * randn(n) # This will be our solution
b = L * xe
b'

x = Vector{Float64}(undef,n)

# Load our triangular solvers
include("../src/trsv.jl")

# x = deepcopy(b)
x .= b
@time trsv_lnn_level0_ijvariant!(L, x)
@assert x == xe

#x = deepcopy(b)
x .= b
@time trsv_lnn_level0_jivariant!(L, x)
@assert x == xe

#x = deepcopy(b)
x .= b
@time trsv_lnn_level0_ijvariant!(L, x)
@assert x == xe

#x = deepcopy(b)
x .= b
@time trsv_lnn_level0_jivariant!(L, x)
@assert x == xe

@time x = LowerTriangular(L)\b
@time x = LowerTriangular(L)\b
@time x = L\b
@assert x == xe

#show opnorm(L,Inf)
@show norm(L*x-b,Inf)/(norm(b,Inf))
@show norm(L*xe-b,Inf)/(norm(b,Inf))
#@show norm(b,Inf)
#@show L
#display(L)
#display(LowerTriangular(L))


#@show x
#x = deepcopy(b)
x .= b
#@show x
#@show xe
@time BLAS.trsv!('L','N','U',L,x)
x = deepcopy(b)
@time BLAS.trsv!('L','N','U',L,x)
#@assert x == xe
@show norm(L*x-b,Inf)
#@show x

;


