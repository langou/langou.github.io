using Printf
using Random
using LinearAlgebra
rng = MersenneTwister(18);

include("/Users/langou/Documents/repositories/darve__numerical_linear_algebra_with_julia.git/src/gees.jl")

n = 4

Λ = diagm(0 => Float32[2.0^(-i) for i=0:n-1]);
#Λ = zeros(n,n); Λ[1:n-2,1:n-2] = diagm(0 => Float32[2.0^(-i) for i=0:n-3]); Λ[n-1:n,n-1:n] = [0 1; -1 0]

X = rand(rng, n, n)

# Λ = diagm(Float32[2.0^(-i) for i=0:n-1])
# A = Λ
# X, = qr(rand(rng, n, n))

A = X * Λ / X

println("Matrix A")
pretty_print(A)
Λ = eigenvalue_sorted(A)

# Applying Householder reflections to make matrix A upper Hessenberg
gehrd!(A)
A0 = copy(A)
#display( triu(A0,-1) )

# Checking that A is upper Hessenberg now
println("Norm for lower part of matrix = ",norm(tril(A,-2)))
println("Matrix A")
display(A)

# Checking that the eigenvalues are the same
D = eigenvalue_sorted(A)
println(norm(D-Λ))

# Sequential version of QR iteration using Givens rotations
for i=1:10
givens_QR_iteration_s!(A)
display(A)
end

# QR iteration with bulge chasing
#A = copy(A0)
#for i=1:100
#givens_QR_iteration!(A)
#display(A)
#end
#;
