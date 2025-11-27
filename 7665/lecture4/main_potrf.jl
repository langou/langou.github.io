
  using LinearAlgebra
  using Printf
  using Random
  rng = MersenneTwister()

# BLAS.set_num_threads(8)

  include("../src/potrf.jl")

  n = 4000

  A = randn(rng,Float64,n,n)
  for i = 1:n, A[i,i]= A[i,i] + n; end
  A0 = copy(A)
  x̂ = rand(rng,Float64,n,1)

  b = Symmetric(A0,:L) * x̂;
# b = Array{Float64, 2}(undef, n, 1)
# BLAS.symm!( 'L', 'L', 1.0, A, x̂, 0.0, b)

  x = Array{Float64, 2}(undef, n, 1)

  x .= b

  t =- time()

# A = cholesky(Symmetric(A,:L))
# x .= A.L \ x
# x .= A.L' \ x

# LAPACK.potrf!( 'L', A )
# BLAS.trsm!( 'L', 'L', 'N', 'N', 1.0, A, x )
# BLAS.trsm!( 'L', 'L', 'T', 'N', 1.0, A, x )

# cholesky_lower_recursive!(A)
# cholesky_lower_leftlooking_level0!(A)
# cholesky_lower_leftlooking_level2!(A)
# cholesky_lower_leftlooking_level3!(A)
  cholesky_lower_rightlooking_level3!(A)
# cholesky_lower_bordered_level3!(A)
  x .= LowerTriangular(A) \ x
  x .= Transpose(LowerTriangular(A)) \ x

  t += time()
  perf = 1.0/3.0*n*n*n/t*1e-9
  @printf("n = %5d, time = %8.4f (sec), perf = %6.2f (GFlops/sec)\n",n,t,perf)

  r = copy(b)
  BLAS.symm!( 'L', 'L', -1.0, A0, x, 1.0, r)

  @printf("forward error  = %6.2e\n",norm(x̂ - x,1)/(norm(x̂,1)))
  @printf("backward error = %6.2e\n",norm(r,1)/(opnorm(Symmetric(A0,:L),1)*norm(x,1)))

# we could also check that the upper part of A is the same as the upper part of A0

