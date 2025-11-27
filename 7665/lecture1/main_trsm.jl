  using Random
  using LinearAlgebra
  using Printf

  include("../src/trsm.jl")

  rng = MersenneTwister(2025);

  n = 4000;
  k = 200;

  A = zeros(Float64,n,n)
  for j=1:n
    A[j,j] = 1
    A[j+1:n,j] = rand(rng, -2:2, n-j) 
  end
# display(A[1:6,1:6])

  x̂ = rand(rng, 0:9, (n,k)) 
  b = A * x̂
  x = deepcopy(b)

# preload some functions before timings
 
  trsm_llnn!(A[1:2,1:2], x[1:2,1:2])
  trsm_llnn_recursive!(A[1:2,1:2], x[1:2,1:2])
  x = LowerTriangular(A[1:2,1:2])\x[1:2,1:2]
  x = A[1:2,1:2]\x[1:2,1:2]
  BLAS.trsm!('L','L','N','U',1.0e+00,A[1:2,1:2],x[1:2,1:2])

# start timings here  

  x = Matrix{Float64}(undef,(n,k))

  x .= b
  t =- time()
  trsm_llnn!(A, x)
  t += time()
  @assert x == x̂
  perf = n*n*k/t*1e-9
  @printf("n = %5d, k = %5d, time = %8.4f (sec), perf = %6.2f (GFlops/sec)\n",n,k,t,perf)

  t =- time()
  x = A \ b
  t += time()
  @assert x == x̂
  perf = n*n*k/t*1e-9
  @printf("n = %5d, k = %5d, time = %8.4f (sec), perf = %6.2f (GFlops/sec)\n",n,k,t,perf)

  t =- time()
  x = UnitLowerTriangular(A) \ b
  t += time()
  @assert x == x̂
  perf = n*n*k/t*1e-9
  @printf("n = %5d, k = %5d, time = %8.4f (sec), perf = %6.2f (GFlops/sec)\n",n,k,t,perf)

  x .= b
  t =- time()
  BLAS.trsm!('L','L','N','U',1.0e+00,A,x)
  t += time()
  @assert x == x̂
  perf = n*n*k/t*1e-9
  @printf("n = %5d, k = %5d, time = %8.4f (sec), perf = %6.2f (GFlops/sec)\n",n,k,t,perf)

  x .= b
  t =- time()
  trsm_llnn_recursive!(A, x)
  t += time()
  @assert x == x̂
  perf = n*n*k/t*1e-9
  @printf("n = %5d, k = %5d, time = %8.4f (sec), perf = %6.2f (GFlops/sec)\n",n,k,t,perf)

