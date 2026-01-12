#

  using LinearAlgebra
  using Printf
  using Random
  rng = MersenneTwister()

  n = 100
  A = randn(rng,Float64,n,n)
  A[3,1:3]   = 1e-6 * randn(rng,Float64,3) + A[1,1:3] + A[2,1:3]
# A[6,1:6]   = 1e-4 * randn(rng,Float64,6) + A[4,1:6] + A[5,1:6]
# A[9,1:9]   = 1e-4 * randn(rng,Float64,9) + A[7,1:9] + A[8,1:9]
# A[12,1:12] = 1e-4 * randn(rng,Float64,12) + A[10,1:12] + A[11,1:12]
# A[15,1:15] = 1e-4 * randn(rng,Float64,15) + A[13,1:15] + A[14,1:15]
# A[18,1:18] = 1e-4 * randn(rng,Float64,18) + A[16,1:18] + A[17,1:18]
# A[21,1:21] = 1e-4 * randn(rng,Float64,21) + A[16,1:21] + A[17,1:21]
# A[24,1:24] = 1e-4 * randn(rng,Float64,24) + A[22,1:24] + A[23,1:24]
  A0 = copy(A)
  x̂ = rand(rng,Float64,n,1)
  b = A * x̂
  x = copy(b)

  @printf("n                     = %8d\n",n)

# 1: LU with no pivoting
  A = copy(A0)
  x = copy(b)
  lu!(A, NoPivot())
  BLAS.trsm!( 'L', 'L', 'N', 'U', 1.0, A, x )
  BLAS.trsm!( 'L', 'U', 'N', 'N', 1.0, A, x )
  fwd1 = norm(x̂ - x,1)/(norm(x̂,1))
  bwd1 = norm(b - A0*x,1)/(opnorm(A0,1)*norm(x,1))
  bwdlu1 = opnorm(UnitLowerTriangular(A)*UpperTriangular(A) - A0 ,1)/opnorm(A0,1)
  bwdlubound1 = opnorm(UnitLowerTriangular(A),1)*opnorm(UpperTriangular(A),1)/opnorm(A0,1)*eps(Float64)

# 2: LU with partial pivoting
  A = copy(A0)
  x = copy(b)
  F = lu(A)
  x = x[F.p]
  x = F.L \ x
  x = F.U \ x
  fwd2 = norm(x̂ - x,1)/(norm(x̂,1))
  bwd2 = norm(b - A0*x,1)/(opnorm(A0,1)*norm(x,1))
  bwdlu2 = opnorm(F.L*F.U - A0[F.p,:],1)/opnorm(A0,1)
  bwdlubound2 = opnorm(F.L,1)*opnorm(F.U,1)/opnorm(A0,1)*eps(Float64)

#
  @printf("\n")
  @printf("condition number of A                         = %6.2e\n",cond(A))
  @printf("\n")
  @printf("[ no pivoting       ] forward error           = %6.2e\n", fwd1)
  @printf("[ no pivoting       ] backward error          = %6.2e\n", bwd1)
  @printf("[ no pivoting       ] backward error LU       = %6.2e\n", bwdlu1)
  @printf("[ no pivoting       ] backward error LU bound = %6.2e\n", bwdlubound1)
  @printf("\n")
  @printf("[ partial pivoting  ] forward error           = %6.2e\n", fwd2)
  @printf("[ partial pivoting  ] backward error          = %6.2e\n", bwd2)
  @printf("[ no pivoting       ] backward error LU       = %6.2e\n", bwdlu2)
  @printf("[ no pivoting       ] backward error LU bound = %6.2e\n", bwdlubound2)


