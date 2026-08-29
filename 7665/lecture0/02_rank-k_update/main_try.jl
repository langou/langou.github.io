
  using LinearAlgebra
  using Random
  using Printf
  using GenericLinearAlgebra


# BLAS.set_num_threads(8)

  m = 9000
  n = 500

  A0 = randn(Float64,m,m)
  B = randn(Float64,m,n)
  A1 = randn(Float64,m,m)
  A2 = randn(Float64,m,m)
  x = randn(Float64,m,1)
  y0 = randn(Float64,m,1)
  y1 = randn(Float64,m,1)

  A2 .= A0; for j=1:m, i=1:j-1, A2[j,i]=A2[i,j]; end
  t =- time(); A2 .= A2 - B * Adjoint(B); t += time();
  y1 .= A2 * x;
  @printf("m = %5d, n = %5d, time = %8.4f (sec)\n",m,n,t)
  
# A1 .= A0
# t =- time(); A1 .= Symmetric(A1) - B * Adjoint(B); t += time();
# y2 = Symmetric(A1) * x; err = norm( y2 - y1 ) / norm( y1 );
# lowerpart = opnorm(tril(A1-A0,-1),1);
# @printf("m = %5d, n = %5d, time = %8.4f (sec), err = %6.2e, lower = %6.2e\n",m,n,t,err,lowerpart)

  A1 .= A0
  t =- time(); BLAS.syrk!( 'U', 'N', -1.0e+00, B, +1.e+00, A1); t += time();
  y2 = Symmetric(A1) * x; err = norm( y2 - y1 ) / norm( y1 );
  lowerpart = opnorm(tril(A1-A0,-1),1);
  @printf("m = %5d, n = %5d, time = %8.4f (sec), err = %6.2e, lower = %6.2e\n",m,n,t,err,lowerpart)

  A1 .= A0
  t =- time(); rankUpdate!(Symmetric(A1), B, -1.0e+00, 1.0e+00); t += time();
  y2 = Symmetric(A1) * x; err = norm( y2 - y1 ) / norm( y1 );
  lowerpart = opnorm(tril(A1-A0,-1),1);
  @printf("m = %5d, n = %5d, time = %8.4f (sec), err = %6.2e, lower = %6.2e\n",m,n,t,err,lowerpart)

  A1 .= A0
  t =- time(); mul!(A1, B, Adjoint(B), -1.0e+00, 1.0e+00); t += time();
  y2 = Symmetric(A1) * x; err = norm( y2 - y1 ) / norm( y1 );
  lowerpart = opnorm(tril(A1-A0,-1),1);
  @printf("m = %5d, n = %5d, time = %8.4f (sec), err = %6.2e, lower = %6.2e\n",m,n,t,err,lowerpart)

