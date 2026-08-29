   using LinearAlgebra
   using Printf
   using Random

   n = 10;

   A = randn(rng,Float64,n,n);
   B = randn(rng,Float64,n,n);
   C = randn(rng,Float64,n,n);

   @show( opnorm( A * ( B * C ) - A * B * C ) )
   @show( opnorm( ( A * B ) * C - A * B * C ) )

   A = randn(rng,Float64,n,n);
   B = randn(rng,Float64,n,n);
   C = randn(rng,Float64,n,1);

   @show( opnorm( A * ( B * C ) - A * B * C ) )
   @show( opnorm( ( A * B ) * C - A * B * C ) )

   A = randn(rng,Float64,n,n);
   B = randn(rng,Float64,n,n);
   C = randn(rng,Float64,n,n);

#  crazy that @show chooses to add parenthesis to make expression clearer in
#  the two last cases! (Note: parenthesis on the left-hand side term.
   @show( opnorm( A \ ( B * C ) - A \ B * C ) )
   @show( opnorm( ( A \ B ) * C - A \ B * C ) )

#  Note that for this last case the two answers are completely different.
   A = randn(rng,Float64,n,n);
   B = randn(rng,Float64,n,n);
   C = randn(rng,Float64,n,n);

   @show( opnorm( A * ( B \ C ) - A * B \ C ) )
   @show( opnorm( ( A * B ) \ C - A * B \ C ) )

;
