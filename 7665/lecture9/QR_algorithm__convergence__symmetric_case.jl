
   using LinearAlgebra
   using Printf
   using Random
   rng = MersenneTwister()

   n = 8 
#  Λ = diagm(0 => Float32[2.0^(-i) for i=0:n-1])
   Λ = diagm(0 => Float32[i for i=1:n])
   X = rand(rng, n, n)
   X,_ = qr(X)
   A = Symmetric(X * Λ * X', :L)

   A0 = copy(A)

   _, A0 = hessenberg(A0)

   residual = [] 
   push!( residual, abs(A0[n,n-1]) )
   display(A0)

   A1 = copy(A0)

## feel welcome to "warm up" the symmetric tridiagonal matrix A so
## so that the smallest eigenvalue is placed in A[n,n] or not

#  for i=1:5
#  global A1, Q, R
#  Q,R = qr(A1)
#  A1 = SymTridiagonal(Symmetric(R*Q))
#  end
#  push!( residual, abs(A1[n,n-1]) )
#  display(A1)

   A2 = copy(A1)
   for i=1:4
   global s, B, Q, R, A2, tau, μ, b, F

## this is a Rayleigh quotient shift 
## it is possible that this shift does not lead to convergence
## see Trefethen and Bau p.222 for a discussion
## see Golub and Van Loan (4th edition) p.461 for a discussion

#  μ = A2[n,n]

## this is the Wilkinson shift
## the Wilkinson shift is defined as that eigenvalue of the 2x2 bottom right block
## that is closer to A[n,n] 
## see Trefethen and Bau p.222 for a discussion
## see Golub and Van Loan (4th edition) p.461 for a discussion

   δ = ( A2[n-1,n-1] - A2[n,n] ) / 2
   μ = A2[n,n] - sign( δ ) * A2[n,n-1]^2 / ( abs( δ ) + sqrt( δ^2 + A2[n,n-1]^2 ) )

## choose how we want to do a step
## 
   case = 3; 

   if ( case == 1 )

##    EITHER case 1 :: a step of Schur iteration step (explicit shift)
##           no optimization done, this is brute force

      B = A2 - μ * Matrix{Float64}(I,n,n)
      Q,R = qr(B)
      A2 = Q' * A2 * Q
      A2 = SymTridiagonal(Symmetric( A2 ))

   elseif ( case == 2 )

##    OR     case 2 :: a step of QR iteration explicit shift 
##           no optimization done, this is brute force

      Q,R = qr(A2 - μ * Matrix{Float64}(I,n,n))
      A2 = R*Q + μ * Matrix{Float64}(I,n,n)
      A2 = SymTridiagonal(Symmetric( A2 ))

   elseif ( case == 3 )

##    OR     case 3 :: a step of QR iteration with implicit Q 
##           no optimization done, this is brute force

      B = A2 - μ * Matrix{Float64}(I,n,n)
      b = B[1:n,1]
      F = qr(b)
      Q = collect(F.Q)
      A2 = Q*A2*Q' 
      A2 = Symmetric( A2, :L )  # it is critical to specify :L here, otherwise the tridiagonal
      _, A2 = hessenberg( A2 )  # reduction algorithm called by "hessenberg" starts at the bottom
      A2 = SymTridiagonal(Symmetric( A2 ))

   end

   push!( residual, abs(A2[n,n-1]) )
   display(A2)

   end

   @printf("\n")
   for i in 1:length(residual)
      @printf("%6.2e ", residual[i] )
   end
   @printf("\n")
;





