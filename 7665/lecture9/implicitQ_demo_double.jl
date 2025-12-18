
   using LinearAlgebra
   using Printf
   using Random
   rng = MersenneTwister()

   include("larfg.jl")

#  the goal is to show that implicit single shift work
#  we are not claiming that this is a good shift!

   n = 8 
#  Λ = diagm(0 => Float32[2.0^(-i) for i=0:n-1])
   Λ = diagm(0 => Float32[i for i=1:n])
   Λ[n-1:n,n-1:n] = 0.25*[1. 2;-2 1]
   X = rand(rng, n, n)
   A = X * Λ / X
   A0 = copy(A)
   A0, tau = LAPACK.gehrd!(A0)
   A0 = triu(A0,-1)

#  let us do a double shift with
#      A^2 - s * A + t * I
   s = -4.
   t = -3.

## step of Schur iteration step (no implicit Q)
   A1 = copy(A0)
   B = A1^2 - s * A1 + t * Matrix{Float64}(I,n,n)
   Q,R = qr(B)
   A1 = Q' * A1 * Q
## Question: why is A1 upper Hessenberg?
## Not so easy to explain!
#  display(A1)
#  A1 = triu(A1,-1)
#  display(A1)
   display(UpperHessenberg(A1))

## step of QR iteration with implicit shift 
   A2 = copy(A0)
   B = A0^2 - s * A0 + t * Matrix{Float64}(I,n,n)
   b = B[1:n,1]
   F = qr(b)
   Q = collect(F.Q)  # you really want the full Q here
   Q2, H2 = hessenberg( Q*A2*Q' )
   display(H2)

## an "optimized" implementation - step 1 - compute only the first column and in a smart way
   A2 = copy(A0)
   v = zeros(n,1);
   v[1:3] = [ A0[1,1]*A0[1,1] + A0[1,2]*A0[2,1] - s*A0[1,1] + t;
         A0[2,1]*(A0[1,1]+A0[2,2]-s);
         A0[2,1]*A0[3,2] ]

   τ = larfg!(v)

   v[1] = 1.0

#  A2 = (I(n) - v * τ * v' ) * A2 * (I(n) - v * τ * v' )

#  A2 = A2 - v * τ * ( v' * A2 )
#  A2 = A2 - ( A2 * v ) * τ * v'

   A2 .-= v * τ * ( v' * A2 )
   A2 .-= ( A2 * v ) * τ * v'

   Q2, H2 = hessenberg( A2 )
   display(H2)


## an "optimized" implementation - step 2 - v is of size 3

   A2 = copy(A0)
   v = [ A0[1,1]*A0[1,1] + A0[1,2]*A0[2,1] - s*A0[1,1] + t;
         A0[2,1]*(A0[1,1]+A0[2,2]-s);
         A0[2,1]*A0[3,2] ]

   τ = larfg!(v)

   v[1] = 1.0

   A2[1:3,1:n].-= v * τ * ( v' * A2[1:3,1:n] )
   A2[1:n,1:3].-= ( A2[1:n,1:3] * v ) * τ * v'

   Q2, H2 = hessenberg( A2 )
   display(H2)

;

## an "optimized" implementation - step 3 - recognize that A2 is Hessenberg with 
## a bulge so we can O(n^2) bulge chasing algorithm to compute H2

#### => see code gees.jl

;


## Question: why is A1 upper Hessenberg?
## Not so easy to explain!

#  Either we say that a polynomial of degree 2 can be factored in two polynomials
#  of degree 1. Or we do as follows

#  B = A1 - 3. I
#  Q,R = qr(B)
#
#  ==> QR = A - 3I
#  ==> R  = Q'( A - 3 I )
#  ==> RQ = Q'( A - 3 I )Q
#  ==> RQ = Q' A Q - 3 I
#  ==> Q' A Q = RQ + 3 I
#  ==> this explains Hessenberg structure of Q' A Q 

#  B = A1^2 + 4 A1 - 3 I 
#  Q,R = qr(B)
#
#  ==> QR = A1^2 + 4 A1 - 3 I
#  ==> R  = Q' A1^2 + 4 Q' A1 - 3 I
#  ==> RQ = ( Q' A1 Q )^2 + 4 Q' A1 Q - 3 I
#  ==> since RQ is 2-Hessenberg, it must be that Q' A1 Q is Hessenberg
#
