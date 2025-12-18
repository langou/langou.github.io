
   using LinearAlgebra
   using Printf
   using Random
   rng = MersenneTwister()

   n = 8 
#  Λ = diagm(0 => Float32[2.0^(-i) for i=0:n-1])
   Λ = diagm(0 => Float32[i for i=1:n])
   Λ[n-1:n,n-1:n] = 0.25*[1. 2;-2 1]
   X = rand(rng, n, n)
   A = X * Λ / X
   A0 = copy(A)
   A0, tau = LAPACK.gehrd!(A0)
   A0 = triu(A0,-1)

#  here we do 10 step of QR algorithm to "prep" the lower 2x2 block to be our complex conjugate pair
#  we could have done the experiment with a :
   A1 = copy(A0)
   for i=1:10
   global A1, Q, R
   Q,R = qr(A1)
   A1 = R*Q
   end

function double_shift_st(A)
    a = A[1,1]
    b = A[1,2]
    c = A[2,1]
    d = A[2,2]
    s = a+d       # sum
    t = a*d - b*c # product
    return s,t
end


## Note difference of B between these two

## this is not a particularly interesting shift
## we observe that B has no particular shape
   B = A0^2 + 4. * A0 - 3. * Matrix{Float64}(I,n,n)
   display(B)

## this is a double shift on the bottom right 2x2 block
## we note that B has a "spike", this is part of the idea in early deflation
   s,t = double_shift_st(A0[n-1:n,n-1:n])
   B = A0^2 - s * A0 + t * Matrix{Float64}(I,n,n)
   display(B)

## this is mostly an observation

## can we prove the shape of B?
## is this useful? how can we exploit this?

;




