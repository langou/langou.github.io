
   n = 8 
#  Λ = diagm(0 => Float32[2.0^(-i) for i=0:n-1])
   Λ = diagm(0 => Float32[i for i=1:n])
   Λ[n-1:n,n-1:n] = 0.25*[1. 2;-2 1]
   X = rand(rng, n, n)
   A = X * Λ / X

   A, tau = LAPACK.gehrd!(A)
   A = triu(A,-1)

   A0 = copy(A)
   B = A0^2 - ( 1. / 2. ) * A0 + ( 5. / 16. ) * Matrix{Float64}(I,n,n)
   Q,R = qr(B)
   display(R)
#  A1 = Q' * A0 * Q
#  display(A1)
   A1 = triu(Q' * A0 * Q,-1)
   display(A1)

;
   
