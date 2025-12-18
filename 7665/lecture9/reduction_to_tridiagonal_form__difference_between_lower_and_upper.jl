
   using LinearAlgebra

   n = 8 
   Λ = diagm(0 => Float32[i for i=1:n])
   X = rand(rng, n, n)
   X,_ = qr(X)

   A = X * Λ * X';

   A = A';

   A_L = Symmetric( A, :L)
   _, T_L = hessenberg(A_L)

   A_U = Symmetric( A )
   _, T_U = hessenberg(A_U)

   @printf("T_L = \n")
   display(T_L)

   @printf("T_U = \n")
   display(T_U)
