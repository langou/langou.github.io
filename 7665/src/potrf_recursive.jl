  using GenericLinearAlgebra

  function cholesky_lower_recursive!(A)

    n = size(A,1)

    if n == 1

        A[1,1] = sqrt(A[1,1]);

    else

        n1 = n ÷ 2

        A11 = view(A, 1:n1, 1:n1)
        A21 = view(A, n1+1:n, 1:n1)
        A22 = view(A, n1+1:n, n1+1:n)
 
        cholesky!(Symmetric(A11,:L))
#       LAPACK.potrf!( 'L', A11 )
#       cholesky_lower_recursive!(A11)
 
#       # A21 = A21 / LowerTriangular(A11)'  # this is not working because we are using "views"
        A21 .= A21 / Transpose(LowerTriangular(A11))
#       rdiv!(A21, LowerTriangular(A11)')
#       BLAS.trsm!( 'R', 'L', 'T', 'N', 1.0e+00, A11, A21)
 
#       # A22 = A22 - A21 * A21' # this is not working because we are using "views"
        A22 .= A22 - A21 * A21'
#       GenericLinearAlgebra.rankUpdate!(Hermitian(A22, :L), A21, -1)
#       BLAS.syrk!( 'L', 'N', -1.0e+00, A21, 1.0e+00, A22)

        cholesky!(Symmetric(A22,:L))
#       LAPACK.potrf!( 'L', A22 )
#       cholesky_lower_recursive!(A22)

    end
    
  end
