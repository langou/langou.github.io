  using GenericLinearAlgebra

  function cholesky_lower_leftlooking_level0!(A)
    n = size(A,1)
    for j=1:n
      for k=1:j-1
        for i=j:n
          A[i,j] -= A[i,k] * A[j,k]
        end
      end
      A[j,j] = sqrt(A[j,j])
      for i=j+1:n
        A[i,j] /= A[j,j]
      end
    end
  end

  function cholesky_lower_leftlooking_level2!(A)
    n = size(A,1)
    for j=1:n
      A[j:n,j:j] .-= A[j:n,1:j-1] * A[j:j,1:j-1]'
      A[j,j] = sqrt(A[j,j])
      A[j+1:n,j] ./= A[j,j]
    end
  end

  function cholesky_lower_leftlooking_level3!(A)

    n = size(A,1)

    nb = 64 

    for j=1:nb:n

#     k is n+1 when n % nb = 0 and so I3 will be ∅ in that case

      k = min( j+nb-1, n) + 1

      I1 = 1:j-1
      I2 = j:k-1
      I3 = k:n

      A21 = view(A,I2,I1)
      A31 = view(A,I3,I1)
      A22 = view(A,I2,I2)
      A32 = view(A,I3,I2)

#     these two steps could be combine in one by adding more views, but that would destroy A22

#     A22 .-= A21 * A21'
#     GenericLinearAlgebra.rankUpdate!(Hermitian(A22, :L), A21, -1)
      BLAS.syrk!( 'L', 'N', -1.0e+00, A21, 1.0e+00, A22)

#     A32 .-= A31 * A21'
#     mul!(A32, A31, Transpose(A21), -1.0e00, 1.0e00)
      BLAS.gemm!( 'N', 'T', -1.0e+00, A31, A21, 1.0e+00, A32)

#     cholesky_lower_leftlooking_level2!(A22)
#     cholesky!(Symmetric(A22,:L))
      LAPACK.potrf!( 'L', A22 )

#     A32 .= A32 / Transpose(LowerTriangular(A22))
#     rdiv!(A32, LowerTriangular(A22)')
      BLAS.trsm!( 'R', 'L', 'T', 'N', 1.0e+00, A22, A32)

    end

  end


  function cholesky_lower_rightlooking_level3!(A)

    n = size(A,1)

    nb = 64

    for j=1:nb:n

      k = min( j+nb-1, n) + 1

      I1 = 1:j-1
      I2 = j:k-1
      I3 = k:n

      A22 = view(A,I2,I2)
      A32 = view(A,I3,I2)
      A33 = view(A,I3,I3)

      cholesky_lower_leftlooking_level2!(A22)
#     cholesky!(Symmetric(A22,:L))
#     LAPACK.potrf!( 'L', A22 )

#     A32 .= A32 / Transpose(LowerTriangular(A22))
#     rdiv!(A32, LowerTriangular(A22)')
      BLAS.trsm!( 'R', 'L', 'T', 'N', 1.0e+00, A22, A32)

#     A33 .= A33 - A32 * A32'
#     GenericLinearAlgebra.rankUpdate!(Hermitian(A33, :L), A32, -1)
      BLAS.syrk!( 'L', 'N', -1.0e+00, A32, 1.0e+00, A33)

    end

  end

  function cholesky_lower_bordered_level3!(A)

    n = size(A,1)

    nb = 11

    for j=1:nb:n

      k = min( j+nb-1, n) + 1

      I1 = 1:j-1
      I2 = j:k-1

      A11 = view(A,I1,I1)
      A21 = view(A,I2,I1)
      A22 = view(A,I2,I2)

      A21 .= A21 / Transpose(LowerTriangular(A11))
#     rdiv!(A21, LowerTriangular(A11)')
#     BLAS.trsm!( 'R', 'L', 'T', 'N', 1.0e+00, A11, A21)

      A22 .= A22 - A21 * A21'
#     GenericLinearAlgebra.rankUpdate!(Hermitian(A22, :L), A21, -1)
#     BLAS.syrk!( 'L', 'N', -1.0e+00, A21, 1.0e+00, A22)

#     cholesky_lower_leftlooking_level2!(A22)
      cholesky!(Symmetric(A22,:L))
#     LAPACK.potrf!( 'L', A22 )

    end

  end

