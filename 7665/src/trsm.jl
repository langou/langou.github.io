
function trsm_llnn!(A, x)

#   trsm for { left, lower, no transpose, non unit } triangular matrix A
#   this is the variant that streams along the columns of A
 
    n,k = size(x)
    for j = 1:n
        x[j,1:k] = x[j,1:k] / A[j,j]
        for i = j+1:n
            x[i,1:k] -= A[i,j] * x[j,1:k]
        end
    end
end


function trsm_llnn_recursive!(A, x)

#   trsm for { left, lower, no transpose, non unit } triangular matrix A
#   this is the recursive variant, we recurse on the dimension of A

#   we stop the recursion at n=1 because this is good enough, we could try to
#   see if we can improve performance by stopping the recursion earlier and
#   backing up on a non-recursive version
 
    n,k = size(x)

    if n == 1

        x[1,1:k] = x[1,1:k] / A[1,1]

    else

        n1 = n ÷ 2

        A11 = view(A, 1:n1, 1:n1)
        A21 = view(A, n1+1:n, 1:n1)
        A22 = view(A, n1+1:n, n1+1:n)
 
        x1 = view(x, 1:n1, 1:k)
        x2 = view(x, n1+1:n, 1:k)
 
        trsm_llnn_recursive!(A11, x1)
#
#       # this uses "broadcast", this makes copy of copies of copy, etc
#       x[n1+1:n, 1:k] = x[n1+1:n, 1:k] - A[n1+1:n, 1:n1] * x[1:n1, 1:k]
#       x[n1+1:n, 1:k] -= A[n1+1:n, 1:n1] * x[1:n1, 1:k]
#
#       # a little better
#       x[n1+1:n, 1:k] .= x[n1+1:n, 1:k] - A[n1+1:n, 1:n1] * x[1:n1, 1:k]
#       x[n1+1:n, 1:k] .-= A[n1+1:n, 1:n1] * x[1:n1, 1:k]
#
#       # a little better using views, (still need to put A21*x1 somewhere)
#       x2 .= x2 - A21 * x1
#
#       BLAS.gemm!( 'N', 'N', -1.0, A21, x1, 1.0, x2 )
#
#       # mul!() uses multiple dispatch, the way to go, can call cublas, dense, sparse, 
#       # other ones to know: rdiv! ldiv! etc.
        mul!( x2, A21, x1, -1.0, 1.0 )
#
        trsm_llnn_recursive!(A22, x2)

    end
    
end





