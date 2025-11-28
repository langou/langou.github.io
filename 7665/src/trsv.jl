function trsv_lnn_level0_jivariant!(A, x)
#   trsv for { lower, no transpose, non unit } 
#   triangular matrix A
#   vector x is inout
#
#   This is the "column" variant because we stream A by column
#   This is fast because the row index is the inner loop
#
    n = length(x)
    for j = 1:n
        x[j] = x[j] / A[j,j]
        for i = j+1:n
            x[i] -= A[i,j] * x[j]
        end
    end
end

function trsv_lnn_level0_ijvariant!(A, x) 
#   trsv for { lower, no transpose, non unit } triangular matrix A
#   triangular matrix A
#   vector x is inout
#
#   This is the "row" variant because we stream A by row
#   This is the slow implementation
#
    n = length(x)
    for i=1:n
        for j = 1:i-1
            x[i] -= A[i,j] * x[j]
        end
        x[i] = x[i] / A[i,i]
    end
end
