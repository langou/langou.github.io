    using LinearAlgebra
    using SparseArrays
    using Plots
    using Random
    using MatrixMarket
    using UnicodePlots
    using Printf


    file_name = "MatrixMarket/bp__1000.mtx"
#   BP 1000: Original Harwell sparse matrix test collection
#   Simplex method basis matrix

#   file_name = "MatrixMarket/trefethen_150.mtx"
#   The entries are zero everywhere except for the prime numbers
#   2, 3, 5, 7, ..., 863 along the main diagonal.
#   Then, we add the number 1 in all the positions A[i,j] with
#   |i-j| = 1, 2, 4, 8, ..., 128.

#   Load the matrix using the matrix market format (this is an optional step).
#   Information data for the matrix
    println(file_name)
    @show MatrixMarket.mmread(file_name,true)
    rows, cols, entries, mat_format, field, symm = MatrixMarket.mmread(file_name,true)

    println("Number of rows    = ", rows)
    println("Number of columns = ", cols)
    println("Number of entries = ", entries)

#   Read the matrix data
    Amm = MatrixMarket.mmread(file_name)

#   Keep only the lower triangular part
    Lmm = tril(Amm)

#   Add nonzero entries on the diagonal
    Lf = Matrix(Lmm)
    for i = 1:size(Lmm)[1]
        Lf[i,i] = 1. + randn()
    end
    Lmm = sparse(Lf)

    display(Lmm)

    include("Sparse.jl")

    L = SparseMatrixCSR(Lmm);

#   Matrix-vector product using the CSR format
    Random.seed!(2026)
    b = Float64[ rand(-9:9) for _=1:L.n ]
    b = L*b

    x = Vector{Float64}(undef,L.m)

### This is the row version of triangular solve using CSR format for L
#   Solving: L x = b
#   L is lower triangular
    x .= b
    for i=1:L.m
        for k=L.rowptr[i]:L.rowptr[i+1]-2
            x[i] -= L.nzval[k] * x[L.colval[k]]
        end
        x[i] /= L.nzval[L.rowptr[i+1]-1]
#
# Note: possible checks
#
#        we could check that 
#            L.colval[ [L.rowptr[i+1]-1] ] == i
#        here we are assuming it and so we are assuming that L.nzval[L.rowptr[i+1]-1] is Aii
#
#        we could check that 
#            L.nzval[L.rowptr[i+1]-1] is not 0
#        before dividing
#
    end    

    @printf("check CSR row variant of TRSV           :: || b - Lf * x ||_1 / || L ||_1 / || x ||_1  = %6.2e\n", norm(b - Lf * x,1) / opnorm(Lf,1) / norm(x,1) )
#   println("Error should be equal to 0: ", norm(b - L * x))
#   @show norm(b - L * x) == 0 ? "PASS" : "FAIL"


### This is the col version of triangular solve using CSC format for L
#   Solving: L x = b
#   L is lower triangular

    x .= b
    for j = 1:Lmm.m
        x[j] = x[j] / Lmm.nzval[Lmm.colptr[j]]
        for k=Lmm.colptr[j]+1:Lmm.colptr[j+1]-1
            x[Lmm.rowval[k]] -= Lmm.nzval[k] * x[j]
        end
    end
    @printf("check CSC col variant of TRSV           :: || b - Lf * x ||_1 / || L ||_1 / || x ||_1  = %6.2e\n", norm(b - Lf * x,1) / opnorm(Lf,1) / norm(x,1) )
#   println("Error should be equal to 0: ", norm(b - L * x))
#   @show norm(b - L * x) == 0 ? "PASS" : "FAIL"

################################# Julien from here on

    Lf = Matrix(Lmm)

    x = Lmm \ b
    @printf("check \\ of julia for CSC                :: || b - Lf * x ||_1 / || L ||_1 / || x ||_1  = %6.2e\n", norm(b - Lf * x,1) / opnorm(Lf,1) / norm(x,1) )

    x = Lf \ b
    @printf("check \\ of julia for dense              :: || b - Lf * x ||_1 / || L ||_1 / || x ||_1  = %6.2e\n", norm(b - Lf * x,1) / opnorm(Lf,1) / norm(x,1) )

    include("trtrs_jl.jl")

    x .= b
    trtrsRow(Lf, x)
    @printf("check homemade trtrsRow for dense-dense :: || b - Lf * x ||_1 / || L ||_1 / || x ||_1  = %6.2e\n", norm(b - Lf * x,1) / opnorm(Lf,1) / norm(x,1) )

    x .= b
    trtrs(Lf, x)
    @printf("check homemade trtrs    for dense-dense :: || b - Lf * x ||_1 / || L ||_1 / || x ||_1  = %6.2e\n", norm(b - Lf * x,1) / opnorm(Lf,1) / norm(x,1) )

;


