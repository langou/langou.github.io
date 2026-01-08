    using LinearAlgebra
    using SparseArrays
    using MatrixMarket
    using Printf

    file_name = "MatrixMarket/bp__1000.mtx"

#   file_name = "MatrixMarket/trefethen_150.mtx"
#   The entries are zero everywhere except for the prime numbers
#   2, 3, 5, 7, ..., 863 along the main diagonal.
#   Then, we add the number 1 in all the positions A[i,j] with
#   |i-j| = 1, 2, 4, 8, ..., 128.

#   Load the matrix using the matrix market format (this is an optional step).
#   Information data for the matrix
#   we set infoonly = true
    rows, cols, entries, mat_format, field, symm = MatrixMarket.mmread(file_name,true)

    println("Number of rows    = ", rows)
    println("Number of columns = ", cols)
    println("Number of entries = ", entries)

#   Read the matrix data
    Amm = MatrixMarket.mmread(file_name)
    @show typeof(Amm)
    @show fieldnames(typeof(Amm))

### Display Amm using julia
#   Amm

### Needs to be executed on its own
#   using UnicodePlots
#   UnicodePlots.spy(Amm)

### Needs to be executed on its own
#   using Plots
#   Plots.spy(Amm)

### Plot the matrix
#   pyplot()
#   n = Amm.m # Size of matrix
#   plot(Plots.spy(Amm), xaxis=((0,n+1), 1:20:n), yaxis=((0,n+1), 1:20:n), markersize = 5, clims = (1,2))

    include("Sparse.jl");

#   Currently, the matrix uses the CSC format (column format).
#   This is often the default, but for this class CSR will be much easier.
    A = SparseMatrixCSR(Amm)

#   println("Matrix A"); println(A)
    println("Indices of the non-zero entries in the 1st row ", A.colval[A.rowptr[1]:A.rowptr[2]-1])
    println("Non-zero entries in the first row              ", A.nzval[A.rowptr[1]:A.rowptr[2]-1])

#   Matrix-vector product using the CSR format
    x = rand(Float64,A.n)
    y = Vector{Float64}(undef,A.m)

#   At the end: y = A * x
    for i=1:A.m
        y[i] = 0.0
        for k=A.rowptr[i]:A.rowptr[i+1]-1
            y[i] += A.nzval[k] * x[A.colval[k]]
        end
    end    

#   println("Error should be equal to 0: ", norm(y - Amm * x))
#   norm(y - Amm * x) == 0 ? "PASS" : "FAIL"

    @printf("check          :: || y - Amm * x ||_1 / || Amm ||_1 / || x ||_1  = %6.2e\n", norm(y - Amm * x,1) / norm(Amm,1) / norm(x,1) )

;

#   check that the overloading of * for our CSR format works too
    @printf("check          :: || A*x - Amm * x ||_1 / || Amm ||_1 / || x ||_1  = %6.2e\n", norm(A*x - Amm * x,1) / norm(Amm,1) / norm(x,1) )

#   Note that "norm(Amm,1)" because Amm is in the julia-supported CSC format
#   for example "norm(A,1)" would not work. We would have to overload norm(.,1) 
#   for our CSR format. 
#   Similarly lu() supported for CSC, triu(), etc.

;
