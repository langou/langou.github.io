    using LinearAlgebra
    using SparseArrays
    using MatrixMarket
    using Printf

#   Read the matrix data
    file_name = "reach_example.mtx"
    L = MatrixMarket.mmread(file_name)
    display(L)
    for i=1:size(L.nzval)[1]
        L.nzval[i] = randn(Float64)
    end

########################################################################
"""
dfs: depth-first search construction of the reach. The ordering
of the nodes is important for correctness of the lower triangular solve.

input j: index of node for which reach is computed
input/output xi: reach set of node
input/output top: pointer that locates where the data
                  is stored in xi.
input/output w: used to mark nodes already visited
input Lp: indexing array of L
input Li: row index array for L
"""
    function dfs(j, xi, top, w, Lp, Li)
#   Mark the current node
        w[j] = true

#   Loop over non-zero entries in column j of L
        for k=Lp[j]:Lp[j+1]-1
            if !w[Li[k]] # Skip node that was marked
                @show j, Li[k]
                @show xi, top, w = dfs(Li[k], xi, top, w, Lp, Li)
#               xi, top, w = dfs(Li[k], xi, top, w, Lp, Li)
                # Call dfs recursively for all nodes in column j
            end
        end

#       Insert j before all nodes inserted so far
        xi[top] = j
#       Move top one position up
        top = top-1

        return xi, top, w
    end
######################################################################

    n = L.m # Size of matrix

#   Starting node
    node = 6
    @assert node>=1
    @assert node<=n

    top = n

    xi = zeros(Int64,n)
    w  = fill(false,n) # Flag whether a node was visited or not
    xi, top, w = dfs(node, xi, top, w, L.colptr, L.rowval)
#   @show xi, top, w = dfs(node, xi, top, w, L.colptr, L.rowval)
    xi = xi[top+1:end] # Keeping only the relevant data

    println("Reach of node ", node)
    println(xi)

######################################################################

    top = n
    xi = zeros(Int64,n)
    w  = fill(false,n) # Flag whether a node was visited or not
    xi, top, w = dfs(4, xi, top, w, L.colptr, L.rowval)
    xi, top, w = dfs(6, xi, top, w, L.colptr, L.rowval)
    xi = xi[top+1:end] # Keeping only the relevant data
    println("Reach of nodes 4 and 6")
    println(xi)

#   Matrix-vector product using the CSR format
    b = zeros(Float64, n, 1) 
    b[4] = randn(Float64)
    b[6] = randn(Float64)

    Lf = Matrix(L)

    x = Lf \ b
    @printf("check          :: || b - Lf * x ||_1 / || L ||_1 / || x ||_1  = %6.2e\n", norm(b - Lf * x,1) / opnorm(Lf,1) / norm(x,1) )

    x .= b
    for j = 1:L.m
        x[j] = x[j] / L.nzval[L.colptr[j]]
        for k=L.colptr[j]+1:L.colptr[j+1]-1
            x[L.rowval[k]] -= L.nzval[k] * x[j]
        end
    end

    @printf("check          :: || b - Lf * x ||_1 / || L ||_1 / || x ||_1  = %6.2e\n", norm(b - Lf * x,1) / opnorm(Lf,1) / norm(x,1) )

    x .= b
    for j = xi
        x[j] = x[j] / L.nzval[L.colptr[j]]
        for k=L.colptr[j]+1:L.colptr[j+1]-1
            x[L.rowval[k]] -= L.nzval[k] * x[j]
        end
    end

    @printf("check          :: || b - Lf * x ||_1 / || L ||_1 / || x ||_1  = %6.2e\n", norm(b - Lf * x,1) / opnorm(Lf,1) / norm(x,1) )

;
