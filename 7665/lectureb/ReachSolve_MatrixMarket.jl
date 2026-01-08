
    using LinearAlgebra
    using SparseArrays
    using MatrixMarket
    using Printf
    using Random

#   Read the matrix data
#   file_name = "reach_example.mtx"
#   L = MatrixMarket.mmread(file_name)
#   display(L)
#   for i=1:size(L.nzval)[1]
#       L.nzval[i] = randn(Float64)
#   end

#   file_name = "MatrixMarket/bp__1000.mtx"
#   L = MatrixMarket.mmread(file_name)
#   for i = 1:L.colptr[L.n+1]-1 L.nzval[i] = randn() end
#   L = tril(L)
#   Lf = Matrix(L)
#   for i = 1:L.n Lf[i,i] = L.n + randn() end
#   L = sparse(Lf)
### BP 1000: Original Harwell sparse matrix test collection
### Simplex method basis matrix

    n = 100
    d = 0.1       # Note: L will be about half as sparse since L will be the lower part
#   d = 5.0  / n  # Note: L will be about half as sparse since L will be the lower part
    L = sparse( Matrix( tril(sprandn(Float64,n,n,d)) + Diagonal(randn(Float64,n)+sqrt(n)*ones(Float64,n)) ) )

### Needs to be executed on its own
#   using UnicodePlots
#   UnicodePlots.spy(L)

### Needs to be executed on its own
#   using Plots
#   Plots.spy(L)

### Darve and Wootters' code, does not work for me
#   pyplot()
#   n = L.m # Size of matrix
#   plot(spy(L), xaxis=((0,n+1), 1:n), yaxis=((0,n+1), 1:n), 
#       markersize = 5, clims = (1,2))

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

### the if statement below is added to Darve and Wooters code and is critical
### to make things work, otherwise, we can crash because too many tops, and we
### can have duplicates, not good
        if !w[j] # Skip node that was marked
        w[j] = true

#   Loop over non-zero entries in column j of L
        for k=Lp[j]:Lp[j+1]-1
            if !w[Li[k]] # Skip node that was marked
#               @show j, Li[k]
#               @show xi, top, w = dfs(Li[k], xi, top, w, Lp, Li)
                xi, top, w = dfs(Li[k], xi, top, w, Lp, Li)
                # Call dfs recursively for all nodes in column j
            end
        end

#       Insert j before all nodes inserted so far
        xi[top] = j
#       Move top one position up
        top = top-1
        end

        return xi, top, w
    end
######################################################################

#   n = L.m # Size of matrix

#   Starting node
#   node = 6
#   @assert node>=1
#   @assert node<=n

#   top = n

#   xi = zeros(Int64,n)
#   w  = fill(false,n) # Flag whether a node was visited or not
#   xi, top, w = dfs(node, xi, top, w, L.colptr, L.rowval)
#   @show xi, top, w = dfs(node, xi, top, w, L.colptr, L.rowval)
#   xi = xi[top+1:end] # Keeping only the relevant data

#   println("Reach of node ", node)
#   println(xi)

######################################################################

#   top = n
#   xi = zeros(Int64,n)
#   w  = fill(false,n) # Flag whether a node was visited or not
#   xi, top, w = dfs(4, xi, top, w, L.colptr, L.rowval)
#   xi, top, w = dfs(6, xi, top, w, L.colptr, L.rowval)
#   xi = xi[top+1:end] # Keeping only the relevant data
#   println("Reach of nodes 4 and 6")
#   println(xi)

#   bf = zeros(Float64, n, 1) 
#   bf[4] = randn(Float64)
#   bf[6] = randn(Float64)

    b_nnz = 10 
    if b_nnz > n/2 b_nnz = ceil(Int,n/2) end
#   b_rowval = sort(Int[ rand(1:n) for _=1:b_nnz ])
    b_rowval = randperm(n)[1:b_nnz]
    b = randn(Float64,b_nnz)

    bf = zeros(Float64, n, 1) 
    for i = 1:b_nnz
        bf[b_rowval[i]] = b[i]
    end

    Lf = Matrix(L)

    xf = Lf \ bf
    @printf("check          :: || bf - Lf * xf ||_1 / || Lf ||_1 / || xf ||_1  = %6.2e\n", norm(bf - Lf * xf,1) / norm(Lf,1) / norm(xf,1) )

    top = n
    xi = zeros(Int64,n)
    w  = fill(false,n)   # flag whether a node was visited or not
    for i = 1:b_nnz
        global xi, top, w
        xi, top, w = dfs(b_rowval[i], xi, top, w, L.colptr, L.rowval)
    end
    xi = xi[top+1:L.n] # Keeping only the relevant data
#   println(xi)

nb_passage = 0
    xf = copy(bf)
    for j = xi
        xf[j] = xf[j] / L.nzval[L.colptr[j]]
        for k=L.colptr[j]+1:L.colptr[j+1]-1
global nb_passage
            nb_passage += 1
            xf[L.rowval[k]] -= L.nzval[k] * xf[j]
        end
    end
    @show nb_passage

### Note: we have what we want: a loop in xi, and then we only have at most the density of a row of L
### storage is O(n)

    @printf("check          :: || bf - Lf * xf ||_1 / || Lf ||_1 / || xf ||_1  = %6.2e\n", norm(bf - Lf * xf,1) / norm(Lf,1) / norm(xf,1) )

#   x = zeros(Float64,L.n-top)
#   w  = fill(false,n)   # flag whether a node was visited or not
#   j = 1
#   for i = 1:L.n-top
#       global j
#       if xi[i] == b_rowval[j]
#           x[i] = b[j]
#           j = j + 1
#       end
#   end
#
##### the issue in the next loop is that we need to work on x[L.rowval[k]]
##### L.rowval[k] refers to the global indexes, so does not work for x
#
#   for i = 1:L.n-top
#       local j
#       j = xi[i]
#       x[i] = x[i] / L.nzval[L.colptr[j]]
#       for k=L.colptr[j]+1:L.colptr[j+1]-1
#           x[L.rowval[k]] -= L.nzval[k] * x[i]
#       end
#   end

#### une autre version qui marche (voir emails)

    xf = copy(bf)

nb_passage = 0

for i = 1:n
    if xf[i] != 0
        xf[i] = xf[i] / L.nzval[L.colptr[i]]
        for j = L.colptr[i]+1:L.colptr[i+1]-1
global nb_passage
            nb_passage += 1
            xf[L.rowval[j]] -= L.nzval[j] * xf[i]
        end
    end
end
    @printf("check          :: || bf - Lf * xf ||_1 / || Lf ||_1 / || xf ||_1  = %6.2e\n", norm(bf - Lf * xf,1) / norm(Lf,1) / norm(xf,1) )
    @show nb_passage




;
