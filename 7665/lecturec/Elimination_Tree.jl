    using MatrixMarket
    using LinearAlgebra
    using SparseArrays
    using Plots
    using Printf

#   file_name = "MatrixMarket/bp__1000.mtx"
#   L = MatrixMarket.mmread(file_name)
#   for i = 1:L.colptr[L.n+1]-1 L.nzval[i] = randn() end
#   L = tril(L)
#   Lf = Matrix(L)
#   for i = 1:L.n Lf[i,i] = L.n + randn() end
#   L = sparse(Lf)
#   A_CSC = L
### BP 1000: Original Harwell sparse matrix test collection
### Simplex method basis matrix

#   n = 100
#   d = 0.1       # Note: L will be about half as sparse since L will be the lower part
#   d = 5.0  / n  # Note: L will be about half as sparse since L will be the lower part
#   L = sparse( Matrix( tril(sprandn(Float64,n,n,d)) + Diagonal(randn(Float64,n)+sqrt(n)*ones(Float64,n)) ) )
#   A_CSC = L


### Read the matrix data
    file_name = "example_A.mtx"
    A_CSC = SparseMatrixCSC{Float64,Int64}( MatrixMarket.mmread(file_name) )

### Graph representing the nonzero entries in A 
### dot -Tpdf A.dot > A.pdf
    open("A.dot", "w") do f
        @printf f "digraph G {\n"
        for i=1:A_CSC.m
            @printf f "%d;\n" i
        end
        for j=1:A_CSC.m
            for k=A_CSC.colptr[j]+1:A_CSC.colptr[j+1]-1
               i = A_CSC.rowval[k]
               @printf f "%d -> %d;\n" j i
            end
        end
        @printf f "label=\"Graph representing the nonzero entries in A\";\n"
        @printf f "labelloc=\"t\";\n"    
        @printf f "}"
    end

#   A is symmetric; we only read the lower triangular part
    for i = 1:A_CSC.m A_CSC[i,i] = 10.0 end
    A = Matrix(A_CSC)
    A = A + A'
    A_CSC = SparseMatrixCSC(A)

### cheat by using dense Cholesky
    A = Matrix(A_CSC)
    C = cholesky(A).L
    C_CSC_dense = SparseMatrixCSC(C)
    display( C_CSC_dense )

### Graph representing the nonzero entries in L 
### dot -Tpdf C_CSC_dense.dot > C_CSC_dense.pdf
    open("C_CSC_dense.dot", "w") do f
        @printf f "digraph G {\n"
        for i=1:C_CSC_dense.m
            @printf f "%d;\n" i
        end
        for j=1:C_CSC_dense.m
            for k=C_CSC_dense.colptr[j]+1:C_CSC_dense.colptr[j+1]-1
               i = C_CSC_dense.rowval[k]
               @printf f "%d -> %d;\n" j i
            end
        end
        @printf f "label=\"Graph representing the nonzero entries in C_CSC_dense\";\n"
        @printf f "labelloc=\"t\";\n"    
        @printf f "}"
    end

### cheat by using sparse Cholesky (CHOLMOD)
    C=cholesky(A_CSC, perm=1:A_CSC.m)  # I am not sure this'll work all the time, but cholmod seems to listen to the user-defined permutation in this case
    C_CSC_cholmode = sparse(C.L)
    display( C_CSC_cholmode )

### only keep the lower part of A from now on
    A_CSC = tril(A_CSC)

    include("Sparse.jl")
    A = SparseMatrixCSR(A_CSC)
    n = A.m # Size of matrix

    display(A_CSC)

###############################################################################
"""Function that calculates the elimination tree given a CSR structure"""
    function etree(rowptr::Vector{Int}, colval::Vector{Int})
        n = length(rowptr) - 1
        parent   = fill(-1,n) # e-tree information
        ancestor = fill(-1,n) # ancestor information to reduce the running time

        # We compute the elimination tree
        for i=1:n
            parent[i]   = -1 # Initialize to -1
            ancestor[i] = -1

            for p = rowptr[i]:rowptr[i+1]-1
                j = colval[p] # column index
                # Traverse row i and stop before the diagonal
                while j != -1 && j < i
                    jnext = ancestor[j] # Search for the root
                    ancestor[j] = i   # Update ancestor for efficiency
                    if jnext == -1    # We have found a root
                        parent[j] = i # Connect to i
                    end
                    j = jnext
                end
            end
            @show i
            @show parent
            @show ancestor
        end
        return parent
    end
###############################################################################
    function row_sparsity(rowptr, colval, parent, i)
        n = length(rowptr) - 1
        s = Vector{Int64}(undef,n)
        w  = fill(false,n) # Used to mark points as visited
        w[i] = true
        len = 1

        for p = rowptr[i]:rowptr[i+1]-1
            j = colval[p] # column index
            # Traverse row i and stop before the diagonal
            while !w[j] && j < i # Stop when marked node is found
                s[len] = j  # Add column j to row i
                w[j] = true # Mark node j
                len += 1
                j = parent[j] # Move to parent in e-tree
            end
        end
        
        len -= 1
        return s[1:len]
    end
###############################################################################
#   Compute the elimination tree
    parent_tree = etree(A.rowptr, A.colval)

    @show parent_tree

# We write the elimination tree to a DOT file.
# Use Graphviz to see the graph of the tree.
# Open the file "etree.dot" using Graphviz.
# http://www.graphviz.org/
# Command line:
# xdot, or
# dot -Tpdf etree.dot > etree.pdf
#
    n = A.m # Size of matrix
    open("etree.dot", "w") do f
        @printf f "digraph G {\n"
        for k=1:n
            @printf f "%d;\n" k
            if parent_tree[k] != -1
                @printf f "%d -> %d;\n" k parent_tree[k]
            end
        end
        @printf f "label=\"Elimination tree of A\";\n"
        @printf f "labelloc=\"t\";\n"    
        @printf f "}"
    end
###############################################################################

n = A.m # Size of matrix

# Select the index of the row subtree

k = 10
k = 8
@assert k>=1
@assert k<=n

# Compute the row sparsity pattern
s = row_sparsity(A.rowptr, A.colval, parent_tree, k)

# We write the row sub-tree to a DOT file.
# Use Graphviz to see the graph of the row sub-tree.
# Open the file "row_subtree.dot".
# dot -Tpdf row_subtree.dot > row_subtree.pdf

open("row_subtree.dot", "w") do f
    @printf f "digraph G {\n"
    @printf f "%d; \n" k
    for i=1:length(s)
        @printf f "%d -> %d;\n" s[i] parent_tree[s[i]]
    end
    @printf f "label=\"Row sub-tree of node %d\";\n" k
    @printf f "labelloc=\"t\";\n"    
    @printf f "}"
end

###############################################################################
# factorization dense

    C = tril(Matrix(A_CSC))

    for k=1:A.m

        for j=1:k-1

            C[k,j] = C[k,j] / C[j,j]

            for i=j+1:k

                C[k,i] = C[k,i] - C[i,j] * C[k,j]

            end

        end

        for j=1:k-1

            C[k,k] = C[k,k] - C[j,k] * C[j,k]

        end

        C[k,k] = sqrt( C[k,k] )

    end

###############################################################################
# factorization using "row_sparsity"

    C = tril(Matrix(A_CSC))

    for k=1:A.m

        global s

#       je ne comprend pas vraiment pourquoi le sort( ) est nécessaire ici
#       mais bon, ça marche avec le sort et ça ne marche pas sans
        s = sort( row_sparsity(A.rowptr, A.colval, parent_tree, k) )

#       for j=1:k-1
#       for j=s
        for jj=1:size(s)[1]

            C[k,s[jj]] = C[k,s[jj]] / C[s[jj],s[jj]]

#           for i=j+1:k
#           for i=s[jj]+1:k
            for ii=jj+1:size(s)[1]

#               C[k,i] = C[k,i] - C[i,s[jj]] * C[k,s[jj]]
                C[k,s[ii]] = C[k,s[ii]] - C[s[ii],s[jj]] * C[k,s[jj]]

            end
            C[k,k] = C[k,k] - C[k,s[jj]] * C[k,s[jj]]

        end

#       for j=1:k-1
        for j=s

            C[k,k] = C[k,k] - C[j,k] * C[j,k]

        end

        C[k,k] = sqrt( C[k,k] )

    end

###############################################################################
# factorization dense

    C_CSR = SparseMatrixCSR(C_CSC_dense)
    C_CSR.nzval = zeros(C_CSR.nzval.size[1])

    A_CSR = SparseMatrixCSR(A_CSC)

### copy A in C
    l = 1
    for j=1:C_CSR.m
        for k=C_CSR.rowptr[j]:C_CSR.rowptr[j+1]-1
            global l
            if ( C_CSR.colval[k] == A_CSR.colval[l] )
                C_CSR.nzval[k] = A_CSR.nzval[l]
                l = l+1
            end
        end
    end

    for k=1:A.m

        for j=1:k-1

            C[k,j] = C[k,j] / C[j,j]

            for i=j+1:k

#               the question is how to access C[i,j] since we store the matrix by row
#               I am not sure how to access C[i,j]
#               note that when we did the triangular solve for sparse matrices, then
#               C was stored in CSC. Now we have C stored in CSR, so what?
                C[k,i] = C[k,i] - C[i,j] * C[k,j]

#               maybe this Cholesky is not the good variant to use for Cholesky
#               we would like a variant that only does row access or column access
#               this variant seems to need row and column access

            end

        end

        for j=1:k-1

            C[k,k] = C[k,k] - C[j,k] * C[j,k]

        end

        C[k,k] = sqrt( C[k,k] )

    end



