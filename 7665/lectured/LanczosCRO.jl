using Printf
using Random
using LinearAlgebra

function compute_Tk(A,q,kmax)
    T = zeros(kmax, kmax)
    r = copy(q[:,1])
    q[:,1] = r / norm(r)
    
    for k = 1:kmax
        if k > 1
            q[:,k] = r / T[k,k-1]
        end

        r = A * q[:,k] # Multiply by A
        for i=1:k
            # Make vector orthogonal
            T[i,k] = dot(q[:,i], r)
            r -= T[i,k] * q[:,i]
        end

        if k<kmax
            T[k+1,k] = norm(r)
        end
    end
    return T
end

function eig_k(T,k)
    Tk = T[1:k,1:k]
    L = eigen(Tk).values # Eigenvalue estimates at step k        
    return L    
end

n = 32
rng = MersenneTwister(18);

Q = qr(rand(n,n)).Q
#D = cos.(LinRange(0,π,n))
#D[div(n,2)] = 2
D = LinRange(-1,1,n).^3

A = Q * diagm(0 => D) * Q'

q = zeros(n,n)
# Random starting vector
q[:,1] = rand(n)
q[:,1] /= norm(q[:,1])

Tk = compute_Tk(A,q,n);

# Testing accuracy
D = eigen(A).values; D1 = eigen(Tk).values
sort!(D); sort!(D1)
@show norm(D-D1) / norm(D)

#include("../load_plot_pkg.jl")
#output = false

L = zeros(n,n)
for k=1:n
    x = eig_k(Tk,k)
    sort!(x)
    for l=1:k
        if l%2 == 1
            L[l,k] = x[(l+1)>>1]
        else
            L[l,k] = x[k-(l-2)>>1]
        end
    end
end

t = Array{PlotlyJS.AbstractTrace,1}()
for k=1:n
    push!(t,scatter(x=k:n,y=L[k,k:n],mode="lines",showlegend=false))
end
push!(t,scatter(x=n*ones(n),y=D,mode="markers",name="Exact"))

p = plot(t,Layout(xaxis_title="Iteration",yaxis_title="Eigenvalue estimates",
        width=500,height=400,margin_l=80))

#if output
    #plotToPDF(p,"lanczos_cro2")
#end
