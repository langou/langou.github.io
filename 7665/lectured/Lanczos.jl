using Printf
using Random
using LinearAlgebra
using PlotlyJS


function compute_Tk(A,q,kmax)
    T = zeros(kmax, kmax)
    r = copy(q[:,1])
    beta = norm(r)
    q1 = zeros(n)

    for k = 1:kmax        
        q0 = copy(q1)
        q1 = r/beta
        r = A * q1
        alpha = dot(q1, r)
        T[k,k] = alpha
        if k > 1
            T[k-1,k] = beta
            T[k,k-1] = beta
        end
        r = r - alpha*q1 - beta*q0
        beta = norm(r)
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

Q, = qr(rand(n,n))

# D = cos.(LinRange(0,π,n))
# D[div(n,2)] = 2
# D = cos.(LinRange(0,π/2,n))
D = LinRange(-1,1,n).^3


x = LinRange(0,1,n)
# D = x
# D = x.^2
# D = x.^2/2 - x.^3/3
# D = x.^5/5 - x.^4/2 + x.^3/3
# D = x.^9/9 - x.^8/2 + 6/7*x.^7 - 2/3*x.^6 + x.^5/5

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
    local x = eig_k(Tk,k)
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
    #plotToPDF(p,"lanczos_cro1")
#end
