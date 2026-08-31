using Printf
using Random
using LinearAlgebra
using Interact
using PlotlyJS

rng = MersenneTwister();

function Arnoldi_Compute_Hk(A,q,kmax)
    H = zeros(kmax, kmax)
    r = copy(q[:,1])
    q[:,1] = r / norm(r)
    
    for k = 1:kmax
        if k > 1
            q[:,k] = r / H[k,k-1]
        end

        r = A * q[:,k] # Multiply by A
        for i=1:k
            # Make vector orthogonal
            H[i,k] = dot(q[:,i], r)
            r -= H[i,k] * q[:,i]
        end

        if k<kmax
            H[k+1,k] = norm(r)
        end
    end
    return H
end

function eig_k(H,k)
    Hk = H[1:k,1:k]
    L = eigen(Hk).values # Eigenvalue estimates at step k        
    return L    
end

n = 128

# Case 1
nh = div(n,2)
@assert 2*nh == n
θ = 2*π*rand(nh)
D = sqrt.(rand(nh)) .* (cos.(θ) + sin.(θ)*im)
D = [D;conj(D)]
D[nh] = 2
D[n] = 2
X = rand(n,nh) + rand(n,nh)*im
X = hcat(X,conj(X))
A = real( X * diagm(0 => D) / X )

# Case 2
#X = rand(n,n)
#D = cos.(LinRange(0,π,n))
#A = real( X * diagm(0 => D) / X )

# Case 3
#A = rand(n,n) / 3

kmax = n # Size of Arnoldi space

q = zeros(n,kmax)
# Random starting vector
q[:,1] = rand(n)
q[:,1] /= norm(q[:,1])

Hk = Arnoldi_Compute_Hk(A,q,kmax);

function lt(z1,z2)
    if abs( abs(z1)-abs(z2) ) > eps(Float32) * (abs(z1) + abs(z2))
        return abs(z1) < abs(z2)
    elseif abs( real(z1-z2) ) > eps(Float32) * (abs(z1) + abs(z2))
        return real(z1) < real(z2)
    else
        return imag(z1) < imag(z2)
    end
end

# Testing accuracy
D = eigen(A).values; D1 = eigen(Hk).values
sort!(D, lt=lt); sort!(D1, lt=lt)
@show norm(D-D1) / norm(D)

#include("../load_plot_pkg.jl")
#output = false

# Pre-computing all the data we are going to plot
L = zeros(ComplexF64,n,kmax)
for k=1:kmax
    L[1:k,k] = eig_k(Hk,k)
end

#@manipulate for k=1:kmax
#    local t1 = scatter(x=real(D),y=imag(D),name="A",mode="markers")
#    local t2 = scatter(x=real(L[1:k,k]),y=imag(L[1:k,k]),name="Arnoldi",mode="markers")
#    local l = PlotlyJS.Layout(xaxis_title="Real part",yaxis_title="Imaginary part",
#        yaxis_scaleanchor="x",yaxis_scaleratio=1,
#        xaxis_range=[-1.1,2.1],yaxis_range=[-1.1,1.1],
#        width=700,height=550,margin_l=80)
#    plot([t1,t2],l)
#end

k = 10
t1 = scatter(x=real(D),y=imag(D),name="A",mode="markers")
t2 = scatter(x=real(L[1:k,k]),y=imag(L[1:k,k]),name="Arnoldi",mode="markers",
    marker=attr(symbol="square-open",line_width=1.3,color="black"))
l = PlotlyJS.Layout(xaxis_title="Real part",yaxis_title="Imaginary part",
    yaxis_scaleanchor="x",yaxis_scaleratio=1,
    xaxis_range=[-1.1,2.1],yaxis_range=[-1.1,1.1],
    width=500,height=300,margin_l=80)
p = plot([t1,t2],l)

#if output
    #plotToPDF(p,"arnoldi")
#end
