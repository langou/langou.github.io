# Set this variable to true if ffmpeg is installed on your computer
#ffmpeg_installed = false
#import Pkg; Pkg.add("PyCall")
using Printf
using Random
using LinearAlgebra
#using PyPlot
#using PyCall
#@pyimport matplotlib.animation as anim

rng = MersenneTwister(18);

# Size of matrix
n = 7
X = rand(rng, n, n)
Λ = [1; 2; 3; 4; 5; 6; 1e-7];
Λ = diagm(Λ);
A = X * Λ / X

Tk = copy(A)

display(Tk)
for k=1:2
    global Tk, Rk
    F = qr(Tk); Uk = F.Q; Rk = F.R;
    Tk = Rk * Uk
    display(Tk)
end

# Print approx. and exact evalues
#n_prt = 3
#Tkdiag = diag(Tk)
#exact = diag(Λ)
#println("Tk    ",Tkdiag[1:n_prt])
#println("Exact ",exact[1:n_prt])
