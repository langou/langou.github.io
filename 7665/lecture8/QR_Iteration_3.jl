
using Printf
using Random
using LinearAlgebra

rng = MersenneTwister(18);

# Size of matrix
n = 7
X = rand(rng, n, n)
Λ = [1e5; 1e3; 5; 4; 3; 2; 1];
Λ = diagm(Λ);
A = X * Λ / X

Tk = copy(A)

display(Tk)
for k=1:6
    global Tk, Rk
    F = qr(Tk); Uk = F.Q; Rk = F.R;
    Tk = Rk * Uk
    display(Tk)
end

