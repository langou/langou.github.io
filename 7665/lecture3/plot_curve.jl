using Plots


N = zeros(Int,1000)
complete = zeros(Float64,1000)
partial = zeros(Float64,1000)

i=0;
for n = 0:1:999
   global i += 1
   N[i] = n
   complete[i] = 1.5 * n ^ ( 3/4 )
end

i=0;
for n = 0:1:999
   global i += 1
   N[i] = n
   partial[i] = 2^n
end

function new_plot()
    plot(xlabel = "x", ylabel = "f(x)",
        xlims = (0,Inf), ylims = (-Inf, 1))
end 

p = new_plot()
plot!(p, x -> x^2, 0, 1)
plot!(p, x -> x^3, 0, 1)

#plot!(p, N, complete, lw=3)
#plot!(p, N, partial, lw=3)

