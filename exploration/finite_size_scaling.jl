cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl");
const global R0 = 1
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

#= This file investigates a potential finite size scaling of the phase space. 
In the same time, I will have to understand the tmax dependence on N =#

# First understand the tmax dependence on N
Ns = round.(Int, logspace(1E2, 1E3, 2, digits=0))
rhos = [1, 2]
T = 0.1
v_sigmas = [(0, 0), (0, 0.1)]#, (0.1, 0), (0.1, 0.1)]
tmax_base100 = 1E2
times = [Int.(logspace(1, tmax_base100 * L, 30, digits=0)) for L in each(Ns)]
R = 8
Ps = Array{Vector{Float64}}(undef, length(Ns), length(rhos), length(v_sigmas), R)
ns = Array{Vector{Float64}}(undef, length(Ns), length(rhos), length(v_sigmas), R)

z = @elapsed for i in each(Ns)
    for p in each(rhos), j in each(v_sigmas)
        Threads.@threads for r in 1:R
            v0, σ = v_sigmas[j]
            N = Ns[i]
            rho = rhos[p]
            L = round(Int, sqrt(N / rho))
            println("σ = $σ, N = $N, v0 = $v0, ρ = $rho, L = $L")
            dt = determine_dt(T, σ, v0, N, rho)
            t = 0.0
            pos, thetas, omegas, psis = initialisation(N, L, L, σ, ["hightemp"])
            token = 1
            ind_neighbours_t0 = get_list_neighbours(pos, N, L, L)
			P_tmp = []
			n_tmp = []
			while t < times[i][end]
                t += dt
                if v0 == 0
                    ind_neighbours = ind_neighbours_t0
                else
                    ind_neighbours = get_list_neighbours(pos, N, L, L)
                end
                pos, thetas = update(pos, thetas, omegas, psis, ind_neighbours, T, v0, N, L, L, dt)
                if t > times[i][token]
                    push!(P_tmp,polarOP(thetas)[1])
                    push!(n_tmp,number_defects(pos, thetas, N, L, L))
                    token += 1
                end
            end
			Ps[i, p, j] = P_tmp
			ns[i, p, j] = n_tmp
        end
    end
end
prinz(z)

# Now plot the results
Ps_avg = mean(Ps, dims=5)[:, :, :, :, 1]  # N, rho, vsig, t, r
ns_avg = mean(ns, dims=5)[:, :, :, :, 1]
plot(xaxis=:log, ylims=(0, 0.2), xlabel=L"1/N", ylabel=L"P", legend=:topleft)
plot!(1 ./ Ns, Ps_avg[:, 1, 1, end])

# Time evolution of the order parameter
p = plot(axis=:log, xlabel=L"t", ylabel=L"P", legend=:topleft);
for n in each(Ns)
    L = round(Int, sqrt(Ns[n] / rho))
    # plot!(times, ns_avg[n, 1,1, :]/L^2)
    plot!(times, Ps_avg[n, 1, :])
end
p