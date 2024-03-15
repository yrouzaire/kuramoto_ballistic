cd("/Users/yrouzaire/Documents/Recherche/GitHub/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl");
# const global R0 = 1
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

## ---------------- Analysis ---------------- ##
# Visualisation scan in phase space
filename = "data/FSS_to_determine_transition.jld2"
@load filename Ntargets v0sigs Ps Cs ns xis params_init aspect_ratio times tmax T comments rho rhoc runtimes R
hrun(runtimes)
Ntargets
v0sigs
v0s = [v for (v,x) in v0sigs]
Ps_avg = nanmean(Ps, 4)[:,:,:,1]
ns_avg = nanmean(ns, 4)[:,:,:,1]
xis_avg = nanmean(xis, 4)[:,:,:,1]

indices = [];
for r in 1:R
    try Cs[:,:,:,r]
		push!(indices, r)
    catch;
    end
end;
indices



Cs_avg = Array{Vector}(undef, length(v0sigs), length(Ntargets), length(times))
for i in 1:length(v0sigs), j in each(Ntargets), k in 1:length(times)
    Cs_avg[i, j, k] = mean([Cs[i, j, k, r] for r in indices])
end

# # ---------------- Plotting ---------------- ##
p=plot(axis=:log, legend=false, xlabel=L"N", ylabel=L"P")
for i in each(v0sigs)
    plot!(Ntargets, Ps_avg[i, :, 30], label="v0sig = $(v0sigs[i])", rib=0, m=true)
end
ylims!(0.03, 1.1)
yticks!([0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], ["0.03", "0.04", "0.05", "", "0.07", "", "", "0.1", "0.2", "0.3", "0.4", "0.5", "", "0.7", "", "", "1"])
p


## ---------------- Investigating the time saturation of the magnetisation ---------------- ##
## ---------------- Investigating the time saturation of the magnetisation ---------------- ##
## ---------------- Investigating the time saturation of the magnetisation ---------------- ##
## ---------------- Investigating the time saturation of the magnetisation ---------------- ##
p=plot(axis=:log, legend=false, xlabel=L"t", ylabel=L"P")
for i in 6#each(v0sigs)
    for j in each(Ntargets)
        plot!(times, Ps_avg[i, j, :], label="v0sig = $(v0sigs[i])", rib=0, m=true)
        N = Ntargets[j]
        scatter!((min(0.2 * N * log(N), tmax), Ps_avg[i, j, end]), label="", m=(5, :black))
    end
end
xticks!([1, 10, 100, 1000, 10000, 100000], [L"10^0", L"10^1", L"10^2", L"10^3", L"10^4", L"10^5"])
p