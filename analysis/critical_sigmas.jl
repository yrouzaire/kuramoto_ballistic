cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl");
const global R0 = 1
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

filename = "data/critical_sigma.jld2"
@load filename runtimes critical_sigmas_fusion times sigmas v0s rhos tmax R
histogram(runtimes / 3600 / 24, bins=20)

critical_sigmas_avg = nanmean(critical_sigmas_fusion, 3)[:, :, 1]
## Plotting
p = plot(uaxis=:log, legend=:topleft)
for i in reverse(each(rhos))
	if rhos[i] ≈ 4.51/pi lab = L"ρ = ρ_c ≈ 1.44"
	else lab = "ρ = $(rhos[i])"
	end
    plot!((v0s), critical_sigmas_avg[:, i], label=lab, rib=0, m=true)
end
p

##
