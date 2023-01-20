cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../../methods.jl");
const global R0 = 1
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()

## This files plots the phase space (v0,sigma) of the Kuramoto model for different values of rho

filename = "data/phase_space_rho_sig_v0_N1E3_tmax2500.jld2"
@load filename Ps Cs ns runtimes Ts Ns v0s rhos sigmas times_log tmax comments R

Ps_avg = nanmean(Ps, 8) # N rho T v0 sigma init t R
Ps_std = nanstd(Ps, 8) # N rho T v0 sigma init t R
ns_avg = nanmean(ns, 8) # N rho T v0 sigma init t R

p1 = heatmap(v0s[2:end], sigmas, Ps_avg[1, 1, 1, 2:end, :, 1, end, 1]',
    xaxis=:log, c=cgrad([:red, :orange, :green]), clims=(0, 1),
    size=(520, 400), xlabel=L"v_0", ylabel="σ", title="ρ = $(rhos[1])",
    colorbartitle="P", colorbar=:right, colorbar_titlefont=font(12))

p2 = heatmap(v0s[2:end], sigmas, Ps_avg[1, end, 1, 2:end, :, 1, end, 1]',
    xaxis=:log, c=cgrad([:red, :orange, :green]), clims=(0, 1),
    size=(520, 400), xlabel=L"v_0", ylabel="σ", title="ρ = $(rhos[end])",
    colorbartitle="P", colorbar=:right, colorbar_titlefont=font(12))


## Find critical velocity at final time for which the order parameter P is == seuil
seuils = [0.4, 0.45, 0.5, 0.55, 0.6]
# seuils = [0.35,0.4,0.45,0.5,0.55,0.6,0.65] # not much better

sigma_crit = zeros(length(rhos), length(v0s), length(seuils))
for i in each(rhos), j in each(v0s), k in each(seuils)
    # Pavg # [N rho T v0 sigma init t R]
    ind = findfirst(Ps_avg[1, i, 1, j, 1:end, end, end, 1] .< seuils[k])
    if ind == nothing
        sigma_crit[i, j, k] = NaN
    else
        sigma_crit[i, j, k] = sigmas[ind]
    end
end

p3 = plot(legend=false, xaxis=:log)
for i in each(rhos)
    plot!(v0s[2:end], nanmean(sigma_crit, 3)[i, 2:end, 1], rib=0)
end
p3

## Final figure
p = plot(p1, p2, p3, layout=(1,3), size=(1400, 400), legend=:topleft)
