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


filename = "data/recover_XYmodel.jld2"
@load filename Ps Cs ns runtimes Ts Ns v0s rhos sigmas times_log tmax comments R
histogram(runtimes / 3600 / 24, bins=20)

Ps_avg = nanmean(Ps, 8) # N rho T v0 sigma t init R
Ps_std = nanstd(Ps, 8) # N rho T v0 sigma t init R
ns_avg = nanmean(ns, 8) # N rho T v0 sigma t init R
# Cs_avg = nanmean(Cs,9)
indices = [];
for r in 1:R
    try
        Cs[r]
        push!(indices, r)
    catch
    end
end
Cs_avg = mean([Cs[i] for i in indices])
xis_avg = zeros(length(rhos), length(v0s), length(sigmas), length(times_log))
for i in each(v0s), j in each(sigmas), t in each(times_log), p in each(rhos)
    L = round(Int, sqrt(Ns[1] / rhos[p]))
    xis_avg[p,i, j, t] = corr_length(Cs_avg[1, p, 1, i, j, 1, t, 1], 0:L)/L
end

## Correlation function at final time for different sigma
ind_rho = 2
p = plot(xlabel="r", ylabel="C(r,∞)", legend=:outerright, size=(550, 400), yaxis=:log, title="ρ = $(rhos[ind_rho])", ylims=(5E-3, 1.2))
for i in 1:length(sigmas)
    L = round(Int, sqrt(Ns[1] / rhos[ind_rho]))
    plot!(0.5:1:L, remove_negative(Cs_avg[1, ind_rho, 1, 1, i, 1, end-5]), label="σ = $(sigmas[i])", rib=0)
end
# savefig("figures\\transition_SR-LR_v0$(v0s[ind_v0]).png")
p

## Order parameter P at final time for different sigma
ind_rho = 2
p = plot(xlabel="t", ylabel="P(r,∞)", legend=:outerright, size=(550, 400), xaxis=:log, title="ρ = $(rhos[ind_rho])", ylims=(0, 1))
for i in 1#:length(sigmas)
    plot!(times_log, Ps_avg[1, ind_rho, 1, 1, i, 1, :], label="σ = $(sigmas[i])", rib=0)
end
# savefig("figures\\transition_SR-LR_v0$(v0s[ind_v0]).png")
p

## Recover XY model 
pyplot(size=(400, 400), legend=nothing, grid=false, markerstrokewidth=0, linewidth=1.3, thickness_scaling=1.5)
p = plot(xlabel="r", ylabel="C(r,∞)", legend=:topright, size=(450, 400), yaxis=:log, ylims=(5E-3, 1.2))
    for i in 1:length(rhos)
        L = round(Int, sqrt(Ns[1] / rhos[i]))
        plot!(1:1:L, remove_negative(Cs_avg[1, i, 1, 1, 1, 1, end]), label="ρ = $(rhos[i])", rib=0)
    end
    # savefig("figures\\transition_SR-LR_v0$(v0s[ind_v0]).png")
    p

##
p=plot(axis=:log)
    for i in each(rhos)
        plot!(times_log, vec(xis_avg[i,:,:,:]))
    end
    p
    plot!(x->1e-2(x)^1,c=:black)
    plot!(x->1e-2(x)^0.5,c=:black,style=:dash)
## Visualise the systems at final time
gr(size=(400, 400), legend=nothing, grid=false, markerstrokewidth=0, linewidth=1.3, thickness_scaling=1.5)
@load filename thetas_saveds pos_saveds
pp = Vector{Any}(undef, length(rhos))
for i in each(rhos)
    N = Ns[1]
    L = round(Int, sqrt(N / rhos[i]))
    pp[i] = plot(pos_saveds[:, :, i, 1, 1, 1, 1, end], thetas_saveds[:, i, 1, 1, 1, 1, end], N, L, L)
end
p = plot(pp..., layout=(1, length(rhos)), size=(400 * length(rhos), 400))

