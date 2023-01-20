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

filename = "data/investigation_v00.jld2"
filename = "data/"
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
xis_avg = zeros(length(v0s), length(sigmas), length(times_log))
for i in each(v0s), j in each(sigmas), t in each(times_log)
    xis_avg[i, j, t] = corr_length(Cs_avg[1, 1, 1, i, j, 1, t, 1], 0:0.5:50)
end

## Phase space
v0s
heatmap(v0s[2:end], sigmas, Ps_avg[1, 1, 1, 2:end, :, 1, end, 1]', xaxis=:log, c=cgrad([:red, :orange, :green]), clims=(0, 1))
rhos

## Correlation function at final time for different v0
ind_sig = 1
p = plot(xlabel="r", ylabel="C(r,∞)", legend=:outerright, size=(550, 400), yaxis=:log, title="σ = $(sigmas[ind_sig])", ylims=(1E-2, 1.2))
for i in 1:length(v0s)
    plot!(0.5:0.5:50, remove_negative(Cs_avg[2:end, 1, 1, 1, i, ind_sig, 1, end, 1]), label="v = $(v0s[i])", rib=0)
end
savefig("figures\\transition_SR-LR_sigma$(sigmas[ind_sig]).png")
p


## Correlation length over time for different v0
ind_sig = 1
p = plot(xlabel="t", ylabel="ξ(t)", legend=:outerright, size=(550, 400), axis=:log, title="σ = $(sigmas[ind_sig])")
for i in 1:length(v0s)
    plot!(times_log, xis_avg[i, ind_sig, :], label="v = $(v0s[i])", rib=0)
end
plot!(x -> sqrt(x), c=:black)
savefig("figures\\xit_SR-LR_sigma$(sigmas[ind_sig]).png")
p

## Correlation functions over time
ind_v0 = 8
ind_sigma = 2
p = plot(xlabel="r", ylabel="C(r,t)", legend=false, axis=:log, ylims=(1E-3, 1.2))
for tt in each(times_log)
    plot!(0.5:0.5:50, remove_negative(Cs_avg[2:end, 1, 1, 1, ind_v0, ind_sigma, 1, tt, 1]))
end
p

## Correlation functions over phase space
styles = [:solid, :dash]
p = plot(xlabel="r", ylabel="C(r,t)", legend=false, axis=:log, ylims=(1E-3, 1.2))
for ind_v0 in each(v0s)
    for ind_sigma in each(sigmas)
        plot!(0.5:0.5:50, remove_negative(Cs_avg[2:end, 1, 1, 1, ind_v0, ind_sigma, 1, end, 1]), c=ind_v0, line=styles[ind_sigma])
    end
end
p

## Correlation length over time
p = plot(xlabel="t", ylabel="ξ(t)", legend=:outerright, axis=:log, size=(600, 400))#,ylims=(1E-3,1.2))
for ind_v0 in 1:length(v0s)
    plot!([NaN, NaN], rib=0, label="v = $(v0s[ind_v0])")
    for ind_sigma in each(sigmas)
        plot!(times_log, xis_avg[ind_v0, ind_sigma, :], c=ind_v0, line=styles[ind_sigma])
    end
end
plot!(x -> sqrt(x / log(x)), c=:black)
plot!(x -> sqrt(x), c=:grey)
hline!([exp(-1) * 100], c=:grey, line=:dot)

## XY Collapse Scaling ?
rr = collect(0.5:0.5:50)
ind_v0 = 1
ind_sigma = 1
p = plot(xlabel="r/ξ(t)", ylabel=L"r^{-η}\,C(r/ξ,t)", legend=false, axis=:log, ylims=(1E-3, 1.2))
for tt in 15:2:length(times_log)
    # plot!(rr ./xis_avg[ind_v0,ind_sigma,tt],rr .^(-0.1/2pi) .* remove_negative(Cs_avg[2:end,1,1,1,ind_v0,ind_sigma,1,tt,1]))
    plot!(rr, remove_negative(Cs_avg[2:end, 1, 1, 1, ind_v0, ind_sigma, 1, tt, 1]))
end
p

## XY Relation Length and number of vortices ?
rr = collect(0.5:0.5:50)
p = plot(xlabel="t", ylabel=L"ξ.\sqrt{n}/L", legend=false, xaxis=:log)
for ind_v0 in 1:length(v0s)
    # plot!([NaN,NaN],rib=0,label="v = $(v0s[ind_v0])")
    for ind_sigma in 1#each(sigmas)
        # plot!(times_log,.(v0s[ind_v0]*times_log).^(-0.01) .* xis_avg[ind_v0,ind_sigma,:] .* sqrt.(ns_avg[1,1,1,ind_v0,ind_sigma,1,:,1])/sqrt(Ns[1]))
        plot!(times_log, xis_avg[ind_v0, ind_sigma, :] .* sqrt.(ns_avg[1, 1, 1, ind_v0, ind_sigma, 1, :, 1]) / sqrt(Ns[1]))
    end
end
p
