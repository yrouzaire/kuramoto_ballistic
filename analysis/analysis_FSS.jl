cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl");
const global R0 = 1
using Plots, ColorSchemes, LaTeXStrings
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

## Load FSS data
filename = "data/FSS.jld2"
@load filename Ntargets v0sigs Ps Cs ns xis params_init aspect_ratio times tmax T comments rho rhoc runtimes R
histogram(runtimes / 3600 / 24, bins=20)
histogram(runtimes / 3600/24 *100, bins=20)

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

times = collect(0:10:100)
Cs_avg = Array{Vector}(undef, length(v0sigs), length(Ntargets), length(times))
for i in 1:length(v0sigs), j in 1:length(Ntargets), k in 1:length(times)
	Cs_avg[i,j,k] = mean([Cs[i,j,k,r] for r in indices])
end

## Plot P vs t
p=plot(axis=:log, legend=:bottomright, xlabel="t", ylabel="P(t)")
for i in 1:length(v0sigs)
	plot!(times[2:end], Ps_avg[i,4,2:end], label="v0 = $(v0sigs[i][1]) , σ = $(v0sigs[i][2])")
end
plot!(t->4E-2sqrt(t/log(10t)), c=:black )
# plot!(t->1E-2sqrt(t), c=:black, line=:dash)
p

## ------------- Check if Steady State reached ------------- 
p=plot(axis=:log, legend=false, xlabel="r", ylabel="C(r)", ylims=(1E-2,1.2))
for nn in 3#:length(Ntargets)
	for i in 2#:length(v0sigs)
		cor = remove_negative(Cs_avg[i,nn,end])
		plot!(1:length(cor), cor, label="v0 = $(v0sigs[i][1]) , σ = $(v0sigs[i][2])")
	end
end
p

## ------------- Impact of N on final quantities ------------- 
## P 
p=plot(axis=:log, legend=:bottomright, xlabel=L"1/N", ylabel=L"P(t_\infty)")
plot!([1/Ntargets[i] for i in 1:length(Ntargets)], Ps_avg[1,:,end],m=true, label="v0 = $(v0sigs[1][1]) , σ = $(v0sigs[1][2])")

## ------------- Regular plots at fixed N ------------- 
## Plot n vs t
p=plot(axis=:log, legend=false, xlabel="t", ylabel="n(t)/L²")
for i in 1:length(v0sigs)
	plot!(times[2:end], ns_avg[i,end,2:end]/(Ntargets[ind_N]), label="v0 = $(v0sigs[i][1]) , σ = $(v0sigs[i][2])")
end
plot!(t->1.5E3(log(10t)/t)/(Ntargets[ind_N]), c=:black )
plot!(t->7E3(1/t)/(Ntargets[ind_N]), c=:black, line=:dash )
p

## Plot ξ vs t
ind_N = 4
p=plot(axis=:log, legend=false, xlabel="t", ylabel="ξ(t)", title="N = $(Ntargets[ind_N])")
for i in 1:length(v0sigs)
	plot!(times[2:end], xis_avg[i,ind_N,2:end]/sqrt(Ntargets[ind_N]), label="v0 = $(v0sigs[i][1]) , σ = $(v0sigs[i][2])")
end
plot!(t->4sqrt(t/log(10t))/sqrt(Ntargets[ind_N]), c=:black )
# plot!(t->1E-2sqrt(t), c=:black, line=:dash)
p


