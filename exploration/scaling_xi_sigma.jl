cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl")
const global R0 = 1
using Plots, ColorSchemes, LaTeXStrings
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

## This file aims at verifying the scaling of the characteristic length as a function of sigma 
# filename = "data/XIN_during_coarsening_v0_sigma_N1E4_tmax2E3.jld2" # only 3 sigmas
# @load filename tmax times_log rhos sigmas R Ts Cs v0s Ns ns Ps runtimes comments pos_saveds thetas_saveds psis_saveds omegas_saveds

## Load data
filename = "data/phase_space_rho_sig_v0_N1E3_tmax2500.jld2" # here 9 sigmas
@load filename tmax times_log rhos sigmas R Ts Cs v0s Ns ns Ps runtimes comments pos_saveds thetas_saveds psis_saveds omegas_saveds

Psavg = mean(Ps, dims=8)
nsavg = mean(ns, dims=8)
#= Cs[1] # N,rho,T,v0,sigma,init?,t  
Ns
Ts
rhos
sigmas
times_log
v0s =#

xis = zeros(length(rhos), length(sigmas), length(times_log), R)
for r in 1:R , j in each(rhos) , k in each(sigmas)
	rho = rhos[j]
	sigma = sigmas[k]
	for t in each(times_log)
		try 
			tmp = Cs[r][1, j, 1, 1, k, 1, t]
			rs = 1:length(tmp)
			xis[j, k, t, r] = corr_length(tmp, rs)
		catch 
			xis[j, k, t, r] = NaN
		end
	end
end
xis_avg = nanmean(xis, 4)

## Plotting ξ over time 
p=plot(axis=:log,legend=false)
	for j in each(rhos)
		rho = rhos[j]
		for k in each(sigmas)
			sigma = sigmas[k]
			plot!(times_log, xis_avg[j,k,:], label="ρ = $rho, σ = $sigma")
		end
	end
	plot!(x->sqrt(x/log(x)),c=:black)
	p

## Plotting ξ against sigma
p=plot(axis=:log,legend=false,xlabel=L"\sigma",ylabel=L"\xi")
	for j in each(rhos)
		rho = rhos[j]
		L = sqrt(Ns[1]/rho)	
		plot!(sigmas[2:end], xis_avg[j,2:end,end]/L, label="ρ = $rho")
	end
	# plot!(x->1/x,c=:black)
	p

