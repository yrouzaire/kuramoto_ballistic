cd("/Users/yrouzaire/Documents/Recherche/GitHub/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl");
using Plots, ColorSchemes, LaTeXStrings
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

## ---------------- Data Generation ---------------- ##
## ---------------- Data Generation ---------------- ##


include("../parameters.jl")
param[:params_init][:init_theta] = "pair"
param[:params_init][:r0] = 32
param[:v0] = 0.6
param[:sigma] = 0.05
param[:Ntarget] = Int(1E4)

R = 5
saving_times = [0, 100, 300, 500]
saving_systems = Array{System}(undef, length(saving_times), R)

z = @elapsed Threads.@threads for r in 1:R
    system = System(param)
    for tt in each(saving_times)
        evolve!(system, saving_times[tt])
        saving_systems[tt,r] = deepcopy(system)
    end
end
prinz(z)

filename = "data/snapshots/snapshot_defects_moving_apart_v$(v0)_σ$(sigma)_ρ$(rho)_time$(saving_times[end]).jld2"
@save filename saving_systems saving_times param R

## ---------------- Plots ---------------- ##
## ---------------- Plots ---------------- ##


# filename = "data/snapshots/snapshot_defects_moving_apart_v0.55_σ0.3_ρ1_time50.jld2"
# @load filename saving_systems saving_times param R

##
for ind_r in 1:R
    spot_def = true
    pp = Array{Any}(undef, length(saving_times))
    for i in each(saving_times)
    pp[i] = plot_thetas(saving_systems[i,ind_r], particles=false, vertical=true, defects=spot_def)
    end    
    p = plot(pp... , layout=(2,2), size=(800, 800))
    display(p)
end


## ---------------- Plots manual pair of defect Appendix ---------------- ##
## ---------------- Plots manual pair of defect Appendix ---------------- ##
## ---------------- Plots manual pair of defect Appendix ---------------- ##
## ---------------- Plots manual pair of defect Appendix ---------------- ##

include("../parameters.jl")
param[:params_init][:init_theta] = "pair"
param[:params_init][:r0] = 32
param[:v0] = 0.6
param[:sigma] = 0.05
param[:Ntarget] = Int(1E4)

R = 5
saving_times = [0, 100, 300, 500]
saving_systems = Array{System}(undef, length(saving_times), R)

system = System(param)
plot_thetas(system, particles=true, vertical=false, defects=true)
savefig("figures_paper/manually_created_pair.png")
