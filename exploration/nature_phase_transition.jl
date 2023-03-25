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

## ---------------- Nature of the Phase Transition ---------------- ##
comments = "The goal of this script is to pass through the transition line, 
in both direction (keeping σ or v0 constant) and to compute correlation functions.
Here for T = 0.1 and ρ = 1."
# Physical Params 
Ntarget = Int(1E3)
aspect_ratio = 1
T = 0.1
R0 = 1
rho = 1
rhoc = 4.51 / π

# Initialisation parameters
init_pos = "random"
init_theta = "hightemp"
r0 = 20.0
q = 1.0
params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)

# Simulation parameters
v0sigs = vcat([(0.25,sigm) for sigm in 0:0.025:0.2],[(v,0.1) for v in logspace(0.03,1,10)])
tmax = 1E2
# times = collect(0:tmax/30:tmax) # linear time
times = logspace(1,tmax,30,digits=1) # log time

P = zeros(length(v0sigs), length(times))
C = Array{Vector{Float64}}(undef, length(v0sigs), length(times))
xi = zeros(length(v0sigs), length(times))
n = zeros(length(v0sigs), length(times))

z = @elapsed for i in each(v0sigs)
    v0, sigma = v0sigs[i]
    println("v0 = $v0, σ = $sigma, $(100i/length(v0sigs))%")
    N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
    dt = determine_dt(T, sigma, v0, N, rho)

    param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
        :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
        :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

    system = System(param)

    t = 0.0
    token = 1

    for tt in eachindex(times)
        evolve(system, times[tt]) # evolves the systems up to times[tt]
        
        P[i,tt]  = polarOP(system)[1]
        corr_tmp    = corr(system)
        C[i,tt]  = corr_tmp
        xi[i,tt] = corr_length(corr_tmp)
        n[i,tt]  = number_defects(system)
    end

end
prinz(z)