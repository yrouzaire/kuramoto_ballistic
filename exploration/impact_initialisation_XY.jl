cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian, LambertW
include("../methods.jl")
const global R0 = 1
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

## ---------------- Impact of init on XY Model ---------------- ##
comments = "Investigates the impact of initialisation for the spatial location of the 
spins in the XY model. Even though no symmetry is broken, it seems that having the spins
on (A) a regular lattice (B) a random lattice (C) a 2D Sobol sequence (D) a RSA lattice changes the behaviour of the defects.
The idea is to investigate the role of the fluctuations in the initialisation of the spins on the defects dynamics."
# Physical Params 
Ntarget = Int(1E4)
aspect_ratio = 1
T = 0.1
R0 = 1
rho = 1
rhoc = 4.51 / π

# Initialisation parameters
init_poss = ["random", "regular", "sobol", "rsa"]
init_theta = "pair"
r0 = 20.0
q = 1.0
params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => NaN, :q => q)

# Simulation parameters
v0sigs = [(0.5,0),(0.5,0.1)]
r0s = 5:10:35
tmax = 1E2
times = 0:5:tmax # linear time

z = @elapsed for i in each(v0sigs), j in each(r0s)
    v0, sigma = v0sigs[i]
    r0 = r0s[j]

    println("v0 = $v0, σ = $sigma, r0 = $r0 ,  $(100i/length(v0sigs))%")
    N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
    dt = determine_dt(T, sigma, v0, N, rho)

    params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
    
    param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
        :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
        :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

    t = 0.0
    system = System(param)
    dft = DefectTracker(system, t)
    dft, pos, thetas, t = track!(dft,system,times)
end
prinz(z)