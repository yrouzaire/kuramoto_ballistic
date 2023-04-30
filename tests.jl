cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian, Sobol
include("methods.jl")
using Plots, ColorSchemes, LaTeXStrings
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()

## ----------------- Initialisation and Visualisation Tests ----------------- ##
include("parameters.jl")
system = System(param)
energy(system)
	# poss = get_pos(system)
# histogram([el[2] for el in poss], bins=50, label="x")
# scatter(poss, markersize=1, legend=false, xlims=(0,Lx), ylims=(0,Ly), xlabel=L"x", ylabel=L"y", title="N = $(length(poss))")

evolve!(system, system.t + 10) ; energy(system)
plot_thetas(system, particles=true, vertical=true, defects=true)
# ndef = number_defects(system)
# plot_thetas(system, particles=false, vertical=true,defects=true,title="n = $ndef")

dft = DefectTracker(system,0)
interdefect_distance(dft.defectsN[1], dft.defectsP[1], Lx, Ly)
times = 5:5:100
track!(dft, system, times, verbose=true)

mmsd = MSD(dft, Lx, Ly)[1]
plot(vcat(0,times), mmsd)


## -------------------- Parameters -------------------- ##
Ntarget = Int(2E3)
aspect_ratio = 1
rho = 1
T = 0.0
sigma = 0.0
v0 = 0.0
R0 = 1.415

N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)

# Initialisation parameters
init_pos = "random"
init_theta = "hightemp"
r0 = 20.0
q = 1.0
params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio, :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0, :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

## ----------------- PBC Tests for R0 non integer ----------------- ##
system = System(param)
evolve!(system, 1E2)
plot_thetas(system, particles=false, vertical=true)
cg(system)


## ------------------ Efficiency Benchmarks ------------------ ##
ind_neighbours = get_list_neighbours(get_pos(system), N, Lx, Ly, R0)
pos = get_pos(system)
@btime get_list_neighbours($pos, N, Lx, Ly, R0)

update_thetas!(system, ind_neighbours)
@btime update_thetas!(system, ind_neighbours)

update_positions!(system)
@btime update_positions!(system)

## ----------------- Other Tests ----------------- ##
corr(system, algo="fast", dr=0.445)
corr(system)
system.Lx

spot_defects(system)
plot_thetas(system, defects=true)



#= 
In the end, for N = 1E4, rho = 2 : 
runtime = 8.3 (find neighbours) + 4 (update thetas) + 3.8 (update positions) = 16.1 ms 
For N = 1E6, rho = 2 : 
runtime = 2.3 (find neighbours) + 0.77 (update thetas) + 0.5 (update positions) = 3.5s 
=#


