include("methods.jl")
using Plots, ColorSchemes, BenchmarkTools
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);

######################## Parameters ########################
Ntarget = Int(1E4)
aspect_ratio = 1
rho = 2
T = 0.0
sigma = 1.0
v0 = 1.0
R0 = 1

N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)

# Initialisation parameters
init_pos = "random"
init_theta = "pair"
r0 = 20.0
q = 1.0
params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio, :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0, :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

######################## Initialisation and Visualisation Tests ########################
system = System(param)
poss = get_pos(system)
thetas = get_theta(system)
cg(poss, thetas, Lx, Ly)

plot_thetas(system)
plot_thetas(system, particles=true)
plot_thetas(system, particles=true, vertical=true)

######################## Efficiency Benchmarks ########################
ind_neighbours = get_list_neighbours(get_pos(system), N, Lx, Ly, R0)
pos = get_pos(system)
@btime get_list_neighbours($pos, N, Lx, Ly, R0)

update_thetas!(system, ind_neighbours)
@btime update_thetas!(system, ind_neighbours)

update_positions!(system)
@btime update_positions!(system)

corr()
#= 
In the end, for N = 1E4, rho = 2 : 
runtime = 8.3 (find neighbours) + 4 (update thetas) + 3.8 (update positions) = 16.1 ms 
For N = 1E6, rho = 2 : 
runtime = 2.3 (find neighbours) + 0.77 (update thetas) + 0.5 (update positions) = 3.5s 
=#