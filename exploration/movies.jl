cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian, Sobol
include("../methods.jl")
const global R0 = 1
using Plots, ColorSchemes, LaTeXStrings
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()

function movies(params, every, tmax, dt; particles=false)
    rho, T, v, Ntarget, Lx,Ly, σ, params_init = params
    N,Lx,Ly = effective_number_particle(Ntarget, rho)
    pos, thetas, omegas, psis = initialisation(N, Lx, Ly, σ, params_init)
    # println("N = $N")
    times = 1:every:round(Int, tmax)
    Q = zeros(length(times))
    P = zeros(length(times))
    anim = @animate for i in 1:length(times)
        println("$(round(i/length(times)*100,digits=2)) %")
        ind_neighbours0 = get_list_neighbours(pos, N, Lx, Ly)
        for j in 1:round(Int, every / dt)
            if v == 0
                ind_neighbours = ind_neighbours0
            else
                ind_neighbours = get_list_neighbours(pos, N, Lx, Ly)
            end
            pos, thetas = update(pos, thetas, omegas, psis, ind_neighbours, T, v, N, Lx, Ly, dt)
        end
        titre = "P=$(round(polarOP(thetas)[1],digits=2))"
        if particles
            p1 = scatter(pos[1, :], pos[2, :], marker_z=mod.(thetas, 2pi), color=cols, clims=(0, 2pi), ms=275 / Lx, size=(512, 512), xlims=(0, Lx), ylims=(0, Ly), title=titre, aspect_ratio=1)
            p2 = heatmap(mod.(cg(pos, thetas, N, Lx, Ly), 2pi)', clims=(0, 2pi), c=cols, aspect_ratio=1, size=(512, 512))
            plot(p1, p2, size=(1024, 512))
        else
            p2 = heatmap(mod.(cg(pos, thetas, N, Lx, Ly), 2pi)', clims=(0, 2pi), c=cols, aspect_ratio=1, size=(512, 512))
        end
    end
    return anim
end

##
N = Int(1E3)
rho = 1
T = 0.1 # temperature for angle diffusion
v = 0.1 # norm of individual velocities
σ = 0

# Other parameters
L = round(Int, sqrt(N / rho))
params_init
params = Any[rho, T, v, N, L, σ, params_init] # any to avoid N being interpreted as a Float
dt = determine_dt(T, σ, v, N, rho)

every = 2;
tmax = 1000;
z = @elapsed anim = movies(params, every, tmax, dt, particles=false)
prinz(z)
mp4(anim, "films/recoverXY_N$(N)_rho$(rho)_T$(T)_v0$(v)_σ$(σ)_tmax$(tmax)_every$(every).mp4", fps=25)

## pre running for efficiency
pos, thetas, psis, omegas = initialisation(N, L, σ)
pos, thetas = update(pos, thetas, psis, omegas, T, v, N, L, dt)


## Impact of initialisation type on the XY model (Sobol, random, regular square lattice)
Ntarget = Int(1E4)
rho = 2
aspect_ratio = 1
N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
T = 0.1 # temperature for angle diffusion
v = 0 ; σ = 0
dt = determine_dt(T, σ, v, N, rho)
every = 2;
tmax = 1000;

params_init_random = ["hightemp","random"]
params_init_sobol  = ["hightemp","sobol"]
params_init_square = ["hightemp","square"]

z = @elapsed for params_init in [params_init_square]
    params = Any[rho, T, v, N, Lx, Ly, σ, params_init] # any to avoid N being interpreted as a Float
    anim = movies(params, every, tmax, dt, particles=false)
    mp4(anim, "films/recoverXY_v00_sigma0/$(params_init[2])__N$(N)_rho$(rho)_T$(T)_v0$(v)_σ$(σ)_tmax$(tmax)_every$(every).mp4", fps=25)
end
prinz(z)