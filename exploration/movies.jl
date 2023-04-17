cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian, Sobol
include("../methods.jl")
using Plots, ColorSchemes, LaTeXStrings
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);

## --------------- Small Wrapper --------------- ##

function movies(param, times; particles=false)   
    system = System(param)
    anim = @animate for tt in each(times)
        println(round(times[tt]*100/times[end],digits=3),"%")
        evolve!(system, times[tt])
        plot_thetas(system,particles=true,size=(512,512))
        title!("t = $(round(times[tt],digits=2))")
    end
    return anim
end

## --------------- Movies --------------- ##
include("../parameters.jl")
inits_pos = ["square_lattice","random","RSA","PDS"]
inits_pos = ["random","RSA","PDS"]
# inits_pos = ["square_lattice"]
tmax = 1E2
times = collect(0:2:tmax) # linear time
# times = logspace(1,tmax,10) # log time
anims = []
z = @elapsed for init_pos in inits_pos
    params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
    param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
        :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
        :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

    anim = movies(param, times, particles=true)
    push!(anims, anim)
    filepath = "impact_init/films/$(init_theta)_$(init_pos)_N$(N)_ρ$(rho)_R0$(round(R0,digits=2))_T$(T).mp4"
    mp4(anim, filepath, fps=30)
end
prinz(z)

