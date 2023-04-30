cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian, Sobol
include("../methods.jl")
using Plots, ColorSchemes, LaTeXStrings
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()

## --------------- Small Wrapper --------------- ##
function movies(param, times; particles=false, verbose=true, nb_neighbours=false, defects=false)   
    system = System(param)
    anim = @animate for tt in each(times)
        if verbose 
            println(round(times[tt]*100/times[end],digits=3),"%")
        end
        evolve!(system, times[tt])
        plot_thetas(system,size=(512,512),
            particles = particles, 
            nb_neighbours = nb_neighbours, 
            defects = defects)
        title!("t = $(round(times[tt],digits=2))")
    end
    return anim
end

## --------------- Warming up --------------- ##
include("../parameters.jl")
system = System(param)
get_pos(system)
evolve!(system, 10)
evolve!(system, 100)
plot_thetas(system,particles=true,nb_neighbours=false,defects=false,vertical=true,size=(512,512))
number_defects(system)
plot_thetas(system,particles=true,nb_neighbours=true,defects=true,vertical=true,size=(512,512))
## --------------- Movies --------------- ##
include("../parameters.jl")
inits_pos = ["square","tri","random","RSA","PDS"]
inits_pos = ["square","tri"]
inits_pos = ["random"]


inits_theta = ["hightemp","pair"]
inits_theta = ["pair"]
tmax = Int(3E2)
times = collect(0:2:tmax) # linear time
# times = logspace(1,tmax,10) # log time
anims = Array{Any}(undef,length(inits_pos),length(inits_theta))
n = 0
z = @elapsed for i in each(inits_pos), j in each(inits_theta) 
    init_pos = inits_pos[i]
    init_theta = inits_theta[j]

    n += 1 ; println("n = $n/$(length(inits_pos)*length(inits_theta))")
    params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
    param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
        :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
        :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

    try 
        anim = movies(param, times, 
                particles=true,
                defects=true,
                nb_neighbours=true,
                verbose=true)

        anims[i,j] = anim
        filepath = "impact_init/films/$(init_theta)/E_$(init_pos)_N$(N)_ρ$(rho)_R0$(round(R0,digits=2))_T$(T).mp4"
        mp4(anim, filepath, fps=30)
    catch e 
        println(e)
    end
end
prinz(z)

