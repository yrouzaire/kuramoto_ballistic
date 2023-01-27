cd("D:/Documents/Research/projects/kuramoto_ballistic")
    using JLD2,StatsBase,Distributions,LinearAlgebra,Parameters,Random,BenchmarkTools
    include("../methods.jl")
    global const R0 = 1
    using Plots,ColorSchemes,LaTeXStrings
    gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    plot()
&

## Reprendre ses marques
N = Int(1E4)
    rho = 1
    v0 = 0.05
    σ  = 0.
    L  = round(Int,sqrt(N/rho))
    T  = 0.1 # to be compared to Tc ~ 1 when \sigma = 0
    tmax = 10 ; dt = determine_dt(T,σ,v0,N,rho)

t = 0.0
params_init = ["hightemp",L/2]
pos,thetas,omegas,psis = initialisation(N,L,σ,params_init)
    # scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/L,size=(512,512),xlims=(0,L),ylims=(0,L))
    # plot(corr_fast(pos,thetas,N,L,0.5))


z = @elapsed while t<300
    t += dt
    global pos,thetas = update(pos,thetas,omegas,psis,T,v0,N,L,dt)
end
    plot(pos,thetas,particles=true)
prinz(z)

spot_defects(pos,thetas,N,L)


## Timing
N = Int(1E4)
    rho = 1
    v = 1
    σ = 0
    L = round(Int,sqrt(N/rho))
    T = 0.1 # to be compared to Tc ~ 1 when \sigma = 0

    tmax = 10 ; dt = determine_dt(T,σ,v,N,rho)

t = 0.0
params_init = ["hightemp",L/2]
pos,vel_angles,thetas,omegas = initialisation(N,L,σ,params_init)

update(pos,vel_angles,thetas,omegas,T,v,N,L,rho,σ,dt)
@btime update(pos,vel_angles,thetas,omegas,T,v,N,L,rho,σ,dt)
@time
#= Runtimes for update, with N = 1E4, rho = 1
Float64 / 32 / 16: 21 / 20.5 / 23.4 ms
sum(sin,...) instead of sum(sin(...)) gain = 0.1 ms
update the same structure instead gain = 0.2 ms
Declaring the type of indices_neighbours made me gain 10 ms out of 20 !!
Having const variables omegas and vel_angles does not seem to change the runtime
=#
