cd("D:/Documents/Research/projects/kuramoto_ballistic")
    using JLD2,StatsBase,Distributions,LinearAlgebra,Parameters,Random,BenchmarkTools,Hungarian
    include("../methods.jl")
    global const R0 = 1
    using Plots,ColorSchemes,LaTeXStrings
    gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    plot()

function movies(params,every,tmax,dt;coarsegrain=false)
    rho,T,v,N,L,σ = params

    pos,thetas,psis,omegas = initialisation(N,L,σ)
    println("N = $N")
    times = 1:every:round(Int,tmax)
    Q = zeros(length(times)) ; P = zeros(length(times))
    anim = @animate for i in 1:length(times)
        println("$(round(i/length(times)*100,digits=2)) %")
        for j in 1:round(Int,every/dt) pos,thetas = update(pos,thetas,psis,omegas,T,v,N,L,dt) end
        titre = "P=$(round(polarOP(thetas)[1],digits=2))"
        if coarsegrain
            p1 = scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/L,size=(512,512),xlims=(0,L),ylims=(0,L),title=titre,aspect_ratio=1)
            p2 = heatmap(mod.(cg(pos,thetas,N,L),2pi)',clims=(0,2pi),c=cols,aspect_ratio=1,size=(512,512))
            plot(p1,p2,size=(1024,512))
        else
            p1 = scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/L,size=(512,512),xlims=(0,L),ylims=(0,L),title=titre,aspect_ratio=1)
        end
    end
    return anim
end

N = Int(1E4)
    rho = 1
    T = 0.1 # temperature for angle diffusion
    v = 1 # norm of individual velocities
    σ = 0.3

    # Other parameters
    L = round(Int,sqrt(N/rho))
    params = Any[rho,T,v,N,L,σ] # any to avoid N being interpreted as a Float
    dt = determine_dt(T,σ,v,N,rho)

every = 1 ; tmax = 1000
z = @elapsed anim = movies(params,every,tmax,dt,coarsegrain=true)
prinz(z)
mp4(anim,"films/KB_N$(N)_rho$(rho)_T$(T)_v0$(v)_σ$(σ)_tmax$(tmax)_every$(every).mp4",fps=25)

## pre running for efficiency
pos,thetas,psis,omegas = initialisation(N,L,σ)
pos,thetas = update(pos,thetas,psis,omegas,T,v,N,L,dt)
