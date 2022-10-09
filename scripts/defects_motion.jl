cd("D:/Documents/Research/projects/kuramoto_ballistic")
    using JLD2,StatsBase,Distributions,LinearAlgebra,Parameters,Random,BenchmarkTools,Hungarian
    include("../methods.jl")
    global const R0 = 1
    using Plots,ColorSchemes,LaTeXStrings
    gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    plot()
&

N = Int(5E3)
    rho = 1
    v0 = 0.5
    σ  = 0.
    L  = round(Int,sqrt(N/rho))
    T  = 0.1 # to be compared to Tc ~ 1 when \sigma = 0
    tmax = 40 ; dt = determine_dt(T,σ,v0,N,rho)

t = 0.0
params_init = ["hightemp",L/2]
pos,thetas,psis,omegas = initialisation(N,L,σ,params_init)
    scatter((pos[1,:],pos[2,:]),marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/L,size=(512,512),xlims=(0,L),ylims=(0,L))

z = @elapsed while t<10
    t += dt
    global pos,thetas = update(pos,thetas,psis,omegas,T,v0,N,L,dt)
end
prinz(z)
# scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/L,size=(512,512),xlims=(0,L),ylims=(0,L))

thetas_cg = cg(pos,thetas,N,L)
    heatmap(mod.(thetas_cg,2pi)',clims=(0,2pi),c=cgrad([:black,:blue,:green,:orange,:red,:black]),size=(400,400))

## Trying to resolve the problem of [-1]
N = Int(1E5)
    rho = 1
    v0 = 0.5
    σ  = 0.
    L  = round(Int,sqrt(N/rho))
    T  = 0.1 # to be compared to Tc ~ 1 when \sigma = 0
    tmax = 40 ; dt = determine_dt(T,σ,v0,N,rho)

# Relax a bit so there is not too many defects
t = 0.
pos,thetas,psis,omegas = initialisation(N,L,σ,["hightemp"])
plot(pos,thetas)

z = @elapsed while t<250
    t += dt
    global pos,thetas = update(pos,thetas,psis,omegas,T,v0,N,L,dt)
end
number_defects(pos,thetas,N,L)
number_defectsP(pos,thetas,N,L)
number_defectsN(pos,thetas,N,L)
plot(pos,thetas,defects=false)
plot(pos,thetas,defects=true)
plot(pos,thetas,defects=false,particles=true)

spot_defects(pos,thetas,N,L)

dft = DefectTracker(pos,thetas,N,L,t)
every = 10; tmax = 500
dft,pos,thetas,t = update_and_track!(dft,pos,thetas,psis,omegas,T,v0,N,L,dt,t,every,tmax)
d=dft.defectsP[1]
plot(rand(10))


msdall,msdp,msdm = MSD(dft,L)
plot()
    plot!(msdall[2:end],axis=:log,m=:circle)
    plot!(msdp[2:end])
    plot!(msdm[2:end])


plot(mean_distance_to_annihilator(dft,L)[1],rib=mean_distance_to_annihilator(dft,L)[2])
