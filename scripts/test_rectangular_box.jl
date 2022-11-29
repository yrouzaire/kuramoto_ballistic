cd("D:/Documents/Research/projects/kuramoto_ballistic")
    using JLD2,StatsBase,Distributions,LinearAlgebra,Parameters,Random,BenchmarkTools,Hungarian
    include("../methods.jl")
    global const R0 = 1
    using Plots,ColorSchemes,LaTeXStrings
    gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    plot()
&

N = Int(1E4)
    rho = 1
    v0 = 1
    σ  = 0.
    aspect_ratio = 1.2# Lx/Ly
    Lx  = round(Int,sqrt(aspect_ratio*N/rho))
    Ly  = round(Int,sqrt(N/rho/aspect_ratio))
    T  = 0.1 # to be compared to Tc ~ 1 when \sigma = 0
    tmax = 10 ; dt = determine_dt(T,σ,v0,N,rho)

## Simple update
t = 0.0
params_init = ["pair",50]
pos,thetas,psis,omegas = initialisation(N,Lx,Ly,σ,params_init)
    scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/Lx,size=(512,512),xlims=(0,Lx),ylims=(0,Ly),aspect_ratio=1/aspect_ratio)
    # plot(corr_fast(pos,thetas,N,L,0.5))
for i in 1:100
    pos,thetas = update(pos,thetas,psis,omegas,T,v0,N,Lx,Ly,dt)
end
scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/Lx,size=(512,512),xlims=(0,Lx),ylims=(0,Ly),aspect_ratio=1/aspect_ratio)
plot(pos,thetas,N,Lx,Ly,particles=false,defects=true)

## Defect Tracker
dft = DefectTracker(pos,thetas,N,Lx,Ly,t)
update_and_track!(dft,pos,thetas,psis,omegas,T,v0,N,Lx,Ly,dt,t,1:5)
dft.defectsN[1].pos
MSD(dft,Lx,Ly)
















&
