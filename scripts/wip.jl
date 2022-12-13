cd("D:/Documents/Research/projects/kuramoto_ballistic")
    using JLD2,StatsBase,Distributions,LinearAlgebra,Parameters,Random
    using Plots,ColorSchemes,LaTeXStrings
    include("../methods.jl")
    pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    plot()
&
cette ligne est mauvaise
## Test specific function for evolution at v0 = 0
N = Int(1E4)
    #rho = 1.3*(4.51/pi)
    rho = 1
    v0 = 0.1
    σ  = 0.3
    aspect_ratio = 1 # Lx/Ly
    Lx  = round(Int,sqrt(N/rho*aspect_ratio))
    Ly  = round(Int,sqrt(N/rho/aspect_ratio))
    T  = 0.1 # to be compared to Tc ~ 1 when \sigma = 0
    tmax = 10 ; dt = determine_dt(T,σ,v0,N,rho)
    global R0 = 1

using BenchmarkTools
t = 0.0
params_init = ["hightemp",L/2]
pos,vel_angles,thetas,omegas = initialisation(N,L,σ,params_init)
    scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/L,size=(512,512),xlims=(0,L),ylims=(0,L))


## Development of vortex localisation (and type identification ?)
N = Int(1E4)
    rho = 1
    v = 5
    σ = 0.5
    L = round(Int,sqrt(N/rho)) ; R0 = 1
    T = 0.1 # to be compared to Tc ~ 1 when \sigma = 0

    tmax = 10 ; dt = determine_dt(T,σ,v,N,rho)

t = 0.0
params_init = ["hightemp",L/2]
gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
pos,vel_angles,thetas,omegas = initialisation(N,L,σ,params_init)
    scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/L,size=(512,512),xlims=(0,L),ylims=(0,L))

z = @elapsed while t<10
    t += dt
    global pos,thetas = update(pos,vel_angles,thetas,omegas,T,v,N,L,rho,R0,σ,dt)
end
prinz(z)
scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/L,size=(512,512),xlims=(0,L),ylims=(0,L))

density, thetas_cg = cg(pos,thetas,N,L,R0)
heatmap(mod.(thetas_cg,2pi)',clims=(0,2pi),c=cgrad([:black,:blue,:green,:orange,:red,:black]),size=(400,400))


## Analysis behaviour with sigma at T = 0
@unpack Ts,Ns,v0s,rhos,tmax,sigmas,times,P,C,runtimes = load("data/OP_T0_scan_sigma.jld")
histogram(runtimes/3600,bin=20)
Pavg = nanmean(P,7)

# Final time
p=plot(legend=:right)
    for i in each(v0s)
        plot!(sigmas,Pavg[1,1,1,i,:,end,1],rib=0,label=L"v_0 = "*string(v0s[i]))
    end
    p

# Over time for v0 = 1
p=plot()
    for i in [1,2,3,5,10,20,30,41]
        plot!(sigmas,Pavg[1,1,1,2,:,i,1],label="t = $(times[i])",rib=0)
    end
    p
