cd("D:/Documents/Research/projects/sync_ballistic")
    using JLD2,StatsBase,Distributions,LinearAlgebra,Parameters,Random
    using Plots,ColorSchemes,LaTeXStrings
    include("methods.jl")
    pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    plot()
&


## Benchmark corr fast
N = 1000
    α = 1
    T  = 0.
    v0 = 1
    σ  = 0.
    R0 = 1
    R  = 10
    tmax = 10 ; times = 0:tmax-0.1:tmax

    dr = R0/2
    L = round(Int,sqrt(pi*N/α)); rho = N/L^2
    param = Any[α,T,v0,N,L,rho,R0,σ]

    pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
    t = 0.0 ; dt = 0.01
while t < 50
        t += dt
        pos,thetas = update(pos,vel_angles,thetas,omegas,param,dt)
    end

C = corr(pos,thetas,param,dr)
    plot(C)
Cfast = corr_fast(pos,thetas,param,dr)
    plot!(Cfast)

plott(pos,thetas)

function plott(pos,thetas,alpha=1)
    N = length(thetas)
    L = sqrt(pi*N/alpha)
    p = scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=500/L,size=(512,512),xlims=(0,L),ylims=(0,L))
    return p
end


## Efficiency of function corr
Ns = Int.(2.0 .^ collect(5:15))
    α = 1
    T  = 0.
    v0 = 1
    σ  = 0.
    R0 = 1
    R  = 10
    tmax = 10 ; times = 0:tmax-0.1:tmax

    dr = R0/2

    runtimes = zeros(length(Ns))
    runtimesfast = zeros(length(Ns))

for n in eachindex(Ns)
    N = Ns[n] ; L = round(Int,sqrt(pi*Ns[n]/α)) ; rho = N/L^2
    println("N = $N")
    params = Any[α,T,v0,Ns[n],L,rho,R0,σ]
    pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
    # runtimes[n] = @elapsed corr(pos,thetas,params,dr)
    runtimesfast[n] = @elapsed corr_fast(pos,thetas,params,dr)
end


# plot(Ns,runtimes,m=true,axis=:log)
plot(Ns,runtimesfast,m=true,axis=:log)
    plot!(Ns,Ns .^2 * 1E-6,m=true,axis=:log)

## Strong Finite Size Effects ? Close to the transitions yes
Ns = Int.(2.0 .^ collect(5:15))
    α = 1
    T  = 0.
    v0 = 1
    σ  = 0.
    R0 = 1
    R  = 1
    tmax = 10 ; times = 0:tmax-0.1:tmax

    P = zeros(length(Ns),R)
    dr = R0/2 ; Crt = Array{Vector{Float64},2}(undef,length(Ns),R)

z = @elapsed for n in eachindex(Ns)
    N = Ns[n]
    R0 = 1
    L = round(Int,sqrt.(π*N/α)*R0)
    rho = N/L^2
    # v0 = v0*L

    params = Any[α,T,v0,N,L,rho,R0,σ]
    dt = determine_dt(T,σ,v0,R0,L)

    println("N = $N started, with L = $L , dt = $dt")
    # for r in 1:R
    Threads.@threads for r in 1:R
        pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
        t = 0.0
        while t < tmax
            t += dt
            pos,thetas = update(pos,vel_angles,thetas,omegas,params,dt)
        end

        P[n,r] = polarOP(thetas)[1]
        Crt[n,r] = corr(pos,thetas,params,dr)
        end
    end # realisations
end # scan parameters
prinz(z)

Pavg = mean(P,dims=2)
Cavg = Array{Vector{Float64},1}(undef,length(Ns))
    for n in eachindex(Ns)
        Cavg[n] = mean([Crt[n,r] for r in 1:R])
    end

p=plot(legend=:bottomright,ylims=(0,1))
    for n in eachindex(Ns)
        rr = collect(0:dr:round(Int,sqrt(R0^2*pi*Ns[n]/α)/2))
        plot!(rr,Cavg[n],label="N = $(Ns[n])")
    end
    p


## Analysis results cluster
@unpack alphas,Ts,v0s,tmax,times,P,runtime,Ns,R0,sigmas = load("data/phasediag_sig_v0_alpha=1_T0.jld")
P
collect(sigmas)
collect(v0s)
Pavg = mean(P,dims=7)
# Movies
phdiag_sig_vo = @animate for i in eachindex(times)
    heatmap(sigmas,v0s,Pavg[1,1,1,:,:,i,1],c=cgrad([:red,:orange, :green]),size=(500,400),clims=(0,1),xlabel="σ",ylabel="v_0/L",colorbar_title="P",title="t=$(times[i])")
end
mp4(phdiag_sig_vo,"films/phase_diagram_sig_v0_alpha=1_T0.mp4",fps=15)


## Analysis of simulations from cluster
@unpack alphas,Ts,v0s,tmax,times,P,runtime,N,L,R0,sigmas = load("data/phasediag_large_sigmas_Ts_alpha=1_v0=0.2.jld")
v0s
Ts
sigmas
P             #  (alphas),(Ts),(v0s),(sigmas),(times),R
Pavg = mean(P,dims=6)
Pavg = nanmean(P,6)

# Movies
anim_varying_time = @animate for i in eachindex(times)
    heatmap(sigmas,Ts,Pavg[1,:,1,:,i,1],c=cgrad([:red,:orange, :green]),size=(500,400),clims=(0,1),xlabel="σ",ylabel="T",colorbar_title="P",title="t=$(round(times[i],digits=1))")
end
mp4(anim_varying_time,"films/phase_diagram_large_sigT_varying_time_vo_0.2.mp4",fps=20)

#
P
