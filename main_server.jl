include("IDrealisation.jl") ;
# using JLD,StatsBase,Distributions,LinearAlgebra,Parameters,Statistics
using JLD2,LinearAlgebra,Statistics
include("methods.jl");
global const R0 = 1

## Phase Diagram
Ns = Int.(1E4)
    rhos = [1]
    Ts     = 0.1
    # v0s    = [0,0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5,1,2,3,5]
    # sigmas = [0,0.1,0.2,0.3,0.5]
    v0s    = [0.01,0.05,0.1,0.2,1,5]
    sigmas = [0,0.15,0.3]

    tmax = 1E4; times_log = logspace(0.1,tmax,31)

    P  = zeros(length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(times_log))
    n  = zeros(length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(times_log))
    @assert length(Ns) == 1
    @assert length(rhos) == 1
    dr = 0.5 ; C = zeros(length(0:dr:round(Int,sqrt(Ns)/2)),length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(times_log))
    m  = 0 ; M = length(Ns)*length(rhos)*length(Ts)*length(v0s)*length(sigmas)

z = @elapsed for a in eachindex(Ns) , i in eachindex(rhos) , j in eachindex(Ts) , k in eachindex(v0s) , q in eachindex(sigmas)
    N = Ns[a] ; rho = rhos[i] ; T = Ts[j] ; v0 = v0s[k] ; σ = sigmas[q]

    L = round(Int,sqrt(N/rho))
    dt = determine_dt(T,σ,v0,N,rho)
    # params = Any[T,v0,N,L,rho,σ]

    # m += 1 ; println("Simulation $m/$M with dt = $(round(dt,digits=2))")

    pos,thetas,psis,omegas = initialisation(N,L,σ)
    t = 0.0 ; token = 1
    while t < tmax
        t += dt
        pos,thetas = update(pos,thetas,psis,omegas,T,v0,N,L,dt)
        if t ≥ times_log[token]
            P[a,i,j,k,q,token]   = polarOP(thetas)[1]
            n[a,i,j,k,q,token]   = number_defects(pos,thetas,N,L)
            C[:,a,i,j,k,q,token] = corr_fast(pos,thetas,N,L,dr)
            token = min(token+1,length(times_log))
        end
    end
end
prinz(z)

comments = ""
filename = "data/scan_v0_sigma_rho1_N1E3_tmax1E4_r$real.jld2"
JLD2.@save filename P n C Ts Ns v0s rhos sigmas times_log tmax runtime=z comments
