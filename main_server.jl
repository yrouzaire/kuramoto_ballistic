include("IDrealisation.jl") ;
# using JLD,StatsBase,Distributions,LinearAlgebra,Parameters,Statistics
using JLD2,LinearAlgebra,Statistics
include("methods.jl");
global const R0 = 1

## Phase Diagram
Ns = Int.(1E3)
    rhos = [1,2]
    Ts     = 0.1
    v0s    = [0.0,0.01,0.016,0.025,0.04,0.063,0.1,0.16,0.25,0.4,0.63,1,1.6,2.5,4,6.3]
    # v0s    = [0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.5,1,2,3,5]
    sigmas = collect(0:0.05:0.5)
    inits  = ["ordered","disordered"]

    tmax = 2E3; times_log = logspace(0.1,tmax,27)

    P  = zeros(length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),length(times_log))
    n  = zeros(length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),length(times_log))
    # @assert length(Ns) == 1
    # @assert length(rhos) == 1
    # dr = 0.5 ; C = zeros(length(0:dr:round(Int,sqrt(Ns)/2)),length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),length(times_log))

    pos_saved  = zeros(Float16,2,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits))
    thetas_saved = zeros(Float16,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits))
    psis_saved   = zeros(Float16,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits))
    omegas_saved = zeros(Float16,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits))

z = @elapsed for a in eachindex(Ns) , i in eachindex(rhos) , j in eachindex(Ts) , k in eachindex(v0s) , q in eachindex(sigmas) , ini in eachindex(inits)
    N = Ns[a] ; rho = rhos[i] ; T = Ts[j] ; v0 = v0s[k] ; σ = sigmas[q] ; init = inits[ini]

    L = round(Int,sqrt(N/rho))
    dt = determine_dt(T,σ,v0,N,rho)
    # params = Any[T,v0,N,L,rho,σ]


    pos,thetas,psis,omegas = initialisation(N,L,σ,[init])
    t = 0.0 ; token = 1
    while t < tmax
        t += dt
        pos,thetas = update(pos,thetas,psis,omegas,T,v0,N,L,dt)
        if t ≥ times_log[token]
            P[a,i,j,k,q,ini,token]   = polarOP(thetas)[1]
            n[a,i,j,k,q,ini,token]   = number_defects(pos,thetas,N,L)
            # C[:,a,i,j,k,q,ini,token] = corr_fast(pos,thetas,N,L,dr)
            token = min(token+1,length(times_log))
        end
    end

    pos_saved[:,:,i,j,k,q,ini]  = pos
    thetas_saved[:,i,j,k,q,ini] = thetas
    psis_saved[:,i,j,k,q,ini]   = psis
    omegas_saved[:,i,j,k,q,ini] = omegas
end
prinz(z)

comments = "No C(r,t) here."
filename = "data/scan_for_discussion_Parisa_v0_sigma_rho_N1E3_tmax1E4_r$real.jld2"
JLD2.@save filename P n Ts Ns v0s rhos inits sigmas times_log tmax runtime=z comments pos_saved thetas_saved psis_saved omegas_saved
# WARNING : "C" n'est plus dans les variables saved, il faut la remettre !
