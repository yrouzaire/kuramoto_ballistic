include("IDrealisation.jl") ;
# using JLD,StatsBase,Distributions,LinearAlgebra,Parameters,Statistics
using JLD2,LinearAlgebra,Statistics
include("methods.jl");
global const R0 = 1


## Phase Diagram
Ns = Int.(1E3)
    rhos = collect(1:0.15:2)
    Ts     = 0.1
    v0s = vcat(0,logspace(0.001,3,15,digits=3))
    # v0s    = [0,0.05]
    sigmas = collect(0:0.05:0.4)
    inits  = ["disordered"]

    tmax = 2E3; times_log = logspace(0.1,tmax,32)

    P  = zeros(length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),length(times_log))
    n  = zeros(length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),length(times_log))
    dr = 0.5 ; C = Array{Vector{Float32},7}(undef,length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),length(times_log))
    pos_saved  = zeros(Float16,2,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits))
    thetas_saved = zeros(Float16,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits))
    psis_saved   = zeros(Float16,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits))
    omegas_saved = zeros(Float16,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits))

z = @elapsed for a in eachindex(Ns) , i in eachindex(rhos) , j in eachindex(Ts) , k in eachindex(v0s) , q in eachindex(sigmas) , ini in eachindex(inits)
    N = Ns[a] ; rho = rhos[i] ; T = Ts[j] ; v0 = v0s[k] ; σ = sigmas[q] ; init = inits[ini]

    Lx = Ly = round(Int,sqrt(N/rho))
    dt = determine_dt(T,σ,v0,N,rho)
    # params = Any[T,v0,N,L,rho,σ]

    pos,thetas,psis,omegas = initialisation(N,Lx,Ly,σ,[init])
    t = 0.0 ; token = 1
    while t < tmax
        t += dt
        pos,thetas = update(pos,thetas,psis,omegas,T,v0,N,Lx,Ly,dt)
        if t ≥ times_log[token]
            P[a,i,j,k,q,ini,token] = polarOP(thetas)[1]
            n[a,i,j,k,q,ini,token] = number_defects(pos,thetas,N,Lx,Ly)
            C[a,i,j,k,q,ini,token] = corr_fast(pos,thetas,N,Lx,Ly,dr)
            token = min(token+1,length(times_log))
        end
    end

    pos_saved[:,:,i,j,k,q,ini]  = pos
    thetas_saved[:,i,j,k,q,ini] = thetas
    psis_saved[:,i,j,k,q,ini]   = psis
    omegas_saved[:,i,j,k,q,ini] = omegas
end
prinz(z)

comments = ""
filename = "data/nature_phase_impact_rho_N1E4_tmax1E4_r$real.jld2"
JLD2.@save filename P C n Ts Ns v0s rhos inits sigmas times_log tmax runtime=z comments pos_saved thetas_saved psis_saved omegas_saved

## FSS
# Ns = round.(Int,logspace(1E2,3E4,15,digits=0))
# rho = 1
# T = 0.1
# v0s    = reverse(logspace(0.01,4,30,digits=2))
# sigmas = 0:0.1:0.5
# tmax = 1E3
#
# vc = zeros(length(Ns),length(sigmas))
# seuil_break = 0.8
# seuil_crit = 0.3
# times_break = range(1, stop=tmax, length=20)
#
# z = @elapsed for j in each(sigmas)
#     σ = sigmas[j]
#     for n in each(Ns)
#         N = Ns[n] ; L = round(Int,sqrt(N/rho))
#         for i in each(v0s)
#             v0 = v0s[i]
#             println("σ = $σ, N = $N, v0 = $v0")
#             dt = determine_dt(T,σ,v0,N,rho)
#             t = 0. ; pos,thetas,psis,omegas = initialisation(N,L,σ,["hightemp"])
#             token = 1
#             while t<tmax
#                 t += dt
#                 pos,thetas = update(pos,thetas,psis,omegas,T,v0,N,L,dt)
#                 if t ≥ times_break[token]
#                     if polarOP(thetas)[1] > seuil_break
#                         println("Simu stopped because P>$(seuil_break) : v0 > vc")
#                         break # stop this simulation, pass to the next v0
#                     end
#                     token = min(token+1,length(times_break))
#                 end
#             end
#             if polarOP(thetas)[1] < seuil_crit
#                 vc[n,j] = v0
#                 println("Critical velocity because P<$(seuil_crit) : v0 < vc")
#                 break # you found the critical velocity
#             end
#         end
#     end
# end
# prinz(z)
#
# comments = "Study the scaling of vc the critical velocity against N, for different sigmas"
# filename = "data/criticalv0_r$real.jld2"
# JLD2.@save filename vc Ns rho seuil_break seuil_crit times_break v0s T sigmas tmax
