include("IDrealisation.jl") ;
# using JLD,StatsBase,Distributions,LinearAlgebra,Parameters,Statistics
using JLD2,LinearAlgebra,Statistics, Hungarian
include("methods.jl");
global const R0 = 1

## Defects Motion
# N = Int(1E4)
#     rho = 2
#     aspect_ratio = 1 # Lx/Ly
#     Lx  = round(Int,sqrt(N/rho*aspect_ratio))
#     Ly  = round(Int,sqrt(N/rho/aspect_ratio))
#     T  = 0.1 # to be compared to Tc ~ 1 when \sigma = 0
#     # In vitro (green region only)
#     v_sigmas = [ (0,0) , (0.01,0) , (0.5,0) , (0.5,0.1)  , (0.5,0.2) , (1,0.1)  , (1,0.2)]
#     # In vivo (whole phase space)
#     # v_sigmas = [ (0,0) , (0.05,0) , (0.1,0) , (0.5,0) , (1,0) ,
#     #              (0,0.1) , (0.05,0.1) , (0.1,0.1) , (0.5,0.1) , (1,0.1) ,
#     #              (0.4,0.2) , (0.65,0.3)]
#     # v_sigmas = [ (0.5,0.1) ]
#     transients = 0 ; tmax = 10 ;
#     # times = logspace(round(Int,transients),tmax,50,digits=1)
#     times = 5:5:tmax
#
# dfts = Array{DefectTracker,1}(undef,length(v_sigmas))
# z = @elapsed for i in each(v_sigmas)
#     v0,σ = v_sigmas[i]
#     dt = determine_dt(T,σ,v0,N,rho)
#     t = 0. ; pos,thetas,omegas,psis = initialisation(N,Lx,Ly,σ,["pair",round(Int,Lx/2)])
#     ind_neighbours = get_list_neighbours(pos,N,Lx,Ly)
#     pos,thetas = update(pos,thetas,omegas,psis,ind_neighbours,T,v0,N,Lx,Ly,dt)
#
#     pos,thetas = evolve(pos,thetas,omegas,psis,T,v0,N,Lx,Ly,dt,transients) # takes care of v0 = 0 if needed
#
#     dft = DefectTracker(pos,thetas,N,Lx,Ly,t)
#     dft,pos,thetas,t = track!(dft,pos,thetas,psis,omegas,T,v0,N,Lx,Ly,dt,t,times)
#     dfts[i] = dft
# end
# prinz(z)
#
# comments = "Study the motion of a pair of defects at various locations in phase space. Goal: understand the role of v and sigma on the MSD. Preliminary results indicate that MSD ~ t even for σ>0. Understand the transition from immobile (MSD ~ t^1.5) to mobile (MSD ~ t)."
# filename = "data/defects_motion_IN_VITRO_v0s_sigmas_N1E4_r$real.jld2"
# JLD2.@save filename dfts T N rho v_sigmas times transients tmax runtime=z comments


## Phase Diagram
Ns = Int.(1E3)
    rhos = [1,2]
    Ts     = 0.1
    # v0s = vcat(0,logspace(0.001,3,15,digits=3))
    v0s    = [0]
    sigmas = collect(0:0.1:0.4)
    inits  = ["disordered"]

    tmax = 1E3; times_log = logspace(0.1,tmax,32)

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

    pos,thetas,omegas,psis = initialisation(N,Lx,Ly,σ,[init])
    ind_neighbours = get_list_neighbours(pos,N,Lx,Ly)
    t = 0.0 ; token = 1
    while t < tmax
        t += dt
        pos,thetas = update(pos,thetas,omegas,psis,ind_neighbours,T,v0,N,Lx,Ly,dt)
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

comments = "Sanity check for v0 = σ = 0, should be XY model since ρ = 2 > 1.44 (critical ρ). For σ, in the green region of phase space, if it exists at v0 = 0, then something contradicts my previous work: let's see what. No Defect Tracking here. For the green zone, the tracking can be found in file defects_motion_IN_VITRO_v0s_sigmas_N1E4.jld2"
filename = "data/investigation_v00_r$real.jld2"
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
#             t = 0. ; pos,thetas,omegas,psis = initialisation(N,L,σ,["hightemp"])
#             token = 1
#             while t<tmax
#                 t += dt
#                 pos,thetas = update(pos,thetas,omegas,psis,T,v0,N,L,dt)
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
