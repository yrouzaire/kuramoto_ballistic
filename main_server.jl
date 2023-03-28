include("IDrealisation.jl");
# using JLD,StatsBase,Distributions,LinearAlgebra,Parameters,Statistics
using JLD2, LinearAlgebra, Statistics, Hungarian
include("methods.jl");

# ## ---------------- Impact of init on XY Model ---------------- ##
# comments = "Investigates the impact of initialisation for the spatial location of the 
# spins in the XY model. Even though no symmetry is broken, it seems that having the spins
# on (A) a regular lattice (B) a random lattice (C) a 2D Sobol sequence (D) a RSA lattice changes the behaviour of the defects.
# The idea is to investigate the role of the fluctuations in the initialisation of the spins on the defects dynamics."
# # Physical Params 
# Ntarget = Int(1E4)
# aspect_ratio = 1
# T = 0.1
# R0 = 1
# rho = 2
# v0 = 0
# sigma = 0
# rhoc = 4.51 / π

# # Initialisation parameters
# inits_pos = ["random", "square_lattice", "Sobol", "RSA"]
# init_theta = "hightemp"
# r0 = 20.0
# q = 1.0
# params_init = Dict(:init_pos => NaN, :init_theta => init_theta, :r0 => r0, :q => q)

# # Simulation parameters
# tmax = 1E2
# times = logspace(1,tmax,10)


# P = zeros(length(inits_pos), length(times))
# C = Array{Vector{Float64}}(undef, length(inits_pos), length(times))
# xi = zeros(length(inits_pos), length(times))
# n = zeros(length(inits_pos), length(times))


# # Impact on the usual quantities
# z = @elapsed for i in each(inits_pos)
# 	init_pos = inits_pos[i]

#     println("Init : $(init_pos)")
#     N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
#     dt = determine_dt(T, sigma, v0, N, rho)

#     params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
    
#     param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
#         :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
#         :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

#     t = 0.0
#     system = System(param)

# 	for tt in eachindex(times)
# 		evolve(system, times[tt]) # evolves the systems up to times[tt]
		
# 		P[i,tt]  = polarOP(system)[1]
# 		corr_tmp = corr(system)
# 		C[i,tt]  = corr_tmp
# 		xi[i,tt] = corr_length(corr_tmp)
# 		n[i,tt]  = number_defects(system)
# 	end
# end
# prinz(z)

# filename = "data/impact_init_XY_r$real.jld2"
# JLD2.@save filename Ntarget v0 sigma inits_pos rho params_init T P C n xi aspect_ratio times tmax comments rhoc runtime = z



# ## ---------------- Tracking a pair of defects ---------------- ##
# comments = "From the defects data, one will be able to infer : \n
# A. the separating distance between the two defects R(t) \n
# B. the MSD and diffusion coeff of an individual defect. "
# # Physical Params 
# Ntarget = Int(1E4)
# aspect_ratio = 1
# T = 0.1
# R0 = 1
# rho = 1
# rhoc = 4.51 / π

# # Initialisation parameters
# init_pos = "random"
# init_theta = "pair"
# r0 = 20.0
# q = 1.0
# params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => NaN, :q => q)

# # Simulation parameters
# v0sigs = [(0.5,0),(0.5,0.1)]
# r0s = 5:10:35
# tmax = 1E2
# times = 0:5:tmax # linear time

# z = @elapsed for i in each(v0sigs), j in each(r0s)
#     v0, sigma = v0sigs[i]
#     r0 = r0s[j]

#     println("v0 = $v0, σ = $sigma, r0 = $r0 ,  $(100i/length(v0sigs))%")
#     N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
#     dt = determine_dt(T, sigma, v0, N, rho)

#     params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
    
#     param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
#         :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
#         :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

#     t = 0.0
#     system = System(param)
#     dft = DefectTracker(system, t)
#     dft, pos, thetas, t = track!(dft,system,times)
# end
# prinz(z)

# filename = "data/DFT_pair_r$real.jld2"
# #JLD2.@save filename Ntarget v0sigs rho params_init T P C n xi aspect_ratio times tmax comments rhoc runtime = z
# #to change


# ## ---------------- Nature of the Phase Transition ---------------- ##
comments = "The goal of this script is to pass through the transition line, 
in both direction (keeping σ or v0 constant) and to compute correlation functions.
Here for T = 0.1 and ρ = 1."
# Physical Params 
Ntarget = Int(1E4)
aspect_ratio = 1
T = 0.1
R0 = 1
rho = 1
rhoc = 4.51 / π

# Initialisation parameters
init_pos = "random"
init_theta = "hightemp"
r0 = 20.0
q = 1.0
params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)

# Simulation parameters
v0sigs = vcat([(0.2,sigm) for sigm in 0:0.025:0.225],[(v,0.1) for v in logspace(0.03,1,10,digits=3)])
tmax = 3E4
# times = collect(0:tmax/30:tmax) # linear time
times = logspace(1,tmax,30,digits=1) # log time

P = zeros(length(v0sigs), length(times))
C = Array{Vector{Float64}}(undef, length(v0sigs), length(times))
xi = zeros(length(v0sigs), length(times))
n = zeros(length(v0sigs), length(times))

z = @elapsed for i in each(v0sigs)
    v0, sigma = v0sigs[i]
    println("v0 = $v0, σ = $sigma, $(100i/length(v0sigs))%")
    N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
    dt = determine_dt(T, sigma, v0, N, rho)

    param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
        :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
        :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

    system = System(param)

    t = 0.0
    token = 1

    for tt in eachindex(times)
        evolve(system, times[tt]) # evolves the systems up to times[tt]
        
        P[i,tt]  = polarOP(system)[1]
        corr_tmp    = corr(system)
        C[i,tt]  = corr_tmp
        xi[i,tt] = corr_length(corr_tmp)
        n[i,tt]  = number_defects(system)
    end

end
prinz(z*19*100)

filename = "data/nature_phase_transition_r$real.jld2"
JLD2.@save filename Ntarget v0sigs rho params_init T P C n xi aspect_ratio times tmax comments rhoc runtime = z


# ## ---------------- No Hysteresis and Nature Phases ---------------- ##
# comments = "The goal of this script is to show that there is no hysteresis 
# and to compute correlation functions at steady state to determine the nature
# of both phases and the transition between them."
# # Physical Params 
# Ntarget = Int(1E3)
# aspect_ratio = 1
# Ts = [0,0.1]
# R0 = 1
# rhos = [1,2]
# rhoc = 4.51 / π

# # Initialisation parameters
# init_pos = "random"
# inits_theta = ["hightemp", "lowtemp"]
# r0 = 20.0
# q = 1.0

# # Simulation parameters
# v0sigs = [(0,0),(0,0.1), (0.05,0), (0.05,0.1), (1,0), (1,0.1)]
# tmax = 1E2
# # times = collect(0:tmax/30:tmax) # linear time
# times = logspace(1,tmax,30,digits=1) # log time

# P = zeros(length(v0sigs), length(Ts), length(rhos),length(inits_theta), length(times))
# C = Array{Vector{Float64}}(undef, length(v0sigs), length(Ts), length(rhos),length(inits_theta), length(times))
# xi = zeros(length(v0sigs), length(Ts), length(rhos),length(inits_theta), length(times))
# n = zeros(length(v0sigs), length(Ts), length(rhos),length(inits_theta), length(times))

# z = @elapsed for i in each(v0sigs), j in each(Ts), k in each(rhos), l in each(inits_theta)
#     v0, sigma = v0sigs[i]
# 	T = Ts[j]
# 	rho = rhos[k]
#     println("v0 = $v0, σ = $sigma, T = $T, ρ = $rho")
#     N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
#     dt = determine_dt(T, sigma, v0, N, rho)

# 	params_init = Dict(:init_pos => init_pos, :init_theta => inits_theta[l], :r0 => r0, :q => q)

#     param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
#         :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
#         :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

#     system = System(param)

#     t = 0.0
#     token = 1

#     for tt in eachindex(times)
#         evolve(system, times[tt]) # evolves the systems up to times[tt]
        
#         P[i,j,k,l,tt]  = polarOP(system)[1]
#         corr_tmp    = corr(system)
#         C[i,j,k,l,tt]  = corr_tmp
#         xi[i,j,k,l,tt] = corr_length(corr_tmp)
#         n[i,j,k,l,tt]  = number_defects(system)
#     end

# end
# prinz(z)

# filename = "data/hysteresis_nature_phases_r$real.jld2"
# JLD2.@save filename Ntarget v0sigs rhos inits_theta Ts P C n xi aspect_ratio times tmax comments rhoc runtime = z


# ## ---------------- Finite Size Scaling ---------------- ##
# comments = "The goal of this script is to evaluate the change in polarisation for given (v0, σ) when L varies."
# # Physical Params 
# aspect_ratio = 1
# T = 0.1
# R0 = 1
# rho = 1.0
# rhoc = 4.51 / π

# # Initialisation parameters
# init_pos = "random"
# init_theta = "hightemp"
# r0 = 20.0
# q = 1.0
# params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)

# # Simulation parameters
# v0sigs = [(1,0), (1,0.1)]
# Ntargets = Int.(logspace(1E2,2E4,10,digits=0))

# tmax = 1E4
# # times = collect(0:tmax/30:tmax) # linear time
# times = logspace(1,tmax,30,digits=1) # log time

# P = zeros(length(v0sigs), length(Ntargets), length(times))
# C = Array{Vector{Float64}}(undef, length(v0sigs), length(Ntargets), length(times))
# xi = zeros(length(v0sigs), length(Ntargets), length(times))
# n = zeros(length(v0sigs), length(Ntargets), length(times))

# z = @elapsed for i in each(v0sigs), nn in each(Ntargets)
#     v0, sigma = v0sigs[i]
#     Ntarget = Ntargets[nn]
#     println("v0 = $v0, σ = $sigma, Ntarget = $Ntarget")
#     N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
#     dt = determine_dt(T, sigma, v0, N, rho)

#     param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
#         :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
#         :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

#     system = System(param)

#     t = 0.0
#     token = 1

#     for tt in eachindex(times)
#         evolve(system, times[tt]) # evolves the systems up to times[tt]
        
#         P[i,nn,tt]  = polarOP(system)[1]
#         corr_tmp    = corr(system)
#         C[i,nn,tt]  = corr_tmp
#         xi[i,nn,tt] = corr_length(corr_tmp)
#         n[i,nn,tt]  = number_defects(system)
#     end

# end
# prinz(z)

# filename = "data/FSS_green_r$real.jld2"
# JLD2.@save filename Ntargets v0sigs P C n xi params_init aspect_ratio rho times tmax T comments rhoc runtime = z



## ---------------- Critical velocity at sigma = 0 for ρ < ρc ---------------- ##
# Physical Params 
# Ntarget = Int(1E3)
# aspect_ratio = 1
# T = 0.1
# sigma = 0.0
# v0 = 1.0
# R0 = 1
# rhoc = 4.51 / π

# # Initialisation parameters
# init_pos = "random"
# init_theta = "hightemp"
# r0 = 20.0
# q = 1.0
# params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)

# # Simulation parameters
# v0s = logspace(1e-3, 0.3, 20, digits=4)
# rhos = collect(1:0.05:rhoc)

# seuil = 0.5 # above P = 0.5, we consider the system to be ordered
# tmax = 3000
# times = 0:tmax/10:tmax

# critical_velocity = v0s[end] * ones(length(rhos))
# z = @elapsed for k in each(rhos)
#     rho = rhos[k]
#     for i in each(v0s)
#         v0 = v0s[i]
#         println("ρ = $rho, v0 = $v0")
#         sigma = 0
#         N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
#         dt = determine_dt(T, sigma, v0, N, rho)

#         param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
#             :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
#             :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

#         system = System(param)

#         t = 0.0
#         token = 1

#         isordered = false
#         while t < tmax
#             t += dt
#             update!(system, N, Lx, Ly, R0)
#             if t ≥ times[token]
#                 token = min(token + 1, length(times))
#                 P = polarOP(system)[1]
#                 if P > seuil
#                     isordered = true
#                     break # gets out of the while loop only 
#                 end
#             end
#         end
#         if isordered
#             critical_velocity[k] = v0
#             P = polarOP(system)[1]
#             println("Critical velocity found for ρ = $rho : $v0 because P = $P")
#             break # gets out of the for loop scanning v0s
#         end
#     end
# end
# prinz(z)

# # comments = "Critical velocity vc against the density ρ. 
# #             From hightemp to make sure the system really has crossed the frontier.
# #             Indeed, here since there is no σ to pertube the dynamics, 
# #             I prefer to start from a disordered state and see whether it can order."
# comments = ""
# filename = "data/critical_velocity_r$real.jld2"
# JLD2.@save filename Ntarget aspect_ratio rhos times tmax critical_velocity T v0s sigmas seuil comments rhoc runtime = z

