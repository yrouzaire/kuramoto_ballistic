include("IDrealisation.jl");
# using JLD,StatsBase,Distributions,LinearAlgebra,Parameters,Statistics
using JLD2, LinearAlgebra, Statistics, Hungarian
include("methods.jl");

## ---------------- Finite Size Scaling ---------------- ##
comments = "The goal of this script is to evaluate the change in polarisation for given (v0, σ) when L varies."
# Physical Params 
aspect_ratio = 1
T = 0.1
R0 = 1
rho = 1.0
rhoc = 4.51 / π

# Initialisation parameters
init_pos = "random"
init_theta = "hightemp"
r0 = 20.0
q = 1.0
params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)

# Simulation parameters
v0sigs = [(1,0), (1,0.1), (0.01,0.1)]
Ntargets = Int.(logspace(1E2,1E4,5,digits=0))
Ntargets = Int(1E5)

tmax = 100
times = 0:tmax/10:tmax

P = zeros(length(v0sigs), length(Ntargets), length(times))
C = Array{Vector{Float64}}(undef, length(v0sigs), length(Ntargets), length(times))
xi = zeros(length(v0sigs), length(Ntargets), length(times))
n = zeros(length(v0sigs), length(Ntargets), length(times))

z = @elapsed for i in each(v0sigs), nn in each(Ntargets)
    v0, sigma = v0sigs[i]
    Ntarget = Ntargets[nn]
    println("v0 = $v0, σ = $sigma, Ntarget = $Ntarget")
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
        
        P[i,nn,tt]  = polarOP(system)[1]
        corr_tmp    = corr(system)
        C[i,nn,tt]  = corr_tmp
        xi[i,nn,tt] = corr_length(corr_tmp)
        n[i,nn,tt]  = number_defects(system)
    end

end
prinz(z)

##
plot(P[1,1,:], label="v0 = $(v0sigs[1][1]), σ = $(v0sigs[1][2])")



##



filename = "data/FSS_r$real.jld2"
JLD2.@save filename Ntargets v0sigs P C n xi params_init aspect_ratio rho times tmax T comments rhoc runtime = z



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

