include("IDrealisation.jl");
# using JLD,StatsBase,Distributions,LinearAlgebra,Parameters,Statistics
using JLD2, LinearAlgebra, Statistics, Hungarian
include("methods.jl");

# ## ---------------- Proba Spinwaves ---------------- ##
# # Fixed important params 
# Ntarget = Int(1E4)
# aspect_ratio = 1
# rho = 1
# T = 0.1
# R0 = 1
# N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
# init_pos = "random"
# init_theta = "hightemp"

# # Useless params
# r0 = round(Int,Lx/2)
# q = 1.0
# phonons = false
# if phonons @assert v0 == 0.0 "Phonons only make sense for immobile particles! " end
# if phonons @assert init_theta ≠ "single" "Phonons only make sense for PBC! " end
# if phonons @assert aspect_ratio == 1 "Phonons only implemented for square box ! (for now) " end
# phonon_amplitude = 1
# phonon_k = 1*(2π/Lx) # wavenumber
# phonon_omega = 0 # "frequency" (up to a factor 2π)

# # Scanned params
# p_threshold = 0.5
# sigmas = [0,0.1,0.2]
# v0s = collect(0.5:0.5:5) 
# # v0s = [5]

# R_per_core = 25
# tmax = 1E3 # max time
# times = tmax/50:tmax/50:tmax

# m = 0 
# M = length(v0s)*length(sigmas)*R_per_core
# nb_detected_spinwave = zeros(Int,length(v0s),length(sigmas))
# times_detected_spinwave = [[] for i in each(v0s), j in each(sigmas)]
# Ps_detected_spinwave = [[] for i in each(v0s), j in each(sigmas)]
# thetas_detected_spinwave = [Vector{Float16}[] for i in each(v0s), j in each(sigmas)]
# pos_detected_spinwave = [Vector{Tuple{Number,Number}}[] for i in each(v0s), j in each(sigmas)]
# systems_detected_spinwave = [System[] for i in each(v0s), j in each(sigmas)]


# z = @elapsed for i in each(v0s)
#     for j in each(sigmas)
#         v0 = v0s[i]
#         sigma = sigmas[j]

#         params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
#         params_phonons = Dict(:phonons => phonons, :phonon_amplitude => phonon_amplitude, :phonon_k => phonon_k, :phonon_omega => phonon_omega)
#         param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
#         :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
#         :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init, :params_phonons => params_phonons)

#         for r in 1:R_per_core
#             global m += 1
#             println("v0 = $v0, σ = $sigma, r = $r : Simu #$m/$M")
#             system = System(param)

#             for tt in eachindex(times)
#                 evolve!(system, times[tt])
#                 n = number_defects(system)
#                 if n == 0
#                     P = round(polarOP(system)[1],digits=2)
#                     if P < p_threshold 
#                         println("v0 = $v0, σ = $sigma, P = $P < $(p_threshold)")
#                         nb_detected_spinwave[i,j] += 1
#                         push!(times_detected_spinwave[i,j],times[tt])
#                         push!(Ps_detected_spinwave[i,j],P)
#                         push!(thetas_detected_spinwave[i,j],get_thetas(system))
#                         push!(pos_detected_spinwave[i,j],get_pos(system))
#                         push!(systems_detected_spinwave[i,j],system)
#                         # p=plot_thetas(system)
#                         # display(p)
#                     end
#                     println("n=0, simulation stopped at t = $(times[tt])")
#                     break # en dehors du if parce que si n == 0 de toute facon continuer la simu n'a aucun sens 
#                 end
#             end
#         end
#     end
# end
# prinz(z)

# # systems_detected_spinwave
# # nb_detected_spinwave
# # pos_detected_spinwave
# # thetas_detected_spinwave
# # Ps_detected_spinwave
# # times_detected_spinwave

# filename = "data/proba_spinwaves_r$real.jld2"
# JLD2.@save filename nb_detected_spinwave times_detected_spinwave systems_detected_spinwave Ps_detected_spinwave thetas_detected_spinwave pos_detected_spinwave R_per_core sigmas v0s tmax times p_threshold init_pos init_theta Ntarget rho T aspect_ratio runtime = z


# ## ---------------- MSD Specifically Tracking a Pair of defects  ---------------- ##
# comments = "From the defect data one can infer the MSD and diffusion coeff of an individual defect. "
# # Physical Params 
# Ntarget = Int(4E3)
# aspect_ratio = 1
# R0 = 1
# rho = 1 
# rhoc = 4.51 / π
# init_theta = "pair"
# init_pos = "random"
# q = 1.0
# r0 = 28
# phonons = false ; phonon_amplitude = 1 ; phonon_k = 1  ; phonon_omega = 0 
# params_phonons = Dict(:phonons => phonons, :phonon_amplitude => phonon_amplitude, :phonon_k => phonon_k, :phonon_omega => phonon_omega)
# params_init = Dict(:init_pos => NaN, :init_theta => init_theta, :r0 => NaN, :q => q)

# R_per_core = 10

# tmax = 2000
# times = collect(0:5:tmax) # linear time

# sigmas = [0,0.05,0.1]
# # sigmas = [0.1]

# Ts = [0.1,0.2,0.3,0.4]
# Ts = [0.1]

# v0s = collect(0.5:0.25:3)
# # v0s = [5]

# xy_pos = Array{Vector{Tuple{Number,Number}}}(undef,length(v0s),length(sigmas),length(Ts),R_per_core)
# xy_neg = Array{Vector{Tuple{Number,Number}}}(undef,length(v0s),length(sigmas),length(Ts),R_per_core)
# rr = Array{Vector{Number}}(undef,length(v0s),length(sigmas),length(Ts),R_per_core)
# times_collision = times[end]*ones(length(v0s),length(sigmas),length(Ts),R_per_core)

# z = @elapsed for i in each(v0s), j in each(sigmas), k in each(Ts), r in 1:R_per_core
#     v0 = v0s[i]
#     sigma = sigmas[j]
# 	T = Ts[k]

#     println("v0 = $v0, σ = $sigma, T = $T")
#     N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
    
#     params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
#     param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
#         :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
#         :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init, :params_phonons => params_phonons)

#     t = 0.0
#     system = System(param)

#     # t = 0
#     defects_pos, defects_neg =  spot_defects(system)
#     xy_pos_tmp = [defects_pos[1][1:2]]
#     xy_neg_tmp = [defects_neg[1][1:2]]
#     r_tmp = [dist(defects_pos[1][1:2],defects_neg[1][1:2],Lx,Ly)]

#     for tt in 2:length(times)
#         evolve!(system, times[tt])

#         defects_pos, defects_neg =  spot_defects(system)
#         @assert length(defects_pos) == length(defects_neg)
#         if length(defects_pos) == 1
#             push!(xy_pos_tmp,defects_pos[1][1:2])
#             push!(xy_neg_tmp,defects_neg[1][1:2])
#             push!(r_tmp,dist(defects_pos[1][1:2],defects_neg[1][1:2],Lx,Ly))
#         elseif length(defects_pos) > 1 
#             #= If there are more than one defect, 
#             consider the closest defect as the most probable. =#
#             distance_pos_tmp = Inf 
#             index_closest_pos_tmp = -1
#             for i in each(defects_pos)
#                 d = dist(defects_pos[i][1:2], xy_pos_tmp[end], Lx, Ly)
#                 if d < distance_pos_tmp
#                     distance_pos_tmp = d
#                     index_closest_pos_tmp = i
#                 end
#             end
#             distance_neg_tmp = Inf 
#             index_closest_neg_tmp = -1
#             for i in each(defects_neg)
#                 d = dist(defects_neg[i][1:2], xy_neg_tmp[end], Lx, Ly)
#                 if d < distance_neg_tmp
#                     distance_neg_tmp = d
#                     index_closest_neg_tmp = i
#                 end
#             end
#             push!(xy_pos_tmp,defects_pos[index_closest_pos_tmp][1:2])
#             push!(xy_neg_tmp,defects_neg[index_closest_neg_tmp][1:2])
#             push!(r_tmp,dist(defects_pos[index_closest_pos_tmp][1:2],defects_neg[index_closest_neg_tmp][1:2],Lx,Ly))
#         elseif length(defects_pos) == 0 
#             println("No defects ! Simulation stopped at t = $(times[tt]).")
#             times_collision[i,j,k,r] = times[tt]
#             break
#         end
#     end
#     xy_pos[i,j,k,r] = xy_pos_tmp
#     xy_neg[i,j,k,r] = xy_neg_tmp
#     rr[i,j,k,r] = r_tmp
# end
# prinz(z) 


# filename = "data/mobility_defects_sigma_v0_r$real.jld2"
# JLD2.@save filename sigmas v0s Ts xy_pos xy_neg rr times_collision R_per_core params_init Ntarget R0 q init_theta init_pos aspect_ratio times tmax comments rhoc runtime = z

## ---------------- MSD Tracking defects  ---------------- ##
comments = "From the defect data one can infer the MSD and diffusion coeff of an individual defect. "
# Physical Params 
Ntarget = Int(4E3)
aspect_ratio = 1
R0 = 1
rho = 1 
rhoc = 4.51 / π
init_theta = "pair"
init_pos = "random"
q = 1.0
r0 = round(Int,sqrt(Ntarget)/rho/2)
phonons = false ; phonon_amplitude = 1 ; phonon_k = 1  ; phonon_omega = 0 
params_phonons = Dict(:phonons => phonons, :phonon_amplitude => phonon_amplitude, :phonon_k => phonon_k, :phonon_omega => phonon_omega)
params_init = Dict(:init_pos => NaN, :init_theta => init_theta, :r0 => NaN, :q => q)

R_per_core = 100

tmax = 1000
times = 0:5:tmax # linear time

# sigmas = [0,0.1]
sigmas = collect(0:0.05:0.3)

# Ts = [0,0.1,0.2,0.3,0.4]
Ts = [0.1]

# v0s = [0.5,0.75,1,1.5,2,3]
v0s = [2]


dfts = Array{DefectTracker}(undef,length(v0s),length(sigmas),length(Ts),R_per_core)

z = @elapsed for i in each(v0s), j in each(sigmas), k in each(Ts), r in 1:R_per_core
    v0 = v0s[i]
    sigma = sigmas[j]
	T = Ts[k]

    println("v0 = $v0, σ = $sigma, T = $T")
    N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
    
    params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
    param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
        :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
        :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init, :params_phonons => params_phonons)

    t = 0.0
    system = System(param)
    dft = DefectTracker(system, t)
    dft, system = track!(dft,system,times,verbose=true)
	dfts[i,j,k,r] = dft
end
prinz(z)

filename = "data/DFT_Rt*_sigmas_r$real.jld2"
JLD2.@save filename Ts dfts R_per_core params_init Ntarget q init_theta sigmas aspect_ratio times tmax comments rhoc runtime = z


# ## ---------------- Tracking a pair of defects for immobile particles ---------------- ##
# comments = "From the defect data one can infer the MSD and diffusion coeff of an individual defect. "
# # Physical Params 
# Ntarget = Int(1E4)
# aspect_ratio = 1
# Ts = [0.1,0.2,0.3,0.4]
# Ts = [0.4]
# R0 = 1
# rho = 1 
# rhoc = 4.51 / π
# v0 = 0 
# sigma = 0
# q = 1.0
# params_init = Dict(:init_pos => NaN, :init_theta => NaN, :r0 => NaN, :q => q)
# tmax = 1E2
# times = 0:5:tmax # linear time

# inits_pos = ["square","rsa","random","pds"]
# R0s = [2,2,1.95]
# init_theta = "single"
# qs = [+1,-1]
# dfts = Array{DefectTracker}(undef, length(inits_pos),length(qs),length(Ts))

# z = @elapsed for i in each(inits_pos), j in each(qs), k in each(Ts)
#     R0 = R0s[i]
#     init_pos = inits_pos[i]
# 	q = qs[j]
# 	T = Ts[k]

#     println("$init_pos, $init_theta R0 = $R0, q = $q, T = $T")
#     N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
    
#     params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
#     param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
#         :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
#         :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

#     t = 0.0
#     system = System(param)
#     dft = DefectTracker(system, t)
#     dft, system = track!(dft,system,times,verbose=false)
# 	dfts[i,j,k] = dft
# end
# prinz(z)

# filename = "data/immobile_DFT_single_r$real.jld2"
# JLD2.@save filename qs R0s Ts inits_pos dfts params_init Ntarget v0 sigma aspect_ratio times tmax comments rhoc runtime = z

# ## ---------------- Tracking a pair of defects for immobile particles ---------------- ##
# comments = "From the defects data, one will be able to infer : \n
# A. the separating distance between the two defects R(t) \n
# B. the MSD and diffusion coeff of an individual defect. "
# # Physical Params 
# Ntarget = Int(1E4)
# aspect_ratio = 1
# Ts = [0.1,0.2,0.4]
# R0 = 1
# rho = 1 
# rhoc = 4.51 / π
# v0 = 0 
# sigma = 0
# q = 1.0
# params_init = Dict(:init_pos => NaN, :init_theta => NaN, :r0 => NaN, :q => q)
# tmax = 1E3
# times = 0:5:tmax # linear time

# inits_pos = ["square","rsa","random"]
# R0s = [2,2,1.95]
# init_theta = "pair"
# r0s = [10,25,50]
# dfts = Array{DefectTracker}(undef, length(inits_pos),length(r0s), length(Ts))

# z = @elapsed for i in each(inits_pos), j in each(r0s), k in each(Ts)
#     R0 = R0s[i]
#     init_pos = inits_pos[i]
# 	r0 = r0s[j]
# 	T = Ts[k]

#     println("$init_pos, $init_theta R0 = $R0, r0 = $r0, T = $T")
#     N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
    
#     params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
#     param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
#         :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
#         :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

#     t = 0.0
#     system = System(param)
#     dft = DefectTracker(system, t)
#     dft, system = track!(dft,system,times,verbose=true)
# 	dfts[i,j,k] = dft
# end
# prinz(z)


# filename = "data/immobile_DFT_pair_r$real.jld2"
# JLD2.@save filename r0s R0s Ts inits_pos dfts params_init Ntarget v0 sigma aspect_ratio times tmax comments rhoc runtime = z


# ## ---------------- Tracking a pair of defects for mobile particles ---------------- ##
# comments = "From the defects data, one will be able to infer : \n
# A. the separating distance between the two defects R(t) \n
# B. the MSD and diffusion coeff of an individual defect. "
# # Physical Params 
# Ntarget = Int(1E4)
# aspect_ratio = 1
# T = 0.1
# R0 = 1 
# rhoc = 4.51 / π
# v0 = 0 
# sigma = 0

# # Initialisation parameters
# rhos = [1,2]
# init_pos = ["random","rsa"]
# init_theta = "pair"
# r0 = 20.0
# q = 1.0
# params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => NaN, :q => q)

# # Simulation parameters
# r0s = 30#5:10:35
# tmax = 1E2
# times = 0:5:tmax # linear time
# params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => NaN, :q => q)

# dfts = Array{DefectTracker}(undef, length(v0sigs), length(r0s))

# z = @elapsed for i in each(v0sigs), j in each(r0s)
#     r0 = r0s[j]
    
#     println("v0 = $v0, σ = $sigma, r0 = $r0")
#     N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
#     dt = determine_dt(T, sigma, v0, N, rho)

    
#     params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
#     param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
#         :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
#         :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

#     t = 0.0
#     system = System(param)
#     dft = DefectTracker(system, t)
#     dft, system = track!(dft,system,times)
# 	dfts[i,j] = dft
# end
# prinz(z)

# filename = "data/immobile_DFT_pair_r$real.jld2"
# JLD2.@save filename Ntarget v0 sigma rhos inits_pos params_init T dfts aspect_ratio times tmax comments rhoc runtime = z


# # ## ---------------- Impact of init on XY Model ---------------- ##
# comments = "Investigates the impact of initialisation for the spatial location of the 
# spins in the XY model. "
# # Physical Params 
# Ntarget = Int(1E4)
# aspect_ratio = 1
# Ts = [0.1,0.2,0.4]
# rho = 1
# v0 = 0
# sigma = 0
# R0c = sqrt(4.51 / π)

# # Initialisation parameters
# inits_pos = ["random", "square", "RSA"]
# R0s = [1.95,2,2]
# init_theta = "hightemp"
# r0 = 20.0
# q = 1.0
# params_init = Dict(:init_pos => NaN, :init_theta => init_theta, :r0 => r0, :q => q)

# # Simulation parameters
# tmax = 1E2
# times = logspace(1,tmax,10)

# P = zeros(length(inits_pos),length(Ts), length(times))
# C = Array{Vector{Float64}}(undef, length(inits_pos),length(Ts), length(times))
# xi = zeros(length(inits_pos),length(Ts), length(times))
# n = zeros(length(inits_pos),length(Ts), length(times))
# E = zeros(length(inits_pos),length(Ts), length(times))

# # Impact on the usual quantities
# z = @elapsed for i in each(inits_pos), j in each(Ts)
# 	init_pos = inits_pos[i]
# 	R0 = R0s[i]
# 	T = Ts[j]

#     println("Init : $(init_pos), $init_theta")
#     N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
#     params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
#     param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
#         :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
#         :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

#     t = 0.0
#     system = System(param)

# 	for tt in eachindex(times)
# 		evolve!(system, times[tt]) # evolves the systems up to times[tt]
		
# 		P[i,j,tt]  = polarOP(system)[1]
# 		corr_tmp = corr(system)
# 		C[i,j,tt]  = corr_tmp
# 		xi[i,j,tt] = corr_length(corr_tmp)
# 		n[i,j,tt]  = number_defects(system)
# 		E[i,j,tt]  = energy(system)
# 	end
# end
# prinz(z)

# filename = "data/impact_init_XY_r$real.jld2"
# JLD2.@save filename inits_pos R0s Ts P C n xi E Ntarget v0 sigma rho params_init aspect_ratio times tmax comments runtime = z

# ## ---------------- Tracking a pair of defects for mobile particles ---------------- ##
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
# v0sigs = [(1,0),(1,0.1)]
# r0s = 5#5:10:35
# tmax = 1E2
# times = 0:5:tmax # linear time
# params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => NaN, :q => q)

# dfts = Array{DefectTracker}(undef, length(v0sigs), length(r0s))

# z = @elapsed for i in each(v0sigs), j in each(r0s)
#     v0, sigma = v0sigs[i]
#     r0 = r0s[j]
    
#     println("v0 = $v0, σ = $sigma, r0 = $r0")
#     N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
#     dt = determine_dt(T, sigma, v0, N, rho)
    
    
#     params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
#     param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
#         :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
#         :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

#     t = 0.0
#     system = System(param)
#     dft = DefectTracker(system, t)
#     dft, system = track!(dft,system,times)
# 	dfts[i,j] = dft
# end
# prinz(z)

# filename = "data/DFT_pair_r$real.jld2"
# JLD2.@save filename Ntarget v0sigs rho params_init T dfts aspect_ratio times tmax comments rhoc runtime = z


# ## ---------------- Nature of the Phase Transition ---------------- ##
# comments = "The goal of this script is to pass through the transition line, 
# in both direction (keeping σ or v0 constant) and to compute correlation functions.
# Here for T = 0.1 and ρ = 1."
# # Physical Params 
# Ntarget = Int(1E4)
# aspect_ratio = 1
# T = 0.1
# R0 = 1
# rho = 1
# rhoc = 4.51 / π

# # Initialisation parameters
# init_pos = "random"
# init_theta = "hightemp"
# r0 = 20.0
# q = 1.0
# params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)

# # Simulation parameters
# v0sigs = [(v,0.1) for v in logspace(0.03,1,10,digits=3)]
# # v0sigs = vcat([(0.2,sigm) for sigm in 0:0.025:0.225],[(v,0.1) for v in logspace(0.03,1,10,digits=3)])
# tmax = 3E1
# # times = collect(0:tmax/30:tmax) # linear time
# times = logspace(1,tmax,30,digits=1) # log time

# P = zeros(length(v0sigs), length(times))
# C = Array{Vector{Float64}}(undef, length(v0sigs), length(times))
# xi = zeros(length(v0sigs), length(times))
# n = zeros(length(v0sigs), length(times))

# z = @elapsed for i in each(v0sigs)
#     v0, sigma = v0sigs[i]
#     println("v0 = $v0, σ = $sigma, $(100i/length(v0sigs))%")
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
        
#         P[i,tt]  = polarOP(system)[1]
#         corr_tmp = corr(system)
#         C[i,tt]  = corr_tmp
#         xi[i,tt] = corr_length(corr_tmp)
#         n[i,tt]  = number_defects(system)
#     end

# end
# prinz(z)

# filename = "data/nature_phase_transition_horizontal _r$real.jld2"
# JLD2.@save filename Ntarget v0sigs rho params_init T P C n xi aspect_ratio times tmax comments rhoc runtime = z


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
# Ntargets = Int.(logspace(50,4E4,14,digits=0))

# tmax = 4E4
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
#         evolve!(system, times[tt]) # evolves the systems up to times[tt]
        
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

