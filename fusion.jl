using JLD2, Parameters
include("methods.jl");

base_filename = "data/immobile_DFT_pair" # look up in main_server.jl
R = 40 # look up into bash_loog.sh
indices = [];
for r in 1:R
    if isfile(base_filename * "_r$r.jld2")
        push!(indices, r)
    end
end;
println("There are $(length(indices))/$R files.")

# # ## Defects Motion (immobile particles)
# @load base_filename*"_r$(indices[1]).jld2" r0s R0s inits_pos inits_thetas dfts params_init Ntarget v0 sigma T aspect_ratio times tmax comments runtime
# runtimes = NaN*zeros(R)
# dfts_fusion_undef = Vector{Array{DefectTracker}}(undef,R)
# dfts_fusion = Array{DefectTracker}[]
# for r in indices
#     println("r = $r")
#     @load base_filename*"_r$r.jld2" dfts runtime
#     dfts_fusion_undef[r] = dfts
#     push!(dfts_fusion,dfts)

#     runtimes[r] = runtime
# end

# @save base_filename*".jld2" R r0s R0s inits_pos inits_thetas dfts_fusion_undef dfts_fusion params_init Ntarget v0 sigma T aspect_ratio times tmax comments runtimes
# println("Fusionned data saved in $(base_filename*".jld2") .")


# ## Defects Motion (mobile particles)
# @load base_filename*"_r$(indices[1]).jld2" dfts T Ntarget v0sigs rho params_init T aspect_ratio times tmax comments rhoc runtime
# runtimes = NaN*zeros(R)
# dfts_fusion_undef = Vector{Array{DefectTracker}}(undef,R)
# dfts_fusion = Array{DefectTracker}[]
# for r in indices
#     println("r = $r")
#     @load base_filename*"_r$r.jld2" dfts runtime
#     dfts_fusion_undef[r] = dfts
#     push!(dfts_fusion,dfts)

#     runtimes[r] = runtime
# end

# @save base_filename*".jld2" dfts_fusion dfts_fusion_undef runtimes T Ntarget v0sigs rho times params_init aspect_ratio tmax comments R
# println("Fusionned data saved in $(base_filename*".jld2") .")


# ## Impact Initialisation on XY
# @load base_filename * "_r$(indices[1]).jld2" Ntarget rho T inits_pos params_init v0 sigma P C n xi aspect_ratio times tmax comments rhoc runtime
# Ps = zeros(length(inits_pos),length(times),R)
# Cs = Array{Vector{Float64}}(undef, length(inits_pos),length(times),R)
# ns = zeros(length(inits_pos),length(times),R)
# xis = zeros(length(inits_pos),length(times),R)
# runtimes = NaN*zeros(R)

# for r in indices
#     println("r = $r")
#     @load base_filename * "_r$r.jld2" runtime P C n xi 
#     Ps[:,:,r] = P
#     Cs[:,:,r] = C
#     ns[:,:,r] = n
#     xis[:,:,r] = xi

#     runtimes[r] = runtime
# end

# @save base_filename * ".jld2" R runtimes inits_pos Ps Cs ns xis rho T v0 sigma Ntarget params_init aspect_ratio times tmax comments rhoc 
# println("Fusionned data saved in $(base_filename*".jld2") .")


## Nature of the phase transition
# @load base_filename * "_r$(indices[1]).jld2" Ntarget rho T params_init v0sigs P C n xi aspect_ratio times tmax comments rhoc runtime
# Ps = zeros(length(v0sigs),length(times),R)
# Cs = Array{Vector{Float64}}(undef, length(v0sigs),length(times),R)
# ns = zeros(length(v0sigs),length(times),R)
# xis = zeros(length(v0sigs),length(times),R)
# runtimes = NaN*zeros(R)

# for r in indices
#     println("r = $r")
#     @load base_filename * "_r$r.jld2" runtime P C n xi 
#     Ps[:,:,r] = P
#     Cs[:,:,r] = C
#     ns[:,:,r] = n
#     xis[:,:,r] = xi

#     runtimes[r] = runtime
# end

# @save base_filename * ".jld2" v0sigs Ps Cs ns xis rho T Ntarget params_init aspect_ratio times tmax comments rhoc runtimes R
# println("Fusionned data saved in $(base_filename*".jld2") .")

## Hysteresis and Nature of Phases
# @load base_filename * "_r$(indices[1]).jld2" Ntarget rhos Ts inits_theta v0sigs P C n xi aspect_ratio times tmax comments rhoc runtime
# Ps = zeros(length(v0sigs), length(Ts), length(rhos),length(inits_theta), length(times),R)
# Cs = Array{Vector{Float64}}(undef, length(v0sigs), length(Ts), length(rhos),length(inits_theta), length(times),R)
# ns = zeros(length(v0sigs), length(Ts), length(rhos),length(inits_theta), length(times),R)
# xis = zeros(length(v0sigs), length(Ts), length(rhos),length(inits_theta), length(times),R)
# runtimes = NaN*zeros(R)

# for r in indices
#     println("r = $r")
#     @load base_filename * "_r$r.jld2" runtime P C n xi 
#     Ps[:,:,:,:, :,r] = P
#     Cs[:,:,:,:, :,r] = C
#     ns[:,:,:,:, :,r] = n
#     xis[:,:,:,:, :,r] = xi

#     runtimes[r] = runtime
# end

# @save base_filename * ".jld2" Ntarget rhos Ts inits_thetas v0sigs Ps Cs ns xis params_init aspect_ratio times tmax comments rhoc runtimes R
# println("Fusionned data saved in $(base_filename*".jld2") .")


## FSS
@load base_filename * "_r$(indices[1]).jld2" Ntargets v0sigs P C n xi params_init aspect_ratio rho times tmax T comments rhoc runtime
Ps = zeros(length(v0sigs), length(Ntargets), length(times),R)
ns = zeros(length(v0sigs), length(Ntargets), length(times),R)
xis = zeros(length(v0sigs), length(Ntargets), length(times),R)
Cs = Array{Vector{Float64}}(undef, length(v0sigs), length(Ntargets), length(times),R)
runtimes = NaN*zeros(R)

for r in indices
    println("r = $r")
    @load base_filename * "_r$r.jld2" runtime P C n xi 
    Ps[:,:,:,r] = P
    Cs[:,:,:,r] = C
    ns[:,:,:,r] = n
    xis[:,:,:,r] = xi

    runtimes[r] = runtime
end
@save base_filename * ".jld2" Ntargets v0sigs Ps Cs ns xis params_init aspect_ratio times tmax T comments rho rhoc runtimes R
println("Fusionned data saved in $(base_filename*".jld2") .")

## Phase Diagram
# @load base_filename*"_r$(indices[1]).jld2" Ts inits Ns v0s rhos sigmas times_log tmax comments

# P_fusion = NaN*zeros(length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),length(times_log),R)
# n_fusion = NaN*zeros(length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),length(times_log),R)
# # C_fusion = Array{Vector{Float32},8}(undef,length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),length(times_log),R)
# C_fusion = Vector{Array{Vector{Float32},7}}(undef,R)

# pos_saveds  = zeros(Float16,2,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),R)
# thetas_saveds = zeros(Float16,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),R)
# psis_saveds   = zeros(Float16,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),R)
# omegas_saveds = zeros(Float16,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),R)
# runtimes = NaN*zeros(R)

# for r in indices
#     println("r = $r")
#     @load base_filename*"_r$r.jld2" P C n pos_saved thetas_saved psis_saved omegas_saved runtime
#     P_fusion[:,:,:,:,:,:,:,r] = P
#     n_fusion[:,:,:,:,:,:,:,r] = n
#     C_fusion[r] = C

#     pos_saveds[:,:,:,:,:,:,:,r]  = pos_saved
#     thetas_saveds[:,:,:,:,:,:,r] = thetas_saved
#     psis_saveds[:,:,:,:,:,:,r]  = psis_saved
#     omegas_saveds[:,:,:,:,:,:,r] = omegas_saved

#     runtimes[r] = runtime
# end

# @save base_filename*".jld2" Ps=P_fusion ns=n_fusion Cs=C_fusion runtimes Ts Ns v0s rhos sigmas times_log tmax comments R pos_saveds thetas_saveds psis_saveds omegas_saveds
# println("Fusionned data saved in $(base_filename*".jld2") .")

# ## FSS
# @load base_filename * "_r$(indices[1]).jld2" Ns rhos times v_sigmas T Ps ns xis

# Ps_fusion = Array{Vector{Float64}}(undef, length(Ns), length(rhos), length(v_sigmas), R)
# ns_fusion = Array{Vector{Float64}}(undef, length(Ns), length(rhos), length(v_sigmas), R)
# xis_fusion = Array{Vector{Float64}}(undef, length(Ns), length(rhos), length(v_sigmas), R)
# runtimes = NaN*zeros(R)
# for r in indices
#     println("r = $r")
#     @load base_filename * "_r$r.jld2" runtime xis ns Ps
#     Ps_fusion[:,:,:,r] = Ps
#     ns_fusion[:,:,:,r] = ns
#     xis_fusion[:,:,:,r] = xis
#     runtimes[r] = runtime
# end
# @save base_filename * ".jld2" Ns rhos times v_sigmas T Ps_fusion ns_fusion xis_fusion runtimes
# println("Fusionned data saved in $(base_filename*".jld2") .")

## Critical sigmas
# @load base_filename * "_r$(indices[1]).jld2" N rhos times tmax critical_sigmas T v0s sigmas seuil vc rhoc runtime comments
# critical_sigmas_fusion = NaN*ones(length(v0s), length(rhos), R)
# runtimes = NaN * zeros(R)
# for r in indices
#     println("r = $r")
#     @load base_filename * "_r$r.jld2" runtime critical_sigmas
#     critical_sigmas_fusion[:, :, r] = critical_sigmas
#     runtimes[r] = runtime
# end
# @save base_filename * ".jld2" N rhos times tmax critical_sigmas_fusion T v0s sigmas seuil vc rhoc runtimes comments R
# println("Fusionned data saved in $(base_filename*".jld2") .")

## Critical velocity
# @load base_filename * "_r$(indices[1]).jld2" N rhos times tmax T v0s seuil rhoc runtime comments
# critical_sigmas_fusion = NaN*ones(length(v0s), length(rhos), R)
# runtimes = NaN * zeros(R)
# for r in indices
#     println("r = $r")
#     @load base_filename * "_r$r.jld2" runtime critical_velocity
#     critical_velocity_fusion[:, :, r] = critical_velocity
#     runtimes[r] = runtime
# end
# @save base_filename * ".jld2" N rhos times tmax critical_velocity_fusion T v0s seuil rhoc runtimes comments R
# println("Fusionned data saved in $(base_filename*".jld2") .")
