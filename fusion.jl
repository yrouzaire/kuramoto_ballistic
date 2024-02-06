using JLD2, Parameters
include("methods.jl");

base_filename = "data/FSS_to_determine_transition" # look up in main_server.jl
R = 40 # look up into bash_loog.sh
indices = [];
for r in 1:R
    if isfile(base_filename * "_r$r.jld2")
        push!(indices, r)
    end
end;
println("There are $(length(indices))/$R files.")

rhoc = 4.51 / pi


## Impact σ on R(t*) 
# @load base_filename*"_r$(indices[1]).jld2" Ts dfts R_per_core params_init Ntarget q init_theta sigmas aspect_ratio times tmax comments rhoc runtime
# runtimes = NaN*zeros(R)
# dfts_fusion = Array{DefectTracker}[] 
# for r in indices
#     println("r = $r")
#     @load base_filename*"_r$r.jld2" dfts runtime
#     push!(dfts_fusion,dfts)

#     runtimes[r] = runtime
# end

# @save base_filename*".jld2" R Ts dfts_fusion R_per_core params_init Ntarget q init_theta sigmas aspect_ratio times tmax comments rhoc runtimes
# println("Fusionned data saved in $(base_filename*".jld2") .")


# ## ---------------- MSD Specifically Tracking a Pair of defects  ---------------- ##
# @load base_filename*"_r$(indices[1]).jld2" sigmas v0s Ts xy_pos xy_neg rr times_collision R_per_core params_init Ntarget R0 q init_theta init_pos aspect_ratio times tmax comments rhoc runtime
# runtimes = NaN*zeros(R)
# Rtot = R_per_core*length(indices)
# all_xy_pos = Array{Vector{Tuple{Number,Number}}}(undef,length(v0s),length(sigmas),length(Ts),Rtot)
# all_xy_neg = Array{Vector{Tuple{Number,Number}}}(undef,length(v0s),length(sigmas),length(Ts),Rtot)
# all_rr = Array{Vector{Number}}(undef,length(v0s),length(sigmas),length(Ts),Rtot)
# all_times_collision = zeros(length(v0s),length(sigmas),length(Ts),Rtot)

# r_ind = 0
# for r in indices
#     global r_ind += 1
#     println("r = $r")
#     @load base_filename*"_r$r.jld2" xy_pos xy_neg rr times_collision runtime
#     all_xy_pos[:,:,:,R_per_core*(r_ind-1)+1:R_per_core*r_ind] = xy_pos
#     all_xy_neg[:,:,:,R_per_core*(r_ind-1)+1:R_per_core*r_ind] = xy_neg
#     all_rr[:,:,:,R_per_core*(r_ind-1)+1:R_per_core*r_ind] = rr
#     all_times_collision[:,:,:,R_per_core*(r_ind-1)+1:R_per_core*r_ind] = times_collision
#     runtimes[r] = runtime
# end

# @save base_filename*".jld2" sigmas v0s Ts all_xy_pos all_xy_neg all_rr all_times_collision R_per_core Rtot params_init Ntarget R0 q init_theta init_pos aspect_ratio times tmax comments rhoc runtimes
# println("Fusionned data saved in $(base_filename*".jld2") .")

## Proba Spinwaves impact v0 and σ
# @load base_filename*"_r$(indices[1]).jld2" nb_detected_spinwave systems_detected_spinwave times_detected_spinwave thetas_detected_spinwave pos_detected_spinwave Ps_detected_spinwave R_per_core sigmas v0s tmax times p_threshold init_pos init_theta Ntarget rho T aspect_ratio runtime
# runtimes = NaN*zeros(R)
# all_nb_detected_spinwave = zeros(length(v0s),length(sigmas))
# all_times_detected_spinwave = [[] for i in each(v0s), j in each(sigmas)]
# all_Ps_detected_spinwave = [[] for i in each(v0s), j in each(sigmas)]
# all_thetas_detected_spinwave = [Vector{Float16}[] for i in each(v0s), j in each(sigmas)]
# all_pos_detected_spinwave = [Vector{Tuple{Number,Number}}[] for i in each(v0s), j in each(sigmas)]
# all_systems_detected_spinwave = [System[] for i in each(v0s), j in each(sigmas)]
# for r in indices
#     println("r = $r")
#     @load base_filename*"_r$r.jld2" nb_detected_spinwave runtime times_detected_spinwave thetas_detected_spinwave pos_detected_spinwave Ps_detected_spinwave systems_detected_spinwave
#     global all_nb_detected_spinwave += nb_detected_spinwave
#     for i in each(v0s), j in each(sigmas)
#         n = length(times_detected_spinwave[i,j])
#         for nn in 1:n
#             push!(all_times_detected_spinwave[i,j], times_detected_spinwave[i,j][nn])
#             push!(all_Ps_detected_spinwave[i,j], Ps_detected_spinwave[i,j][nn])
#             push!(all_thetas_detected_spinwave[i,j], thetas_detected_spinwave[i,j][nn])
#             push!(all_pos_detected_spinwave[i,j], pos_detected_spinwave[i,j][nn])
#             # push!(all_systems_detected_spinwave[i,j], systems_detected_spinwave[i,j][nn])
#         end
#     end
#     runtimes[r] = runtime
# end
# Rtot = R_per_core*length(indices)

# @save base_filename * ".jld2" R_per_core Rtot R all_nb_detected_spinwave all_times_detected_spinwave all_Ps_detected_spinwave all_thetas_detected_spinwave all_pos_detected_spinwave sigmas v0s tmax times p_threshold init_pos init_theta Ntarget rho T aspect_ratio runtimes
# # @save base_filename * ".jld2" R_per_core Rtot R all_systems_detected_spinwave all_nb_detected_spinwave all_times_detected_spinwave all_Ps_detected_spinwave all_thetas_detected_spinwave all_pos_detected_spinwave sigmas v0s tmax times p_threshold init_pos init_theta Ntarget rho T aspect_ratio runtimes
# println("Fusionned data saved in $(base_filename*".jld2") .")


# ## Impact R0 on MSD (Square Lattice)
# @load base_filename*"_r$(indices[1]).jld2" R0s Ts inits_pos dfts R_per_core params_init Ntarget v0 q sigma aspect_ratio times tmax comments rhoc runtime
# runtimes = NaN*zeros(R)
# dfts_fusion = Array{DefectTracker}[] 
# for r in indices
#     println("r = $r")
#     @load base_filename*"_r$r.jld2" dfts runtime
#     push!(dfts_fusion,dfts)

#     runtimes[r] = runtime
# end

# @save base_filename*".jld2" R R0s Ts inits_pos dfts_fusion R_per_core params_init Ntarget v0 q sigma aspect_ratio times tmax comments rhoc runtimes
# println("Fusionned data saved in $(base_filename*".jld2") .")

# ## Single Defects  Motion (immobile particles)
# @load base_filename*"_r$(indices[1]).jld2" qs R0s Ts inits_pos dfts params_init Ntarget v0 sigma aspect_ratio times tmax comments rhoc runtime
# runtimes = NaN*zeros(R)
# # dfts_fusion_undef = Vector{Array{DefectTracker}}(undef,R) 
# dfts_fusion = Array{DefectTracker}[] # easier analysis if ≥ 1 simulation is missing
# for r in indices
#     println("r = $r")
#     @load base_filename*"_r$r.jld2" dfts runtime
#     # dfts_fusion_undef[r] = dfts
#     push!(dfts_fusion,dfts)

#     runtimes[r] = runtime
# end

# @save base_filename*".jld2" R qs R0s Ts inits_pos dfts_fusion params_init Ntarget v0 sigma aspect_ratio times tmax comments rhoc runtimes
# println("Fusionned data saved in $(base_filename*".jld2") .")

# ## Pair Defects  Motion (immobile particles)
# @load base_filename*"_r$(indices[1]).jld2" r0s R0s Ts inits_pos dfts params_init Ntarget v0 sigma aspect_ratio times tmax comments rhoc runtime
# runtimes = NaN*zeros(R)
# # dfts_fusion_undef = Vector{Array{DefectTracker}}(undef,R) 
# dfts_fusion = Array{DefectTracker}[] # easier analysis if ≥ 1 simulation is missing
# for r in indices
#     println("r = $r")
#     @load base_filename*"_r$r.jld2" dfts runtime
#     # dfts_fusion_undef[r] = dfts
#     push!(dfts_fusion,dfts)

#     runtimes[r] = runtime
# end

# @save base_filename*".jld2" R r0s R0s Ts inits_pos dfts_fusion params_init Ntarget v0 sigma aspect_ratio times tmax comments rhoc runtimes
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
# @load base_filename * "_r$(indices[1]).jld2" inits_pos R0s Ts P C n xi E Ntarget v0 sigma rho params_init aspect_ratio times tmax comments runtime
# Ps = zeros(length(inits_pos),length(Ts),length(times),R)
# Cs = Array{Vector{Float64}}(undef, length(inits_pos),length(Ts),length(times),R)
# ns = zeros(length(inits_pos),length(Ts),length(times),R)
# xis = zeros(length(inits_pos),length(Ts),length(times),R)
# Es = zeros(length(inits_pos),length(Ts),length(times),R)
# runtimes = NaN*zeros(R)

# for r in indices
#     println("r = $r")
#     @load base_filename * "_r$r.jld2" runtime P C n xi E
#     Ps[:,:,:,r] = P
#     Cs[:,:,:,r] = C
#     ns[:,:,:,r] = n
#     xis[:,:,:,r] = xi
#     Es[:,:,:,r] = E

#     runtimes[r] = runtime
# end

# @save base_filename * ".jld2" R inits_pos R0s Ts Ps Cs ns xis Es Ntarget v0 sigma rho params_init aspect_ratio times tmax comments runtimes
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


# ## FSS
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

# Critical sigmas
# @load base_filename * "_r$(indices[1]).jld2" Ntarget aspect_ratio rhos sigmas v0s times tmax critical_sigmas T seuil comments rhoc runtime
# critical_sigmas_fusion = NaN * ones(length(rhos), length(v0s), R)
# runtimes = NaN * zeros(R)
# for r in indices
#     println("r = $r")
#     @load base_filename * "_r$r.jld2" runtime critical_sigmas
#     critical_sigmas_fusion[:, :, r] = critical_sigmas
#     runtimes[r] = runtime
# end
# @save base_filename * ".jld2" Ntarget aspect_ratio rhos sigmas times tmax critical_sigmas_fusion T seuil comments rhoc runtimes R
# println("Fusionned data saved in $(base_filename*".jld2") .")

# ## Critical velocity
# @load base_filename * "_r$(indices[1]).jld2" Ntarget rhos times tmax T v0s seuil rhoc runtime comments
# critical_velocity_fusion = NaN*ones(length(rhos), R)
# runtimes = NaN * zeros(R)
# for r in indices
#     println("r = $r")
#     @load base_filename * "_r$r.jld2" critical_velocity runtime
#     critical_velocity_fusion[:, r] = critical_velocity
#     runtimes[r] = runtime
# end
# @save base_filename * ".jld2" Ntarget rhos times tmax critical_velocity_fusion T v0s seuil rhoc runtimes comments R Reff=length(indices)
# println("Fusionned data saved in $(base_filename*".jld2") .")


