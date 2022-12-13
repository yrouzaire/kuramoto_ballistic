using JLD2,Parameters

base_filename = "data/defects_motion_IN_VITRO_v0s_sigmas_N1E4" # look up in main_server.jl
R = 40 # look up into bash_loog.sh
indices = [] ; for r in 1:R  if isfile(base_filename*"_r$r.jld2") push!(indices,r) end end
println("There are $(length(indices))/$R files.")


## Defects Motion
@load base_filename*"_r$(indices[1]).jld2" dfts T N rho v_sigmas times transients tmax runtime comments
runtimes = NaN*zeros(R)
dfts_fusion_undef = Vector{DefectTracker}(undef,R)
dfts_fusion = DefectTracker[]
for r in indices
    println("r = $r")
    @load base_filename*"_r$r.jld2" dfts runtime
    dfts_fusion_undef[r] = dfts
    push!(dfts_fusion,dfts)

    runtimes[r] = runtime
end

@save base_filename*".jld2" runtimes Ts Ns v0s rhos sigmas times_log tmax comments R pos_saveds thetas_saveds psis_saveds omegas_saveds
println("Fusionned data saved in $(base_filename*".jld2") .")

## Phase Diagram
# @load base_filename*"_r$(indices[1]).jld2" Ts inits Ns v0s rhos sigmas times_log tmax comments
#
# P_fusion = NaN*zeros(length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),length(times_log),R)
# n_fusion = NaN*zeros(length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),length(times_log),R)
# # C_fusion = Array{Vector{Float32},8}(undef,length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),length(times_log),R)
# C_fusion = Vector{Array{Vector{Float32},7}}(undef,R)
#
# pos_saveds  = zeros(Float16,2,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),R)
# thetas_saveds = zeros(Float16,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),R)
# psis_saveds   = zeros(Float16,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),R)
# omegas_saveds = zeros(Float16,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),R)
# runtimes = NaN*zeros(R)
#
# for r in indices
#     println("r = $r")
#     @load base_filename*"_r$r.jld2" P C n pos_saved thetas_saved psis_saved omegas_saved runtime
#     P_fusion[:,:,:,:,:,:,:,r] = P
#     n_fusion[:,:,:,:,:,:,:,r] = n
#     C_fusion[r] = C
#
#     pos_saveds[:,:,:,:,:,:,:,r]  = pos_saved
#     thetas_saveds[:,:,:,:,:,:,r] = thetas_saved
#     psis_saveds[:,:,:,:,:,:,r]  = psis_saved
#     omegas_saveds[:,:,:,:,:,:,r] = omegas_saved
#
#     runtimes[r] = runtime
# end
#
# @save base_filename*".jld2" Ps=P_fusion ns=n_fusion Cs=C_fusion runtimes Ts Ns v0s rhos sigmas times_log tmax comments R pos_saveds thetas_saveds psis_saveds omegas_saveds
# println("Fusionned data saved in $(base_filename*".jld2") .")

# ## FSS
# @load base_filename*"_r$(indices[1]).jld2" Ns rho seuil_break seuil_crit times_break v0s T sigmas tmax
#
# vc_fusion = NaN*zeros(length(Ns),length(sigmas),R)
#
# for r in indices
#     println("r = $r")
#     @load base_filename*"_r$r.jld2" vc
#     vc_fusion[:,:,r]  = vc
# end
# @save base_filename*".jld2" Ns rho seuil_break seuil_crit times_break v0s T sigmas tmax vcs=vc_fusion
# println("Fusionned data saved in $(base_filename*".jld2") .")
