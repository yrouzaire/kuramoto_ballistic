using JLD2,Parameters

base_filename = "data/bigscan_v0_sigma_N1E3_rho1_tmax1E4" # look up in main_server.jl
R = 40 # look up into bash_loog.sh
indices = [] ; for r in 1:R  if isfile(base_filename*"_r$r.jld2") push!(indices,r) end end
println("There are $(length(indices))/$R files.")

## Phase Diagram
@load base_filename*"_r$(indices[1]).jld2" Ts Ns v0s rhos sigmas times_log tmax comments

P_fusion = NaN*zeros(length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(times_log),R)
n_fusion = NaN*zeros(length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(times_log),R)
C_fusion = NaN*zeros(length(0:0.5:round(Int,sqrt(Ns)/2)),length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(times_log),R)
runtimes = NaN*zeros(R)

for r in indices
    println("r = $r")
    @load base_filename*"_r$r.jld2" P n C runtime
    P_fusion[:,:,:,:,:,:,r] = P
    n_fusion[:,:,:,:,:,:,r] = n
    C_fusion[:,:,:,:,:,:,:,r] = C
    runtimes[r] = runtime
end

@save base_filename*".jld2" Ps=P_fusion ns=n_fusion Cs=C_fusion runtimes Ts Ns v0s rhos sigmas times_log tmax comments R
println("Fusionned data saved in $(base_filename*".jld2") .")
