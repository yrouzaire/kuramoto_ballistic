cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian, LambertW
include("../methods.jl")
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

## ---------------- Analysis DFT Pair Mobile Particles ---------------- ##
filename = "data/DFT_pair.jld2"
@load filename dfts_fusion dfts_fusion_undef runtimes T Ntarget v0sigs rho times params_init aspect_ratio tmax comments R
hrun(runtimes)
times = 0:5:tmax
r0s = 5:10:35
Reff = length(dfts_fusion)

T
Ntarget
rho
v0sigs

# R(t)
rts = zeros(length(v0sigs), length(r0s), length(times),Reff)
rts_reverse = NaN*zeros(length(v0sigs), length(r0s), length(times),Reff)
L = sqrt(Ntarget / aspect_ratio / rho)
for r in each(dfts_fusion), i in each(v0sigs), j in each(r0s), k in each(times)
    dft = dfts_fusion[r][i,j]
    tmp = interdefect_distance(dft.defectsP[1], dft.defectsN[1], L, L)
    ll = min(length(times),length(tmp))
    rts[i,j,1:ll,r] = tmp[1:ll]
    rts_reverse[i,j,1:ll,r] = reverse(tmp[1:ll])
end
rts_avg = mean(rts, dims=4)[:,:,:,1]
rts_reverse_avg = nanmean(rts_reverse,4)[:,:,:,1]

# MSD(t) & G(x,t) of individual defects
dx_P = [[] for i in each(v0sigs), j in each(r0s), k in each(times)]
dy_P = [[] for i in each(v0sigs), j in each(r0s), k in each(times)]
SD_P = [[] for i in each(v0sigs), j in each(r0s), k in each(times)]
dx_M = [[] for i in each(v0sigs), j in each(r0s), k in each(times)]
dy_M = [[] for i in each(v0sigs), j in each(r0s), k in each(times)]
SD_M = [[] for i in each(v0sigs), j in each(r0s), k in each(times)]


for r in each(dfts_fusion), i in each(v0sigs), j in each(r0s)
    dft = dfts_fusion[r][i,j]
    for tt in each(times)
        try
            d = dft.defectsP[1].pos[tt] .- dft.defectsP[1].pos[1]
            push!(dx_P[i,j,tt], d[1])
            push!(dy_P[i,j,tt], d[2])
            push!(SD_P[i,j,tt], dist2(dft.defectsP[1].pos[tt],dft.defectsP[1].pos[1],L,L))
        catch ;
        end

        try
            d = dft.defectsN[1].pos[tt] .- dft.defectsN[1].pos[1]
            push!(dx_M[i,j,tt], d[1])
            push!(dy_M[i,j,tt], d[2])
            push!(SD_M[i,j,tt], dist2(dft.defectsN[1].pos[tt],dft.defectsN[1].pos[1],L,L))
        catch ;
        end
    end
end

MSD_P = zeros(length(v0sigs), length(r0s), length(times))
for i in each(v0sigs), j in each(r0s), k in each(times)
    if !isempty(SD_P[i,j,k])
        MSD_P[i,j,k] = mean(SD_P[i,j,k])
    end
end

MSD_M = zeros(length(v0sigs), length(r0s), length(times))
for i in each(v0sigs), j in each(r0s), k in each(times)
    if !isempty(SD_M[i,j,k])
        MSD_M[i,j,k] = mean(SD_M[i,j,k])
    end
end
MSD_all = (MSD_P + MSD_M) / 2

## ---------------- Plots R(t) ---------------- ##
p1=plot(ylabel="R(t)", xlabel="t", legend=true,uaxis=:log)
for i in 2#each(v0sigs)
    for j in each(r0s)
        plot!(times[2:end], remove_negative(rts_avg[i,j,2:end]), label=L"r_0 = "*string(r0s[j]),rib=0)
        # plot!(times, (rts_avg[i,j,:]), label="v0 = $(v0sigs[i][1]), r0 = $(r0s[j])")
    end
end
p1

p2=plot(ylabel=L"R(t^*)", xlabel=L"t^*", legend=:bottomright,axis=:log)
for i in 2#each(v0sigs)
    for j in each(r0s)
        plot!(times[2:end], remove_negative(rts_reverse_avg[i,j,2:end]))
    end
end
using LambertW
plot!(times[2:end],x->exp(lambertw(2π/0.5*x)/2),c=:black,line=:dash,label="XY prediction")
plot!(times[2:end],x->3E-1x^0.5,c=:black,label=L"\sqrt{t^*}")
p2

plot(p1,p2,layout=(1,2),size=(800,400))
# savefig("figures/two_defects/Rt_v0_1_sigma_0.1.png")


## ---------------- Histogram Displacements ---------------- ##
ind_r = 4
tt = 100#length(times) 
histogram(abs.(vcat(dx_P[1,ind_r,tt],dy_P[1,ind_r,tt],dx_M[1,ind_r,tt],dy_M[1,ind_r,tt])), bins=30)

## ---------------- Plot MSD(t) ---------------- ##
p=plot(xlabel="t", ylabel="MSD", legend=:topleft, axis=:log)
for i in each(v0sigs), j in 4#each(r0s)
    plot!(times[2:end], remove_negative(MSD_all[i,j,2:end]), label=L"(v_0,σ) = "*string(v0sigs[i]),rib=0)
end
plot!(times[2:end],x->3E-1x^1,c=:black,line=:dash,label="MSD ∼ t")
p
# savefig("figures/two_defects/MSD_mobile.png")


## ---------------- Analysis DFT Pair Immobile Particles ---------------- ##
filename = "data/DFT_pair_immobile.jld2"
filename = "data/immobile_DFT_pair.jld2"
@load filename dfts_fusion dfts_fusion_undef inits_pos inits_thetas runtimes T Ntarget v0 sigma rhos times params_init aspect_ratio tmax comments R
hrun(runtimes)
times = 0:5:tmax
r0s = 5:10:35
Reff = length(dfts_fusion)
times 
rhos
inits_pos 
inits_thetas
r0s



# R(t) for inits_thetas = "pair"
rts = zeros(length(rhos), length(inits_pos), length(r0s), length(times),Reff)
rts_reverse = NaN*zeros(length(rhos), length(inits_pos), length(r0s), length(times),Reff)
for r in each(dfts_fusion), i in each(rhos), j in each(inits_pos), l in each(r0s), m in each(times)
    dft = dfts_fusion[r][i,j,1,l] # 1 = inits_thetas = "pair"
    L = sqrt(Ntarget / aspect_ratio / rhos[i])
    tmp = interdefect_distance(dft.defectsP[1], dft.defectsN[1], L, L)
    ll = min(length(times),length(tmp))
    rts[i,j,l,1:ll,r] = tmp[1:ll]
    rts_reverse[i,j,l,1:ll,r] = reverse(tmp[1:ll])
end
rts_avg = mean(rts, dims=5)[:,:,:,:,1]
rts_reverse_avg = nanmean(rts_reverse,5)[:,:,:,:,1]

# MSD(t) & G(x,t) of individual defects, for inits_thetas = "pair" and "single"
dx_P = [[] for i in each(rhos), j in each(inits_pos), k in each(inits_thetas), l in each(r0s), m in each(times)]
dy_P = [[] for i in each(rhos), j in each(inits_pos), k in each(inits_thetas), l in each(r0s), m in each(times)]
dx_M = [[] for i in each(rhos), j in each(inits_pos), k in each(inits_thetas), l in each(r0s), m in each(times)]
dy_M = [[] for i in each(rhos), j in each(inits_pos), k in each(inits_thetas), l in each(r0s), m in each(times)]
SD_P = [[] for i in each(rhos), j in each(inits_pos), k in each(inits_thetas), l in each(r0s), m in each(times)]
SD_M = [[] for i in each(rhos), j in each(inits_pos), k in each(inits_thetas), l in each(r0s), m in each(times)]

for r in each(dfts_fusion), i in each(rhos), j in each(inits_pos), k in each(inits_thetas), l in each(r0s)
    dft = dfts_fusion[r][i,j,k,l]
    L = sqrt(Ntarget / aspect_ratio / rhos[i])
    for tt in min(length(times),length(dft.defectsP[1].pos))
        d = dft.defectsP[1].pos[tt] .- dft.defectsP[1].pos[1]
        push!(dx_P[i,j,k,l,tt], d[1])
        push!(dy_P[i,j,k,l,tt], d[2])
        push!(SD_P[i,j,k,l,tt], dist2(dft.defectsP[1].pos[tt],dft.defectsP[1].pos[1],L,L))
    end

    for tt in min(length(times),length(dft.defectsN[1].pos))
        d = dft.defectsN[1].pos[tt] .- dft.defectsN[1].pos[1]
        push!(dx_M[i,j,k,l,tt], d[1])
        push!(dy_M[i,j,k,l,tt], d[2])
        push!(SD_M[i,j,k,l,tt], dist2(dft.defectsN[1].pos[tt],dft.defectsN[1].pos[1],L,L))
    end
end

MSD_P = zeros(length(rhos), length(inits_pos), length(inits_thetas), length(r0s), length(times))
for i in each(rhos), j in each(inits_pos), k in each(inits_thetas), l in each(r0s), m in each(times)
    if !isempty(SD_P[i,j,k,l,m])
        MSD_P[i,j,k,l,m] = mean(SD_P[i,j,k,l,m])
    end
end
MSD_M = zeros(length(rhos), length(inits_pos), length(inits_thetas), length(r0s), length(times))
for i in each(rhos), j in each(inits_pos), k in each(inits_thetas), l in each(r0s), m in each(times)
    if !isempty(SD_M[i,j,k,l,m])
        MSD_M[i,j,k,l,m] = mean(SD_M[i,j,k,l,m])
    end
end
MSD_all = (MSD_P + MSD_M) / 2

## ---------------- Plots R(t) ---------------- ##
p1=plot(ylabel="R(t)", xlabel="t", legend=false,uaxis=:log)
for i in [2] # rhos
    for j in [1,2,3] # inits_pos
        for l in [3,4] # 
            plot!(times[2:end], remove_negative(rts_avg[i,j,l,2:end]), label=L"r_0 = "*string(r0s[j]),rib=0)
        end
    end
end
p1

##
p2=plot(ylabel=L"R(t^*)", xlabel=L"t^*", legend=false,axis=:log)
for i in 2#each(rhos)
    for j in [1,2,3]
        for l in 3#each(r0s)
            plot!(times[2:end], remove_negative(rts_reverse_avg[i,j,l,2:end]), label=L"r_0 = "*string(r0s[j]),rib=0)
        end
    end
end
using LambertW
plot!(times[2:end],x->exp(lambertw(2π/1*x)/2),c=:black,line=:dash,label="XY prediction")
plot!(times[2:end],x->3E-1x^0.5,c=:black,label=L"\sqrt{t^*}")
p2

plot(p1,p2,layout=(1,2),size=(800,400))
# savefig("figures/two_defects/Rt_v0_1_sigma_0.1.png")


## ---------------- Histogram Displacements ---------------- ##
tt = length(times) ; histogram(vcat(dx_P[1,1,tt],dy_P[1,1,tt],dx_M[1,1,tt],dy_M[1,1,tt]), bins=50)

## ---------------- Plot MSD(t) ---------------- ##
p=plot(xlabel="t", ylabel="MSD", legend=false, axis=:log)
for i in [2] # rhos
    for j in [1,2,3] # inits_pos
        for k in 1 # inits_thetas
            for l in 4 # r0s
                plot!(times[2:end], remove_negative(MSD_P[i,j,k,l,2:end]), label=L"(v_0,σ) = "*string(v0sigs[i]),rib=0)
            end
        end
    end
end
plot!(times[2:end],x->3E-1x^1,c=:black,line=:dash,label="MSD ∼ t")
p
# savefig("figures/two_defects/MSD.png")


# ## ---------------- Monitor the defects in the green area  ---------------- ##
# comments = "From the defects data, one will be able to infer : \n
# A. the separating distance between the two defects R(t) \n
# B. the MSD and diffusion coeff of an individual defect. "
# # Physical Params 
# Ntarget = Int(1E3)
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


## Compute R(t) from dft;  Linear time
# Rts = Array{Vector{Float64}}(undef, length(r0s), R)
# for i in each(r0s), r in 1:R
#     dft = dfts[i, r]
#     Rts[i, r] = interdefect_distance(dft.defectsP[1], dft.defectsN[1], L, L)
# end

# Rts_avg = Vector{Vector{Float64}}(undef, length(r0s))
# for i in each(r0s)
#     tmp = vector_of_vector2matrix(Rts[i, :])
#     Rts_avg[i] = nanmean(tmp, 2)[:, 1]
# end

# p = plot(legend=false, xlabel="t", ylabel="R(t)")
# for i in each(r0s)
#     rt = Rts_avg[i]
#     plot!(rt, label="r0 = $(r0s[i])")
# end
# p

## Compute R(t*) from dft; Reversed time = time to annihilation
# Rts_reverse = Array{Vector{Float64}}(undef, length(r0s), length(v0sigs), R)
# for i in each(r0s), r in 1:R
#     dft = dfts[i, r]
#     Rts_reverse[i, r] = reverse(interdefect_distance(dft.defectsP[1], dft.defectsN[1], L, L))
# end

# Rts_reverse_avg = Vector{Vector{Float64}}(undef, length(r0s))
# for i in each(r0s)
#     tmp = vector_of_vector2matrix(Rts_reverse[i, :])
#     Rts_reverse_avg[i] = nanmean(tmp, 2)[:, 1]
# end

# p = plot(legend=:bottomright, xlabel="t*", ylabel="R(t*)", axis=:log, title="ρ = $(round(rho,digits=2))")
# for i in each(r0s)
#     rt = Rts_reverse_avg[i]
#     plot!(rt[2:end], label="r0 = $(r0s[i])", m=true)
# end
# # savefig("figures/two_defects/Rt_several_r0s_rho$(round(rho,digits=2))_v0$(v0)_sigma$(sigma)_tmax$(tmax)_R$(R)_N$(Ntarget)..png")
# p
# plot!(x -> 3x^0.33, c=:black, line=:dot)
# plot!(x -> exp(lambertw(2π / 2 * x) / 2), c=:black, line=:dash)
# ## B. Study MSD with two defects far apart
