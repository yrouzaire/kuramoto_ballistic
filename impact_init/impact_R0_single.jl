cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian, LambertW
include("../methods.jl");
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

## ---------------- Analysis DFT Single Immobile Particles ---------------- ##
filename = "data/immobile_DFT_single_impactR0.jld2"
@load filename R R0s Ts inits_pos R_per_core dfts_fusion params_init Ntarget v0 sigma aspect_ratio times tmax comments rhoc runtimes
hrun(runtimes)
times = 0:5:tmax
Reff = length(dfts_fusion)
Rtot = R_per_core * Reff
inits_pos 
R0s

dfts_fusion[1]
# MSD(t) & G(x,t) of individual defects
dx_P = [[] for i in each(inits_pos), j in each(R0s), k in each(Ts), l in each(times)]
dy_P = [[] for i in each(inits_pos), j in each(R0s), k in each(Ts), l in each(times)]
SD_P = [[] for i in each(inits_pos), j in each(R0s), k in each(Ts), l in each(times)]

for r in 1:Reff, i in each(inits_pos), j in each(R0s), k in each(Ts)
    for r_per_core in 1:R_per_core
        dft = dfts_fusion[r][i,j,k,r_per_core]
        L = sqrt(Ntarget / aspect_ratio / param[:rho])
        for tt in 1:min(length(times),length(dft.defectsP[1].pos))
            d = dft.defectsP[1].pos[tt] .- dft.defectsP[1].pos[1]
            push!(dx_P[i,j,k,tt], d[1])
            push!(dy_P[i,j,k,tt], d[2])
            push!(SD_P[i,j,k,tt], dist2(dft.defectsP[1].pos[tt],dft.defectsP[1].pos[1],L,L))
        end
    end
end

MSD_P = zeros(length(inits_pos), length(R0s), length(Ts), length(times))
std_P = zeros(length(inits_pos), length(R0s), length(Ts), length(times))
for j in each(inits_pos), k in each(R0s), l in each(Ts), m in each(times)
    if !isempty(SD_P[j,k,l,m])
        MSD_P[j,k,l,m] = mean(SD_P[j,k,l,m])
        std_P[j,k,l,m] = std(SD_P[j,k,l,m])
    end
end

## ---------------- Histogram Displacements ---------------- ##
tt = length(times) ; histogram(vcat(dx_P[1,1,tt],dy_P[1,1,tt],dx_M[1,1,tt],dy_M[1,1,tt]), bins=50)

## ---------------- Plot MSD(t) ---------------- ##
R0s_real = [1, sqrt(2), 2, sqrt(5)]
p=plot(xlabel="t", ylabel="MSD", legend=:topleft, uaxis=:log)
for i in [1] #(inits_pos)
    for j in 1:length(R0s)
        for k in 1:length(Ts)
            plot!(times[2:end], remove_negative(MSD_P[i,j,k,2:end]),rib=0std_P[i,j,k,2:end],c=j)
        end
    end
end
p
plot!([NaN,NaN],c=1,rib = 0,label=L"R_0 = 1 (4NN)")
plot!([NaN,NaN],c=2,rib = 0,label=L"R_0 = \sqrt{2} (8NN)")
plot!([NaN,NaN],c=3,rib = 0,label=L"R_0 = 2 (12NN)")
plot!([NaN,NaN],c=4,rib = 0,label=L"R_0 = \sqrt{5} (20NN)")
plot!(times[2:end],x->1E-1x^1,c=:black,line=:dash,label=L"t")
ylims!(0,20)
p
# savefig("impact_init/figures/impactR0_MSD.png")
