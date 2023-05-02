cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian, LambertW
include("../methods.jl");
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

## ---------------- Analysis DFT Single Immobile Particles ---------------- ##
filename = "data/immobile_DFT_single.jld2"
@load filename R qs R0s Ts inits_pos dfts_fusion params_init Ntarget v0 sigma aspect_ratio times tmax comments rhoc runtimes
hrun(runtimes)
# times = 0:5:tmax
Reff = length(dfts_fusion)
inits_pos 
R0s
qs

dfts_fusion[1]
# MSD(t) & G(x,t) of individual defects
dx_P = [[] for i in each(inits_pos), j in each(qs), k in each(Ts), l in each(times), r in each(dfts_fusion)]
dy_P = [[] for i in each(inits_pos), j in each(qs), k in each(Ts), l in each(times), r in each(dfts_fusion)]
SD_P = [[] for i in each(inits_pos), j in each(qs), k in each(Ts), l in each(times), r in each(dfts_fusion)]
dx_M = [[] for i in each(inits_pos), j in each(qs), k in each(Ts), l in each(times), r in each(dfts_fusion)]
dy_M = [[] for i in each(inits_pos), j in each(qs), k in each(Ts), l in each(times), r in each(dfts_fusion)]
SD_M = [[] for i in each(inits_pos), j in each(qs), k in each(Ts), l in each(times), r in each(dfts_fusion)]

for r in each(dfts_fusion), i in each(inits_pos), j in each(Ts)
    dft = dfts_fusion[r][i,1,j]
    L = sqrt(Ntarget / aspect_ratio / param[:rho])
    for tt in 1:min(length(times),length(dft.defectsP[1].pos))
        d = dft.defectsP[1].pos[tt] .- dft.defectsP[1].pos[1]
        push!(dx_P[i,1,j,tt], d[1])
        push!(dy_P[i,1,j,tt], d[2])
        push!(SD_P[i,1,j,tt], dist2(dft.defectsP[1].pos[tt],dft.defectsP[1].pos[1],L,L))
    end
    
    dft = dfts_fusion[r][i,2,j]
    for tt in 1:min(length(times),length(dft.defectsN[1].pos))
        d = dft.defectsN[1].pos[tt] .- dft.defectsN[1].pos[1]
        push!(dx_M[i,2,j,tt], d[1])
        push!(dy_M[i,2,j,tt], d[2])
        push!(SD_M[i,2,j,tt], dist2(dft.defectsN[1].pos[tt],dft.defectsN[1].pos[1],L,L))
    end
end

MSD_P = zeros(length(inits_pos), length(qs), length(Ts), length(times))
for j in each(inits_pos), k in each(qs), l in each(Ts), m in each(times)
    if !isempty(SD_P[j,k,l,m])
        MSD_P[j,k,l,m] = mean(SD_P[j,k,l,m])
    end
end
MSD_M = zeros(length(inits_pos), length(qs), length(Ts), length(times))
for j in each(inits_pos), k in each(qs), l in each(Ts), m in each(times)
    if !isempty(SD_M[j,k,l,m])
        MSD_M[j,k,l,m] = mean(SD_M[j,k,l,m])
    end
end

MSD_all = (MSD_P + MSD_M) / 2

## ---------------- Histogram Displacements ---------------- ##
tt = length(times) ; histogram(vcat(dx_P[1,1,tt],dy_P[1,1,tt],dx_M[1,1,tt],dy_M[1,1,tt]), bins=50)

## ---------------- Plot MSD(t) ---------------- ##
p=plot(xlabel="t", ylabel="MSD", legend=:topleft, axis=:log)
for i in [1,2,3] #(inits_pos)
    for j in 1# 1+ 2-
        for k in [1] #(Ts)
            if j == 1 
                plot!(times[2:end], remove_negative(MSD_P[i,j,k,2:end]),rib=0)
            else
                plot!(times[2:end], remove_negative(MSD_M[i,j,k,2:end]),rib=0)
            end
        end
    end
end
plot!([NaN,NaN],c=1,rib = 0,label="Square")
plot!([NaN,NaN],c=2,rib = 0,label="RSA")
plot!([NaN,NaN],c=3,rib = 0,label="Random")
plot!(times[2:end],x->3E-1x^1,c=:black,line=:dash)
p
# savefig("figures/two_defects/MSD.png")

