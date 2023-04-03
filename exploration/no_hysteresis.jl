cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl");
const global R0 = 1
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

## ---------------- Analysis ---------------- ##
filename = "data/hysteresis_nature_phases.jld2"
@load filename Ntarget rhos Ts inits_theta v0sigs Ps Cs ns xis aspect_ratio times tmax comments rhoc runtimes R
histogram(runtimes / 3600 /24, bins=20)
Ps_avg = nanmean(Ps, 6)[:,:,:,:,:,1]
ns_avg = nanmean(ns, 6)[:,:,:,:,:,1]
xis_avg = nanmean(xis, 6)[:,:,:,:,:,1]

indices = [];
for r in 1:R
    try Cs[:,:,:,:,:,r]
		push!(indices, r)
    catch;
    end
end;
indices

Cs_avg = Array{Vector}(undef, length(v0sigs),length(Ts), length(rhos),length(inits_theta), length(times))
for i in 1:length(v0sigs), j in 1:length(Ts), k in 1:length(rhos), l in 1:length(inits_theta), m in 1:length(times)
	Cs_avg[i,j,k,l,m] = mean([Cs[i,j,k,l,m,r] for r in indices])
end

## Les 6 courbes de P(t)
ind_rho = 2
ind_T = 2
L = sqrt(Ntarget/rhos[ind_rho])
pa = plot(xlabel=L"t", ylabel=L"P(t)",axis=:log, legend=false)#:topleft)
for i in each(v0sigs)
    plot!(times, Ps_avg[i,ind_T,ind_rho,1,:], label=v0sigs[i], c=i, rib=0, line=:solid)
    plot!(times, Ps_avg[i,ind_T,ind_rho,2,:], label=v0sigs[i], c=i, rib=0, line=:dash)
end
plot!(times, x->3.2E-2sqrt(x/log(8x)),c=:black, label=L"\sqrt{t/\log(t)}")

pb = plot(xlabel=L"t", ylabel=L"n(t)/L^2",axis=:log, legend=:outerright,size=(700,400))#:topleft)
for i in each(v0sigs)
    plot!(times, remove_negative(ns_avg[i,ind_T,ind_rho,1,:]/L^2), label=L"v_0"*" = $(v0sigs[i][1]) , σ = $(v0sigs[i][2])", c=i, rib=0, line=:solid)
    plot!(times, remove_negative(ns_avg[i,ind_T,ind_rho,2,:]/L^2), c=i, rib=0, line=:dash)
end
plot!(times, x->3.E-2log(8x)/x,c=:black, label="XY model")
plot!([NaN,NaN],[NaN,NaN],c=:grey, label="Init = disordered")
plot!([NaN,NaN],[NaN,NaN],line=:dash,c=:grey, label="Init = ordered")

# pc = plot(xlabel=L"t", ylabel=L"ξ(t)",uaxis=:log, legend=false)#:topleft)
# for i in each(v0sigs)
#     plot!(times, xis_avg[i,1,1,1,:]/L, label=v0sigs[i], c=i, rib=0, line=:solid)
#     plot!(times, xis_avg[i,1,1,2,:]/L, label=v0sigs[i], c=i, rib=0, line=:dash)
# end
# plot!(times, x->3.2E-2sqrt(x/log(8x)),c=:black, label=L"\sqrt{t/\log(t)}")

layy = @layout [a{0.27w} b{0.63w}]
plot(pa, pb, layout=layy, size=(1100,400))
title!("T = $(Ts[ind_T]), ρ = $(rhos[ind_rho])")
# savefig("figures/no_hysteresis/hysteresis_nature_phases_T$(Ts[ind_T])_rho$(rhos[ind_rho]).png")

## Summary figure
ind_rho = 1
L = sqrt(Ntarget/rhos[ind_rho])
ind_T = 2
pa = plot(xlabel=L"t", ylabel=L"P(t)",axis=:log, legend=false)#:topleft)
for i in each(Ts)
    plot!(times, Ps_avg[i,ind_T,ind_rho,1,:], label=v0sigs[i], c=i, rib=0, line=:solid)
    plot!(times, Ps_avg[i,ind_T,ind_rho,2,:], label=v0sigs[i], c=i, rib=0, line=:dash)
end
plot!(times, x->3.2E-2sqrt(x/log(8x)),c=:black, label=L"\sqrt{t/\log(t)}")


## ---------------- No Hysteresis and Nature Phases ---------------- ##
comments = "The goal of this script is to show that there is no hysteresis 
and to compute correlation functions at steady state to determine the nature
of both phases and the transition between them."
# Physical Params 
Ntarget = Int(1E3)
aspect_ratio = 1
Ts = [0,0.1]
R0 = 1
rhos = [1,2]
rhoc = 4.51 / π

# Initialisation parameters
init_pos = "random"
inits_theta = ["hightemp", "lowtemp"]
r0 = 20.0
q = 1.0
params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)

# Simulation parameters
v0sigs = [(0,0),(0,0.1), (0.05,0), (0.05,0.1), (1,0), (1,0.1)]
tmax = 1E2
# times = collect(0:tmax/30:tmax) # linear time
times = logspace(1,tmax,30,digits=1) # log time

P = zeros(length(v0sigs), length(Ts), length(rhos),length(inits_theta), length(times))
C = Array{Vector{Float64}}(undef, length(v0sigs), length(Ts), length(rhos),length(inits_theta), length(times))
xi = zeros(length(v0sigs), length(Ts), length(rhos),length(inits_theta), length(times))
n = zeros(length(v0sigs), length(Ts), length(rhos),length(inits_theta), length(times))

z = @elapsed for i in each(v0sigs), j in each(Ts), k in each(rhos), l in each(inits_theta)
    v0, sigma = v0sigs[i]
	T = Ts[j]
	rho = rhos[k]
    println("v0 = $v0, σ = $sigma, T = $T, ρ = $rho")
    N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
    dt = determine_dt(T, sigma, v0, N, rho)

	params_init = Dict(:init_pos => init_pos, :init_theta => inits_theta[l], :r0 => r0, :q => q)

    param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
        :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
        :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

    system = System(param)

    t = 0.0
    token = 1

    for tt in eachindex(times)
        evolve(system, times[tt]) # evolves the systems up to times[tt]
        
        P[i,j,k,l,tt]  = polarOP(system)[1]
        corr_tmp    = corr(system)
        C[i,j,k,l,tt]  = corr_tmp
        xi[i,j,k,l,tt] = corr_length(corr_tmp)
        n[i,j,k,l,tt]  = number_defects(system)
    end

end
prinz(z)