cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian, LambertW
include("../methods.jl")
const global R0 = 1
using Plots, ColorSchemes, LaTeXStrings
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

## ---------------- Analysis ---------------- ##
filename = "data/impact_init_XY.jld2"
@load filename R runtimes inits_pos Ps Cs ns xis rho T v0 sigma Ntarget params_init aspect_ratio times tmax comments rhoc 
L = round(Int,sqrt(Ntarget/rho))
histogram(runtimes / 3600 , bins=20)
times
Ps_avg = nanmean(Ps, 3)[:,:,1]
ns_avg = nanmean(ns, 3)[:,:,1]
xis_avg = nanmean(xis, 3)[:,:,1]

indices = [];
for r in 1:R
    try Cs[:,:,:,r]
		push!(indices, r)
    catch;
    end
end;

Cs_avg = Array{Vector}(undef, length(inits_pos), length(times))
for i in 1:length(inits_pos), k in 1:length(times)
	Cs_avg[i,k] = mean([Cs[i,k,r] for r in indices])
end

xis_avg = zeros(length(inits_pos), length(times))
for i in each(inits_pos)
	for t in each(times)
		xis_avg[i,t] = corr_length(Cs_avg[i,t])
	end
end
## 
p1 = plot(xlabel=L"t", ylabel=L"P", xscale=:log10, yscale=:log10, legend=:topleft)
for i in 1:length(inits_pos)
	plot!(times, Ps_avg[i,:], label=inits_pos[i], c=i, rib=0)
end
plot!(times, x->3.2E-2sqrt(x/log(10x)),line=:dash,c=:black, label=L"\sqrt{t/\log(t)}")
p1

##
p2 = plot(xlabel=L"t", ylabel=L"n+1", xscale=:log10, yscale=:log10, legend=false)#:bottomleft)
for i in 1:length(inits_pos)
	plot!(times, remove_negative(ns_avg[i,:]/L^2), label=inits_pos[i], c=i, rib=0)
end
plot!(times, x->1.2E-2log(10x)/x,line=:dash,c=:black, label=L"\log(t)/t}")
# plot!(times, x->1E3/x,line=:dot,c=:black, label=L"1/t}")
p2

##
rr = 0:round(Int,L/2)
p3 = plot(xlabel=L"r", ylabel=L"C(r,t_∞)", axis=:log, legend=:bottomleft)
for i in 1:length(inits_pos)
	plot!(rr[2:end],remove_negative(Cs_avg[i,40])[2:end], label=inits_pos[i], c=i, rib=0)
end
plot!(rr[2:end], r->1E0 * r^(-T/2π),line=:dash,c=:black, label=L"r^{-T/2\pi}")
p3

##
p4 = plot(xlabel=L"t", ylabel=L"ξ", axis=:log, legend=:topleft)
for i in 1:length(inits_pos)
	plot!(times, xis_avg[i,:]/L, label=inits_pos[i], c=i, rib=0)
end
plot!(times, x->3.2E-2sqrt(x/log(10x)),line=:dash,c=:black, label=L"\sqrt{t/\log(t)}")
p4

##
p5 = plot(xlabel=L"t", ylabel=L"ξ\,\sqrt{n}", axis=:log, legend=:bottomleft)#:topright)
for i in 1:length(inits_pos)
	plot!(times, xis_avg[i,:].*sqrt.(ns_avg[i,:]), label=inits_pos[i], c=i, rib=0)
end
p5

##
plot(p1,p2,p3,p4, layout=(2,2), size=(800,800))


## ---------------- Impact of init on XY Model ---------------- ##
comments = "Investigates the impact of initialisation for the spatial location of the 
spins in the XY model. Even though no symmetry is broken, it seems that having the spins
on (A) a regular lattice (B) a random lattice (C) a 2D Sobol sequence (D) a RSA lattice changes the behaviour of the defects.
The idea is to investigate the role of the fluctuations in the initialisation of the spins on the defects dynamics."
# Physical Params 
Ntarget = Int(1E4)
aspect_ratio = 1
T = 0.1
R0 = 1
rho = 2
v0 = 0
sigma = 0
rhoc = 4.51 / π

# Initialisation parameters
inits_pos = ["random", "square_lattice", "Sobol", "RSA"]
init_theta = "hightemp"
r0 = 20.0
q = 1.0
params_init = Dict(:init_pos => NaN, :init_theta => init_theta, :r0 => r0, :q => q)

# Simulation parameters
tmax = 1E2
times = logspace(1,tmax,10)


P = zeros(length(inits_pos), length(times))
C = Array{Vector{Float32}}(undef, length(inits_pos), length(times))
xi = zeros(length(inits_pos), length(times))
n = zeros(length(inits_pos), length(times))


# Impact on the usual quantities
z = @elapsed for i in each(inits_pos)
	init_pos = inits_pos[i]

    println("Init : $(init_pos)")
    N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
    dt = determine_dt(T, sigma, v0, N, rho)

    params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
    
    param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
        :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
        :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

    t = 0.0
    system = System(param)

	for tt in eachindex(times)
		evolve(system, times[tt]) # evolves the systems up to times[tt]
		
		P[i,tt]  = polarOP(system)[1]
		corr_tmp = corr(system)
		C[i,tt]  = corr_tmp
		xi[i,tt] = corr_length(corr_tmp)
		n[i,tt]  = number_defects(system)
	end
end
prinz(z)

## ------------ Visually compare effective number of neighbours ------------ ##
Ntarget = Int(1E5)
aspect_ratio = 1
T = 0.1
R0 = 1
rho = 1
v0 = 0
sigma = 0
rhoc = 4.51 / π

# Initialisation parameters
inits_pos = ["square_lattice", "random", "Sobol", "RSA"]
init_theta = "hightemp"
r0 = 20.0
q = 1.0
params_init = Dict(:init_pos => NaN, :init_theta => init_theta, :r0 => r0, :q => q)

pp = Vector{Any}(undef,length(inits_pos))
hh = Vector{Any}(undef,length(inits_pos))

for i in each(inits_pos)
	init_pos = inits_pos[i]
    N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
    dt = determine_dt(T, sigma, v0, N, rho)

    params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
    
    param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
        :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
        :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

    t = 0.0
    system = System(param)
	nnns = get_number_neighbours(system)
	hh[i] = nnns
	# p=plot()
	# scatter!(get_pos(system),zcolor=nnns,markersize=2,
	# markerstrokecolor=:black,markerstrokewidth=0.,
	# c=cgrad([:blue,:green,:red]),aspect_ratio=1)#,clims=(0,10))
	# title!(string(init_pos)*" $(round(mean(nnns),digits=1)) ± $(round(std(nnns),digits=1))")
	# pp[i] = p 
end
plot(pp..., layout=(2,2), size=(800,800))

## Comparing the Distributions of the number of neighbours 
using StatsBase
p=plot(yaxis=:log)
for i in [2,3,4]
	h = fit(Histogram, hh[i], nbins=10)
	plot!(h.edges[1][1:end-1],h.weights/length(hh[i]), label=inits_pos[i], c=i, rib=0)
	# histogram!(hh[i],normalize=true,bins=0:10,alpha=0.99)#, label=inits_pos[i], c=i, rib=0)
	annotate!(2.5,1/10^i,text(string(round.(mean(hh[i]),digits=2))*"±"*string(round.(std(hh[i]),digits=2)),8,ColorSchemes.tab10.colors[i]))
end
p
