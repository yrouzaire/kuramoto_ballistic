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
filename = "data/nature_phase_transition.jld2"
@load filename v0sigs Ps Cs ns xis rho T Ntarget params_init aspect_ratio times tmax comments rhoc runtimes R
histogram(runtimes / 3600 /24, bins=20)

Ps_avg = nanmean(Ps, 3)[:,:,1]
ns_avg = nanmean(ns, 3)[:,:,1]
xis_avg = nanmean(xis, 3)[:,:,1]

indices = [];
for r in 1:R
    try Cs[:,:,r]
		push!(indices, r)
    catch;
    end
end;
indices

Cs_avg = Array{Vector}(undef, length(v0sigs), length(times))
for i in 1:length(v0sigs), k in 1:length(times)
	Cs_avg[i,k] = mean([Cs[i,k,r] for r in indices])
end


in_green = [1,2,3,4,5]#,13,14,15,16]
in_red = [6,7,8,9]#,10,11,12]
in_all = 1:16
## 
p1 = plot(xlabel=L"t", ylabel=L"P", xscale=:log10, yscale=:log10, legend=false)#:topleft)
for i in in_green
	plot!(times, Ps_avg[i,:], label=v0sigs[i], c=i, rib=0)
end
plot!(times, x->3.2E-2sqrt(x/log(8x)),line=:dash,c=:black, label=L"\sqrt{t/\log(t)}")
p1

##
L = sqrt(Ntarget/rho)
p2 = plot(xlabel=L"t", ylabel=L"n+1",axis=:log, legend=false)#:bottomleft)
for i in in_green
	plot!(times, remove_negative(ns_avg[i,:]/L^2), label=v0sigs[i], c=i, rib=0)
end
plot!(times, x->7E-2log(3x)/x,line=:dash,c=:black, label=L"\log(t)/t}")
# plot!(times, x->3E-2/x,line=:dot,c=:black, label=L"1/t}")
p2

##
rr = 0:round(Int,sqrt(Ntarget/rho)/2)
p3 = plot(xlabel=L"r", ylabel=L"C(r,t_∞)", yaxis=:log, ylims=(1E-2,1.1), legend=false)#:topright)
for i in in_green
	plot!(rr[2:end],remove_negative(Cs_avg[i,end])[2:end], label=v0sigs[i], c=i, rib=0)
end
plot!(rr[2:end], r->0.9 * r^(-0.15),line=:dash,c=:black, label=L"r^{-T/2\pi}")
# plot!(rr[2:end], r->0.9 * r^(-T/2π),line=:dash,c=:black, label=L"r^{-T/2\pi}")
p3


##
p4 = plot(xlabel=L"r", ylabel=L"C(r,t_∞) - P^2",xaxis=:log, legend=false)#:topright)
for i in in_green
    tt = length(times)
	plot!(rr[2:end-1],remove_negative(Cs_avg[i,tt][2:end-1] .- Ps_avg[i,tt].^2), label=v0sigs[i], c=i, rib=0)
	# plot!(rr[2:end-1],abs.(remove_negative(Cs_avg[i,end])[2:end-1] .- Cs_avg[i,end][end]), label=v0sigs[i], c=i, rib=0)
end
# plot!(rr, r->1E0 * r^(-T/2π),line=:dash,c=:black, label=L"r^{-T/2\pi}")
p4


##
p5 = plot(xlabel=L"t", ylabel=L"ξ\,\sqrt{n}", xaxis=:log, legend=false)#:topright)
for i in 1:length(v0sigs)
	plot!(times, xis_avg[i,:].*sqrt.(ns_avg[i,:]), label=v0sigs[i], c=i, rib=0)
end
# plot!(rr, r->1E0 * r^(-T/2π),line=:dash,c=:black, label=L"r^{-T/2\pi}")
p5

##
plot(p1,p2,p3,p4, layout=(2,2), size=(800,800))


## ---------------- Nature of the Phase Transition ---------------- ##
comments = "The goal of this script is to pass through the transition line, 
in both direction (keeping σ or v0 constant) and to compute correlation functions.
Here for T = 0.1 and ρ = 1."
# Physical Params 
Ntarget = Int(1E3)
aspect_ratio = 1
T = 0.1
R0 = 1
rho = 1
rhoc = 4.51 / π

# Initialisation parameters
init_pos = "random"
init_theta = "hightemp"
r0 = 20.0
q = 1.0
params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)

# Simulation parameters
v0sigs = vcat([(0.25,sigm) for sigm in 0:0.025:0.2],[(v,0.1) for v in logspace(0.03,1,10)])
tmax = 1E2
# times = collect(0:tmax/30:tmax) # linear time
times = logspace(1,tmax,30,digits=1) # log time

P = zeros(length(v0sigs), length(times))
C = Array{Vector{Float64}}(undef, length(v0sigs), length(times))
xi = zeros(length(v0sigs), length(times))
n = zeros(length(v0sigs), length(times))

z = @elapsed for i in each(v0sigs)
    v0, sigma = v0sigs[i]
    println("v0 = $v0, σ = $sigma, $(100i/length(v0sigs))%")
    N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
    dt = determine_dt(T, sigma, v0, N, rho)

    param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
        :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
        :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

    system = System(param)

    t = 0.0
    token = 1

    for tt in eachindex(times)
        evolve(system, times[tt]) # evolves the systems up to times[tt]
        
        P[i,tt]  = polarOP(system)[1]
        corr_tmp    = corr(system)
        C[i,tt]  = corr_tmp
        xi[i,tt] = corr_length(corr_tmp)
        n[i,tt]  = number_defects(system)
    end

end
prinz(z)