cd("/Users/yrouzaire/Documents/Recherche/GitHub/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl");
# const global R0 = 1
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

## ---------------- Analysis ---------------- ##
# Visualisation scan in phase space
filename = "data/phase_space_rho_sig_v0_N1E3_tmax2500.jld2"
@load filename Ps Cs ns runtimes Ts Ns v0s rhos sigmas times_log tmax comments R

Ps_avg = nanmean(Ps, 8) # N rho T v0 sigma init t R
Ps_std = nanstd(Ps, 8) # N rho T v0 sigma init t R
ns_avg = nanmean(ns, 8) # N rho T v0 sigma init t R

phsp = heatmap(v0s[2:end], sigmas, Ps_avg[1, 1, 1, 2:end, :, 1, end, 1]',
    xaxis=:log, c=cgrad([:red, :orange, :green]), clims=(0, 1),
    size=(520, 400), xlabel=L"v_0", ylabel="σ",
    colorbartitle="P", colorbar=:right, colorbar_titlefont=font(12),
    xticks=([1E-3,1E-2,1E-1,1],[L"10^{-3}",L"10^{-2}",L"10^{-1}",L"10^{0}"]))

# Horizontal
filename = "data/nature_phase_transition_horizontal.jld2"
@load filename v0sigs Ps Cs ns xis rho T Ntarget params_init aspect_ratio times tmax comments rhoc runtimes R
# histogram(runtimes / 3600 /24, bins=20)
v0sigs
v0sigs_horizontal = v0sigs[1:end]
scatter!(phsp,v0sigs_horizontal[1:1:end],c=:black,m=:circle,ms=3.5)
Ps_avg_horizontal = nanmean(Ps, 3)[:,:,1]
ns_avg_horizontal = nanmean(ns, 3)[:,:,1]
xis_avg_horizontal = nanmean(xis, 3)[:,:,1]

indices = [];
for r in 1:R
    try Cs[:,:,r]
		push!(indices, r)
    catch;
    end
end;
indices

Cs_avg_horizontal = Array{Vector}(undef, length(v0sigs_horizontal), length(times))
for i in 1:length(v0sigs_horizontal), k in 1:length(times)
	Cs_avg_horizontal[i,k] = mean([Cs[i,k,r] for r in indices])
end

# Vertical
filename = "data/nature_phase_transition_vertical.jld2"
@load filename v0sigs Ps Cs ns xis rho T Ntarget params_init aspect_ratio times tmax comments rhoc runtimes R
# histogram(runtimes / 3600 /24, bins=20)
v0sigs_vertical = v0sigs
scatter!(phsp,v0sigs_vertical,c=:black,m=:utriangle,ms=4)
ylims!(0,0.4)
xlims!(1E-3,3)
# savefig(phsp,"figures/nature_phase_transition/vizu_phase_space_scan_transition.png")

Ps_avg_vertical = nanmean(Ps, 3)[:,:,1]
ns_avg_vertical = nanmean(ns, 3)[:,:,1]
xis_avg_vertical = nanmean(xis, 3)[:,:,1]

indices = [];
for r in 1:R
    try Cs[:,:,r]
		push!(indices, r)
    catch;
    end
end;
indices

Cs_avg_vertical = Array{Vector}(undef, length(v0sigs_vertical), length(times))
for i in 1:length(v0sigs_vertical), k in 1:length(times)
	Cs_avg_vertical[i,k] = mean([Cs[i,k,r] for r in indices])
end

## ------------------ Plots Horizontal ------------------ ##
p1 = plot(xlabel=L"t", ylabel=L"P", xscale=:log10, yscale=:log10,
     legend=false,yticks=([1E-2,1E-1,1],[L"10^{-2}",L"10^{-1}",L"10^{0}"]),
     xticks=([1,10,100,1000,1E4],[L"10^{0}",L"10^{1}",L"10^{2}",L"10^{3}",L"10^{4}"]))
for i in each(v0sigs_horizontal)
	plot!(times, Ps_avg_horizontal[i,:], c=i, rib=0,m=:circle,
    ms=3,line=true, label="σ = $(round(v0sigs_horizontal[i][1],digits=2))")
end
plot!(times, x->3.2E-2sqrt(x/log(8x)),line=:dash,c=:black, label=L"\sqrt{t/\log(t)}")
annotate!(1.4,0.012,text("(a)",12))
p1


L = sqrt(Ntarget/rho)
p2 = plot(xlabel=L"t", ylabel=L"n/L^2", xscale=:log10, yscale=:log10,legend=false,
    yticks=([1E-5,1E-4,1E-3,1E-2],[L"10^{-5}",L"10^{-4}",L"10^{-3}",L"10^{-2}"]),
    xticks=([1,10,100,1000,1E4],[L"10^{0}",L"10^{1}",L"10^{2}",L"10^{3}",L"10^{4}"]))
for i in each(v0sigs_horizontal)
    plot!(times, remove_negative(ns_avg_horizontal[i,:]/L^2), label=v0sigs_horizontal[i], c=i, rib=0,m=:circle,ms=3,line=true)
end

plot!(times[9:end], x->4E-2log(8x)/x,line=:dash,c=:black, label=L"\sqrt{t/\log(t)}")
annotate!(1.4,6E-6,text("(b)",12))
p2


rr = 0:round(Int,L/2)
p3 = plot(xlabel=L"r", ylabel=L"C(r,t_∞)", axis=:log, ylims=(1E-1,1.1), legend=true,
     yticks = ([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],[L"10^{-1}","","","","","","","","",L"10^{0}"]))
for i in 1:1:length(v0sigs_horizontal) 
    v = v0sigs_horizontal[i][1]
    # plot!(rr[2:end], r->r^(-sqrt(max(0,0.005/(v-0.22)))),line=:dot,c=i, label=L"r^{-T/2\pi}")
	plot!(rr[2:end],remove_negative(Cs_avg_horizontal[i,end])[2:end], label=L"v_0 = "*string(v0sigs_horizontal[i][1]), c=i, rib=0,m=:circle,ms=3)
end
plot!(rr[2:end], r->r^(-T/2π),line=:dot,c=:black, label=L"r^{-T/2\pi}")
plot!(rr[2:end], r->r^(-0.25),line=:dash,c=:black, label=L"r^{-1/4}")
annotate!(1.25,0.12,text("(c)",12))
p3

p4 = plot(xlabel=L"r", ylabel=L"C(r,t_∞) - P^2",yaxis=:log, legend=false,
ylims=(1E-2,1.1), yticks = ([0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],[L"10^{-2}","","","","","","","","",L"10^{-1}","","","","","","","","",L"10^{0}"]),
xticks = ([1,2,3,4,5,6,7,8,9,10,20,30,40,50],[L"10^{0}","","","","","","","","",L"10^{1}","","","",""]))
tt = length(times)
for i in 1:length(v0sigs_horizontal) 
	plot!(rr[2:end],remove_negative(Cs_avg_horizontal[i,tt] .- Ps_avg_horizontal[i,tt].^2)[2:end], label=v0sigs_horizontal[i], c=i, rib=0,m=:circle,ms=3)
end
plot!(rr[2:end], r->0.25*r^(-0.25),line=:dash,c=:black, label=L"r^{-T/2\pi}")
annotate!(1.2,0.013,text("(d)",12))
p4

p5 = plot(xlabel=L"t", ylabel=L"ξ\,\sqrt{n}", xaxis=:log, legend=false)#:topright)
for i in 1:length(v0sigs_horizontal) 
	plot!(times[2:end],remove_negative(xis_avg_horizontal[i,2:end].*sqrt.(ns_avg_horizontal[i,2:end])), label=v0sigs_horizontal[i], c=i, rib=0,m=:circle,ms=3)
end
# plot!(rr, r->1E0 * r^(-T/2π),line=:dash,c=:black, label=L"r^{-T/2\pi}")
p5

plot(p1,p2,p3,p4, layout=(2,2), size=(800,800))
# savefig("figures/nature_phase_transition/PnCC_horizontal.png")

#
p4_loglin = plot(xlabel=L"r", ylabel=L"C(r,t_∞) - P^2",yaxis=:log, legend=false,
ylims=(1E-2,1.1)) 
tt = length(times)
for i in 1:length(v0sigs_horizontal) 
	plot!(rr[2:end],remove_negative(Cs_avg_horizontal[i,tt] .- Ps_avg_horizontal[i,tt].^2)[2:end], label=v0sigs_horizontal[i], c=i, rib=0)
end
plot!(rr[2:end], r->0.25*r^(-0.25),line=:dash,c=:black, label=L"r^{-T/2\pi}")
annotate!(4,0.014,text("(d)",12))
p4_loglin
# savefig("figures/nature_phase_transition/loglinCC_horizontal.png")

## ------------------ Plots Vertical ------------------ ##
p1 = plot(xlabel=L"t", ylabel=L"P", xscale=:log10, yscale=:log10,
     legend=false,yticks=([1E-2,1E-1,1],[L"10^{-2}",L"10^{-1}",L"10^{0}"]),
     xticks=([1,10,100,1000,1E4],[L"10^{0}",L"10^{1}",L"10^{2}",L"10^{3}",L"10^{4}"]))
for i in each(v0sigs_vertical)
	plot!(times, Ps_avg_vertical[i,:], c=i, rib=0,m=:utriangle,
    ms=3,line=true, label="σ = $(round(v0sigs_vertical[i][2],digits=2))")
end
plot!(times, x->3.2E-2sqrt(x/log(8x)),line=:dash,c=:black, label=L"\sqrt{t/\log(t)}")
annotate!(1.4,0.012,text("(a)",12))
p1

L = sqrt(Ntarget/rho)
p2 = plot(xlabel=L"t", ylabel=L"n/L^2", xscale=:log10, yscale=:log10,legend=false,
    yticks=([1E-5,1E-4,1E-3,1E-2],[L"10^{-5}",L"10^{-4}",L"10^{-3}",L"10^{-2}"]),
    xticks=([1,10,100,1000,1E4],[L"10^{0}",L"10^{1}",L"10^{2}",L"10^{3}",L"10^{4}"]))
for i in each(v0sigs_vertical)
    plot!(times, remove_negative(ns_avg_vertical[i,:]/L^2), label=v0sigs_vertical[i], c=i, rib=0,m=:utriangle,ms=3,line=true)
end

plot!(times[9:end], x->4E-2log(8x)/x,line=:dash,c=:black, label=L"\sqrt{t/\log(t)}")
annotate!(1.4,6E-6,text("(b)",12))
p2


rr = 0:round(Int,L/2)
p3 = plot(xlabel=L"r", ylabel=L"C(r,t_∞)", axis=:log, ylims=(1E-1,1.1), legend=false,
     yticks = ([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],[L"10^{-1}","","","","","","","","",L"10^{0}"]))
for i in 1:1:length(v0sigs_vertical) 
    v = v0sigs_vertical[i][1]
    # plot!(rr[2:end], r->r^(-sqrt(max(0,0.005/(v-0.22)))),line=:dot,c=i, label=L"r^{-T/2\pi}")
	plot!(rr[2:end],remove_negative(Cs_avg_vertical[i,end])[2:end], label=L"σ = "*string(v0sigs_vertical[i][2]), c=i, rib=0,m=:utriangle,ms=3)
end
plot!(rr[2:end], r->r^(-T/2π),line=:dot,c=:black, label=L"r^{-T/2\pi}")
plot!(rr[2:end], r->r^(-0.25),line=:dash,c=:black, label=L"r^{-1/4}")
annotate!(1.25,0.12,text("(c)",12))
p3

p4 = plot(xlabel=L"r", ylabel=L"C(r,t_∞) - P^2",axis=:log, legend=false,
ylims=(1E-2,1.1), yticks = ([0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],[L"10^{-2}","","","","","","","","",L"10^{-1}","","","","","","","","",L"10^{0}"]),
xticks = ([1,2,3,4,5,6,7,8,9,10,20,30,40,50],[L"10^{0}","","","","","","","","",L"10^{1}","","","",""]))
tt = length(times)
for i in 1:length(v0sigs_vertical) 
	plot!(rr[2:end],remove_negative(Cs_avg_vertical[i,tt] .- Ps_avg_vertical[i,tt].^2)[2:end], label=v0sigs_vertical[i], c=i, rib=0,m=:utriangle,ms=3)
end
plot!(rr[2:end], r->0.25*r^(-0.25),line=:dash,c=:black)
annotate!(1.2,0.013,text("(d)",12))
p4

p5 = plot(xlabel=L"t", ylabel=L"ξ\,\sqrt{n}", xaxis=:log, legend=false)#:topright)
for i in 1:length(v0sigs_vertical) 
	plot!(times[2:end],remove_negative(xis_avg_vertical[i,2:end].*sqrt.(ns_avg_vertical[i,2:end])), label=v0sigs_vertical[i], c=i, rib=0,m=:utriangle,ms=3)
end
# plot!(rr, r->1E0 * r^(-T/2π),line=:dash,c=:black, label=L"r^{-T/2\pi}")
p5

#
p4_loglin = plot(xlabel=L"r", ylabel=L"C(r,t_∞) - P^2",yaxis=:log, legend=false,
ylims=(1E-2,1.1)) 
tt = length(times)
for i in 1:length(v0sigs_vertical) 
	plot!(rr[2:end],remove_negative(Cs_avg_vertical[i,tt] .- Ps_avg_horizontal[i,tt].^2)[2:end], label=v0sigs_horizontal[i], c=i, rib=0)
end
plot!(rr[2:end], r->0.25*r^(-0.25),line=:dash,c=:black, label=L"r^{-T/2\pi}")
annotate!(4,0.014,text("(d)",12))
p4_loglin
# savefig("figures/nature_phase_transition/loglin_CC_vertical.png")
#

plot(p1,p2,p3,p4, layout=(2,2), size=(800,800))
# savefig("figures/nature_phase_transition/PnCC_vertical.png")

# ## ---------------- Old Analysis ---------------- ##
# filename = "data/nature_phase_transition.jld2"
# @load filename v0sigs Ps Cs ns xis rho T Ntarget params_init aspect_ratio times tmax comments rhoc runtimes R
# histogram(runtimes / 3600 /24, bins=20)

# Ps_avg = nanmean(Ps, 3)[:,:,1]
# ns_avg = nanmean(ns, 3)[:,:,1]
# xis_avg = nanmean(xis, 3)[:,:,1]

# indices = [];
# for r in 1:R
#     try Cs[:,:,r]
# 		push!(indices, r)
#     catch;
#     end
# end;
# indices

# Cs_avg = Array{Vector}(undef, length(v0sigs), length(times))
# for i in 1:length(v0sigs), k in 1:length(times)
# 	Cs_avg[i,k] = mean([Cs[i,k,r] for r in indices])
# end


# in_green = [1,2,3,4,5]#,13,14,15,16]
# in_red = [6,7,8,9]#,10,11,12]
# in_all = 1:16
# ## 
# p1 = plot(xlabel=L"t", ylabel=L"P", xscale=:log10, yscale=:log10, legend=false)#:topleft)
# for i in in_green
# 	plot!(times, Ps_avg[i,:], label=v0sigs[i], c=i, rib=0)
# end
# plot!(times, x->3.2E-2sqrt(x/log(8x)),line=:dash,c=:black, label=L"\sqrt{t/\log(t)}")
# p1

# ##
# L = sqrt(Ntarget/rho)
# p2 = plot(xlabel=L"t", ylabel=L"n+1",axis=:log, legend=false)#:bottomleft)
# for i in in_green
# 	plot!(times, remove_negative(ns_avg[i,:]/L^2), label=v0sigs[i], c=i, rib=0)
# end
# plot!(times, x->7E-2log(3x)/x,line=:dash,c=:black, label=L"\log(t)/t}")
# # plot!(times, x->3E-2/x,line=:dot,c=:black, label=L"1/t}")
# p2

# ##
# rr = 0:round(Int,sqrt(Ntarget/rho)/2)
# p3 = plot(xlabel=L"r", ylabel=L"C(r,t_∞)", yaxis=:log, ylims=(1E-2,1.1), legend=false)#:topright)
# for i in in_green
# 	plot!(rr[2:end],remove_negative(Cs_avg[i,end])[2:end], label=v0sigs[i], c=i, rib=0)
# end
# plot!(rr[2:end], r->0.9 * r^(-0.15),line=:dash,c=:black, label=L"r^{-T/2\pi}")
# # plot!(rr[2:end], r->0.9 * r^(-T/2π),line=:dash,c=:black, label=L"r^{-T/2\pi}")
# p3


# ##
# p4 = plot(xlabel=L"r", ylabel=L"C(r,t_∞) - P^2",xaxis=:log, legend=false)#:topright)
# for i in in_green
#     tt = length(times)
# 	plot!(rr[2:end-1],remove_negative(Cs_avg[i,tt][2:end-1] .- Ps_avg[i,tt].^2), label=v0sigs[i], c=i, rib=0)
# 	# plot!(rr[2:end-1],abs.(remove_negative(Cs_avg[i,end])[2:end-1] .- Cs_avg[i,end][end]), label=v0sigs[i], c=i, rib=0)
# end
# # plot!(rr, r->1E0 * r^(-T/2π),line=:dash,c=:black, label=L"r^{-T/2\pi}")
# p4


# ##
# p5 = plot(xlabel=L"t", ylabel=L"ξ\,\sqrt{n}", xaxis=:log, legend=false)#:topright)
# for i in 1:length(v0sigs)
# 	plot!(times, xis_avg[i,:].*sqrt.(ns_avg[i,:]), label=v0sigs[i], c=i, rib=0)
# end
# # plot!(rr, r->1E0 * r^(-T/2π),line=:dash,c=:black, label=L"r^{-T/2\pi}")
# p5

# ##
# plot(p1,p2,p3,p4, layout=(2,2), size=(800,800))


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