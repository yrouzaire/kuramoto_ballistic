cd("/Users/yrouzaire/Documents/Recherche/GitHub/kuramoto_ballistic/figures_paper")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl");
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

## ------------------- for Kuramoto ------------------- ##
## ------------------- for Kuramoto ------------------- ##
## ------------------- for Kuramoto ------------------- ##
## ------------------- for Kuramoto ------------------- ##


## Load FSS data
filename = "data/FSS_green.jld2" 
@load filename Ntargets v0sigs Ps Cs ns xis params_init aspect_ratio times tmax T comments rho rhoc runtimes R
hrun(runtimes) 
Ps_avg = nanmean(Ps, 4)[:,:,:,1]
ns_avg = nanmean(ns, 4)[:,:,:,1]
xis_avg = nanmean(xis, 4)[:,:,:,1]
Ntargets

indices = [];
for r in 1:R
    try Cs[:,:,:,r]
		push!(indices, r)
    catch;
    end
end;

Cs_avg = Array{Vector}(undef, length(v0sigs), length(Ntargets), length(times))
for i in 1:length(v0sigs), j in 1:length(Ntargets), k in 1:length(times)
	Cs_avg[i,j,k] = mean([Cs[i,j,k,r] for r in indices])
end

## Plot P vs t
p=plot(xaxis=:log, legend=:bottomright, xlabel="t", ylabel="P(t)")
for i in 1:length(v0sigs)
	plot!(times[2:end], Ps_avg[i,10,2:end], m=true,label="v0 = $(v0sigs[i][1]) , σ = $(v0sigs[i][2])")
end
plot!(t->4E-2sqrt(t/log(10t)), c=:black )
# plot!(t->1E-2sqrt(t), c=:black, line=:dash)
p

## Time to order for several system sizes N
plot(xlabel="N", ylabel="Relaxation Time τ", legend=:topleft)
plot!([NaN, NaN], [NaN, NaN], c=1, rib=0, label=L"v_0"*"=$(v0sigs[1][1]), σ=$(v0sigs[1][2])")
plot!([NaN, NaN], [NaN, NaN], c=2, rib=0, label=L"v_0"*"=$(v0sigs[2][1]), σ=$(v0sigs[2][2])")
# σ = 0
finalPs = Ps_avg[1,:,30]
finalts = [findfirst(x->x>0.95finalPs[i], Ps_avg[1,i,2:end]) for i in 1:length(Ntargets)]
plot!(Ntargets, times[finalts],m=true,axis=:log,c=1)
# σ = 0.1
finalPs = Ps_avg[2,:,30]
finalts = [findfirst(x->x>0.95finalPs[i], Ps_avg[2,i,2:end]) for i in 1:length(Ntargets)]
plot!(Ntargets, times[finalts],m=true,axis=:log,c=2)
# fits
plot!(Ntargets, 2.5E-2Ntargets.*log.(Ntargets),c=:black,axis=:log, label=L"N\,\log(N)")
plot!(Ntargets, 1.5E-1Ntargets,c=:black,axis=:log, label=L"N", line=:dash,lw=0.7)
# savefig("figures/FSS_relaxation_time/relaxation_time.png")

## ------------- Check if Steady State reached ------------- 
p=plot(axis=:log, legend=false, xlabel="r", ylabel="C(r)", ylims=(1E-1,1.2))
for nn in 1:length(Ntargets)
	for i in 2#:length(v0sigs)
		cor = remove_negative(Cs_avg[i,nn,end])
		plot!(1:length(cor), cor, label="v0 = $(v0sigs[i][1]) , σ = $(v0sigs[i][2])")
	end
end
p

p2=plot(axis=:log, legend=false, xlabel="r", ylabel=L"C(r) - C(r_{max}) ", ylims=(1E-1,1.2))
for nn in 1:length(Ntargets)
	for i in 2#:length(v0sigs)
		cor = remove_negative(Cs_avg[i,nn,end])
		plot!(1:length(cor)-1, cor[1:end-1] .- Ps_avg[i,nn,end]^2, label="v0 = $(v0sigs[i][1]) , σ = $(v0sigs[i][2])")
	end
end
p2

plot(p,p2,layout=(1,2),size=(800,400))

## ------------- Impact of N on final quantities ------------- 
## P 
times_to_plot = [10,15,20,23,25,26,28,29,30]
p1=plot(xaxis=:log, legend=false, xlabel="1/N", ylabel="P",size=(400,400))
for t in times_to_plot
	plot!(1 ./Ntargets, Ps_avg[2,:,t], rib=0,label="t = $(round(Int,times[t]))",m=true)
end
plot!(1 ./Ntargets, 6 * (Ntargets).^-0.5, c=:black,line=:dash, label=L"1/\sqrt{N}")
plot!(1 ./Ntargets, 8.6E-1 * (Ntargets).^-0.015, c=:black,line=:dot, label=L"1/\sqrt{N}")
# plot!(1 ./Ntargets, 8E-0 ./ log.(Ntargets), c=:black,line=:dot, label=L"1/\sqrt{N}")
p1


p2=plot(axis=:log, legend=:outerright, xlabel="1/N", ylabel="P",size=(600,400))
for t in times_to_plot
	plot!(1 ./Ntargets, Ps_avg[2,:,t], rib=0,label="t = $(round(Int,times[t]))",m=true)
end
plot!(1 ./Ntargets, 6 * (Ntargets).^-0.5, c=:black,line=:dash, label=L"1/\sqrt{N}")
plot!(1 ./Ntargets, 8.6E-1 *  (Ntargets).^-0.015, c=:black,line=:dot, label=L"N^{-0.015}")
p2

plot(p1,p2,layout=(1,2),size=(1000,400))
# savefig("figures/FSS_relaxation_time/FSS.svg")
## ------------- Regular plots at fixed N ------------- 
## Plot n vs t
p=plot(axis=:log, legend=false, xlabel="t", ylabel="n(t)/L²")
for i in 1:length(v0sigs)
	plot!(times[2:end], remove_negative(ns_avg[i,end,2:end])/(Ntargets[ind_N]), label="v0 = $(v0sigs[i][1]) , σ = $(v0sigs[i][2])")
end
plot!(t->1.5E3(log(10t)/t)/(Ntargets[ind_N]), c=:black )
plot!(t->7E3(1/t)/(Ntargets[ind_N]), c=:black, line=:dash)
p

## Plot ξ vs t
ind_N = 4
p=plot(uaxis=:log, legend=false, xlabel="t", ylabel="ξ(t)", title="N = $(Ntargets[ind_N])")
for i in 1:length(v0sigs)
	plot!(times[2:end], xis_avg[i,ind_N,2:end]/sqrt(Ntargets[ind_N]), label="v0 = $(v0sigs[i][1]) , σ = $(v0sigs[i][2])")
end
# plot!(t->4sqrt(t/log(10t))/sqrt(Ntargets[ind_N]), c=:black )
# plot!(t->1E-2sqrt(t), c=:black, line=:dash)
p


## ------------------- for XY ------------------- ##
## ------------------- for XY ------------------- ##
## ------------------- for XY ------------------- ##
## ------------------- for XY ------------------- ##

filename = "data/FSS_XY_Squ_Tri.jld2"
@load filename Ns Ls lattices Ps_avg times params comments runtimes

times
Ns
Ps_avg2 = nanmean(Ps_avg, 4)[:,:,:,1]
p = plot(legend=false, axis=:log)
couleur = 1
for tt in length(times)
    plot!(1 ./ Ns, Ps_avg2[1, :, tt], rib=0, c=couleur, label="t = $(round(Int,times[tt]))", m=true)
    plot!(1 ./ Ns, Ps_avg2[2, :, tt], rib=0, c=couleur, line=:dash, m=true)
    couleur += 1
end
plot!(1 ./ Ns, 8.4E-1 * (Ns) .^ -0.0025, c=:black, line=:dot)
p