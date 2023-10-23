cd("/Users/yrouzaire/Documents/Recherche/GitHub/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl");
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&


## ---------------- Phase Spaces rho = 1, 1.9 ---------------- ##
## ---------------- Phase Spaces rho = 1, 1.9 ---------------- ##
## ---------------- Phase Spaces rho = 1, 1.9 ---------------- ##
## ---------------- Phase Spaces rho = 1, 1.9 ---------------- ##


filename = "data/phase_space_rho_sig_v0_N1E3_tmax2500.jld2"
@load filename Ps Cs ns runtimes Ts Ns v0s rhos sigmas times_log tmax comments R

Ps_avg = nanmean(Ps, 8) # N rho T v0 sigma init t R
Ps_std = nanstd(Ps, 8) # N rho T v0 sigma init t R
ns_avg = nanmean(ns, 8) # N rho T v0 sigma init t R

p_phase_space_rho1 = heatmap(v0s[2:end], sigmas, Ps_avg[1, 1, 1, 2:end, :, 1, end, 1]',
    xaxis=:log, c=cgrad([:red, :orange, :green]), clims=(0, 1),
    size=(400, 400), xlabel=L"v_0", ylabel="σ",
    colorbartitle="P", colorbar=false, colorbar_titlefont=font(12), colorbar_titlefontrotation=90)
xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])


p_phase_space_rho19 = heatmap(v0s[2:end], sigmas, Ps_avg[1, end, 1, 2:end, :, 1, end, 1]',
    xaxis=:log, c=cgrad([:red, :orange, :green]), clims=(0, 1),
    size=(470, 400), xlabel=L"v_0", yticks=false,
    colorbartitle="P", colorbar=:right, colorbar_titlefont=font(12))
xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])


## ----------------  Critical Sigmas ---------------- ##
## ----------------  Critical Sigmas ---------------- ##
## ----------------  Critical Sigmas ---------------- ##
## ----------------  Critical Sigmas ---------------- ##

filename = "data/critical_sigma.jld2"
@load filename critical_sigmas_fusion times sigmas v0s rhos tmax R #runtimes
hrun(runtimes)

p_critical_sigmas = plot(uaxis=:log, legend=(0.0, 0.32))#, xlabel=L"\sqrt{v_0}", ylabel=L"\sigma_c")
for i in each(rhos)
    if rhos[i] ≈ 4.51 / π
        lab = L"ρ = ρ_c"
    else
        lab = L"ρ = "*string(rhos[i])
    end
    plot!(sqrt.(v0s), critical_sigmas_avg[:, i], rib=0, label=lab, m=true, ms=3)
end
# plot!(x -> 0.5sqrt(1.8) * (0 + (x + 0.0)), c=:black, l=:solid)
# plot!(sqrt.(v0s[11:end]), x -> 0.5max(x - 0.3, 0), c=:black, l=:dash, label="Slope 1")
p_critical_sigmas
annotate!((0.08, 0.9), text(L"\sigma_c", 15, :center, :bottom, :black))
annotate!((0.89, 0.03), text(L"\sqrt{v_0}", 15, :center, :bottom, :black))

plot(p_phase_space_rho1, p_phase_space_rho19, p_critical_sigmas, layout=(1, 3), size=(1300, 400))
savefig("figures_paper/phase_space_and_critical_sigmas.svg")

## ---------------- 1/N scaling of P ---------------- ##
## ---------------- 1/N scaling of P ---------------- ##
## ---------------- 1/N scaling of P ---------------- ##
## ---------------- 1/N scaling of P ---------------- ##

filename = "data/FSS_green.jld2"
@load filename Ntargets v0sigs Ps Cs ns xis params_init aspect_ratio times tmax T comments rho rhoc runtimes R
hrun(runtimes)
Ps_avg
v0sigs

times_to_plot = [10, 15, 20, 23, 25, 26, 28, 29, 30]
p2 = plot(axis=:log, legend=:topright, size=(400, 400))
for t in times_to_plot
    plot!(1 ./ Ntargets, Ps_avg[2, :, t], rib=0, m=true)# label="t = $(round(Int,times[t]))")
end
plot!(1 ./ Ntargets, 6 * (Ntargets) .^ -0.5, c=:black, line=:dash, label=L"1/\sqrt{N}")
plot!(1 ./ Ntargets, 8.6E-1 * (Ntargets) .^ -0.015, c=:black, line=:dot, label=L"N^{-0.015}")
ylims!(0.028, 1.25)
yticks!([0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
    ["0.03", "", "0.05", "", "0.07", "", "", "0.1", "0.2", "0.3", "", "0.5", "", "0.7", "", "", "1"])
annotate!((0.06, 0.89), text(L"P", 15, :center, :bottom, :black))
annotate!((0.9, 0.03), text(L"1/N", 15, :center, :bottom, :black))
# xticks!([1E-5, 1E-4, 1E-3, 1E-2], [L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}"])
xticks!([3E-5, 4E-5, 5E-5, 6E-5, 7E-5, 8E-5, 9E-5, 1E-4, 2E-4, 3E-4, 4E-4, 5E-4, 6E-4, 7E-4, 8E-4, 9E-4, 1E-3, 2E-3, 3E-3, 4E-3, 5E-3, 6E-3, 7E-3, 8E-3, 9E-3, 1E-2], 
    ["", "", "", "", "", "", "", L"10^{-4}", "", "", "", "", "", "", "", "", L"10^{-3}", "", "", "", "", "", "", "", "", L"10^{-2}"])
p2
