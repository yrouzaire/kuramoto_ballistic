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
filename = "data/nature_phase_transition_horizontal.jld2"
@load filename v0sigs Ps Cs ns xis rho T Ntarget params_init aspect_ratio times tmax comments rhoc runtimes R
# histogram(runtimes / 3600 /24, bins=20)
v0sigs_horizontal = v0sigs[1:end]


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
scatter!(v0sigs_horizontal[1:1:end], c=:black, m=:circle, ms=3)
xlims!(minimum(v0s[2:end]), maximum(v0s[2:end]))
ylims!(minimum(sigmas), maximum(sigmas))
scatter!((0.2, 0.3), c=:black, m=:square, ms=5)
scatter!((2, 0.3), c=:black, m=:star5, ms=9)


p_phase_space_rho19 = heatmap(v0s[2:end], sigmas, Ps_avg[1, end, 1, 2:end, :, 1, end, 1]',
    xaxis=:log, c=cgrad([:red, :orange, :green]), clims=(0, 1),
    size=(470, 400), xlabel=L"v_0", yticks=false,
    colorbartitle="P", colorbar=:right, colorbar_titlefont=font(12))
xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])

plot(p_phase_space_rho1, p_phase_space_rho19, layout=(1, 2), size=(870, 400))
# savefig("figures_paper/phase_spaces.svg")

p_phase_space_n_rho1 = heatmap(v0s[2:end], sigmas, log10.(ns_avg[1, 1, 1, 2:end, :, 1, end, 1]' .+1),
    xaxis=:log, c=reverse(cgrad([:red, :orange, :green])),
    size=(400, 400), xlabel=L"v_0", ylabel="σ", clims=(0, 1.5),
    colorbartitle="n", colorbar=false, colorbar_titlefont=font(12), colorbar_titlefontrotation=90)
xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])


p_phase_space_n_rho19 = heatmap(v0s[2:end], sigmas, log10.(ns_avg[1, end, 1, 2:end, :, 1, end, 1]' .+1),
    xaxis=:log, c=reverse(cgrad([:red, :orange, :green])),
    size=(470, 400), xlabel=L"v_0", yticks=false, clims=(0, 1.5),
    colorbartitle=L"\log_{10}(n+1)", colorbar=:right, colorbar_titlefont=font(12))
xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])

plot(p_phase_space_n_rho1, p_phase_space_n_rho19, layout=(1, 2), size=(870, 400))
# savefig("figures_paper/phase_spaces_N.svg")

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
# savefig(p_critical_sigmas,"figures_paper/critical_sigmas.svg")

## ---------------- 1/N scaling of P ---------------- ##
## ---------------- 1/N scaling of P ---------------- ##
## ---------------- 1/N scaling of P ---------------- ##
## ---------------- 1/N scaling of P ---------------- ##

filename = "data/FSS_green.jld2"
@load filename Ntargets v0sigs Ps Cs ns xis params_init aspect_ratio times tmax T comments rho rhoc runtimes R
# hrun(runtimes)
Ps_avg = nanmean(Ps, 4)[:,:,:,1]
v0sigs

times_to_plot = [10, 15, 20, 23, 25, 26, 28, 29, 30]
p2 = plot(axis=:log, legend=:top, legend_column=2, size=(400, 400))
for t in times_to_plot
    plot!(1 ./ Ntargets, Ps_avg[2, :, t], rib=0, m=true)# label="t = $(round(Int,times[t]))")
end
plot!(1 ./ Ntargets, 8.6E-1 * (Ntargets) .^ -0.015, c=:black, line=:dot, label=L"N^{-0.015}")
plot!(1 ./ Ntargets, 6 * (Ntargets) .^ -0.5, c=:black, line=:dash, label=L"1/\sqrt{N}")
ylims!(0.028, 1.25)
yticks!([0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
    ["0.03", "", "0.05", "", "0.07", "", "", "0.1", "0.2", "0.3", "", "0.5", "", "0.7", "", "", "1"])
annotate!((0.06, 0.89), text(L"P", 15, :center, :bottom, :black))
annotate!((0.9, 0.03), text(L"1/N", 15, :center, :bottom, :black))
# xticks!([1E-5, 1E-4, 1E-3, 1E-2], [L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}"])
xticks!([3E-5, 4E-5, 5E-5, 6E-5, 7E-5, 8E-5, 9E-5, 1E-4, 2E-4, 3E-4, 4E-4, 5E-4, 6E-4, 7E-4, 8E-4, 9E-4, 1E-3, 2E-3, 3E-3, 4E-3, 5E-3, 6E-3, 7E-3, 8E-3, 9E-3, 1E-2], 
    ["", "", "", "", "", "", "", L"10^{-4}", "", "", "", "", "", "", "", "", L"10^{-3}", "", "", "", "", "", "", "", "", L"10^{-2}"])
p2

# savefig(p2, "figures_paper/FSS_1N.svg")

## ---------------- Through the transition ---------------- ##
## ---------------- Through the transition ---------------- ##
## ---------------- Through the transition ---------------- ##
## ---------------- Through the transition ---------------- ##
filename = "data/phase_space_rho_sig_v0_N1E3_tmax2500.jld2"
@load filename Ps Cs ns runtimes Ts Ns v0s rhos sigmas times_log tmax comments R
v0ss = v0s
Ps_avg_rho1_sigma01 = nanmean(Ps, 8)[1,1,1,:,3,1,end,1] # for rho = 1, and sigma = 0.1

filename = "data/nature_phase_transition_horizontal.jld2"
@load filename v0sigs Ps Cs ns xis rho T Ntarget params_init aspect_ratio times tmax comments rhoc runtimes R
indices = [];
for r in 1:R
    try
        Cs[:, :, r]
        push!(indices, r)
    catch
    end
end;
indices

Cs_avg_horizontal = Array{Vector}(undef, length(v0sigs_horizontal), length(times))
for i in 1:length(v0sigs_horizontal), k in 1:length(times)
    Cs_avg_horizontal[i, k] = mean([Cs[i, k, r] for r in indices])
end

cols_P = cgrad([:green, :orange, :red]);

##
p1 = plot(xlabel=L"t", ylabel=L"P", xscale=:log10, yscale=:log10,
    legend=false, yticks=([1E-2, 1E-1, 1], [L"10^{-2}", L"10^{-1}", L"10^{0}"]),
    xticks=([1, 10, 100, 1000, 1E4], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}"]))
for i in each(v0sigs_horizontal)
    v0 = v0sigs_horizontal[i][1]
    ind = findfirst(x -> x > v0, v0ss)
    couleur = cols_P[Ps_avg_rho1_sigma01[ind]]
    plot!(times, Ps_avg_horizontal[i, :], c=couleur, rib=0, m=:circle,
        ms=3, line=true, label="σ = $(round(v0sigs_horizontal[i][1],digits=2))")
end
plot!(times, x -> 3.2E-2sqrt(x / log(8x)), line=:dash, c=:black, label=L"\sqrt{t/\log(t)}")
# annotate!((0.1, -0.1), text("(a)", 12, :black))
p1

##
L = sqrt(Ntarget / rho)
p2 = plot(xlabel=L"t", ylabel=L"n/L^2", xscale=:log10, yscale=:log10, legend=false,
    yticks=([1E-5, 1E-4, 1E-3, 1E-2], [L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}"]),
    xticks=([1, 10, 100, 1000, 1E4], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}"]))
for i in each(v0sigs_horizontal)
    v0 = v0sigs_horizontal[i][1]
    ind = findfirst(x -> x > v0, v0ss)
    couleur = cols_P[Ps_avg_rho1_sigma01[ind]]
    plot!(times, remove_negative(ns_avg_horizontal[i, :] / L^2), label=v0sigs_horizontal[i], c=couleur, rib=0, m=:circle, ms=3, line=true)
end

plot!(times[9:end], x -> 4E-2log(8x) / x, line=:dash, c=:black, label=L"\sqrt{t/\log(t)}")
annotate!(1.4, 6E-6, text("(b)", 12))
p2

##
rr = 0:round(Int, L / 2)
p3 = plot(xlabel=L"r", ylabel=L"C(r,t_∞)", axis=:log, ylims=(1E-1, 1.1), legend=false)
for i in 1:1:length(v0sigs_horizontal)
    v0 = v0sigs_horizontal[i][1]
    ind = findfirst(x -> x > v0, v0ss)
    couleur = cols_P[Ps_avg_rho1_sigma01[ind]]
    plot!(rr[2:end], remove_negative(Cs_avg_horizontal[i, end])[2:end], label=L"v_0 = " * string(v0sigs_horizontal[i][1]), c=couleur, rib=0, m=:circle, ms=3)
end
plot!(rr[2:end], r -> r^(-T / 2π), line=:dot, c=:black, label=L"r^{-T/2\pi}")
plot!(rr[2:end], r -> r^(-0.25), line=:dash, c=:black, label=L"r^{-1/4}")
annotate!(1.25, 0.12, text("(c)", 12))
yticks!([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], [L"10^{-1}", "", "", "", "", "", "", "", "", L"10^{0}"])
xticks!([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50], [L"10^{0}", "", "", "", "", "", "", "", "", L"10^{1}", "", "", "", ""])
p3

##
p4 = plot(xlabel=L"t", ylabel=L"ξ\,\sqrt{n}", xaxis=:log, legend=false)#:topright)
for i in length(v0sigs_horizontal):-1:1
    v0 = v0sigs_horizontal[i][1]
    ind = findfirst(x -> x > v0, v0ss)
    couleur = cols_P[Ps_avg_rho1_sigma01[ind]]
    plot!(times[2:end], remove_negative(xis_avg_horizontal[i, 2:end] .* sqrt.(ns_avg_horizontal[i, 2:end])), label=v0sigs_horizontal[i], c=couleur, rib=0, m=:circle, ms=3)
end
xticks!([1, 10, 100, 1000, 1E4], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}"])
p4


plot(p1, p2, p3, p4, layout=(2, 2), size=(800, 800))