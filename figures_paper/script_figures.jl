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
    colorbartitle="P", colorbar=true, colorbar_titlefont=font(12), colorbar_titlefontrotation=90)
xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])
scatter!(v0sigs_horizontal[1:1:end], c=:black, m=:circle, ms=3)
xlims!(minimum(v0s[2:end]), maximum(v0s[2:end]))
ylims!(minimum(sigmas), maximum(sigmas))
scatter!((0.2, 0.3), c=:black, m=:square, ms=5)
scatter!((2, 0.3), c=:black, m=:star5, ms=9)

plot(p_phase_space_rho1, size=(500, 400))


# p_phase_space_rho19 = heatmap(v0s[2:end], sigmas, Ps_avg[1, end, 1, 2:end, :, 1, end, 1]',
#     xaxis=:log, c=cgrad([:red, :orange, :green]), clims=(0, 1),
#     size=(470, 400), xlabel=L"v_0", yticks=false,
#     colorbartitle="P", colorbar=:right, colorbar_titlefont=font(12))
# xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])

# plot(p_phase_space_rho1, p_phase_space_rho19, layout=(1, 2), size=(870, 400))
# savefig("figures_paper/phase_spaces_P.svg")

# p_phase_space_n_rho1 = heatmap(v0s[2:end], sigmas, log10.(ns_avg[1, 1, 1, 2:end, :, 1, end, 1]' .+1),
#     xaxis=:log, c=reverse(cgrad([:red, :orange, :green])),
#     size=(400, 400), xlabel=L"v_0", ylabel="σ", clims=(0, 1.5),
#     colorbartitle="n", colorbar=false, colorbar_titlefont=font(12), colorbar_titlefontrotation=90)
# xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])


# p_phase_space_n_rho19 = heatmap(v0s[2:end], sigmas, log10.(ns_avg[1, end, 1, 2:end, :, 1, end, 1]' .+1),
#     xaxis=:log, c=reverse(cgrad([:red, :orange, :green])),
#     size=(470, 400), xlabel=L"v_0", yticks=false, clims=(0, 1.5),
#     colorbartitle=L"\log_{10}(n+1)", colorbar=:right, colorbar_titlefont=font(12))
# xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])

# plot(p_phase_space_n_rho1, p_phase_space_n_rho19, layout=(1, 2), size=(870, 400))
# savefig("figures_paper/phase_spaces_N.svg")

## ----------------  Critical Sigmas ---------------- ##
## ----------------  Critical Sigmas ---------------- ##
## ----------------  Critical Sigmas ---------------- ##
## ----------------  Critical Sigmas ---------------- ##

filename = "data/critical_sigma.jld2"
@load filename critical_sigmas_fusion times sigmas v0s rhos tmax R #runtimes
# hrun(runtimes)
critical_sigmas_avg = nanmean(critical_sigmas_fusion, 3)[:, :, 1]
p_critical_sigmas = plot(uaxis=:log, legend=(0.0, 0.32))#, xlabel=L"\sqrt{v_0}", ylabel=L"\sigma_c")
for i in each(rhos)
    if rhos[i] ≈ 4.51 / π
        lab = L"ρ_{perco} = 1.44"
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


## ---------------- Through the transition ---------------- ##
## ---------------- Through the transition ---------------- ##
## ---------------- Through the transition ---------------- ##
## ---------------- Through the transition ---------------- ##
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);

filename = "data/phase_space_rho_sig_v0_N1E3_tmax2500.jld2"
@load filename Ps Cs ns runtimes Ts Ns v0s rhos sigmas times_log tmax comments R
v0ss = v0s
Ps_avg_rho1_sigma01 = nanmean(Ps, 8)[1,1,1,:,3,1,end,1] # for rho = 1, and sigma = 0.1

filename = "data/nature_phase_transition_horizontal.jld2"
@load filename v0sigs Ps Cs ns xis rho T Ntarget params_init aspect_ratio times tmax comments rhoc runtimes R
Ntarget
v0sigs_horizontal = v0sigs
Ps_avg_horizontal = nanmean(Ps, 3)[:, :, 1]
ns_avg_horizontal = nanmean(ns, 3)[:, :, 1]
xis_avg_horizontal = nanmean(xis, 3)[:, :, 1]

# hrun(runtimes)

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

for i in each(v0sigs_horizontal)
    for k in each(times)
        xis_avg_horizontal[i,k] = corr_length(Cs_avg_horizontal[i,k])
    end
end

cols_P = cgrad([:red, :orange, :green])

##
# p1 = plot(axis=:log,legend=false)
# for i in each(v0sigs_horizontal)
#     couleur = cols_P[Ps_avg_horizontal[i,end]]
#     plot!(times, Ps_avg_horizontal[i, :], c=couleur, rib=0, m=:circle,
#         ms=3, line=true, label="σ = $(round(v0sigs_horizontal[i][1],digits=2))")
# end
# plot!(times[5:end-8], x ->7E-2sqrt(x / log(8x)), line=:dash, c=:black, label=L"\sqrt{t/\log(t)}")
# annotate!((0.2, 0.05), text("(d)", 12))
# yticks!([0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], ["0.02", "0.03", "", "0.05", "", "0.07", "", "", "0.1", "0.2", "0.3", "", "0.5", "", "0.7", "", "", "1"])
# xticks!([1, 10, 100, 1000, 1E4, 1E5], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}", L"10^{5}"])
# annotate!((0.07,0.92), text(L"P", 15, :center, :black))
# annotate!((0.94, 0.1), text(L"t", 15, :top, :black))
# xlims!(0.9, 1E5)
# p1

##
L = sqrt(Ntarget / rho)
p2 = plot(xscale=:log10, yscale=:log10, legend=:top,
    yticks=([1E-5, 1E-4, 1E-3, 1E-2], [L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}"]),
    xticks=([1, 10, 100, 1000, 1E4], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}"]))
for i in each(v0sigs_horizontal)
        couleur = cols_P[Ps_avg_horizontal[i,end]]
    v00 = sqrt(v0sigs_horizontal[i][1])
    # plot!(times, remove_negative(ns_avg_horizontal[i, :] / L^2)*v00, c=couleur, rib=0, m=:circle, ms=3, line=true)
    plot!(times, remove_negative(ns_avg_horizontal[i, :] / L^2), c=couleur, rib=0, m=:circle, ms=3, line=true)
end
ylims!(1E-5, 1E-1)
annotate!((0.15, 0.93), text(L"n/L^2", 15, :center, :black))
annotate!((0.94, 0.1), text(L"t", 15, :top, :black))
plot!(times[5:end-8], x -> 1.3E-2log(10x) / x, line=:dash, c=:black, label=L"\log(t)/t")
annotate!((0.92, 0.93), text("(b)", 15))
p2

##
rr = 0:round(Int, L / 2)
p3 = plot(axis=:log, ylims=(1E-1, 1.3), legend=:top)
for i in 1:1:length(v0sigs_horizontal)
        couleur = cols_P[Ps_avg_horizontal[i,end]]
    plot!(rr[2:end], remove_negative(Cs_avg_horizontal[i, end])[2:end], c=couleur, rib=0, m=:circle, ms=3)
end
# plot!(rr[2:end], r -> r^(-T / 2π), line=:dot, c=:black, label=L"r^{-T/2\pi}")
plot!(rr[2:end], r -> 0.96r^(-0.25), line=:dash, c=:black, label=L"r^{-1/4}")
yticks!([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], ["0.1", "0.2", "0.3", "", "0.5", "", "0.7", "", "", "1"])
xticks!([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50], [L"10^{0}", "", "", "", "", "", "", "", "", L"10^{1}", "", "", "", ""])
annotate!((0.1, 0.945), text(L"C(r)", 15, :center, :black))
annotate!((0.94, 0.1), text(L"r", 15, :top, :black))
annotate!((0.92, 0.93), text("(c)", 15))

p3

## #
p4 = plot(xaxis=:log, legend=false, size=(290,180))#:topright)
for i in each(v0sigs_horizontal)
        couleur = cols_P[Ps_avg_horizontal[i,end]]
    data = remove_negative(xis_avg_horizontal[i, 2:end-4] .* sqrt.(ns_avg_horizontal[i, 2:end-4])) / L
    plot!(times[2:end-4], data, c=couleur, m=:circle, ms=2, lw=0.8)
end
ylims!(0.38,0.62)
yticks!([0.4, 0.5, 0.6])
annotate!((0.2, 0.89), text(L"ξ\,\sqrt{n}/L", 12, :center))
annotate!((0.93, 0.08), text(L"t", 12, :center))
xticks!([1, 10, 100, 1000, 1E4], [L"10^{0}", "", L"10^{2}", "", L"10^{4}"])
p4
# savefig(p4, "figures_paper/inset_xi.svg")

##
p5 = plot(axis=:log,legend=false)
for i in each(v0sigs_horizontal)
    couleur = cols_P[Ps_avg_horizontal[i, end]]
    v00 = sqrt(v0sigs_horizontal[i][1])
    # plot!(times, remove_negative(xis_avg_horizontal[i, :]) / L / v00, c=couleur, rib=0, m=:circle, ms=3, line=true)
    plot!(times, remove_negative(xis_avg_horizontal[i, :])/L, c=couleur, rib=0, m=:circle, ms=3, line=true)
end
p5
yticks!([ 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], ["0.03", "", "0.05", "", "0.07", "", "", "0.1", "0.2", "0.3", "", "0.5", "", "0.7"])
plot!(times[5:end-11], x -> 4.5 / L * sqrt(x / log(8x)), line=:dash, c=:black, label=L"\sqrt{t/\log(t)}")
annotate!((0.09, 1.06), text(L"\xi/L", 15, :center, :black))
annotate!((0.94, 0.1), text(L"t", 15, :top, :black))
annotate!((0.065, 0.05), text("(b)", 12))
xlims!(0.2, 4E4)


plot(p5, p2, p3, p3, layout=(1,4), size=(1600, 400))

## ---------------- 1/N scaling of P ---------------- ##
## ---------------- 1/N scaling of P ---------------- ##
## ---------------- 1/N scaling of P ---------------- ##
## ---------------- 1/N scaling of P ---------------- ##

filename = "data/FSS_green.jld2"
@load filename Ntargets v0sigs Ps Cs ns xis params_init aspect_ratio times tmax T comments rho rhoc runtimes R
# hrun(runtimes)
# times[times_to_plot]

Ps_avg = nanmean(Ps, 4)[:, :, :, 1]
v0sigs

times_to_plot = [10, 15, 20, 23, 26, 30]

inset_FSS_time = plot(axis=:log, legend=:bottomleft, size=(250, 250), box=false)
for t in times_to_plot
    plot!(Ntargets, Ps_avg[2, :, t], rib=0, m=true)# label="t = $(round(Int,times[t]))")
    # plot!(1 ./ Ntargets, Ps_avg[2, :, t], rib=0, m=true)# label="t = $(round(Int,times[t]))")
end
plot!(Ntargets, 8.6E-1 * (Ntargets) .^ -0.015, c=:black, line=:dot, label=L"N^{-0.015}")
plot!(Ntargets[2:end-2], 6 * (Ntargets[2:end-2]) .^ -0.5, c=:black, line=:dash, label=L"1/\sqrt{N}")
# plot!(1 ./ Ntargets, 8.6E-1 * (Ntargets) .^ -0.015, c=:black, line=:dot, label=L"N^{-0.015}")
# plot!(1 ./ Ntargets, 6 * (Ntargets) .^ -0.5, c=:black, line=:dash, label=L"1/\sqrt{N}")
ylims!(0.028, 1.25)
yticks!([0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
    ["0.05", "", "0.07", "", "", "0.1", "0.2", "0.3", "", "0.5", "", "0.7", "", "", "1"])
annotate!((0.07, 0.905), text(L"P", 15, :center, :bottom, :black))
annotate!((0.87, -0.006), text(L"N", 15, :center, :bottom, :black))
# xticks!([1E-5, 1E-4, 1E-3, 1E-2], [L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}"])
# xticks!([3E-5, 4E-5, 5E-5, 6E-5, 7E-5, 8E-5, 9E-5, 1E-4, 2E-4, 3E-4, 4E-4, 5E-4, 6E-4, 7E-4, 8E-4, 9E-4, 1E-3, 2E-3, 3E-3, 4E-3, 5E-3, 6E-3, 7E-3, 8E-3, 9E-3, 1E-2],
    # ["", "", "", "", "", "", "", L"10^{-4}", "", "", "", "", "", "", "", "", L"10^{-3}", "", "", "", "", "", "", "", "", L"10^{-2}"])
# annotate!((0.19, 0.05), text("(d)", 12))

ylims!(0.04, 1.1)
inset_FSS_time
# savefig(inset_FSS_time, "figures_paper/inset_FSS_time.svg")

##
filename = "data/FSS_to_determine_transition.jld2"
@load filename Ntargets v0sigs Ps Cs ns xis params_init aspect_ratio times tmax T comments rho rhoc runtimes R
hrun(runtimes)
Ntargets
v0sigs
v0s = [v for (v, x) in v0sigs]
Ps_avg = nanmean(Ps, 4)[:, :, :, 1]
Ps_std = nanstd(Ps, 4)[:, :, :, 1]
ns_avg = nanmean(ns, 4)[:, :, :, 1]
xis_avg = nanmean(xis, 4)[:, :, :, 1]

p = plot(axis=:log, legend=false, xlabel=L"N", ylabel=L"P")
for i in each(v0sigs)
    couleur = cols_P[Ps_avg[i, end, end]]
    # couleur = cols_P[Ps_avg_horizontal[i, end]]
    plot!(Ntargets, Ps_avg[i, :, end], label="v0sig = $(v0sigs[i])", rib=0Ps_std[i, :, end], m=true, c=couleur)
end
ylims!(0.03, 1.15)
yticks!([0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], ["0.03", "0.04", "0.05", "", "0.07", "", "", "0.1", "0.2", "0.3", "0.4", "0.5", "", "0.7", "", "", "1"])
plot!(Ntargets[2:end-1], x-> 3x ^ -0.5, c=:black, line=:dot, label=L"N^{-1/2}")
p


##

plot(p2, p5, p3, pFSS, layout=(2, 2), size=(800, 800))
# savefig("figures_paper/figure_through_transition.svg")

plot(p2, p5, p3,p1, pFSS, layout=(1,5), size=(2000, 400))
# savefig("figures_paper/figure_through_transition_horizontal.svg")
# savefig("figures_paper/figure_through_transition_horizontal.pdf")


## ---------------- R(t) defects  ---------------- ##
## ---------------- R(t) defects  ---------------- ##
## ---------------- R(t) defects  ---------------- ##
## ---------------- R(t) defects  ---------------- ##

filename = "data/mobility_defects_sigma_v0.jld2"
@load filename sigmas v0s Ts all_xy_pos all_xy_neg all_rr all_times_collision R_per_core Rtot params_init Ntarget R0 q init_theta init_pos aspect_ratio times tmax comments rhoc runtimes
sigmas
Rtot
rho

L = round(Int, sqrt(Ntarget))
all_rr_reverse = NaN * zeros(length(v0s), length(sigmas), length(Ts), length(times), Rtot)
all_rr_ = NaN * zeros(length(v0s), length(sigmas), length(Ts), length(times), Rtot)
for i in each(v0s)
    for j in each(sigmas)
        for k in each(Ts)
            for r in 1:Rtot
                data = all_rr[i, j, k, r]
                ll = length(data)
                all_rr_reverse[i, j, k, 1:ll, r] = reverse(data)
                all_rr_[i, j, k, 1:ll, r] = (data)
                all_rr_[i, j, k, 1+ll:end, r] .= 0
            end
        end
    end
end
all_rr_reverse_avg = nanmean(all_rr_reverse, 5)[:, :, :, :, 1]
all_rr_avg = nanmean(all_rr_, 5)[:, :, :, :, 1]

## R(t) 
ind_T = 2

phistogram = plot(size=(230, 350), legend=(0.43, 0.4), legend_title=L"v_0")
histogram!(log10.(all_times_collision[1, 1, ind_T, :]), bins=30, c=1, lw=0.2, label=L"0.5")
histogram!(log10.(all_times_collision[5, 1, ind_T, :]), bins=15, c=5, lw=0.2, label=L"2.5")
histogram!(log10.(all_times_collision[end, 1, ind_T, :]), bins=8, c=10, lw=0.2, label=L"5")
ylims!(0, 93)
xticks!(1:3, [L"10^{1}", L"10^{2}", L"10^{3}"])
xlims!(1, 4)
annotate!((0.5, 0.9), text("Distribution of " * L"\tau", 11, :center, :bottom, :black))
annotate!((0.93, 0.02), text(L"\tau", 13, :center, :bottom, :black))

##
p = plot(xaxis=:log, legend=:bottomleft, legend_title=L"v_0")
for i in each(v0s)
    for j in each(sigmas)
        for k in ind_T#each(Ts)
            plot!(times[2:end], remove_negative(all_rr_avg[i, j, k, 2:end]), label=string(v0s[i]), rib=0)
            # plot!(times[2:end], all_rr_avg[i, j, k, 2:end])
        end
    end
end
p
xticks!([1, 10, 100, 1000, 1E4], [L"10^0", L"10^1", L"10^2", L"10^3", L"10^4"])
ylims!(-0., 37)
xlims!(1, 1.6E4)
annotate!((0.12, 0.885), text(L"R(t)", 15, :center, :bottom, :black))
annotate!((0.93, 0.02), text(L"t", 15, :center, :bottom, :black))
using LambertW
mean_annihilation_time = mean(all_times_collision[1, 1, ind_T, :])
normalisation = exp(0.5 * (1 + lambertw((mean_annihilation_time - 1) * 2 / exp(1))))
# plot!(0:1:mean_annihilation_time, x -> 31.8/normalisation*exp(0.5 * (1 + lambertw((mean_annihilation_time - x) * 2 / exp(1)))), c=:black, line=:dash)
# plot!(0:1:mean_annihilation_time,x->31.8*sqrt((-x+mean_annihilation_time)/mean_annihilation_time),c=:black)
p

# savefig(p, "figures_paper/Rt.svg")
# savefig(phistogram, "figures_paper/Rt_inset.svg")

## R(t*) 
p = plot(xaxis=:log, background_color=:transparent)
for i in 2:length(v0s)
    for j in 1#each(sigmas)
        for k in ind_T
            data = all_rr_reverse_avg[i, j, k, :]
            ll = round(Int, length(data) * 0.1)
            plot!(times[2:ll], data[2:ll])
        end
    end
end
p
xlims!(0.5, 600)
xticks!([1, 10, 100], [L"10^0", L"10^1", L"10^2"])
annotate!((0.15, 0.885), text(L"R(t\!^*)", 15, :center, :bottom, :black))
annotate!((0.96, 0.02), text(L"t\!^*", 15, :center, :bottom, :black))
##

# pcollapse = plot(xlabel=L"\sqrt{v_0}t^*", ylabel=L"R(t^*)/\sqrt{v_0}", xaxis=:log, legend=false, size=(250, 250))
pcollapse = plot(xaxis=:log, legend=false, size=(250, 250))
for i in 2:length(v0s)
    for j in 1#each(sigmas)
        for k in ind_T
            data = all_rr_reverse_avg[i, j, k, :]
            ll = round(Int, length(data) * 0.1)
            plot!(sqrt(v0s[i]) * times[2:ll], 1 / sqrt(v0s[i]) * data[2:ll],
                label=L"v_0 = " * "$(v0s[i])", rib=0, c=i)
        end
    end
end
mu = 1 / 2
# plot!(times[2:700], x -> exp(0.5 * lambertw(2π * x / mu)), c=:black)
xlims!(0.9, 1E3)
ylims!(0, 30)
xticks!([1, 10, 100, 1000], [L"10^0", L"10^1", L"10^2", L"10^3"])
annotate!((0.32, 0.85), text(L"R(t\!^*)/\sqrt{v_0}", 12, :center, :bottom, :black))
annotate!((0.83, 0.04), text(L"\sqrt{v_0}t\!^*", 12, :center, :bottom, :black))
pcollapse
##
# savefig(p, "figures_paper/Rtstar.svg")
# savefig(pcollapse, "figures_paper/Rtstar_inset.svg")

## ---------------- Plots spinwaves ---------------- ##
## ---------------- Plots spinwaves ---------------- ##
## ---------------- Plots spinwaves ---------------- ##
## ---------------- Plots spinwaves ---------------- ##

filename = "data/proba_spinwaves.jld2"
@load filename R_per_core Rtot R all_nb_detected_spinwave all_times_detected_spinwave all_Ps_detected_spinwave all_thetas_detected_spinwave all_pos_detected_spinwave sigmas v0s tmax times p_threshold init_pos init_theta Ntarget rho T aspect_ratio runtimes

proba_spinwave = all_nb_detected_spinwave / Rtot
all_nb_detected_spinwave
all_times_detected_spinwave
all_Ps_detected_spinwave
all_thetas_detected_spinwave
all_pos_detected_spinwave

ind_v0, ind_sig, rr = 4, 3, 8
scatter(all_pos_detected_spinwave[ind_v0, ind_sig][rr],
    marker_z=mod.(all_thetas_detected_spinwave[ind_v0, ind_sig][rr], 2π),
    c=cols, markersize=2, aspect_ratio=1, legend=false, axis=false)
# savefig("figures_paper/spinwave_r$(rr)_v0$(v0s[ind_v0])_sig$(sigmas[ind_sig]).svg")


ind_v0, ind_sig, rr = 2,3,2
scatter(all_pos_detected_spinwave[ind_v0, ind_sig][rr],
    marker_z=mod.(all_thetas_detected_spinwave[ind_v0, ind_sig][rr], 2π),
    c=cols, markersize=2, aspect_ratio=1, legend=false, axis=false)
# savefig("figures_paper/spinwave_r$(rr)_v0$(v0s[ind_v0])_sig$(sigmas[ind_sig]).svg")


