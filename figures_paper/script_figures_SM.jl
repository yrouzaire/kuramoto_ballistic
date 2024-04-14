cd("/Users/yrouzaire/Documents/Recherche/GitHub/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl");
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&


## ---------------- Phase Spaces P and n ---------------- ##
## ---------------- Phase Spaces P and n ---------------- ##
## ---------------- Phase Spaces P and n ---------------- ##
## ---------------- Phase Spaces P and n ---------------- ##

filename = "data/phase_space_rho_sig_v0_N1E3_tmax2500.jld2"
@load filename Ps Cs ns runtimes Ts Ns v0s rhos sigmas times_log tmax comments R

Ps_avg = nanmean(Ps, 8) # N rho T v0 sigma init t R
Ps_std = nanstd(Ps, 8) # N rho T v0 sigma init t R
ns_avg = nanmean(ns, 8) # N rho T v0 sigma init t R
ns_std = nanstd(ns, 8) # N rho T v0 sigma init t R

##

p_phase_space_rho1 = heatmap(v0s[2:end], sigmas, Ps_avg[1, 1, 1, 2:end, :, 1, end, 1]',
    xaxis=:log, c=cgrad([:red, :orange, :green]), clims=(0, 1),
    size=(420, 420), xlabel=L"v_0", ylabel="σ",
    colorbartitle="P", colorbar=false, colorbar_titlefont=font(12), colorbar_titlefontrotation=90)
xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])
# scatter!(v0sigs_horizontal[1:1:end], c=:black, m=:circle, ms=3)
xlims!(minimum(v0s[2:end]), maximum(v0s[2:end]))
ylims!(minimum(sigmas), maximum(sigmas))
# scatter!((0.2, 0.3), c=:black, m=:square, ms=5)
# scatter!((2, 0.3), c=:black, m=:star5, ms=9)
title!("ρ = 1")
annotate!((0.1, 0.9), text("(a)", 15))

# #

p_phase_space_rho19 = heatmap(v0s[2:end], sigmas, Ps_avg[1, end, 1, 2:end, :, 1, end, 1]',
    xaxis=:log, c=cgrad([:red, :orange, :green]), clims=(0, 1),
    size=(500, 420), xlabel=L"v_0", yticks=false,
    colorbartitle="P", colorbar=:right, colorbar_titlefont=font(12))
xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])
title!("ρ = 1.9")
annotate!((0.1, 0.9), text("(b)", 15))

# #

p_phase_space_n_rho1 = heatmap(v0s[2:end], sigmas, log10.(ns_avg[1, 1, 1, 2:end, :, 1, end, 1]' .+ 1),
    xaxis=:log, c=reverse(cgrad([:red, :orange, :green])),
    size=(420, 400), xlabel=L"v_0", ylabel="σ", clims=(0, 1.5),
    colorbartitle=L"\log_{10}\, n", colorbar=false, colorbar_titlefont=font(12), colorbar_titlefontrotation=90)
xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])
p_phase_space_n_rho1
annotate!((0.1, 0.9), text("(c)", 15))


# #
p_phase_space_n_rho19 = heatmap(v0s[2:end], sigmas, log10.(ns_avg[1, end, 1, 2:end, :, 1, end, 1]' .+ 1),
    xaxis=:log, c=reverse(cgrad([:red, :orange, :green])),
    size=(500, 400), xlabel=L"v_0", clims=(0, 1.5), yticks=false,
    colorbartitle=L"\log_{10}\, n", colorbar=true, colorbar_titlefont=font(12), colorbar_titlefontrotation=90)
xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])
p_phase_space_n_rho19
annotate!((0.1, 0.9), text("(d)", 15))



# #

data = (ns_std[1, 1, 1, 2:end, :, 1, end, 1] ./ (ns_avg[1, 1, 1, 2:end, :, 1, end, 1] .+ 0.001))'
p_phase_space_stdn_rho1 = heatmap(v0s[2:end], sigmas, log.(data),
    xaxis=:log, c=reverse(cgrad([:red, :orange, :green])),
    size=(420, 400), xlabel=L"v_0", ylabel="σ", clims=(-1.5, 1.5),
    colorbartitle="log(std(n) / avg(n))", colorbar=false, colorbar_titlefont=font(12), colorbar_titlefontrotation=90)
xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])
p_phase_space_stdn_rho1
annotate!((0.1, 0.9), text("(e)", 15))


# #

data = (ns_std[1, end, 1, 2:end, :, 1, end, 1] ./ (ns_avg[1, end, 1, 2:end, :, 1, end, 1] .+ 0.01))'
p_phase_space_stdn_rho19 = heatmap(v0s[2:end], sigmas, log.(data),
    xaxis=:log, c=reverse(cgrad([:red, :orange, :green])),
    size=(500, 400), xlabel=L"v_0", clims=(-1.5, 1.5),yticks=false,
    colorbartitle="log(std(n) / avg(n))", colorbar=true, colorbar_titlefont=font(12), colorbar_titlefontrotation=90)
xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])
p_phase_space_stdn_rho19
annotate!((0.1, 0.9), text("(f)", 15))


# #

plot(p_phase_space_rho1, p_phase_space_rho19, p_phase_space_n_rho1, p_phase_space_n_rho19, p_phase_space_stdn_rho1, p_phase_space_stdn_rho19, layout=(3,2), size=(1000, 1220))
# savefig("figures_paper/SM/phase_space_P_n_std.svg")
# savefig("figures_paper/SM/phase_space_P_n_std.pdf")

## ---------------- Through the transition Rescaling ---------------- ##
## ---------------- Through the transition Rescaling ---------------- ##
## ---------------- Through the transition Rescaling ---------------- ##
## ---------------- Through the transition Rescaling ---------------- ##

filename = "data/phase_space_rho_sig_v0_N1E3_tmax2500.jld2"
@load filename Ps Cs ns runtimes Ts Ns v0s rhos sigmas times_log tmax comments R
v0ss = v0s
Ps_avg_rho1_sigma01 = nanmean(Ps, 8)[1, 1, 1, :, 3, 1, end, 1] # for rho = 1, and sigma = 0.1

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
        xis_avg_horizontal[i, k] = corr_length(Cs_avg_horizontal[i, k])
    end
end

cols_P = cgrad([:red, :orange, :green])

##
L = sqrt(Ntarget / rho)
p1 = plot(xscale=:log10, yscale=:log10, legend=:bottomleft,
    yticks=([1E-5, 1E-4, 1E-3, 1E-2], [L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}"]),
    xticks=([1, 10, 100, 1000, 1E4], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}"]))
for i in each(v0sigs_horizontal)
    couleur = cols_P[Ps_avg_horizontal[i, end]]
    v00 = sqrt(v0sigs_horizontal[i][1])
    plot!(times * v00, remove_negative(ns_avg_horizontal[i, :] / L^2) * v00^2, c=couleur, rib=0, m=false, ms=3, lw=1.5)
    # plot!(times, remove_negative(ns_avg_horizontal[i, :] / L^2) , c=couleur, rib=0, m=false, ms=3, lw=1.5)
end
ylims!(1E-5, 1E-1)
annotate!((0.19, 0.93), text(L"v_0 \,n/L^2", 15, :center, :black))
annotate!((0.94, 0.1), text(L"t\sqrt{v_0}", 15, :top, :black))
plot!(times[5:end-8], x -> 1.3E-2log(10x) / x, line=:dash, c=:black, label=L"\log(t)/t")
annotate!((0.92, 0.93), text("(a)", 15))
xticks!([1, 10, 100, 1000, 1E4, 1E5], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}", L"10^{5}"])
p1

##
p2 = plot(axis=:log, legend=:top)
for i in each(v0sigs_horizontal)
    couleur = cols_P[Ps_avg_horizontal[i, end]]
    v00 = sqrt(v0sigs_horizontal[i][1])
    # plot!(times * v00, remove_negative(xis_avg_horizontal[i, :]) / L, c=couleur, rib=0, m=:false, ms=3, lw=1.5)
    plot!(times * v00, remove_negative(xis_avg_horizontal[i, :]) / L, c=couleur, rib=0, m=:false, ms=3, lw=1.5)
    # plot!(times * v00, remove_negative(xis_avg_horizontal[i, :]) / v00 / L, c=couleur, rib=0, m=:false, ms=3, lw=1.5)
end
p2
yticks!([0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], ["0.03", "", "0.05", "", "0.07", "", "", "0.1", "0.2", "0.3", "", "0.5", "", "0.7"])
plot!(times[5:end-13], x -> 4.3 / L * sqrt(x / log(8x)), line=:dash, c=:black, label=L"\sqrt{t/\log(t)}")
annotate!((0.09, 0.93), text(L"\xi/L", 15, :center, :black))
annotate!((0.874, 0.12), text(L"t\sqrt{v_0}", 15, :top, :black))
annotate!((0.92, 0.93), text("(b)", 15))
xlims!(0.2, 14E4)
ylims!(0.022, 0.7)
xticks!([1, 10, 100, 1000, 1E4, 1E5], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}", L"10^{5}"])

## 
plot(p1, p2, layout=(1, 2), size=(800, 400))
# savefig("figures_paper/SM/phase_space_n_xi_rescaled.svg")




## ---------------- Through the transition Vertical / Horizontal ---------------- ##
## ---------------- Through the transition Vertical / Horizontal ---------------- ##
## ---------------- Through the transition Vertical / Horizontal ---------------- ##
## ---------------- Through the transition Vertical / Horizontal ---------------- ##
cols_P = cgrad([:red, :orange, :green])
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);

filename = "data/phase_space_rho_sig_v0_N1E3_tmax2500.jld2"
@load filename Ps Cs ns runtimes Ts Ns v0s rhos sigmas times_log tmax comments R
v0ss = v0s
Ps_avg_rho1_sigma01 = nanmean(Ps, 8)[1, 1, 1, :, 3, 1, end, 1] # for rho = 1, and sigma = 0.1

filename = "data/nature_phase_transition_horizontal.jld2"
@load filename v0sigs Ps Cs ns xis rho T Ntarget params_init aspect_ratio times tmax comments rhoc runtimes R
Ntarget
v0sigs_horizontal = v0sigs
Ps_avg_horizontal = nanmean(Ps, 3)[:, :, 1]
ns_avg_horizontal = nanmean(ns, 3)[:, :, 1]
xis_avg_horizontal = nanmean(xis, 3)[:, :, 1]



phsp = heatmap(v0s[2:end], sigmas, Ps_avg[1, 1, 1, 2:end, :, 1, end, 1]',
    xaxis=:log, c=cgrad([:red, :orange, :green]), clims=(0, 1),
    size=(520, 400), xlabel=L"v_0", ylabel="σ",
    colorbartitle="P", colorbar=:right, colorbar_titlefont=font(12),
    xticks=([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"]))

scatter!(phsp, v0sigs_horizontal[1:1:end], c=:black, m=:circle, ms=4)

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
        xis_avg_horizontal[i, k] = corr_length(Cs_avg_horizontal[i, k])
    end
end

# Vertical
filename = "data/nature_phase_transition_vertical.jld2"
@load filename v0sigs Ps Cs ns xis rho T Ntarget params_init aspect_ratio times tmax comments rhoc runtimes R
v0sigs_vertical = v0sigs
scatter!(phsp, v0sigs_vertical, c=:black, m=:utriangle, ms=5.5)
display(phsp)
ylims!(0, 0.4)
xlims!(1E-3, 3)
# savefig(phsp, "figures_paper/SM/phase_space_P_vertical_horizontal.svg")
# savefig(phsp, "figures_paper/SM/phase_space_P_vertical_horizontal.pdf")

Ps_avg_vertical = nanmean(Ps, 3)[:, :, 1]
ns_avg_vertical = nanmean(ns, 3)[:, :, 1]
xis_avg_vertical = nanmean(xis, 3)[:, :, 1]

indices = [];
for r in 1:R
    try
        Cs[:, :, r]
        push!(indices, r)
    catch
    end
end;
indices

Cs_avg_vertical = Array{Vector}(undef, length(v0sigs_vertical), length(times))
for i in 1:length(v0sigs_vertical), k in 1:length(times)
    Cs_avg_vertical[i, k] = mean([Cs[i, k, r] for r in indices])
end


##
L = sqrt(Ntarget / rho)
p1_horizontal = plot(xscale=:log10, yscale=:log10, legend=:bottomleft,
    yticks=([1E-5, 1E-4, 1E-3, 1E-2], [L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}"]),
    xticks=([1, 10, 100, 1000, 1E4], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}"]))
for i in each(v0sigs_horizontal)
    couleur = cols_P[Ps_avg_horizontal[i, end]]
    plot!(times, remove_negative(ns_avg_horizontal[i, :] / L^2), c=couleur, rib=0, m=false, ms=3, lw=1.5)
end
ylims!(1E-5, 1E-1)
annotate!((0.19, 0.93), text(L"n/L^2", 15, :center, :black))
annotate!((0.94, 0.1), text(L"t", 15, :top, :black))
plot!(times[5:end-8], x -> 1.3E-2log(10x) / x, line=:dash, c=:black, label=L"\log(t)/t")
annotate!((0.8, 0.93), text(" ●", 10))
annotate!((0.92, 0.93), text("(a)", 15))
xticks!([1, 10, 100, 1000, 1E4, 1E5], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}", L"10^{5}"])
p1_horizontal

##
p1_vertical = plot(xscale=:log10, yscale=:log10, legend=:bottomleft,
    yticks=([1E-5, 1E-4, 1E-3, 1E-2], [L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}"]),
    xticks=([1, 10, 100, 1000, 1E4], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}"]))
for i in each(v0sigs_vertical)
    couleur = cols_P[Ps_avg_vertical[i, end]]
    plot!(times, remove_negative(ns_avg_vertical[i, :] / L^2), c=couleur, rib=0, m=false, ms=3, lw=1.5)
end
ylims!(1E-5, 1E-1)
annotate!((0.19, 0.93), text(L"n/L^2", 15, :center, :black))
annotate!((0.98, 0.1), text(L"t", 15, :top, :black))
plot!(times[7:end-3], x -> 3E-2log(10x) / x, line=:dash, c=:black, label=L"\log(t)/t")
annotate!((0.92, 0.93), text(L"\blacktriangle \, (d)", 15))
xticks!([1, 10, 100, 1000, 1E4, 1E5], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}", L"10^{5}"])
p1_vertical

##
p2_horizontal = plot(axis=:log, legend=:top)
for i in each(v0sigs_horizontal)
    couleur = cols_P[Ps_avg_horizontal[i, end]]
    plot!(times, remove_negative(xis_avg_horizontal[i, :]) / L, c=couleur, rib=0, m=:false, ms=3, lw=1.5)
end
yticks!([0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], ["0.03", "", "0.05", "", "0.07", "", "", "0.1", "0.2", "0.3", "", "0.5", "", "0.7"])
plot!(times[5:end-13], x -> 4.3 / L * sqrt(x / log(8x)), line=:dash, c=:black, label=L"\sqrt{t/\log(t)}")
annotate!((0.09, 0.93), text(L"\xi/L", 15, :center, :black))
annotate!((0.95, 0.12), text(L"t", 15, :top, :black))
annotate!((0.8, 0.93), text(" ●", 10))
annotate!((0.92, 0.93), text("(b)", 15))
xlims!(0.2, 14E4)
ylims!(0.022, 0.7)
xticks!([1, 10, 100, 1000, 1E4, 1E5], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}", L"10^{5}"])
p2_horizontal

##
p2_vertical = plot(axis=:log, legend=:top)
for i in each(v0sigs_vertical)
    couleur = cols_P[Ps_avg_vertical[i, end]]
    plot!(times, remove_negative(xis_avg_vertical[i, :]) / L, c=couleur, rib=0, m=:false, ms=3, lw=1.5)
end
yticks!([0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], ["0.03", "", "0.05", "", "0.07", "", "", "0.1", "0.2", "0.3", "", "0.5", "", "0.7"])
plot!(times[11:end-8], x -> 2.7 / L * sqrt(x / log(8x)), line=:dash, c=:black, label=L"\sqrt{t/\log(t)}")
annotate!((0.09, 0.93), text(L"\xi/L", 15, :center, :black))
annotate!((0.95, 0.12), text(L"t", 15, :top, :black))
annotate!((0.89, 0.93), text(L"\blacktriangle \, (e)", 15))
xlims!(0.2, 14E4)
ylims!(0.022, 0.7)
xticks!([1, 10, 100, 1000, 1E4, 1E5], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}", L"10^{5}"])
p2_vertical


##
rr = 0:round(Int, L / 2)
p3_horizontal = plot(axis=:log, ylims=(1E-1, 1.3), legend=:top)
for i in 1:1:length(v0sigs_horizontal)
    couleur = cols_P[Ps_avg_horizontal[i, end]]
    plot!(rr[2:end], remove_negative(Cs_avg_horizontal[i, end])[2:end], c=couleur, rib=0, m=false, lw=1.5, ms=3)
end
# plot!(rr[2:end], r -> r^(-T / 2π), line=:dot, c=:black, label=L"r^{-T/2\pi}")
plot!(rr[2:end], r -> 0.96r^(-0.25), line=:dash, c=:black, label=L"r^{-1/4}")
yticks!([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], ["0.1", "0.2", "0.3", "", "0.5", "", "0.7", "", "", "1"])
xticks!([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50], [L"10^{0}", "", "", "", "", "", "", "", "", L"10^{1}", "", "", "", ""])
annotate!((0.1, 0.945), text(L"C(r)", 15, :center, :black))
annotate!((0.94, 0.1), text(L"r", 15, :top, :black))
annotate!((0.91, 0.93), text("(c)", 15))
annotate!((0.8, 0.93), text(" ●", 10))
p3_horizontal
## 

rr = 0:round(Int, L / 2)
p3_vertical = plot(axis=:log, ylims=(1E-1, 1.3), legend=:top)
for i in 1:1:length(v0sigs_vertical)
    couleur = cols_P[Ps_avg_vertical[i, end]]
    plot!(rr[2:end], remove_negative(Cs_avg_vertical[i, end])[2:end], c=couleur, rib=0, m=false, lw=1.5, ms=3)
end
# plot!(rr[2:end], r -> r^(-T / 2π), line=:dot, c=:black, label=L"r^{-T/2\pi}")
plot!(rr[2:end], r -> 0.96r^(-0.25), line=:dash, c=:black, label=L"r^{-1/4}")
yticks!([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], ["0.1", "0.2", "0.3", "", "0.5", "", "0.7", "", "", "1"])
xticks!([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50], [L"10^{0}", "", "", "", "", "", "", "", "", L"10^{1}", "", "", "", ""])
annotate!((0.1, 0.945), text(L"C(r)", 15, :center, :black))
annotate!((0.94, 0.1), text(L"r", 15, :top, :black))
annotate!((0.9, 0.93), text(L"\blacktriangle \, (f)", 15))
p3_vertical
## 
p_vertical_horizontal = plot(p1_horizontal, p2_horizontal, p3_horizontal, p1_vertical, p2_vertical, p3_vertical, layout=(2, 3), size=(1200, 800));
# savefig(p_vertical_horizontal, "figures_paper/SM/phase_space_vertical_horizontal.svg")


## ---------------- 1/N scaling of P ---------------- ##
## ---------------- 1/N scaling of P ---------------- ##
## ---------------- 1/N scaling of P ---------------- ##
## ---------------- 1/N scaling of P ---------------- ##

filename = "data/FSS_green.jld2"
@load filename Ntargets v0sigs Ps Cs ns xis params_init aspect_ratio times tmax T comments rho rhoc runtimes R
hrun(runtimes)

Ps_avg = nanmean(Ps, 4)[:, :, :, 1]
v0sigs

times_to_plot = [10, 15, 20, 23, 26, 30]

FSS_time = plot(axis=:log, legend=:bottomleft, size=(400, 400), box=false)
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
FSS_time


## Time to order for several system sizes N
relaxtime = plot(xlabel=L"N", ylabel="Relaxation Time", legend=:topleft)
plot!([NaN, NaN], [NaN, NaN], c=1, rib=0, label=L"v_0" * "=$(v0sigs[1][1]), σ=$(v0sigs[1][2])")
plot!([NaN, NaN], [NaN, NaN], c=2, rib=0, label=L"v_0" * "=$(v0sigs[2][1]), σ=$(v0sigs[2][2])")
# σ = 0
finalPs = Ps_avg[1, :, 30]
finalts = [findfirst(x -> x > 0.95finalPs[i], Ps_avg[1, i, 2:end]) for i in 1:length(Ntargets)]
plot!(Ntargets, times[finalts], m=true, axis=:log, c=1)
# σ = 0.1
finalPs = Ps_avg[2, :, 30]
finalts = [findfirst(x -> x > 0.95finalPs[i], Ps_avg[2, i, 2:end]) for i in 1:length(Ntargets)]
plot!(Ntargets, times[finalts], m=true, axis=:log, c=2)
# fits
plot!(Ntargets, 2.5E-2Ntargets .* log.(Ntargets), c=:black, axis=:log, label=L"N\,\log(N)")
plot!(Ntargets, 1.5E-1Ntargets, c=:black, axis=:log, label=L"N", line=:dash, lw=0.7)

FSS_relax = plot(FSS_time, relaxtime, layout=(1, 2), size=(800, 400))
# savefig(FSS_relax, "figures_paper/SM/FSS_relaxtime.svg")
# savefig(FSS_relax, "figures_paper/SM/FSS_relaxtime.pdf")


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
plot!(times[2:700], x -> exp(0.5 * lambertw(2π * x / mu)), c=:black)
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


