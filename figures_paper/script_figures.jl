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
hrun(runtimes)



filename = "data/phase_space_rho_sig_v0_N1E3_tmax2500.jld2"
@load filename Ps Cs ns runtimes Ts Ns v0s rhos sigmas times_log tmax comments R

Ps_avg = nanmean(Ps, 8) # N rho T v0 sigma init t R
Ps_std = nanstd(Ps, 8) # N rho T v0 sigma init t R
ns_avg = nanmean(ns, 8) # N rho T v0 sigma init t R
ns_std = nanstd(ns, 8) # N rho T v0 sigma init t R


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

p_phase_space_n_rho1 = heatmap(v0s[2:end], sigmas, log10.(ns_avg[1, 1, 1, 2:end, :, 1, end, 1]' .+ 1),
    xaxis=:log, c=reverse(cgrad([:red, :orange, :green])),
    size=(500, 400), xlabel=L"v_0", ylabel="σ", clims=(0, 1.5),
    colorbartitle="n", colorbar=true, colorbar_titlefont=font(12), colorbar_titlefontrotation=90)
xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])
p_phase_space_n_rho1
vline!([1], c=:black, l=:dash, lw=0.8)

##
p_phase_space_n_rho1 = heatmap(v0s[2:end], sigmas, log10.(ns_avg[1, 1, 1, 2:end, :, 1, end, 1]' .+ 1),
    xaxis=:log, c=reverse(cgrad([:red, :orange, :green])),
    size=(400, 400), xlabel=L"v_0", ylabel="σ", clims=(0, 1.5),
    colorbartitle="n", colorbar=false, colorbar_titlefont=font(12), colorbar_titlefontrotation=90)
xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])
p_phase_space_n_rho1

p_phase_space_nstdavg_rho1 = heatmap(v0s[2:end], sigmas, (ns_std[1, 1, 1, 2:end, :, 1, end, 1] ./ (ns_avg[1, 1, 1, 2:end, :, 1, end, 1] .+ 1))',
    xaxis=:log,
    size=(400, 400), xlabel=L"v_0", ylabel="σ", 
    colorbartitle="n", colorbar=true, colorbar_titlefont=font(12), colorbar_titlefontrotation=90)
xticks!([1E-3, 1E-2, 1E-1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])
p_phase_space_nstdavg_rho1

plot(p_phase_space_n_rho1, p_phase_space_nstdavg_rho1, layout=(1, 2), size=(800, 400))


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



## ----------------  Critical Velocities ---------------- ##
## ----------------  Critical Velocities ---------------- ##
## ----------------  Critical Velocities ---------------- ##
## ----------------  Critical Velocities ---------------- ##


filename = "data/critical_velocity_N1E4_extended2.jld2"
@load filename Ntarget critical_velocity_fusion rhos v0s times tmax T seuil comments rhoc runtimes
critical_velocity_fusion_avg_sigma0 = nanmean(critical_velocity_fusion, 2)[:, 1]
critical_velocity_fusion_std_sigma0 = nanstd(critical_velocity_fusion, 2)[:, 1]
rhos0 = rhos

filename = "data/critical_velocity_sigma0.1_rhos_N1E4.jld2"
@load filename Ntarget critical_velocity_fusion critical_velocity_fusion_avg aspect_ratio rhos v0s sigmas times tmax T seuil comments rhoc runtimes
critical_velocity_fusion_avg_sigma01 = critical_velocity_fusion_avg
critical_velocity_fusion_std_sigma01 = nanstd(critical_velocity_fusion, 3)[:, 1, 1]
rhos01 = rhos

filename = "data/critical_velocity_sigma0.1_complement_rhos_N1E4.jld2"
@load filename Ntarget critical_velocity_fusion critical_velocity_fusion_avg aspect_ratio rhos v0s sigmas times tmax T seuil comments rhoc runtimes
critical_velocity_fusion_avg_sigma01_complement = critical_velocity_fusion_avg
critical_velocity_fusion_std_sigma01_complement = nanstd(critical_velocity_fusion, 3)[:, 1, 1]
rhos01_complement = rhos

filename = "data/critical_velocity_sigma0.2_rhos_N1E4.jld2"
@load filename Ntarget critical_velocity_fusion critical_velocity_fusion_avg aspect_ratio rhos v0s sigmas times tmax T seuil comments rhoc runtimes
critical_velocity_fusion_avg_sigma02 = critical_velocity_fusion_avg
critical_velocity_fusion_std_sigma02 = nanstd(critical_velocity_fusion, 3)[:, 1, 1]
rhos02 = rhos

p = plot(legend_title=L"σ", legend=:topright, size=(400, 400))
plot!(rhos0, critical_velocity_fusion_avg_sigma0, rib=critical_velocity_fusion_std_sigma0, m=true, c=1)
plot!(rhos01, critical_velocity_fusion_avg_sigma01, rib=critical_velocity_fusion_std_sigma01, m=true, c=2)
plot!(rhos01_complement, critical_velocity_fusion_avg_sigma01_complement, rib=critical_velocity_fusion_std_sigma01_complement, m=true, c=2)
plot!([rhos01[end], rhos01_complement[1]], [critical_velocity_fusion_avg_sigma01[end], critical_velocity_fusion_avg_sigma01_complement[1]], c=2, rib=critical_velocity_fusion_std_sigma01_complement)
plot!(rhos02, critical_velocity_fusion_avg_sigma02, rib=critical_velocity_fusion_std_sigma02, m=true, c=3)
hline!([0], c=:black, l=:solid, lw=0.7)
vline!([1.435], c=:grey, l=:dash, lw=0.7)
annotate!(1.35, 0.58, text(L"ρ_{perco}", 10, :center, :center, 90.0, :grey))
annotate!((0.95, 0.07), text(L"ρ", 17, :right, :bottom))
annotate!((0.26, 0.88), text(L"v_c", 17, :right, :bottom))
plot!([NaN,NaN], [NaN,NaN], c=1, l=:solid, rib=0, label="0")
plot!([NaN,NaN], [NaN,NaN], c=2, l=:solid, rib=0, label="0.1")
plot!([NaN,NaN], [NaN,NaN], c=3, l=:solid, rib=0, label="0.2")

xxx = minimum(rhos02):0.01:maximum(rhos02)
plot!(xxx, x -> max(0, 0.2(1.435 - x)) / x, c=:black, l=:dash)
plot!(xxx, x -> max(0, 0.2(2 - x)) / x, c=:black, l=:dash)
plot!(xxx, x -> max(0, 0.3(2.4 - x)) / x, c=:black, l=:dash)
# savefig(p, "figures_paper/critical_velocities.svg")
p


## Inset with the critical velocity dependance on sigma
filename = "data/critical_sigma.jld2"
@load filename critical_sigmas_fusion times sigmas v0s rhos tmax R #runtimes
critical_sigmas_avg = nanmean(critical_sigmas_fusion, 3)[:, :, 1]
rhos
v0s
sigmas_to_consider = [0, 0.1, 0.2]
crit_vel_sigma = zeros(length(rhos), length(sigmas_to_consider))
for i in each(rhos)
    for j in each(sigmas_to_consider)
        println(i, " ", j)
        ind = findfirst(x -> x > sigmas_to_consider[j], critical_sigmas_avg[:, j])
        x1 = v0s[ind]
        y1 = critical_sigmas_avg[i, j]
        if ind < length(v0s)
            x2 = v0s[ind + 1]
            y2 = critical_sigmas_avg[i, ind + 1]

            tmp = mean([x1, x2])
            # tmp = mean([sqrt(a), sqrt(b)]) ^2
            # tmp = mean([sqrt(a), sqrt(b)]) ^2
            # a = (y2 - y1) / (x2 - x1)
            # b = y1 - a * x1
            # crit_vel_sigma[i, j] = (sigmas[j] - b) / a
        else
            x2 = v0s[ind-1]
            y2 = critical_sigmas_avg[i, ind-1]

            tmp = mean([x1, x2])
            # tmp = mean([sqrt(a), sqrt(b)]) ^2
            # tmp = mean([sqrt(a), sqrt(b)]) ^2
            # a = (y2 - y1) / (x2 - x1)
            # b = y1 - a * x1
            # crit_vel_sigma[i, j] = (sigmas[j] - b) / a
        end

        crit_vel_sigma[i, j] = tmp
    end
end


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
L = sqrt(Ntarget / rho)
p1 = plot(xscale=:log10, yscale=:log10, legend=:top,
    yticks=([1E-5, 1E-4, 1E-3, 1E-2], [L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}"]),
    xticks=([1, 10, 100, 1000, 1E4], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}"]))
for i in each(v0sigs_horizontal)
        couleur = cols_P[Ps_avg_horizontal[i,end]]
    v00 = sqrt(v0sigs_horizontal[i][1])
    # plot!(times, remove_negative(ns_avg_horizontal[i, :] / L^2)*v00, c=couleur, rib=0, m=:circle, ms=3, line=true)
    plot!(times, remove_negative(ns_avg_horizontal[i, :] / L^2), c=couleur, rib=0, m=false, ms=3, lw=1.5)
end
times
ylims!(1E-5, 1E-1)
xlims!(2E-1, 1E5)
annotate!((0.15, 0.93), text(L"n/L^2", 15, :center, :black))
annotate!((0.94, 0.1), text(L"t", 15, :top, :black))
plot!(times[5:end-8], x -> 1.3E-2log(10x) / x, line=:dash, c=:black, label=L"\log(t)/t")
annotate!((0.92, 0.93), text("(b)", 15))
xticks!([1, 10, 100, 1000, 1E4, 1E5], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}", L"10^{5}"])
p1

##
rr = 0:round(Int, L / 2)
p3 = plot(axis=:log, ylims=(1E-1, 1.3), legend=:top)
for i in 1:1:length(v0sigs_horizontal)
        couleur = cols_P[Ps_avg_horizontal[i,end]]
    plot!(rr[2:end], remove_negative(Cs_avg_horizontal[i, end])[2:end], c=couleur, rib=0, m=false, lw=1.5, ms=3)
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
inset_xi = plot(xaxis=:log, legend=false, size=(250,200))#:topright)
hline!([0.5], c=:black, lw=0.8)
for i in each(v0sigs_horizontal)
        couleur = cols_P[Ps_avg_horizontal[i,end]]
    data = remove_negative(xis_avg_horizontal[i, 2:end-4] .* sqrt.(ns_avg_horizontal[i, 2:end-4])) / L
    plot!(times[2:end-4], data, c=couleur, m=false, ms=2, lw=0.99)
end
ylims!(0.36,0.64)
yticks!([0.4, 0.5, 0.6])
annotate!((0.25, 0.89), text(L"ξ\,\sqrt{n}/L", 12, :center))
annotate!((0.93, 0.08), text(L"t", 12, :center))
xticks!([1, 10, 100, 1000, 1E4], [L"10^{0}", "", L"10^{2}", "", L"10^{4}"])
inset_xi
# savefig(inset_xi, "figures_paper/fig2/inset_xi.svg")

##
p2 = plot(axis=:log,legend=:top)
for i in each(v0sigs_horizontal)
    couleur = cols_P[Ps_avg_horizontal[i, end]]
    v00 = sqrt(v0sigs_horizontal[i][1])
    # plot!(times, remove_negative(xis_avg_horizontal[i, :]) / L / v00, c=couleur, rib=0, m=:circle, ms=3, line=true)
    plot!(times, remove_negative(xis_avg_horizontal[i, :])/L, c=couleur, rib=0, m=:false, ms=3, lw=1.5)
end
p2
yticks!([0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], ["0.03", "", "0.05", "", "0.07", "", "", "0.1", "0.2", "0.3", "", "0.5", "", "0.7"])
plot!(times[5:end-13], x -> 4.3 / L * sqrt(x / log(8x)), line=:dash, c=:black, label=L"\sqrt{t/\log(t)}")
annotate!((0.09, 0.93), text(L"\xi/L", 15, :center, :black))
annotate!((0.94, 0.1), text(L"t", 15, :top, :black))
annotate!((.92, .93), text("(a)", 15))
xlims!(0.2, 14E4)
ylims!(0.022, 0.7)
xticks!([1, 10, 100, 1000, 1E4, 1E5], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}", L"10^{5}"])


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

p4 = plot(axis=:log, legend=:top)
for i in each(v0sigs)
    couleur = cols_P[Ps_avg[i, end, end]]
    # couleur = cols_P[Ps_avg_horizontal[i, end]]
    plot!(Ntargets, Ps_avg[i, :, end], rib=0Ps_std[i, :, end], m=false, c=couleur, lw=1.5, ms=3)
end
ylims!(0.03, 1.4)
yticks!([0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], ["0.03", "0.04", "0.05", "", "0.07", "", "", "0.1", "0.2", "0.3", "0.4", "0.5", "", "0.7", "", "", "1"])
plot!(Ntargets[2:end-1], x -> 3.7x^-0.5, c=:black, line=:dash)
plot!(Ntargets[1:end], x -> 1.09x^-0.025, c=:black, line=:dot)
p4
annotate!((0.09, 0.06), text("(d)", 15))
annotate!((0.345, 0.35), text(L"\sim N^{-1/2}", 10))
annotate!((0.65, 0.935), text(L"\sim N^{-0.025}", 10))

annotate!((0.08, 0.9), text(L"P", 15, :center, :bottom, :black))
annotate!((0.9, 0), text(L"N", 15, :center, :bottom, :black))


##
fig2=plot(p2, p1, p3, p4, layout=(1, 4), size=(1600, 400));
# savefig(fig2, "figures_paper/fig2/fig2.svg")

## Critical Density


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

phistogram = plot(size=(200, 300), legend=(0.43, 0.4), legend_title=L"v_0")
histogram!(log10.(all_times_collision[1, 1, ind_T, :]), bins=20, c=1, lw=0.2, label=L"0.5")
histogram!(log10.(all_times_collision[5, 1, ind_T, :]), bins=10, c=5, lw=0.2, label=L"2.5")
histogram!(log10.(all_times_collision[end, 1, ind_T, :]) , bins=5, c=10, lw=0.2, label=L"5")

# ## Histograms mean_annihilation_time
# p = plot(legend=:outerright)
# for i in 1:1:length(v0s)
#     v0 = v0s[i]
#     lab = "v0 = $v0"
#     data = log10.(all_times_collision[i, 1, ind_T, :]) * sqrt(v0)
#     # histogram!(data, bins=20, lw=0.2, label=L"5")
    
#     # fit the data
#     h = fit(Histogram, data, nbins=10)
#     h = normalize(h, mode=:density)
#     plot!(h.edges, h.weights, label = lab, rib=0)


# end
# p




ylims!(0, 150)
xticks!(1:4, [L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}"])
xlims!(1, 4)
annotate!((0.25, 0.88), text(L"\mathbb{P}(\tau)", 13, :center, :bottom, :black))
annotate!((0.9, 0.02), text(L"\tau", 13, :center, :bottom, :black))

savefig(phistogram, "figures_paper/Rt_inset.svg")
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

## ---------------- Plots Proba spinwaves ---------------- ##
## ---------------- Plots Proba spinwaves ---------------- ##
## ---------------- Plots Proba spinwaves ---------------- ##
## ---------------- Plots Proba spinwaves ---------------- ##

filename = "data/proba_spinwaves_scan_phase_space_N4000.jld2"
# filename = "data/proba_spinwaves.jld2"
@load filename R_per_core Rtot R all_nb_detected_spinwave all_times_detected_spinwave all_Ps_detected_spinwave all_thetas_detected_spinwave all_pos_detected_spinwave sigmas v0s tmax times p_threshold init_pos init_theta Ntarget rho T aspect_ratio runtimes
proba_spinwave = all_nb_detected_spinwave / Rtot
hrun(runtimes)

Ntarget
Rtot
tmax
sigmas
v0s

#
proba_spinwave = all_nb_detected_spinwave / Rtot
proba_spinwave = all_nb_detected_spinwave / Rtot + all_nb_detected_spinwave_complement / Rtot_complement


colss = cgrad([:black, :red, :orange, :gold])
plot(xaxis=:log, size=(470, 400))
heatmap!(v0s, sigmas, 100 * proba_spinwave', clims=(0, 10),
    c=colss, colorbartitle=L"\mathbb{P}\," * "(TPS) [%]")
plot!(v0s, x -> 1 / 2 * max(0, x - (0.23)^2), c=:white, lw=0.8)
ylims!(-0.001, 0.4)
xlims!(minimum(v0s), 1.18maximum(v0s))
xticks!([1E-2, 1E-1, 1E-0], [L"10^{-2}", L"10^{-1}", L"10^{0}"])
# xticks!([1E-2, 2E-2, 3E-2, 4E-2, 5E-2, 6E-2, 7E-2, 8E-2, 9E-2, 1E-1, 2E-1, 3E-1, 4E-1, 5E-1, 6E-1, 7E-1, 8E-1, 9E-1, 1E-0, 2, 3], 
#     [L"10^{-2}", "", "", "", "", "", "", "", "", L"10^{-1}", "", "", "", "", "", "", "", "", L"10^{0}", "", ""])
# critere : is_green_region = sigma < 1 / 2 * max(0, sqrt(v0) - 0.25) , at rho=1
annotate!((0.05, 0.98), text(L"\sigma", 17, :left, :top, :white))
annotate!((0.96, 0.03), text(L"v_0", 17, :right, :bottom, :white))