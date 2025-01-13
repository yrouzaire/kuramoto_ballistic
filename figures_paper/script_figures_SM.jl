cd("/Users/yrouzaire/Documents/Recherche/GitHub/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl");
using CairoMakie, ColorSchemes, LaTeXStrings

CairoMakie.activate!()

custom_theme = Theme(
    Axis=(
        xgridvisible=false, ygridvisible=false,
        xtickalign=1, ytickalign=1,
        xticklabelpad=2, yticklabelpad=4
    ),
    fontsize=24,
    # Text=(font=:bold_italic, space=:relative, fontsize=24),
)
set_theme!(custom_theme)
cols_py = cgrad(ColorSchemes.viridis.colors[1:end])

&

## ---------------- Comparison T=0 or T>0 ---------------- ##
## ---------------- Comparison T=0 or T>0 ---------------- ##
## ---------------- Comparison T=0 or T>0 ---------------- ##
## ---------------- Comparison T=0 or T>0 ---------------- ##
cols_P = cgrad([:red, :orange, :green])

fig = Figure(size=(950, 800))
ax_Rt = Axis(fig[1, 1], xlabel=L"t", ylabel=L"R/R_0", xscale=log10,
    # xticks=([1, 10, 100, 1000, 1E4], [L"10^0", L"10^1", L"10^2", L"10^3", L"10^4"]),
    # yticks=(collect(0:0.2:1), latexstring.(0:0.2:1)),
    # limits=(1, 1E4, 0, 1)
)
ax_xi = Axis(fig[1, 2], xlabel=L"t", ylabel=L"\xi/L", xscale=log10, yscale=log10,
    # xticks=([1, 10, 100, 1000, 1E4], [L"10^0", L"10^1", L"10^2", L"10^3", L"10^4"]),
    # yticks=(collect(0:0.2:1), latexstring.(0:0.2:1)),
    # limits=(1, 1E4, 0, 1)
)
ax_n = Axis(fig[2,1], xlabel=L"t", ylabel=L"n/L^2", xscale=log10, yscale=log10,
    # xticks=([1, 10, 100, 1000, 1E4], [L"10^0", L"10^1", L"10^2", L"10^3", L"10^4"]),
    # yticks=(collect(0:0.2:1), latexstring.(0:0.2:1)),
    # limits=(1, 1E4, 0, 1)
)
ax_Cr = Axis(fig[2,2], xlabel=L"r", ylabel=L"C(r)", xscale=log10, yscale=log10,
limits=(nothing, nothing, 1E-1, 1.15), 
    xticks=(vcat(1:10, 20:10:50), string.([1, 2, 3, "", 5, "", "", "", "", 10, 20, 30, "", 50])),
)


# R(t)
filename = "data/mobility_defects_sigma_v0_distribution_sigmas_gaussian_T0.jld2"
@load filename sigmas v0s Ts all_xy_pos all_xy_neg all_rr all_times_collision R_per_core Rtot params_init Ntarget R0 q init_theta init_pos aspect_ratio times tmax comments rhoc runtimes

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


for i in 1:length(v0s)
    data_T0 = remove_negative(all_rr_avg[i, 1, 1, 2:end])/28
    lines!(ax_Rt, times[2:end], data_T0, color = cols_py[i/length(v0s)])

    data_Tfinite = remove_negative(all_rr_avg[i, 1, 2, 2:end])/28
    lines!(ax_Rt, times[2:end], data_Tfinite, color = cols_py[i/length(v0s)], linestyle=:dash)
end
arrows!(ax_Rt, [1300], [5/28], [-1280], [0], linecolor=:black, arrowcolor=:black)
text!(ax_Rt, 13, 5/28, text=L"v_0", color=:black,
    fontsize=22, font=:bold_italic,
    align=(:center, :center)) # in the coordinate system of the axis)


# Coarsening 
filename = "data/nature_phase_transition_horizontal_T0.jld2"
@load filename Ntarget v0sigs rho params_init T Ps Cs ns xis aspect_ratio times tmax comments rhoc runtimes 

L = 100
Ps_avg = mean(Ps, dims=3)[:,:,1]
Cs_avg = mean(Cs, dims=3)[:,:,1]
ns_avg = mean(ns, dims=3)[:,:,1]
xis_avg = mean(xis, dims=3)[:,:,1]

for i in 1:length(v0sigs)
    data = xis_avg[i, :]/L
    lines!(ax_xi, times, data, color=cols_P[Ps_avg[i, end]])
end
lines!(ax_xi, times[5:end-5], x -> 2.5e-2 * sqrt(x/log(x)), color=:black, linestyle=:dash)

for i in 1:length(v0sigs)
    data = ns_avg[i, :]/L^2
    lines!(ax_n, times, data, color=cols_P[Ps_avg[i, end]])
end
lines!(ax_n, times[5:end], x->2.5e-2*log(x)/x, color=:black, linestyle=:dash)

rr = 1:51
for i in 1:length(v0sigs)
    # data = remove_negative(Cs_avg[i, end] .- Ps_avg[i, end]^2)
    data = remove_negative(Cs_avg[i, end])
    lines!(ax_Cr, data, color=cols_P[Ps_avg[i, end]])
end
lines!(ax_Cr, 1:51, x->x^(-1/4), color=:black, linestyle=:dash)

subgrid = GridLayout(fig[:, end+1], tellheight=false)
Label(subgrid[1, 1],L"P_{t\to\infty}", fontsize=40)

Colorbar(subgrid[2,1],#, label=L"P_{SS}", #label=L"P_{t\to\infty}",
    ticklabelsize=20, labelrotation=0, labelsize=40,
    tickalign=1,
    labelpadding=-10, # distance between the colorbar ticks and the label
    # ticks=([0, 0.25, 0.5, 0.75, 1], [L"0", "", L"0.5", "", L"1"]),
    ticks=(0:0.1:1, ["0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1"]),
    colormap=cols_P
)


text!(ax_Rt, 0.1, 0.1, text="(a)",
    fontsize=24, font=:bold_italic,
    align=(:right, :top), offset=(8, 4), space=:relative) # in the coordinate system of the axis)
text!(ax_xi, 0.9, 0.1, text="(b)",
    fontsize=24, font=:bold_italic,
    align=(:left, :top), offset=(-8, 4), space=:relative) # in the coordinate system of the axis)
text!(ax_n, 0.1, 0.1, text="(c)",
    fontsize=24, font=:bold_italic,
    align=(:right, :top), offset=(8, 4), space=:relative) # in the coordinate system of the axis)
text!(ax_Cr, 0.1, 0.1, text="(d)",
    fontsize=24, font=:bold_italic,
    align=(:right, :top), offset=(8, 4), space=:relative) # in the coordinate system of the axis)


resize_to_layout!(fig)
fig

# save("figures_paper/comparison_T0.pdf", fig)
# save("figures_paper/comparison_T0.png", fig)





## ---------------- R(t) Comparison different distributions ---------------- ##
## ---------------- R(t) Comparison different distributions ---------------- ##
## ---------------- R(t) Comparison different distributions ---------------- ##
## ---------------- R(t) Comparison different distributions ---------------- ##

distributions_type = ["gaussian", "uniform", "laplace", "truncated_cauchy"]

L = round(Int, sqrt(Ntarget))
all_rr_reverse = NaN * zeros(length(distributions_type), 3, 5, 1, 601, 400) # 3 v0s, 5 sigmas, 1 Temperature, 601 times, 400 realisations
all_rr_ = NaN * zeros(length(distributions_type), 3, 5, 1, 601, 400)


filename = "data/mobility_defects_sigma_v0_distribution_sigmas_gaussian.jld2"
@load filename sigmas v0s Ts all_xy_pos all_xy_neg all_rr all_times_collision R_per_core Rtot params_init Ntarget R0 q init_theta init_pos aspect_ratio times tmax comments rhoc runtimes

for (ind_distrib, distrib_type) in enumerate(distributions_type)
    filename = "data/mobility_defects_sigma_v0_distribution_sigmas_" * distrib_type * ".jld2"
    @load filename sigmas v0s Ts all_xy_pos all_xy_neg all_rr all_times_collision R_per_core Rtot params_init Ntarget R0 q init_theta init_pos aspect_ratio times tmax comments rhoc runtimes



    for i in each(v0s)
        for j in each(sigmas)
            for k in 1#each(Ts)
                for r in 1:Rtot
                    data = all_rr[i, j, k, r]
                    ll = length(data)
                    all_rr_reverse[ind_distrib, i, j, k, 1:ll, r] = reverse(data)
                    all_rr_[ind_distrib, i, j, k, 1:ll, r] = (data)
                    all_rr_[ind_distrib, i, j, k, 1+ll:end, r] .= 0
                end
            end
        end
    end
end
all_rr_reverse_avg = nanmean(all_rr_reverse, 6)[:,:, :, :, :, 1]
all_rr_avg = nanmean(all_rr_, 6)[:,:, :, :, :, 1]
 
# # R(t) 
cols = cgrad([:blue, :orange, :red, :green,]);
R0 = 28
lss = [:solid, :dash, :dot]
using LambertW


fig=Figure(size=(1300, 400))
ax = Axis(fig[1, 1], xlabel=L"t", title=L"R(t)/R_0", xscale=log10,
    xticks=([1, 10, 100, 1000, 1E4], [L"10^0", L"10^1", L"10^2", L"10^3", L"10^4"]),
    yticks=(collect(0:0.2:1), latexstring.(0:0.2:1)), 
    limits=(nothing, nothing, 0, 1)
)
ax_star = Axis(fig[1, 2], xlabel=L"\sqrt{v_0} t^*", title=L"R(t^*)/(R_0 \, \sqrt{v_0}) ", xscale=log10,
    xticks=([1, 10, 100, 1000, 1E4], [L"10^0", L"10^1", L"10^2", L"10^3", L"10^4"]),
    yticks = (collect(0:0.2:1), latexstring.(0:0.2:1)),
    limits=(nothing, nothing, nothing, 1.2)
)

ax_Cr = Axis(fig[1, 3], xlabel=L"r", title=L"C(r)", xscale=log10, yscale=log10,
    limits=(nothing, nothing, 1E-1, 1.15),
    xticks=(vcat(1:10, 20:10:50), string.([1, 2, 3, "", 5, "", "", "", "", 10, 20, 30, "", 50])),
)

for ind_distrib in 1:length(distributions_type)
    for i in 1:length(v0s)
        for j in 3#length(sigmas)
            for k in 1#:length(Ts)
                if k == 1
                    lab = string(v0s[i])
                else
                    lab = ""
                end
                data = remove_negative(all_rr_avg[ind_distrib, i, j, k, 2:end])/R0
                lines!(ax, times[2:end], data, label=lab, linestyle=lss[k], color = cols[ind_distrib])

                data = 1 / sqrt(v0s[i]) * remove_negative(all_rr_reverse_avg[ind_distrib, i, j, k, 2:end]) / R0
                lines!(ax_star, sqrt(v0s[i]) * times[2:end], data, label=lab, linestyle=lss[k], color = cols[ind_distrib])
                # lines!(times[2:end], all_rr_avg[i, j, k, 2:end])
            end
        end
    end
end
lines!([NaN, NaN], color=:black, linestyle=:solid, label="T=0")
lines!([NaN, NaN], color=:black, linestyle=:dash, label="T=0.1")
mu = 1 / 2
lines!(ax_star, times[2:140], x -> exp(0.5 * lambertw(2π * x / mu))/R0, color=:black, linewidth=3)


text!(ax, 1, 1, text="(a)",
    fontsize=24, font=:bold_italic,
    align=(:right, :top), offset=(-8, -4), space=:relative) # in the coordinate system of the axis)

text!(ax_star, 0, 1, text="(b)",
    fontsize=24, font=:bold_italic,
    align=(:left, :top), offset=(+8, -4), space=:relative) # in the coordinate system of the axis)

text!(ax_Cr, 1, 1, text="(c)",
    fontsize=24, font=:bold_italic,
    align=(:right, :top), offset=(-8, -4), space=:relative) # in the coordinate system of the axis)


text!(ax, 0.07, 0.65, text="Gaussian", color=cols[1],
    fontsize=24, font=:bold_italic,
    align=(:left, :top), offset=(-8, -4), space=:relative) # in the coordinate system of the axis)

text!(ax, 0.07, 0.5, text="Uniform", color=cols[2],
    fontsize=24, font=:bold_italic,
    align=(:left, :top), offset=(-8, -4), space=:relative) # in the coordinate system of the axis)

text!(ax, 0.07, 0.35, text="Laplace", color=cols[3],
    fontsize=24, font=:bold_italic,
    align=(:left, :top), offset=(-8, -4), space=:relative) # in the coordinate system of the axis)

text!(ax, 0.07, 0.2, text="Cauchy (trunc.)", color=cols[4],
    fontsize=24, font=:bold_italic,
    align=(:left, :top), offset=(-8, -4), space=:relative) # in the coordinate system of the axis)


arrows!(ax, [480], [0.35], [-380], [0], linecolor=:black, arrowcolor=:black)
text!(ax, 72, 0.35, text=L"v_0", color=:black,
    fontsize=22, font=:bold_italic,
    align=(:center, :center)) # in the coordinate system of the axis)


arrows!(ax_star, [300], [1.03], [0], [-0.5], linecolor=:black, arrowcolor=:black)
text!(ax_star, 300, 0.47, text=L"v_0", color=:black,
    fontsize=22, font=:bold_italic,
    align=(:center, :center)) # in the coordinate system of the axis)


# Coarsening for Uniform distribution
filename = "data/nature_phase_transition_horizontal_uniform.jld2"
@load filename Ntarget v0sigs rho params_init T Ps Cs ns xis aspect_ratio times tmax comments rhoc runtimes

Cs_avg = mean(Cs, dims=3)[:, :, 1]

rr = 1:51
for i in 1:length(v0sigs)
    # data = remove_negative(Cs_avg[i, end] .- Ps_avg[i, end]^2)
    data = remove_negative(Cs_avg[i, end])
    lines!(ax_Cr, data, color=cols_P[Ps_avg[i, end]])
end
lines!(ax_Cr, 1:51, x -> x^(-1 / 4), color=:black, linestyle=:dash)

subgrid = GridLayout(fig[:, end+1], tellheight=false)
Label(subgrid[1, 1], L"P_{t\to\infty}", fontsize=40)

Colorbar(subgrid[2, 1],#, label=L"P_{SS}", #label=L"P_{t\to\infty}",
    ticklabelsize=20, labelrotation=0, labelsize=40,
    tickalign=1,
    labelpadding=-10, # distance between the colorbar ticks and the label
    # ticks=([0, 0.25, 0.5, 0.75, 1], [L"0", "", L"0.5", "", L"1"]),
    ticks=(0:0.1:1, ["0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1"]),
    colormap=cols_P
)


fig


save("figures_paper/fig_Rt_coarsening_distributions.pdf", fig)
save("figures_paper/fig_Rt_coarsening_distributions.png", fig)
