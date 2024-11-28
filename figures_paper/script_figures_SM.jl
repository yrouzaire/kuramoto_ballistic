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
    size=(450, 400),
    # Text=(font=:bold_italic, space=:relative, fontsize=24),
)
set_theme!(custom_theme)


&



## ---------------- R(t) Comparison different distributions ---------------- ##
## ---------------- R(t) Comparison different distributions ---------------- ##
## ---------------- R(t) Comparison different distributions ---------------- ##
## ---------------- R(t) Comparison different distributions ---------------- ##


distributions_type = ["gaussian", "uniform", "laplace", "truncated_cauchy"]

L = round(Int, sqrt(Ntarget))
all_rr_reverse = NaN * zeros(length(distributions_type), 3, 5, 1, 601, 400) # 3 v0s, 5 sigmas, 1 Temperature, 601 times, 400 realisations
all_rr_ = NaN * zeros(length(distributions_type), 3, 5, 1, 601, 400)

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

## R(t) 
cols = cgrad([:blue, :orange, :red, :green,]);
R0 = 28
lss = [:solid, :dash, :dot]
using LambertW

fig=Figure(size=(800, 400))
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
lines!(times[2:140], x -> exp(0.5 * lambertw(2π * x / mu))/R0, color=:black, linewidth=3)


text!(ax, 1, 1, text="(a)",
    fontsize=24, font=:bold_italic,
    align=(:right, :top), offset=(-8, -4), space=:relative) # in the coordinate system of the axis)

text!(ax_star, 0, 1, text="(b)",
    fontsize=24, font=:bold_italic,
    align=(:left, :top), offset=(+8, -4), space=:relative) # in the coordinate system of the axis)


distributions_type
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



fig

save("figures_paper/fig_Rt_distributions.pdf", fig)
save("figures_paper/fig_Rt_distributions.png", fig)

