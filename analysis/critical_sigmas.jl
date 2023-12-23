cd("/Users/yrouzaire/Documents/Recherche/GitHub/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl");
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&
## Critical Sigmas
filename = "data/critical_sigma.jld2"
@load filename critical_sigmas_fusion times sigmas v0s rhos tmax R #runtimes
hrun(runtimes)
critical_sigmas_avg = nanmean(critical_sigmas_fusion, 3)[:, :, 1]
## Plotting
p1 = plot(xaxis=:log, legend=:topleft, xlabel=L"v_0", ylabel=L"\sigma_c")
for i in each(rhos)
    if rhos[i] ≈ 4.51 / π
        lab = L"ρ_c ≈ 1.44"
    else
        lab = "ρ = $(rhos[i])"
    end
    plot!((v0s), critical_sigmas_avg[:, i], label=lab, rib=0, m=true)
end
p1
##
p2 = plot(uaxis=:log, legend=:topleft, xlabel=L"\sqrt{v_0}", ylabel=L"\sigma_c")
for i in each(rhos)
    plot!(sqrt.(v0s), critical_sigmas_avg[:, i], rib=0, m=true, ms=3)
end
# plot!(x -> 0.5sqrt(1.8) * (0 + (x + 0.0)), c=:black, l=:solid)
# plot!(sqrt.(v0s[11:end]), x -> 0.5max(x - 0.3, 0), c=:black, l=:dash, label="Slope 1")
p2
# xlims!(0, 0.08)
# ylims!(0, 0.08)

##

plot(p1, p2, layout=(1, 2), size=(800, 400))
# savefig("figures/critical_sigmas_different_rhos.png")

## Critical velocity
filename = "data/critical_velocity.jld2"
@load filename critical_velocity_fusion times v0s rhos tmax R #runtimes
# histogram(runtimes / 3600, bins=20)
critical_velocity_avg = nanmean(critical_velocity_fusion, 2)[:, 1]
critical_velocity_std = nanstd(critical_velocity_fusion, 2)[:, 1]

rhoc = 4.51 / π
nc = 4.51
p = plot(uaxis=:log, legend=false, xlabel=L"ρ", ylabel=L"v_c")
plot!((rhos), (critical_velocity_avg) .- 0(v0s[1]), m=true, rib=critical_velocity_std)
# R0 = 1
# plot!(rhos, x -> 3E-2 / π / R0 / x^2 * (nc^2 - π * R0^2 * x^2), c=:black, l=:solid)
plot!(rhos, x -> 0.08(rhoc - x) / x, c=:black, l=:solid)
plot!(rhos, x -> 0.6(rhoc - x)^3 / x^3, c=:black, l=:dash)
# plot!(rhos, x -> 5E-2/x^2, c=:black, l=:dot)
# plot!(rhos, x->1E-1/x^2)
scatter!([1,1.1,1.2], [0.06,0.025,0.018], m=true, ms=3, c=:black,label="From crit. σ simu")
# 0.23^2
# 0.16^2
# 0.135^2

##
filename = "data/critical_velocity_N1E4.jld2"
@load filename critical_velocities_fusion times v0s rhos tmax #runtimes
critical_velocities_avg = nanmean(critical_velocities_fusion, 2)[:, 1]
critical_velocities_std = nanstd(critical_velocities_fusion, 2)[:, 1]

# p = plot(uaxis=:log, legend=false, xlabel=L"ρ", ylabel=L"v_c")
plot!((rhos), (critical_velocities_avg) .- 0(v0s[1]), m=true, rib=critical_velocities_std)
plot!(rhos, x -> 0.3(rhoc - x) / x, c=:black, l=:solid)
xlims!(0.96,1.44)
# plot!(rhos, x -> 0.6(rhoc - x)^3 / x^3, c=:black, l=:dash)
# scatter!([1, 1.1, 1.2], [0.06, 0.025, 0.018], m=true, ms=3, c=:black, label="From crit. σ simu")

