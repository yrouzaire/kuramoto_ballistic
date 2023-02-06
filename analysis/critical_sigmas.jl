cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl");
const global R0 = 1
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&
## Critical Sigmas
filename = "data/critical_sigma.jld2"
@load filename critical_sigmas_fusion times sigmas v0s rhos tmax R #runtimes
histogram(runtimes / 3600 / 24, bins=20)

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

p2 = plot(uaxis=:log, legend=:topleft, xlabel=L"\sqrt{v_0}", ylabel=L"\sigma_c")
for i in each(rhos)
    plot!(sqrt.(v0s), critical_sigmas_avg[:, i], rib=0, m=true)
end
plot!(sqrt.(v0s), x -> 0.1sqrt(1) * sqrt(1 + (x - 0.24)))
plot!(sqrt.(v0s[11:end]), x -> 0.5max(x - 0.3, 0), c=:black, l=:dash, label="Slope 1")
p2

##

plot(p1, p2, layout=(1, 2), size=(800, 400))
# savefig("figures/critical_sigmas_different_rhos.png")

## Critical velocity
filename = "data/critical_velocity.jld2"
@load filename runtimes critical_velocity_fusion times v0s rhos tmax R
histogram(runtimes / 3600, bins=20)
critical_velocity_avg = nanmean(critical_velocity_fusion, 2)[:, 1] 

rhoc = 4.51 / π
nc = 4.51
p = plot(uaxis=:log, legend=:topleft, xlabel=L"ρ", ylabel=L"v_c")
plot!((rhos), (critical_velocity_avg) .- (v0s[1]), m=true)
plot!(rhos, x -> 1E-1 / π / R0 / x * (nc - π * R0^2 * x), c=:black, l=:solid)
# plot!(rhos, x -> 0.08(rhoc - x) / x, c=:black, l=:solid)
# plot!(rhos, x -> 0.15(rhoc - x)^2 / x^2, c=:black, l=:dash)
# plot!(rhos, x -> 800exp(-10x), c=:black, l=:dot)
