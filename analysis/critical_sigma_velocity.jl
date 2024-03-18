cd("/Users/yrouzaire/Documents/Recherche/GitHub/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl");
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&


rhoc = 4.51 / π


## ---------------------- Critical Velocity N=1E3 ---------------------- ##
## ---------------------- Critical Velocity N=1E3 ---------------------- ##
## ---------------------- Critical Velocity N=1E3 ---------------------- ##
## ---------------------- Critical Velocity N=1E3 ---------------------- ##
p = plot(uaxis=:log, legend=:top)

filename = "data/critical_velocity_N1E3_extended2.jld2"
@load filename critical_velocity_fusion times v0s rhos tmax R runtimes
# hrun(runtimes)
critical_velocity_avg = nanmean(critical_velocity_fusion, 2)[:, 1]
critical_velocity_std = nanstd(critical_velocity_fusion, 2)[:, 1]
# plot!((rhos), (critical_velocity_fusion[:, 1]), m=true, rib=0critical_velocity_std, label=L"N = 10^3")
hline!([0], c=:black, l=:solid, lw=0.7)
plot!((rhos), (critical_velocity_avg), m=true, rib=critical_velocity_std, label=L"N = 10^3", c=1)
# plot!(0.5:0.01:1.5, x -> max(0, 0.2(1.435 - x) / x), c=:black, l=:dash)
# plot!(0.5:0.01:1.5, x -> max(0, 0.2(1.16 - x) / x), c=:black, l=:solid)
# ylims!(0, 0.02)

critical_velocity_avg

ylims!(-0.02, 0.5)
annotate!((0.95, 0.07), text(L"ρ", 17, :right, :bottom))
annotate!((0.26, 0.88), text(L"v_c", 17, :right, :bottom))

filename = "data/critical_velocity_N1E4_extended2.jld2"
@load filename critical_velocity_fusion times v0s rhos tmax R runtimes comments
# hrun(runtimes)
comments
tmax
critical_velocity_avg = nanmean(critical_velocity_fusion, 2)[:, 1]
critical_velocity_std = nanstd(critical_velocity_fusion, 2)[:, 1]
# plot!((rhos), (critical_velocity_fusion[:, 1]), m=true, rib=0critical_velocity_std, label=L"N = 10^3")
plot!((rhos), (critical_velocity_avg), m=true, rib=critical_velocity_std, label=L"N = 10^4", c=2)
plot!(0.5:0.01:1.6, x -> max(0, 0.2(1.435 - x) / x), c=:black, l=:solid)#, label=L"v_c(ρ)"*" (theory)")
# ylims!(0, 0.02)
xlims!(0.45, 1.65)
vline!([1.435], c=:grey, l=:dash, lw=0.7)
annotate!(1.5, 0.23, text(L"ρ_{perco}=1.44", 10, :center, :center, 90.0, :grey))

# savefig("figures_paper/fig1/critical_velocity_rho.svg")

## -------------- Critical Sigmas Output Reconstitution -------------- ##
## -------------- Critical Sigmas Output Reconstitution -------------- ##
## -------------- Critical Sigmas Output Reconstitution -------------- ##
## -------------- Critical Sigmas Output Reconstitution -------------- ##

# v0s = round.(collect(0.02:0.02:0.5) .^ 2, digits=3)
# rhos = [1, 1.1, 1.2, 1.3, 1.435, 1.5, 1.6, 1.8, 2]
# Déclaration des valeurs de rho et v0
rhos = [1, 1.1, 1.2, 1.3, 1.435, 1.5, 1.6, 1.8, 2]
v0s = [0.0, 0.002, 0.004, 0.006, 0.01, 0.014, 0.02, 0.026, 0.032, 0.04, 0.048, 0.058, 0.068, 0.078, 0.09, 0.102, 0.116, 0.13, 0.144, 0.16, 0.176, 0.194, 0.212, 0.23, 0.25]
R = 40

# Création de la matrice critical_sigmas_fusion
critical_sigmas_fusion = 0.3 * ones(length(rhos), length(v0s), R)
runtimes = zeros(R)

# Chemin du fichier à lire
fichier_root = pwd() * "/data/datacrit/crit_sig.o6848339"

# Fonction pour extraire la valeur numérique après " : " dans une ligne
function extraire_valeur_ligne(ligne::String)
    if occursin("Critical sigma", ligne)
        matchh = match(r" : (0.\d+)", ligne)
        if matchh !== nothing
            return parse(Float32, split(matchh.match, " : ")[end])
        end
    end
    return nothing
end

function extraire_runtime(ligne::String)
    if occursin("Runtime :", ligne)
        matchh = match(r" : (\d+.\d+)", ligne) # will find the runtime in hours
        if matchh !== nothing
            return 3600 * parse(Float32, split(matchh.match, " : ")[end]) # now in seconds
        end
    end
    return nothing
end



for rr in 1:R
    fichier = pwd() * "/data/datacrit/crit_sig.o68483" * string(38 + rr)

    # Lire le fichier ligne par ligne
    open(fichier, "r") do file
        for (i, ligne) in enumerate(eachline(file))
            valeur_xxx = extraire_valeur_ligne(ligne)
            if !isnothing(valeur_xxx)
                # Trouver les indices rho et v0 correspondants
                rho_index = findfirst(x -> occursin("ρ = $x", ligne), rhos)
                v0_index = findlast(x -> occursin("v0 = $x", ligne), v0s)

                if !isnothing(rho_index) && !isnothing(v0_index)
                    critical_sigmas_fusion[rho_index, v0_index, rr] = valeur_xxx
                end
            end

            valeur_runtime = extraire_runtime(ligne)
            if !isnothing(valeur_runtime)
                runtimes[rr] = valeur_runtime
            end
        end
    end
end
critical_sigmas_fusion = round.(critical_sigmas_fusion, digits=3)
critical_sigmas_fusion_avg = nanmean(critical_sigmas_fusion, 3)[:, :, 1]
critical_sigmas_fusion_std = nanstd(critical_sigmas_fusion, 3)[:, :, 1]
hrun(runtimes)

v0s
critical_sigmas_fusion[:, 2, 1]

# Plots of the critical sigmas
p = plot(uaxis=:log, legend=false)#, xlabel=L"\sqrt{v_0}", ylabel=L"\sigma_c")
for i in each(rhos)
    plot!(v0s[2:end] .^ 0.5, critical_sigmas_fusion_avg[i, 2:end],
        label="ρ = $(rhos[i])", rib=(i == 0) * critical_sigmas_fusion_std[i, 2:end], m=:circle, ms=2)
end
p
annotate!((0.05, 0.99), text(L"\sigma_c", 17, :left, :top))
annotate!((0.98, 0.01), text(L"\sqrt{v_0}", 17, :right, :bottom))
# savefig("figures/critical_sigmas.svg")
# savefig("figures/critical_sigmas.png")

