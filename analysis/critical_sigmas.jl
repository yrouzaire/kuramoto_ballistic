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
p = plot(uaxis=:log, legend=:topright)#, xlabel=L"ρ", ylabel=L"v_c")
plot!((rhos), (critical_velocity_avg) .- 0(v0s[1]), m=true, rib=critical_velocity_std, label=L"N = 10^3")
# R0 = 1
# plot!(rhos, x -> 3E-2 / π / R0 / x^2 * (nc^2 - π * R0^2 * x^2), c=:black, l=:solid)
plot!(1:0.01:1.5, x -> max(0,0.08(rhoc - x) / x), c=:black, l=:solid)
# plot!(rhos, x -> 0.6(rhoc - x)^3 / x^3, c=:black, l=:dash)
# plot!(rhos, x -> 5E-2/x^2, c=:black, l=:dot)
# plot!(rhos, x->1E-1/x^2)
# scatter!([1,1.1,1.2], [0.06,0.025,0.018], m=true, ms=3, c=:black,label="From crit. σ simu")
# 0.23^2
# 0.16^2
# 0.135^2

#
filename = "data/critical_velocity_N1E4.jld2"
@load filename critical_velocities_fusion times v0s rhos tmax #runtimes
critical_velocities_avg = nanmean(critical_velocities_fusion, 2)[:, 1]
critical_velocities_std = nanstd(critical_velocities_fusion, 2)[:, 1]

# p = plot(uaxis=:log, legend=false, xlabel=L"ρ", ylabel=L"v_c")
plot!((rhos), (critical_velocities_avg) .- 0(v0s[1]), m=true, rib=critical_velocities_std, label=L"N = 10^4")
plot!(1:0.01:1.5, x -> 0.3*max((rhoc - x)/x,0) , c=:black, l=:solid)
# xlims!(0.96,1.5)
# plot!(rhos, x -> 0.6(rhoc - x)^3 / x^3, c=:black, l=:dash)
# scatter!([1, 1.1, 1.2], [0.06, 0.025, 0.018], m=true, ms=3, c=:black, label="From crit. σ simu")
annotate!((0.05, 0.99), text(L"v_c", 17, :left, :top))
annotate!((0.97, 0.05), text(L"\rho", 17, :right, :bottom))
# savefig("figures/critical_velocity.svg")
# savefig("figures/critical_velocity.png")




## -------------- Critical Sigmas Output Reconstitution with ChatGPT -------------- ##
## -------------- Critical Sigmas Output Reconstitution with ChatGPT -------------- ##
## -------------- Critical Sigmas Output Reconstitution with ChatGPT -------------- ##
## -------------- Critical Sigmas Output Reconstitution with ChatGPT -------------- ##
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
            return parse(Float32,split(matchh.match, " : ")[end])
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
    fichier = pwd() * "/data/datacrit/crit_sig.o68483"*string(38+rr)

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
critical_sigmas_fusion = round.(critical_sigmas_fusion,digits=3)
critical_sigmas_fusion_avg = nanmean(critical_sigmas_fusion, 3)[:, :, 1]
critical_sigmas_fusion_std = nanstd(critical_sigmas_fusion, 3)[:, :, 1]
hrun(runtimes)

v0s
critical_sigmas_fusion[:,2,1]

# Plots of the critical sigmas
p=plot(uaxis=:log, legend=false)#, xlabel=L"\sqrt{v_0}", ylabel=L"\sigma_c")
for i in each(rhos)
    plot!(v0s[2:end].^0.5, critical_sigmas_fusion_avg[i, 2:end], 
        label="ρ = $(rhos[i])", rib = (i==0) * critical_sigmas_fusion_std[i, 2:end], m=:circle, ms=2)
end
p
annotate!((0.05, 0.99), text(L"\sigma_c", 17, :left, :top))
annotate!((0.98, 0.01), text(L"\sqrt{v_0}", 17, :right, :bottom))
# savefig("figures/critical_sigmas.svg")
# savefig("figures/critical_sigmas.png")

## Scatter the critical velocities
filename = "data/critical_velocity_N1E4.jld2"
@load filename critical_velocities_fusion times v0s rhos tmax #runtimes
critical_velocities_avg = nanmean(critical_velocities_fusion, 2)[:, 1]
for i in 1:2:length(rhos)
    scatter!((sqrt(critical_velocities_avg[i]), 0), c=Int((i+1)/2), m=:circle, ms=5, label="ρ = $(rhos[i])")
end
p

