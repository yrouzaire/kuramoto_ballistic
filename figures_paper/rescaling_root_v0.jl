cd("/Users/yrouzaire/Documents/Recherche/GitHub/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl");
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

filename = "data/nature_phase_transition_horizontal.jld2"
@load filename v0sigs Ps Cs ns xis rho T Ntarget params_init aspect_ratio times tmax comments rhoc runtimes R
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

## ---------------- Plotting for n/L^2 defect density ---------------- ##
## No rescaling
L = sqrt(Ntarget / rho)
p1 = plot(xscale=:log10, yscale=:log10, legend=true,
    yticks=([1E-5, 1E-4, 1E-3, 1E-2], [L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}"]),
    xticks=([1, 10, 100, 1000, 1E4], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}"]))
for i in each(v0sigs_horizontal)
    couleur = cols_P[Ps_avg_horizontal[i, end]]
    v00 = sqrt(v0sigs_horizontal[i][1])
    plot!(times, remove_negative(ns_avg_horizontal[i, :] / L^2), c=couleur, rib=0, m=:circle, ms=3, line=true)
end
ylims!(1E-5, 1E-1)
annotate!((0.15, 0.93), text(L"n/L^2", 15, :center, :black))
annotate!((094, 0.1), text(L"t", 15, :top, :black))

plot!(times[5:end-8], x -> 1E-2log(10x) / x, line=:dash, c=:black, label=L"\log(t)/t")
annotate!((0.07, 0.05), text("(a)", 12))
p1

## Rescaling y-axis
L = sqrt(Ntarget / rho)
p2 = plot(xscale=:log10, yscale=:log10, legend=true,
    yticks=([1E-5, 1E-4, 1E-3, 1E-2], [L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}"]),
    xticks=([1, 10, 100, 1000, 1E4], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}"]))
for i in each(v0sigs_horizontal)
    couleur = cols_P[Ps_avg_horizontal[i, end]]
    v00 = sqrt(v0sigs_horizontal[i][1])
    plot!(times, remove_negative(ns_avg_horizontal[i, :] / L^2) * v00, c=couleur, rib=0, m=:circle, ms=3, line=true)
end
ylims!(1E-5, 1E-1)
annotate!((0.22, 0.93), text(L"n\,\sqrt{v_0}/L^2", 15, :center, :black))
annotate!((0.94, 0.1), text(L"t", 15, :top, :black))
plot!(times[5:end-8], x -> 1E-2log(10x) / x, line=:dash, c=:black, label=L"\log(t)/t")
annotate!((0.07, 0.05), text("(a)", 12))
p2

## Rescaling xy-axis
L = sqrt(Ntarget / rho)
p3 = plot(xscale=:log10, yscale=:log10, legend=true,
    yticks=([1E-5, 1E-4, 1E-3, 1E-2], [L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}"]),
    xticks=([1, 10, 100, 1000, 1E4], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}"]))
for i in each(v0sigs_horizontal)
    couleur = cols_P[Ps_avg_horizontal[i, end]]
    v00 = sqrt(v0sigs_horizontal[i][1])
    plot!(times * v00, remove_negative(ns_avg_horizontal[i, :] / L^2) * v00, c=couleur, rib=0, m=:circle, ms=3, line=true)
end
ylims!(1E-5, 1E-1)
annotate!((0.22, 0.93), text(L"n\,\sqrt{v_0}/L^2", 15, :center, :black))
annotate!((0.85, 0.1), text(L"t\,\sqrt{v_0}", 15, :top, :black))
plot!(times[5:end-8], x -> 1E-2log(10x) / x, line=:dash, c=:black, label=L"\log(t)/t")
annotate!((0.07, 0.05), text("(a)", 12))
p3

L = sqrt(Ntarget / rho)
p3bis = plot(xscale=:log10, yscale=:log10, legend=true,
    yticks=([1E-5, 1E-4, 1E-3, 1E-2], [L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}"]),
    xticks=([1, 10, 100, 1000, 1E4], [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}"]))
for i in each(v0sigs_horizontal)
    couleur = cols_P[Ps_avg_horizontal[i, end]]
    v00 = sqrt(v0sigs_horizontal[i][1])
    plot!(times, remove_negative(ns_avg_horizontal[i, :] / L^2) * v00^2, c=couleur, rib=0, m=:circle, ms=3, line=true)
end
ylims!(1E-5, 1E-1)
annotate!((0.22, 0.93), text(L"n\,v_0/L^2", 15, :center, :black))
annotate!((0.95, 0.1), text(L"t", 15, :top, :black))
plot!(times[5:end-8], x -> 1E-2log(10x) / x, line=:dash, c=:black, label=L"\log(t)/t")
annotate!((0.07, 0.05), text("(a)", 12))
p3bis



## ---------------- Plotting for xi correlation length ---------------- ##
## ---------------- Plotting for xi correlation length ---------------- ##
## ---------------- Plotting for xi correlation length ---------------- ##
## ---------------- Plotting for xi correlation length ---------------- ##
## No rescaling 
p4 = plot(axis=:log, legend=false)
for i in each(v0sigs_horizontal)
    couleur = cols_P[Ps_avg_horizontal[i, end]]
    v00 = sqrt(v0sigs_horizontal[i][1])
    plot!(times, remove_negative(xis_avg_horizontal[i, :]) / L, c=couleur, rib=0, m=:circle, ms=3, line=true)
end
p4
yticks!([0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], ["0.03", "", "0.05", "", "0.07", "", "", "0.1", "0.2", "0.3", "", "0.5", "", "0.7"])
plot!(times[5:end-11], x -> 4.5 / L * sqrt(x / log(8x)), line=:dash, c=:black, label=L"\sqrt{t/\log(t)}")
annotate!((0.1, 1.06), text(L"\xi/L", 15, :center, :black))
annotate!((0.94, 0.1), text(L"t", 15, :top, :black))
annotate!((0.065, 0.05), text("(b)", 12))
xlims!(0.2, 4E4)

## y axis rescaling 
p5 = plot(axis=:log, legend=false)
for i in each(v0sigs_horizontal)
    couleur = cols_P[Ps_avg_horizontal[i, end]]
    v00 = sqrt(v0sigs_horizontal[i][1])
    plot!(times, remove_negative(xis_avg_horizontal[i, :]) / L / v00, c=couleur, rib=0, m=:circle, ms=3, line=true)
end
p5
yticks!([0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], ["0.03", "", "0.05", "", "0.07", "", "", "0.1", "0.2", "0.3", "", "0.5", "", "0.7"])
plot!(times[5:end-11], x -> 4.5 / L * sqrt(x / log(8x)), line=:dash, c=:black, label=L"\sqrt{t/\log(t)}")
annotate!((0.2, 0.93), text(L"\xi/L/\sqrt{v_0}", 15, :center, :black))
annotate!((0.9, 0.1), text(L"t", 15, :top, :black))
annotate!((0.065, 0.05), text("(b)", 12))
xlims!(0.2, 4E4)

## xy axes rescaling 
p6 = plot(axis=:log, legend=false)
for i in each(v0sigs_horizontal)
    couleur = cols_P[Ps_avg_horizontal[i, end]]
    v00 = sqrt(v0sigs_horizontal[i][1])
    plot!(times * v00, remove_negative(xis_avg_horizontal[i, :]) / L / v00, c=couleur, rib=0, m=:circle, ms=3, line=true)
end
p6
yticks!([0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], ["0.03", "", "0.05", "", "0.07", "", "", "0.1", "0.2", "0.3", "", "0.5", "", "0.7"])
plot!(times[5:end-11], x -> 4.5 / L * sqrt(x / log(8x)), line=:dash, c=:black, label=L"\sqrt{t/\log(t)}")
annotate!((0.2, 0.93), text(L"\xi/L/\sqrt{v_0}", 15, :center, :black))
annotate!((0.9, 0.1), text(L"t\sqrt{v_0}", 15, :top, :black))
annotate!((0.065, 0.05), text("(b)", 12))
xlims!(0.2, 4E4)

## xi root n constant 
## xi root n constant 
## xi root n constant 
## xi root n constant 
pcst = plot(xaxis=:log, legend=false)#, size=(290, 180))#:topright)
for i in each(v0sigs_horizontal)
    couleur = cols_P[Ps_avg_horizontal[i, end]]
    v00 = (v0sigs_horizontal[i][1])
    data = remove_negative(xis_avg_horizontal[i, 2:end-4] .* sqrt.(ns_avg_horizontal[i, 2:end-4])) / L
    plot!(times[2:end-4] * v00, data * v00^(-0.25), c=couleur, m=:circle, ms=3.5, lw=1.3)
    # plot!(times[2:end-4], data, c=couleur, m=:circle, ms=2, lw=0.8)
end
ylims!(0.3, 0.69)
yticks!([0.4, 0.5, 0.6])
annotate!((0.2, 0.89), text(L"ξ\,\sqrt{n}/L", 15, :center))
annotate!((0.93, 0.08), text(L"t", 15, :center))
xticks!([1, 10, 100, 1000, 1E4], [L"10^{0}", "", L"10^{2}", "", L"10^{4}"])
pcst


##

plot(p1, p2, p3, p3bis, p4, p5, p6, pcst, layout=(2, 4), size=(1600, 800))



## Correlation Functions 
## Correlation Functions 
## Correlation Functions 
## Correlation Functions 
rr = 0:round(Int, L / 2)
p3 = plot(axis=:log, legend=false)
for i in 1:1:length(v0sigs_horizontal)
    v00 = sqrt(v0sigs_horizontal[i][1])
    couleur = cols_P[Ps_avg_horizontal[i, end]]
    plot!(rr[2:end]/v00, (v00)^(-0.5) * remove_negative(Cs_avg_horizontal[i, end])[2:end], c=couleur, rib=0, m=:circle, ms=3)
end
# plot!(rr[2:end], r -> r^(-T / 2π), line=:dot, c=:black, label=L"r^{-T/2\pi}")
plot!(rr[2:end], r -> 0.96r^(-0.25), line=:dash, c=:black, label=L"r^{-1/4}")
annotate!(1.25, 0.12, text("(c)", 12))
yticks!([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], ["0.1", "0.2", "0.3", "", "0.5", "", "0.7", "", "", "1"])
xticks!([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50], [L"10^{0}", "", "", "", "", "", "", "", "", L"10^{1}", "", "", "", ""])
annotate!((0.1, 0.945), text(L"C(r)", 15, :center, :black))
annotate!((0.94, 0.1), text(L"r", 15, :top, :black))
p3