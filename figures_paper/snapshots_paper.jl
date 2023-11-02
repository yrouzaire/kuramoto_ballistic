cd("/Users/yrouzaire/Documents/Recherche/GitHub/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl");
using Plots, ColorSchemes, LaTeXStrings
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

## ---------------- Data Generation ---------------- ##
## ---------------- Data Generation ---------------- ##


include("../parameters.jl")
system = System(param)
tmax = 1000
@time evolve!(system, tmax)

gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot_thetas(system, particles=true, vertical=true, defects=false)
filename = "data/snapshots/snapshot_v$(v0)_σ$(sigma)_ρ$(rho)_time$(tmax).jld2"
@save filename system tmax param


## ---------------- Plots green region ---------------- ##
## ---------------- Plots green region ---------------- ##

filename = "data/snapshots/snapshot_v2_σ0.3_ρ1_time50.jld2"
@load filename system tmax param
plot_thetas(system, particles=true, vertical=true, defects=true)
savefig("figures_paper/snapshot_v$(v0)_σ$(sigma)_ρ$(rho)_particles_time$(tmax).svg")

pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);

plot_thetas(system, particles=false, vertical=false, defects=true)
savefig("figures_paper/snapshot_v$(v0)_σ$(sigma)_ρ$(rho)_time$(tmax).svg")

filename = "data/snapshots/snapshot_v2_σ0.3_ρ1_time500.jld2"
@load filename system tmax param
plot_thetas(system, particles=false, vertical=false, defects=true)
savefig("figures_paper/snapshot_v$(v0)_σ$(sigma)_ρ$(rho)_time$(tmax).svg")


filename = "data/snapshots/snapshot_v2_σ0.3_ρ1_time1000.jld2"
@load filename system tmax param
plot_thetas(system, particles=false, vertical=false, defects=true)
savefig("figures_paper/snapshot_v$(v0)_σ$(sigma)_ρ$(rho)_time$(tmax).svg")

gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);


## ---------------- Plots red region ---------------- ##
## ---------------- Plots red region ---------------- ##

filename = "data/snapshots/snapshot_v0.2_σ0.3_ρ1_time50.jld2"
@load filename system tmax param
plot_thetas(system, particles=true, vertical=true, defects=true)
savefig("figures_paper/snapshot_v$(v0)_σ$(sigma)_ρ$(rho)_particles_time$(tmax).svg")

pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);

plot_thetas(system, particles=false, vertical=false, defects=true)
savefig("figures_paper/snapshot_v$(v0)_σ$(sigma)_ρ$(rho)_time$(tmax).svg")

filename = "data/snapshots/snapshot_v0.2_σ0.3_ρ1_time500.jld2"
@load filename system tmax param
plot_thetas(system, particles=false, vertical=false, defects=true)
savefig("figures_paper/snapshot_v$(v0)_σ$(sigma)_ρ$(rho)_time$(tmax).svg")


filename = "data/snapshots/snapshot_v0.2_σ0.3_ρ1_time1000.jld2"
@load filename system tmax param
plot_thetas(system, particles=false, vertical=false, defects=true)
savefig("figures_paper/snapshot_v0.2_σ0.3_ρ1_time1000.svg")

gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);

## ---------------- Plots yellow region ---------------- ##
## ---------------- Plots yellow region ---------------- ##
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
filename = "data/snapshots/snapshot_v0.55_σ0.3_ρ1_time50.jld2"
@load filename system tmax param
plot_thetas(system, particles=true, vertical=true, defects=true)
savefig("figures_paper/snapshot_v$(v0)_σ$(sigma)_ρ$(rho)_particles_time$(tmax).svg")

pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);

plot_thetas(system, particles=false, vertical=false, defects=true)
savefig("figures_paper/snapshot_v$(v0)_σ$(sigma)_ρ$(rho)_time$(tmax).svg")

filename = "data/snapshots/snapshot_v0.55_σ0.3_ρ1_time500.jld2"
@load filename system tmax param
plot_thetas(system, particles=false, vertical=false, defects=true)
savefig("figures_paper/snapshot_v$(v0)_σ$(sigma)_ρ$(rho)_time$(tmax).svg")


filename = "data/snapshots/snapshot_v0.55_σ0.3_ρ1_time1000.jld2"
@load filename system tmax param
plot_thetas(system, particles=false, vertical=false, defects=true)
savefig("figures_paper/snapshot_v$(v0)_σ$(sigma)_ρ$(rho)_time$(tmax).svg")

gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);


