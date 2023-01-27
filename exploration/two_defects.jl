cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian
include("../methods.jl")
const global R0 = 1
using Plots, ColorSchemes, LaTeXStrings
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

#= This file investigates the behaviour of two defects. 
In this project, impossible to isolate one defect, 
because Free BC are not compatible with the motion of individual particles. 
Plan : 
A. Study R(t)
B. If unbounded, study MSD with two defects far apart.=#

## A. Study R(t) 
Ntarget = Int(1E4)
rho = 4.51 / π
aspect_ratio = 1
N, L, L = effective_number_particle(Ntarget, rho, aspect_ratio)
T = 0.1 # temperature for angle diffusion
v0 = 0.25
sigma = 0.05
every = 2;
tmax = 30;
r0s = 10:10:50
dt = determine_dt(T, sigma, v0, N, rho)
every = 2;
tmax = 100;
times = every:every:tmax;

dfts = Vector{DefectTracker}(undef, length(r0s))
z = @elapsed for i in each(r0s)
    r0 = r0s[i]
    params_init = ["random", "pair", r0]
    pos, thetas, omegas, psis = initialisation(N, L, L, sigma, params_init)
    t = 0
    dft = DefectTracker(pos, thetas, N, L, L, t)
    dft, pos, thetas, t = track!(dft, pos, thetas, omegas, psis, T, v0, N, L, L, dt, t, times)
    dfts[i] = dft
end
prinz(z)

## Compute R(t) from dft
Rts = Vector{Vector{Float64}}(undef, length(r0s))
for i in each(r0s)
    dft = dfts[i]
    Rts[i] = interdefect_distance(dft.defectsP[1], dft.defectsN[1], L, L)
end
plot()
for i in each(r0s)
    rt = Rts[i]
    plot!(rt, label="r0 = $(r0s[i])")
end

## B. Study MSD with two defects far apart



