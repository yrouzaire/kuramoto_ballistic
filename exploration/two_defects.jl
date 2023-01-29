cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian, LambertW
include("../methods.jl")
const global R0 = 1
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
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
sigma = 0.0
every = 2;
r0s = 5:5:25
dt = determine_dt(T, sigma, v0, N, rho)
every = 2;
tmax = 400;
times = every:every:tmax;
R = 50

dfts = Array{Union{DefectTracker,Nothing}}(undef, length(r0s), R)
z = @elapsed for i in each(r0s)
    Threads.@threads for r in 1:R
        r0 = r0s[i]
        params_init = ["random", "pair", r0]
        pos, thetas, omegas, psis = initialisation(N, L, L, sigma, params_init)
        t = 0
        dft = DefectTracker(pos, thetas, N, L, L, t)
        # try 
        dft, pos, thetas, t = track!(dft, pos, thetas, omegas, psis, T, v0, N, L, L, dt, t, times)
        dfts[i, r] = dft
        # catch e
        #     println("Error at r0 = $(r0) and r = $(r)")
        #     println(e)
        #     dfts[i, r] = nothing
        # end    
    end
end
prinz(z)

filename = "data/two_defects_r0s_rho$(rho)_v0$(v0)_sigma$(sigma)_tmax$(tmax)_R$(R)_N$(Ntarget).jld2"
@save filename dfts tmax times v0 sigma dt r0s Ntarget N L R T rho aspect_ratio runtime = z
# @load filename dfts tmax times v0 sigma dt r0s Ntarget N L R T rho aspect_ratio runtime 

## Compute R(t) from dft;  Linear time
# Rts = Array{Vector{Float64}}(undef, length(r0s), R)
# for i in each(r0s), r in 1:R
#     dft = dfts[i, r]
#     Rts[i, r] = interdefect_distance(dft.defectsP[1], dft.defectsN[1], L, L)
# end

# Rts_avg = Vector{Vector{Float64}}(undef, length(r0s))
# for i in each(r0s)
#     tmp = vector_of_vector2matrix(Rts[i, :])
#     Rts_avg[i] = nanmean(tmp, 2)[:, 1]
# end

# p = plot(legend=false, xlabel="t", ylabel="R(t)")
# for i in each(r0s)
#     rt = Rts_avg[i]
#     plot!(rt, label="r0 = $(r0s[i])")
# end
# p

## Compute R(t*) from dft; Reversed time = time to annihilation
Rts_reverse = Array{Vector{Float64}}(undef, length(r0s), R)
for i in each(r0s), r in 1:R
    dft = dfts[i, r]
    Rts_reverse[i, r] = reverse(interdefect_distance(dft.defectsP[1], dft.defectsN[1], L, L))
end

Rts_reverse_avg = Vector{Vector{Float64}}(undef, length(r0s))
for i in each(r0s)
    tmp = vector_of_vector2matrix(Rts_reverse[i, :])
    Rts_reverse_avg[i] = nanmean(tmp, 2)[:, 1]
end

p = plot(legend=:bottomright, xlabel="t*", ylabel="R(t*)", axis=:log, title="ρ = $(round(rho,digits=2))")
for i in each(r0s)
    rt = Rts_reverse_avg[i]
    plot!(rt[2:end], label="r0 = $(r0s[i])", m=true)
end
# savefig("figures/two_defects/Rt_several_r0s_rho$(round(rho,digits=2))_v0$(v0)_sigma$(sigma)_tmax$(tmax)_R$(R)_N$(Ntarget)..png")
p
plot!(x -> 3x^0.33, c=:black, line=:dot)
plot!(x -> exp(lambertw(2π / 2 * x) / 2), c=:black, line=:dash)
## B. Study MSD with two defects far apart
