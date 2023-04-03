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

## ---------------- Monitor the defects in the green area  ---------------- ##
comments = "From the defects data, one will be able to infer : \n
A. the separating distance between the two defects R(t) \n
B. the MSD and diffusion coeff of an individual defect. "
# Physical Params 
Ntarget = Int(1E3)
aspect_ratio = 1
T = 0.1
R0 = 1
rho = 1
rhoc = 4.51 / π

# Initialisation parameters
init_pos = "random"
init_theta = "pair"
r0 = 20.0
q = 1.0
params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => NaN, :q => q)

# Simulation parameters
v0sigs = [(0.5,0),(0.5,0.1)]
r0s = 5:10:35
tmax = 1E2
times = 0:5:tmax # linear time

z = @elapsed for i in each(v0sigs), j in each(r0s)
    v0, sigma = v0sigs[i]
    r0 = r0s[j]

    println("v0 = $v0, σ = $sigma, r0 = $r0 ,  $(100i/length(v0sigs))%")
    N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
    dt = determine_dt(T, sigma, v0, N, rho)

    params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
    
    param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
        :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
        :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)

    t = 0.0
    system = System(param)
    dft = DefectTracker(system, t)
    dft, pos, thetas, t = track!(dft,system,times)
end
prinz(z)


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
Rts_reverse = Array{Vector{Float64}}(undef, length(r0s), length(v0sigs), R)
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
