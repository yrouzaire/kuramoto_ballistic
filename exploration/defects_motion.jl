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

## Cluster results analysis
filename = "data/defects_motion_IN_VITRO_v0s_sigmas_N1E4.jld2"
@load filename runtimes T N rho v_sigmas times tmax comments R dfts_fusion dfts_fusion_undef
histogram(runtimes / 3600 / 24, bins=20)
L = round(Int, sqrt(N / rho))


times
dfts_fusion[4]
dfts_fusion_undef

# for each v_sigma, gather al DefectTrackers from all realisations
dfts = [DefectTracker[] for i in each(v_sigmas)]
for i in each(v_sigmas)
    for r in each(dfts_fusion)
        push!(dfts[i], dfts_fusion[r][i])
    end
end
plot(MSD(dfts[1],L,L)[3][2:end],axis=:log)
plot(MSD(dfts[1][3],L,L)[1])
plot(MSD(dfts[1][4],L,L)[1])
plot(MSD(dfts[1][5],L,L)[1])
plot(MSD(dfts[1][6],L,L)[1])
plot(MSD(dfts[1][7],L,L)[1])





# MSD
mySD = zeros(length(v_sigmas), length(times) + 1, R)
for i in each(v_sigmas)
    for r in 1:R
        positions = dfts_fusion[r][i].defectsP[1].pos
        pos0 = positions[1]
        mySD[i, :, r] = [sum(abs2, pos .- pos0) for pos in positions]
    end
end
myMSD = mean(mySD, dims=3)[:, :, 1]
p = plot(axis=:log)
for i in each(v_sigmas)
    plot!(times, myMSD[i, 2:end])
end
# plot!(x->x)
p

# R(t)
Rt = zeros(length(v_sigmas), length(times) + 1, R)
for i in each(v_sigmas)
    for r in 1:R
        positions = dfts_fusion[r][i].defectsP[1].pos
        def1 = dfts_fusion[r][i].defectsP[1]
        def2 = dfts_fusion[r][i].defectsN[1]
        Rt[i, :, r] = interdefect_distance(def1, def2, L, L)
    end
end
meanRt = mean(Rt, dims=3)[:, :, 1]
p = plot(axis=:log)
for i in each(v_sigmas)
    plot!(times, meanRt[i, 2:end])
end
# plot!(x->x)
p

#= This file investigates the motion of defects.
1. Departure from XY model as v varies at fixed T
2. Departure from XY model as T varies at fixed v
3. Departure from AXY model as v varies at fixed T, sigma
4. Departure from AXY model as sigma varies at fixed v,T

Basically consists in scan on two parameters, v and sigma
=#

## For cluster
N = Int(1E4)
rho = 2
aspect_ratio = 1 # Lx/Ly
Lx = round(Int, sqrt(N / rho * aspect_ratio))
Ly = round(Int, sqrt(N / rho / aspect_ratio))
T = 0.1 # to be compared to Tc ~ 1 when \sigma = 0
# In vitro (green region only)
v_sigmas = [(0, 0), (0.05, 0), (0.1, 0), (0.5, 0), (1, 0),
    (0, 0.1), (0.05, 0.1), (0.1, 0.1), (0.5, 0.1), (1, 0.1),
    (0.4, 0.2), (0.65, 0.3)]
# In vivo (whole phase space)
# v_sigmas = [ (0,0) , (0.05,0) , (0.1,0) , (0.5,0) , (1,0) ,
#              (0,0.1) , (0.05,0.1) , (0.1,0.1) , (0.5,0.1) , (1,0.1) ,
#              (0.4,0.2) , (0.65,0.3)]
v_sigmas = [(0.5, 0.1)]
transients = 2;
tmax = 100;
times = logspace(round(Int, transients), tmax, 50, digits=1)
times = 5:5:tmax

dfts = Array{DefectTracker,1}(undef, length(v_sigmas))
z = @elapsed for i in each(v_sigmas)
    v0, σ = v_sigmas[i]
    dt = determine_dt(T, σ, v0, N, rho)
    t = 0.0
    pos, thetas, omegas, psis = initialisation(N, Lx, Ly, σ, ["pair", round(Int, Lx / 2)])
    ind_neighbours = get_list_neighbours(pos, N, Lx, Ly)
    pos, thetas = update(pos, thetas, omegas, psis, ind_neighbours, T, v0, N, Lx, Ly, dt)

    pos, thetas = evolve(pos, thetas, omegas, psis, T, v0, N, Lx, Ly, dt, transients)

    dft = DefectTracker(pos, thetas, N, Lx, Ly, t)
    dft, pos, thetas, t = track!(dft, pos, thetas, psis, omegas, T, v0, N, Lx, Ly, dt, t, times)
    dfts[i] = dft
end
prinz(z)

plot(MSD(dfts, Lx, Ly)[1])
plot(interdefect_distance(dfts[1].defectsP[1], dfts[1].defectsN[1], Lx, Ly))

comments = "Study the motion of defects at various locations in phase space. Goal: understand the role of v and sigma on the MSD. Preliminary results indicate that MSD ~ t even for σ>0. Understand the transition from immobile (MSD ~ t^1.5) to mobile (MSD ~ t)."
filename = "data/defects_motion_IN_VITRO_v0s_sigmas_N1E4_r$real.jld2"
JLD2.@save filename dfts Ts Ns rhos inits v_sigmas times tmax runtime = z comments

## For this PC
N = Int(1E4)
rho = 1
v0s = [1]
σs = [0.2]
L = round(Int, sqrt(N / rho))
T = 0.1 # to be compared to Tc ~ 1 when \sigma = 0
R = 2
transients = 1;
tmax = 2000;
times = logspace(round(Int, transients), tmax, 50, digits=1)

dfts = Array{DefectTracker,3}(undef, length(v0s), length(σs), R)
z = @elapsed for i in each(v0s)
    for j in each(σs)
        σ = σs[j]
        v0 = v0s[i]
        for r in 1:R
            dt = determine_dt(T, σ, v0, N, rho)
            t = 0.0
            pos, thetas, psis, omegas = initialisation(N, L, σ, ["pair", round(Int, L / 2)])
            while t < transients
                t += dt
                pos, thetas = update(pos, thetas, psis, omegas, T, v0, N, L, dt)
            end
            dft = DefectTracker(pos, thetas, N, L, t)
            dft, pos, thetas, t = update_and_track!(dft, pos, thetas, psis, omegas, T, v0, N, L, dt, t, times)
            dfts[i, j, r] = dft
        end
    end
end

msda = Array{Number,3}(undef, length(v0s), length(σs), length(times) + 1)
msdp = Array{Number,3}(undef, length(v0s), length(σs), length(times) + 1)
msdn = Array{Number,3}(undef, length(v0s), length(σs), length(times) + 1)
for i in each(v0s)
    for j in each(σs)
        msda[i, j, :], msdp[i, j, :], msdn[i, j, :] = MSD(dfts[i, j, :], L)
    end
end

# findfirst(!iszero,msdm[1,1,:])
# ind_start = max(findfirst(!iszero,msdp[1,1,:]),findfirst(!iszero,msdm[1,1,:]))

p = plot(legend=false, axis=:log)
for i in each(v0s)
    for j in each(σs)
        ind_start = max(findfirst(!iszero, msdp[i, j, :]), findfirst(!iszero, msdn[i, j, :]))
        plot!(times[ind_start:end], msda[i, j, ind_start:end-1] .+ 1)
        plot!(times[ind_start:end], msdp[i, j, ind_start:end-1] .+ 1)
        plot!(times[ind_start:end], msdn[i, j, ind_start:end-1] .+ 1)
    end
end
p
plot!(x -> x / 5)

plot(mean_distance_to_annihilator(dfts[1, 1, :], L)[1], rib=0)

# Old version, without different i,j,r
t = 0.0;
pos, thetas, psis, omegas = initialisation(N, L, σ, ["pair", round(Int, L / 3)]);
# t = 0. ; pos,thetas,psis,omegas = initialisation(N,L,σ,["single",+1])
plot(pos, thetas, N, L)

z = @elapsed while t < 1
    t += dt
    global pos, thetas = update(pos, thetas, psis, omegas, T, v0, N, L, dt)
end
plot(pos, thetas, N, L)
number_defects(pos, thetas, N, L)
# number_defectsP(pos,thetas,N,L)
# number_defectsN(pos,thetas,N,L)
# plot(pos,thetas,defects=true)
# plot(pos,thetas,defects=false,particles=true)

dft = DefectTracker(pos, thetas, N, L, t)
tmax = 1000
# every = 1 ; times = round(Int,t):every:tmax
times = logspace(round(Int, t), tmax, 40, digits=1)
dft, pos, thetas, t = update_and_track!(dft, pos, thetas, psis, omegas, T, v0, N, L, dt, t, times)
# d=dft.defectsP[1]


msdall, msdp, msdm = MSD(dft, L)
plot()
ind_start = max(findfirst(!iszero, msdp), findfirst(!iszero, msdm))
plot!(msdall[ind_start:end] .+ 1, yaxis=:log, m=:circle)
plot!(msdp[ind_start:end] .+ 1)
plot!(msdm[ind_start:end] .+ 1)
plot!(x -> 1E-1x^2)


plot(mean_distance_to_annihilator(dft, L)[1], rib=mean_distance_to_annihilator(dft, L)[2])

## Archives
# N = Int(1E4)
#     # rho = 1.3*(4.51/pi)
#     rho = 2
#     v0 = 0
#     σ  = 0.3
#     aspect_ratio = 1 # Lx/Ly
#     Lx  = round(Int,sqrt(N/rho*aspect_ratio))
#     Ly  = round(Int,sqrt(N/rho/aspect_ratio))
#     T  = 0.1 # to be compared to Tc ~ 1 when \sigma = 0
#     tmax = 10 ; dt = determine_dt(T,σ,v0,N,rho)
#     global R0 = 1
#
# t = 0.0
#     params_init = ["hightemp",Lx/2]
#     pos,thetas,omegas,psis = initialisation(N,Lx,Ly,σ,params_init)
#     plot(pos,thetas,N,Lx,Ly,particles=false,defects=false,title="t = $(round(Int,t))")
#
# ind_neighbours0 = get_list_neighbours(pos,N,Lx,Ly)
# z = @elapsed while t < tmax
#     t += dt
#     pos,thetas = update(pos,thetas,omegas,psis,ind_neighbours0,T,v0,N,Lx,Ly,dt)
# end
# prinz(z)
#
# z = @elapsed while t < tmax
#     t += dt
#     ind_neighbours = get_list_neighbours(pos,N,Lx,Ly)
#     pos,thetas = update(pos,thetas,omegas,psis,ind_neighbours,T,v0,N,Lx,Ly,dt)
# end
# prinz(z)
# plot(pos,thetas,N,Lx,Ly,particles=false,defects=false,title="t = $(round(Int,t))")
