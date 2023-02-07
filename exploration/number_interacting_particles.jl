cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian, Sobol
include("../methods.jl")
const global R0 = 1
using Plots, ColorSchemes, LaTeXStrings
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot()
&

#= The file aims at computing the number of different particles
a test particles interacts with over time. Underlying idea: understanding
the relation between moving particles and MF approximation.
=#

##
include("../parameters.jl");
# scanned parameters 
rhos = collect(1:0.25:2)
v0s = [0.1, 0.5, 1, 2]
v0s = [1]
R0s = [1]
Niter = Int(1E4)
R = 5

number_interactions = zeros(Ntarget, Niter + 1, length(rhos), length(v0s), length(R0s), R)

z = @elapsed for j in each(rhos), k in each(v0s), m in each(R0s)
    Threads.@threads for r in 1:R
        rho = rhos[j]
        v0 = v0s[k]
        R0 = R0s[m]
        N, L, L = effective_number_particle(Ntarget, rho)

        param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio, :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0, :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)
        
        system = System(param)
        ID_interactions = [Int[] for i in 1:N]

        dt = 0.1
        for t in 0:Niter
            ind_neighbours = get_list_neighbours(get_pos(system), N, L, L, R0)
            for n in 1:min(Ntarget, N)
                push!(ID_interactions[n], ind_neighbours[n]...)
                unique!(ID_interactions[n])
                number_interactions[n, t+1, j, k, m, r] = length(ID_interactions[n])
            end
            update_positions!(system)
        end
    end
end
nintavg = mean(number_interactions, dims=(1, 6))[1, :, :, :, :, 1]
prinz(z)

##
cst_atan = 0.2
cst_exp = 1.9
p = plot(legend=false, uaxis=:log)
for j in each(rhos)
    for k in 1#each(v0s)
        for m in 1#each(R0s)
            rho = rhos[j]
            v0 = v0s[k]
            R0 = R0s[m]
            times = collect(0:Niter)
            plot!(times, nintavg[1:Niter+1, j, k, m] / N, label="ρ = $rho , v0 = $v0, R0 = $R0", c=j)

            N0 = π * rho * R0^2
            # plot!(times, t -> N0/N + (N - N0)/N * (1 - exp(-cst_exp * t/N)), c=:black)
            # plot!(t -> N0 / N + 2 / pi * (N - N0) / N * atan(cst_atan * rho * v0 / N * t), c=j, ls=:dash)
        end
    end
end
# plot!(x->0.00071(x)^1)
p

##
nintavg[1, :, :] - pi * (rhos)


## Effective density
N = 10000
rhos = collect(1:0.2:2)
Ls = round.(Int, sqrt.(N ./ rhos))
N ./ (Ls .^ 2) - rhos

## Visualizing the area of seen particules over time
for v0 in [0.1, 0.5, 1, 2, 5]
    N = Int(1E4)
    rho = 5
    L = round(Int, sqrt(N / rho))
    T = 0.1 # to be compared to Tc ~ 1 when \sigma = 0
    println("v0 = $v0")
    sigma = 0

    global pos, thetas, omegas, psis = initialisation(N, L, L, sigma)

    dt = 0.1
    Niter = 500
    seen = zeros(Int, N)
    z = @elapsed anim = @animate for t in 1:Niter

        for n in eachindex(psis)
            pos[:, n] += v0 * dt * [cos(psis[n]), sin(psis[n])]
        end
        pos = mod.(pos, L)

        nb_cells_1D = Int(div(L, R0)) + 1
        head = -ones(Int, nb_cells_1D, nb_cells_1D) # head[i,j] contains the index of the first particle in cell (i,j). -1 if empty
        list = -ones(Int, N) # list[n] contains the index of the particle to which particle n points. -1 if it points to no one
        for n in 1:N
            cellx, celly = Int(div(pos[1, n], R0)) + 1, Int(div(pos[2, n], R0)) + 1 # cell to which particle n belongs
            list[n] = head[cellx, celly]
            head[cellx, celly] = n
        end

        for n in 1:N
            poscur = pos[:, n]

            cellx, celly = Int(div(pos[1, n], R0)) + 1, Int(div(pos[2, n], R0)) + 1 # cell to which particle n belongs
            should_take_mod = (cellx == 1) || (cellx == nb_cells_1D) || (celly == 1) || (celly == nb_cells_1D)
            if should_take_mod
                neighbouring_cells = Vector{Int}[[cellx, celly], [cellx, mod1(celly + 1, nb_cells_1D)], [mod1(cellx + 1, nb_cells_1D), celly], [cellx, mod1(celly - 1, nb_cells_1D)], [mod1(cellx - 1, nb_cells_1D), celly], [mod1(cellx + 1, nb_cells_1D), mod1(celly + 1, nb_cells_1D)], [mod1(cellx - 1, nb_cells_1D), mod1(celly - 1, nb_cells_1D)], [mod1(cellx - 1, nb_cells_1D), mod1(celly + 1, nb_cells_1D)], [mod1(cellx + 1, nb_cells_1D), mod1(celly - 1, nb_cells_1D)]]
            else
                neighbouring_cells = Vector{Int}[[cellx, celly], [cellx, celly + 1], [cellx + 1, celly], [cellx, celly - 1], [cellx - 1, celly], [cellx + 1, celly + 1], [cellx - 1, celly - 1], [cellx - 1, celly + 1], [cellx + 1, celly - 1]]
            end

            ind_neighbours = Int[]
            for (i, j) in neighbouring_cells
                next = head[i, j]
                if next ≠ -1
                    if 0 < dist(poscur, pos[:, next], L, L) < R0
                        push!(ind_neighbours, next)
                    end
                    while list[next] ≠ -1
                        if 0 < dist(poscur, pos[:, list[next]], L, L) < R0
                            push!(ind_neighbours, list[next])
                        end
                        next = list[next]
                    end
                end
            end
            if n == 1
                seen[ind_neighbours] .= 1
            end
        end
        p = scatter(pos[1, :], pos[2, :], marker_z=seen[:, end], c=cgrad([:green, :red]), m=2, clims=(0, 1), ms=275 / L, size=(512, 512), xlims=(0, L), ylims=(0, L))
        scatter!((pos[1, 1], pos[2, 1]), c=:black, m=5)# display(p)
    end
    prinz(z)
    mp4(anim, "films/dispersion_seen_v0$(v0).mp4")
end

## Visualizing analytic predictions for the phase space separation
rhoc = 1.44;
cst = 3.9;
vc(rho) = (rhoc - rho) / rho / cst * π^2 * R0 / 2
vv = logspace(1e-2, 5, 100, digits=3)
p = plot(xaxis=:log, legend=false)
for rho in 1:0.1:2
    plot!(vv, (max.(0, (vv) .- vc(rho))) .^ 1, label="rho = $rho")
end
ylims!(0, 0.5)
p

## Some Stuff 
seuil = 1
xx = 0.01:0.01:5
function f1(v, rho)
    if v > vc(rho)
        return (v - vc(rho)) * sqrt(rho * π)
    else
        return 0
    end
end

function f2(x)
    if x > seuil
        return x / (0.01 + sqrt(x))
    else
        return 0
    end
end
p = plot(xaxis=:log)
for rho in 1:0.2:2
    plot!(xx, f1.(xx, rho), label="rho = $rho")
end
p


## Investigate the behaviour of critical sigma σc(v0) both for ρ ≤ ρc and ρ > ρc
rhoc = 1.44
vc(rho) = (rhoc - rho) / rho / cst * π^2 * R0 / 2
v0s = logspace(1e-2, 1, 10, digits=3)
sigmas = collect(0:0.05:0.2)
rhos = [1, 1.44, 2]
R = 1
N = Int(1E3)
T = 0.1
seuil = 0.5
tmax = 100
times = 0:tmax/10:tmax

critical_sigmas = NaN * zeros(length(v0s), length(rhos))
z = @elapsed for i in each(v0s), k in each(rhos)
    for r in 1:R
        for j in each(sigmas)
            v0 = v0s[i]
            sigma = sigmas[j]
            rho = rhos[k]
            println("v0 = $v0, σ = $sigma, ρ = $rho")
            L = round(Int, sqrt(N / rho))

            dt = determine_dt(T, sigma, v0, N, rho)
            pos, thetas, omegas, psis = initialisation(N, L, L, sigma, ["lowtemp"])
            t = 0.0
            token = 1

            already_broken_at_time = -1
            while t < tmax
                t += dt
                ind_neighbours = get_list_neighbours(pos, N, L, L)
                pos, thetas = update(pos, thetas, omegas, psis, ind_neighbours, T, v0, N, L, L, dt)
                if t ≥ times[token]
                    token += 1
                    P = polarOP(thetas)[1]
                    if P < seuil
                        already_broken_at_time = t
                        println("Broken at t = $already_broken_at_time")
                        break # gets out of the while loop only 
                    end
                end
            end

            P = polarOP(thetas)[1]
            if (already_broken_at_time > 0) || P < seuil
                critical_sigmas[i, k, r] = sigma
                println("σc = $sigma for v0 = $v0 and rho = $rho, at time = $already_broken_at_time")
                break # gets out of the sigma for loop
            end
        end
    end
end
prinz(z)




## Solving the dynamics of the fully connected Kuramoto model
# N = Int(1E4)
# function complexOP(thetas::Vector{Float64})
#     z = mean(exp.(im*thetas))
#     return abs(z),angle(z)
# end

# dt = 1E-1 ; tmax = 500
# times = logspace(dt,tmax,40,digits=2)
# rs = zeros(length(times))
# thetas = 2pi*rand(N)
# t = 0 ; token = 1
# z = @elapsed while t < tmax
#     t += dt
#     r,psi = complexOP(thetas)
#     thetas += dt*0.05r*sin.(psi .- thetas)
#     if t ≥ times[token]
#         rs[token] = r
#         token = min(token+1,length(times))
#     end
# end
# prinz(z)
# plot(times,rs,uaxis=:log,m=true)
#     plot!(x->1-exp(-x/10000),uaxis=:log)
#     plot!(x->2/pi*atan(x/1000))
#     histogram(mod.(thetas,2pi),bins=50,normalize=true)

