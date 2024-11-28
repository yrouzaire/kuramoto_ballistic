using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian, Sobol, Dates
include("../methods.jl")
include("../parameters.jl")
using CairoMakie,GLMakie, ColorSchemes, LaTeXStrings
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);

## --------------- Data Generation --------------- ##
## --------------- Data Generation --------------- ##
## --------------- Data Generation --------------- ##
## --------------- Data Generation --------------- ##
R = 4
for r in 1:R 
    param[:T] = 0.1
    param[:v0] = 1.0
    param[:sigma] = 0.1
    param[:distribution_type] = "gaussian"
    system = System(param)
    times = collect(0:5:2000)
    xss = Array{Float32}(undef, system.N, length(times))
    yss = Array{Float32}(undef, system.N, length(times))
    thetass = Array{Float32}(undef, system.N, length(times))
    thetas_cgss = Array{Float32}(undef, system.Lx, system.Ly, length(times))

    z = @elapsed for tt in each(times)
        evolve!(system, times[tt])
        xss[:,tt] = [get_pos(system)[i][1] for i in 1:system.N]
        yss[:,tt] = [get_pos(system)[i][2] for i in 1:system.N]
        thetass[:,tt] = mod.(get_thetas(system), 2pi)
        thetas_cgss[:,:,tt] = mod.(cg(system), 2pi)
        println("t = $(system.t), ($tt/$(length(times)))")
    end
    prinz(z)

    ## --------------- Movie Generation --------------- ##
    ## --------------- Movie Generation --------------- ##
    ## --------------- Movie Generation --------------- ##
    ## --------------- Movie Generation --------------- ##

    function initial_figure(xs, ys, thetas, thetas_cg, Lx, Ly)
        fig = Figure()

        ax = Axis(fig[1, 1:2], title=L"t = %$0")
        hidedecorations!(ax)
        ax1 = Axis(fig[1, 1], aspect=1,
            xticks=([1, Lx]), yticks=([1, Ly]),
            limits=(1, Lx, 1, Ly))
        CairoMakie.scatter!(ax1, xs, ys, markersize=6, color=thetas, colormap=cols, colorrange=(0, 2pi))

        ax2 = Axis(fig[1, 2], aspect=1,
            xticks=([1, Lx]), yticks=([1, Ly]))
        p_cg = CairoMakie.heatmap!(ax2, thetas_cg, colormap=cols, colorrange=(0, 2pi))

        Colorbar(fig[1, 3], p_cg,
            label=L"\theta",
            ticklabelsize=20, labelrotation=0, labelsize=40,
            tickalign=1,
            labelpadding=-10, # distance between the colorbar ticks and the label
            ticks=([0, pi / 2, pi, 3pi / 2, 2pi], [L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"]),
        )
        colgap!(fig.layout, 1, 10)
        colgap!(fig.layout, 2, 10)
        rowsize!(fig.layout, 1, Aspect(1, 1))
        resize_to_layout!(fig)
        fig
    end

    function name_directory_distribution(distribution_type)
        distribution_type = lowercase(distribution_type)
        if distribution_type in ["gaussian", "normal"]
            return "gaussian"
        elseif distribution_type in ["uniform"]
            return "uniform"
        elseif distribution_type in ["laplace", "exponential"]
            return "laplace"
        elseif distribution_type in ["lorentzian", "cauchy"]
            return "cauchy"
        elseif distribution_type in ["trunc_lorentzian", "trunc_cauchy", "truncated_lorentzian", "truncated_cauchy"]
            return "cauchy_trunc"
        else
            error("Unknown distribution type")
        end
    end


    using GLMakie
    GLMakie.activate!()

    xs = Observable(xss[:,1])
    ys = Observable(yss[:,1])
    thetas = Observable(thetass[:,1])
    thetas_cg = Observable(thetas_cgss[:,:,1])


    fig = initial_figure(xs, ys, thetas, thetas_cg, system.Lx, system.Ly)


    path_to_movie = "/Users/yrouzaire/Documents/Recherche/GitHub/kuramoto_ballistic/movies/$(name_directory_distribution(param[:distribution_type]))/"
    filename = "Kuramoto_T$(param[:T])_v0$(param[:v0])_σ$(param[:sigma])_tmax$(tmax)_r$(r).mp4"

    GLMakie.record(fig, path_to_movie*filename, 1:length(times)) do tt
        println("Frame $(100tt / length(times)) %")
        
        xs[] = xss[:, tt]
        ys[] = yss[:, tt]
        thetas[] = thetass[:, tt]
        thetas_cg[] = thetas_cgss[:, :, tt]

        ax.title = L"t = %$(times[tt])"
    end


end