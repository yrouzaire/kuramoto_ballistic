#"Copyright (c) 2022 Y.Rouzaire All Rights Reserved."

#################### Structure Definitions ####################
using BenchmarkTools, Sobol, Parameters
using Plots, ColorSchemes
mutable struct Agent{F<:AbstractFloat}
    pos::Tuple{F,F}
    pos0::Tuple{F,F} # position at t=0, used to implement phonons 
    theta::F
    omega::F
    psi::F
end

mutable struct System{F<:AbstractFloat}
    agents::Vector{Agent{F}}
    Lx::Int
    Ly::Int
    rho::Number
    Ntarget::Int
    N::Int
    T::Number
    sigma::Number
    v0::Number
    R0::Number
    dt::Number
    t::Number
    params::Dict
    periodic::Bool
    phonons::Bool
    phonon_amplitude::Number
    phonon_k::Number
    phonon_omega::Number
end

function System(params; float_type=Float32)
    @unpack Ntarget, aspect_ratio, rho, T, sigma, v0, R0, params_init, params_phonons = params
    @unpack init_pos, init_theta, r0, q = params_init
    @unpack phonons, phonon_amplitude, phonon_k, phonon_omega = params_phonons
    Ntarget2, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
    
    pos,N  = initialisation_positions(Ntarget2, Lx, Ly, init_pos)
    thetas = initialisation_thetas(N, init_theta, r0, q, pos=pos, Lx=Lx, Ly=Ly)
    omegas = initialisation_omegas(N, sigma)
    psis   = initialisation_psis(N)
    vec_agents = Vector{Agent{float_type}}(undef, N)
    for n in 1:N
        vec_agents[n] = Agent{float_type}(Tuple(pos[:, n]), Tuple(pos[:, n]), thetas[n], omegas[n], psis[n])
    end
    if init_theta == "single" 
        @assert v0 == 0 "Error : v0 should be 0 if PBC not satisfied."
        periodic = false
        @assert phonons == false "Error : phonons cannot exist if PBC not satisfied."
    else 
        periodic = true
    end
    dt = determine_dt(T, sigma, v0, N, rho)
    return System(vec_agents, Lx, Ly, rho, Ntarget, N, T, sigma, v0, R0, dt, 0, params, periodic,phonons, phonon_amplitude, phonon_k, phonon_omega)
end

## -----------------  Get Methods ----------------- ##
get_pos(syst::System) = [syst.agents[n].pos for n in eachindex(syst.agents)]
get_theta(syst::System) = [syst.agents[n].theta for n in eachindex(syst.agents)]
get_omega(syst::System) = [syst.agents[n].omega for n in eachindex(syst.agents)]
get_psi(syst::System) = [syst.agents[n].psi for n in eachindex(syst.agents)]
# and for convenience :
get_thetas = get_theta
get_omegas = get_omega
get_psis = get_psi


## ----------------- Initialisation functions ----------------- ##
initialisation_psis(N) = 2π * rand(N)
initialisation_omegas(N, σ) = σ * randn(N)

function initialisation_positions(N, Lx, Ly, init_pos)
    if init_pos in ["random", "rand"]
        pos = rand(2, N)
        pos[1, :] *= Lx
        pos[2, :] *= Ly
    elseif init_pos in ["square_lattice", "square", "squ"]
        pos = zeros(2, N)
        for n in 1:N
            pos[1, n] = mod(n, Lx)
            pos[2, n] = Ly/Lx*mod(floor(n / Ly), Ly)
        end
    # elseif init_pos in ["triangular_lattice", "triangular","tri"]
    #     pos = zeros(2, N)
    #     for n in 1:N
    #         pos[1, n] = mod(n, Lx)
    #         pos[2, n] = mod(floor(n / Ly), Ly)
    #         if pos[2, n] % 2 == 1
    #             pos[1, n] += 0.5
    #         end
    #     end
    elseif init_pos in ["sobol", "Sobol"]
        pos = sobol(N, Lx, Ly)
    elseif init_pos in ["RSA", "rsa"] # Random Sequential Adsorption
        pos = rsa(N, Lx, Ly)
    elseif init_pos in ["PDS", "pds", "Poisson", "PoissonDisc", "PoissonDisk"] # Poisson Disc Sampling
        # Since the number of points generated is not always N, we try 10 times to generate as close to N points
        attempts = 1
        pos, generated_points = pds(N, Lx, Ly)
        for r in 1:attempts-1
            pos_tmp, generated_points_tmp = pds(N, Lx, Ly)
            closer_to_N_than_before = abs(N - generated_points) ≥ abs(N - generated_points_tmp)
            if closer_to_N_than_before
                pos = pos_tmp
                generated_points = generated_points_tmp
                if generated_points == N
                    break
                end
            end
        end
    else
        println("Error : init_pos should be \"random\" or \"square\" or \"sobol\" or \"RSA\".")
    end
    Neff = size(pos, 2)
    return pos,Neff
end

function pds(N, Lx, Ly)
    #= Poisson Disc Sampling
    A set of 2D points is said to be Poisson-disc if no two points 
    are closer than a given distance r.
    To acheive this rapidly, mainly to avoid computing N² distances everytime, 
    one maintains 2 datastrcutures :
    - a list of active points, from which we try to generate new points. 
    Note that a point is active if it has not yet generated k failed points. 
    - a grid of size Lx x Ly, where each cell contains the index of the points contained in it.
    The grid is used to check if a point is too close to another one.
    =#
    
    # Ugly but benchmarked so that for N ≤ 1E5 at ρ = 1, one is sure to obtain N points and no big void in the system
    if N > 1E4 constant = 0.791
    elseif N > 1E3 constant = 0.7925
    elseif N > 1E2 constant = 0.793
    else constant = 0.794
    end 
    r = constant*sqrt(Lx * Ly / N)
    k = 30

    pos = zeros(2, 2N) 
    #= Preallocate and take margin for safety 
        because the final number of generated points 
        is a random variable. In the end, only returns 
        actually generated points. =#
    generated_points = 0
    active = []
    nb_cells_x = ceil(Int, Lx / r)
    nb_cells_y = ceil(Int, Ly / r)
    grid = [Int[] for i in 1:nb_cells_x, j in 1:nb_cells_y]

    # initial point
    generated_points += 1 
    push!(active,generated_points) # list of points from which we try to generate new points
    x,y = rand(2) .* [Lx, Ly]
    pos[:, 1] = [x,y]
    push!(grid[ceil(Int, x / r), ceil(Int, y / r)],generated_points)

    # while generated_points < N && !isempty(active)
    while !isempty(active)
        i = rand(active)
        for j in 1:k
            if j == k  # Garde-fou contre les loops infinies
                deleteat!(active, findfirst(x -> x == i, active)) 
                # println("Point $i deleted from active list, which now contains $(length(active)) points.")
            end

            # Generate a new point
            theta = 2π * rand()
            rj = r * (1 + rand())
            x = mod(pos[1, i] + rj * cos(theta),Lx)
            y = mod(pos[2, i] + rj * sin(theta),Ly)

            # Get the cell of the new point
            cellx = ceil(Int, x / r)
            celly = ceil(Int, y / r)

            # Get neighbouring cells
            cells = get_neighbouring_cells(cellx, celly, nb_cells_x, nb_cells_y)

            # Get points ids in the neighbouring cells
            id_points = [grid[cell...] for cell in cells]
            id_points = vcat(id_points...)
            
            # Check if the new point is too close to another one
            too_close = false
            for p in id_points
                if dist(pos[:, p] , [x,y], Lx,Ly) < r
                    too_close = true
                    # println("too close")
                    break # exits this for loop
                end
            end

            # If the new point is not too close to another one, add it to the list of active points
            if !too_close
                generated_points += 1
                pos[:, generated_points] = [x,y]
                push!(active, generated_points)
                push!(grid[cellx, celly], generated_points)
                # println("Point $generated_points added to active list, which now contains $(length(active)) points.")
                break # exits the "for j in 1:k" loop, pick another active point and try to generate a new point from it
            end
        end
    end
    println("PDS generated $generated_points points.")
    return pos[:,1:generated_points],generated_points
end


function sobol(N, Lx, Ly)
    @assert Lx == Ly "Error : Lx and Ly should be equal for Sobol initialisation"
    dimm = 10
    good_sobol_axes = [(1, 3), (1, 4), (1, 8), (2, 3), (2, 4), (2, 8), (3, 1), (3, 2), (3, 5), (3, 8), (3, 7), (4, 1), (4, 2), (4, 5), (4, 7), (5, 3), (6, 5), (6, 9), (7, 3), (7, 4), (7, 8), (7, 9), (8, 1), (8, 2), (8, 3), (8, 7), (9, 6), (9, 7), (5, 10), (7, 10), (10, 5), (10, 7)]
    # this 32-element list of good sobol axis (for a total dimension of dimm = 10) was obtained by visual inspection. We wanted to have a good spread of points in the 2D plane, with low discrepancy.
    sobol_axis = rand(good_sobol_axes)
    sobol = SobolSeq(dimm) # points in [0,1]^2 (therefore we multiply by Lx and Ly and Lx should be equal to Ly)
    p = reduce(hcat, next!(sobol) for i = 1:N)
    pos = zeros(2, N)
    pos[1, :] = p[sobol_axis[1], :] * Lx
    pos[2, :] = p[sobol_axis[2], :] * Ly

    return pos
end

function rsa(N, Lx, Ly)
    #= The algorithm goes as follows. 
    1. Draw an initial random point in the box and add it to the final list.
    2. Define a radius (imagine points as the centers of physical coins)
    radius = sqrt(Lx * Ly / N / π) * constant
    sqrt(Lx * Ly / N / π) is the radius of a circle that would fit N points in the box
    constant is a parameter that we can tune to make the radius smaller or bigger
    3. While the number of rejections is smaller than a maximum number of rejections
        3.1. Draw a random point in the box
        3.2. If the point is at a distance of at least 2 * radius from all the points in the final list, add it to the final list
        3.3. If the final list has N points, stop
        3.4. If the number of rejections is bigger than the maximum number of rejections, decrease the radius by 10% and loop again until the final list has N points.
    =#
    pos = Tuple{Float32,Float32}[]
    max_number_rejections = 1000
    
    number_rejections = 0
    constant = 1
    while number_rejections < max_number_rejections
        radius = sqrt(Lx * Ly / N / π) * constant
        x = Float32(Lx*rand())
        y = Float32(Ly*rand())
        if all(dist((x,y),position,Lx,Ly) > 2radius for position in pos)
            push!(pos, (x, y))
            if length(pos) == N
                break
            end
        else 
            number_rejections += 1
            if number_rejections ≥ max_number_rejections
                constant = 0.9 * constant # skrink all the disks and try again for a maximum of max_number_rejections rejections
                number_rejections = 0
            end
        end
    end
    pos = vecTuple2matrix(pos)
    return pos
end

function initialisation_thetas(N, init_theta, r0, q, nb_horizontal_waves=1, nb_vertical_waves=1; pos, Lx, Ly)
    if init_theta in ["disordered", "hightemp", "random"]
        thetas = 2π * rand(N)
    elseif init_theta in ["ordered", "lowtemp"]
        thetas = zeros(N)
    elseif init_theta in ["isolated", "single"]
        thetas = zeros(N)
        for n in 1:N
            x, y = pos[:, n]
            thetas[n] = q * atan(Ly / 2 - y, Lx / 2 - x)
        end
    elseif init_theta == "pair"
        thetas = zeros(N)
        for n in 1:N
            x, y = pos[:, n]
            thetas[n] = q * (atan(Ly / 2 - y, Lx / 2 - x + r0 / 2) - atan(Ly / 2 - y, Lx / 2 - x - r0 / 2)) # +1 defect by convention
            # smooth the transition between top and bottom
            if (y < 0.05 * Ly) || (y > 0.95 * Ly)
                thetas[n] = 0
            end
        end
    elseif init_theta in ["1D_horizontal_wave", "1D_hor", "hor_wave", "wave_hor", "horizontal_wave"]
        thetas = zeros(N)
        for n in 1:N
            thetas[n] = mod(nb_horizontal_waves * 2pi * pos[1, n] / Lx, 2pi)
        end
    elseif init_theta in ["1D_vertical_wave", "1D_vert", "vert_wave", "wave_vert", "vertical_wave"]
        for n in 1:N
            thetas[n] = mod(nb_vertical_waves * 2pi * pos[2, n] / Ly, 2pi)
        end
    elseif init_theta == "2Dwave"
        thetas = zeros(N)
        for n in 1:N
            thetas[n] = mod(π .+ π / 2 * (sin(nb_vertical_waves * pos[1, n] * 2pi / Lx) + sin(nb_horizontal_waves * pos[2, n] * 2pi / Lx)), 2pi)
        end
    else
        error("Error : init_theta should be \"disordered\", \"ordered\", \"isolated\", \"pair\", \"1D_horizontal_wave\", \"1D_vertical_wave\" or \"2Dwave\"")
    end
    return thetas
end

function determine_dt(T, σ, v0, N, rho)
    estimation_max = sqrt(2log(N)) * σ # from math.stackexchange.com/questions/89030, simple answer of Sivaraman
    arbitrary_coeff = π / 20
    return minimum(
            [1 / 5 / v0,
            arbitrary_coeff / estimation_max,
            arbitrary_coeff / π / rho,
            arbitrary_coeff^2 * π / 4 / T])
end

function effective_number_particle(Ntarget, rho, aspect_ratio=1)
    Lx = round(Int, sqrt(aspect_ratio * Ntarget / rho))
    Ly = round(Int, sqrt(Ntarget / rho / aspect_ratio))
    N_effective = round(Int, Lx * Ly * rho)
    return N_effective, Lx, Ly
end
effective_number_particles = effective_number_particle

#################### Visualisation methods ####################
function plot_thetas(system; particles=false, vertical=false, 
    size=(512, 512), defects=false, title="", nb_neighbours = false, 
    cols = cgrad([:black, :blue, :green, :orange, :red, :black]))

    pos, thetas, N, Lx, Ly = get_pos(system), get_thetas(system), system.N, system.Lx, system.Ly
    
    if particles 
        p_particles = scatter(pos, marker_z=mod.(thetas, 2pi), color=cols, clims=(0, 2pi), ms=275 / Lx, size=size, aspect_ratio=Ly / Lx, xlims=(0, Lx), ylims=(0, Ly))
    end

    if nb_neighbours
        msss = 1
        list_nnn = length.(get_list_neighbours(system))    
        p_nb_neighbours = plot(aspect_ratio=1)
        scatter!(pos, marker_z=list_nnn, 
        markerstrokewidth=0, markersize=msss,size=size,
        c=cgrad([:black, :blue, :green, :orange, :red]))
        if defects 
            defects_p, defects_m =  spot_defects(system)
            for defp in defects_p
                scatter!(defp[1:2].*system.R0,marker=:circle,ms=5)
            end
            for defm in defects_m
                scatter!(defm[1:2].*system.R0,marker=:circle,ms=5)
            end
        end
    end
    thetas_cg = cg(system)
    p_cg = heatmap(mod.(thetas_cg, 2pi)', clims=(0, 2pi), c=cols, size=size, aspect_ratio=Ly / Lx, xlims=(0, Lx/system.rho/system.R0), ylims=(0, Ly/system.rho/system.R0))
    if defects 
        highlight_defects!(p_cg, system)
    end

    if vertical 
        if particles && nb_neighbours
            final_plot = plot(p_particles,p_cg, p_nb_neighbours, layout=(3, 1), size=(size[1], size[2] * 3), title=title)
        elseif particles
            final_plot = plot(p_particles, p_cg, layout=(2, 1), size=(size[1], size[2] * 2), title=title)
        elseif nb_neighbours
            final_plot = plot(p_cg, p_nb_neighbours, layout=(2, 1), size=(size[1], size[2] * 2), title=title)
        else
            final_plot = plot(p_cg, layout=(1, 1), size=size, title=title)
        end
    else 
        if particles && nb_neighbours
            final_plot = plot(p_particles, p_cg, p_nb_neighbours, layout=(1, 3), size=(size[1] * 3, size[2]), title=title)
        elseif particles
            final_plot = plot(p_particles, p_cg, layout=(1, 2), size=(size[1] * 2, size[2]), title=title)
        elseif nb_neighbours
            final_plot = plot(p_cg, p_nb_neighbours, layout=(1, 2), size=(size[1] * 2, size[2]), title=title)
        else
            final_plot = plot(p_cg, layout=(1, 1), size=size, title=title)
        end
    end 
    # if particles
    #     if vertical
    #         p1 = scatter(pos, marker_z=mod.(thetas, 2pi), color=cols, clims=(0, 2pi), ms=275 / Lx, size=size, aspect_ratio=Ly / Lx, xlims=(0, Lx), ylims=(0, Ly))
    #         final_plot = plot(p1, p2, layout=(2, 1), size=(size[1], size[2] * 2), title=title)
    #     else
    #         p1 = scatter(pos, marker_z=mod.(thetas, 2pi), color=cols, clims=(0, 2pi), ms=275 / Lx, size=size, aspect_ratio=Ly / Lx, xlims=(0, Lx), ylims=(0, Ly))
    #         thetas_cg = cg(system)
    #         p2 = heatmap(mod.(thetas_cg, 2pi)', clims=(0, 2pi), c=cols, size=size, aspect_ratio=Ly / Lx, xlims=(0, Lx/system.rho/system.R0), ylims=(0, Ly/system.rho/system.R0))
    #         final_plot = plot(p1, p2, layout=(1, 2), size=(size[1] * 2, size[2]), title=title)
    #     end
    # else
    #     thetas_cg = cg(system)
    #     final_plot = heatmap(mod.(thetas_cg, 2pi)', clims=(0, 2pi), c=cols, size=size, aspect_ratio=Ly / Lx, xlims=(0, Lx/system.rho/system.R0), ylims=(0, Ly/system.rho/system.R0), title=title)
    #     if defects 
    #         highlight_defects!(final_plot, system)
    #     end
    # end
    return final_plot
end

function highlight_defects!(p, system, symbP=:circle, symbM=:utriangle)
    defects_p, defects_m = spot_defects(system)
    for defect in defects_p
        # scatter!((defect[2],Ly-defect[1]), m=(1.5, 6., symbP, :transparent, stroke(1.2, :grey85)))
        scatter!((defect[1:2] ./ system.R0), m=(1.5, 6., symbP, :transparent, stroke(1.2, :grey85)))
    end
    for defect in defects_m
        scatter!((defect[1:2] ./ system.R0), m=(1.5, 6., symbM, :transparent, stroke(1.2, :grey85)))
    end
    # xlims!((1, Lx))
    # ylims!((1, Ly))
    return p
end

function mod1_2D(xx::Tuple{T,T}, Lx::Int, Ly::Int) where {T<:Number}
    return (mod1(xx[1], Lx), mod1(xx[2], Ly))
end


function cg(system::System{T}) where {T<:AbstractFloat}
    Lx, Ly = system.Lx, system.Ly
    N = system.N
    R0 = system.R0
    pos, thetas = get_pos(system), get_thetas(system)
    mesh_size = R0
    cutoff = 4R0 # for contributions

    list, head = construct_cell_list(system) 
    
    nb_cells_x = ceil(Int,Lx / mesh_size)
    nb_cells_y = ceil(Int,Ly / mesh_size)
    fine_grid = NaN * zeros(nb_cells_x, nb_cells_y)
    fine_grid_density = NaN * zeros(nb_cells_x, nb_cells_y)
    for i in 1:nb_cells_x, j in 1:nb_cells_y # scan fine_grid cells
        center_finegrid_cell = T.([(i - 0.5) * mesh_size, (j - 0.5) * mesh_size])

        # find cell from coarse mesh correponding to cell i,j belonging to fine_grid
        cellx, celly = Int(div(center_finegrid_cell[1], mesh_size)) + 1, Int(div(center_finegrid_cell[2], mesh_size)) + 1
        a = round(Int, cutoff / mesh_size)
        neighbouring_cells = vec([(ii, jj) for ii in -a:a, jj in -a:a])
        neighbouring_cells = [mod1_2D(el .+ (cellx, celly), nb_cells_x, nb_cells_y) for el in neighbouring_cells]

        # get indices of particles within those cells (cells belonging to the coarse mesh)
        ind_neighbours = Int[]
        for (i, j) in neighbouring_cells
            next = head[i, j]
            if next ≠ -1
                if dist(center_finegrid_cell, pos[next], Lx, Ly) < cutoff
                    push!(ind_neighbours, next)
                end
                while list[next] ≠ -1
                    if dist(center_finegrid_cell, pos[list[next]], Lx, Ly) < cutoff
                        push!(ind_neighbours, list[next])
                    end
                    next = list[next]
                end
            end
        end

        # compute contributions from those "neighbouring" particles
        tmp = ComplexF32[]
        sizehint!(tmp, length(ind_neighbours))
        tmp_density = Float64[]
        for m in ind_neighbours
            r = dist(center_finegrid_cell, pos[m], Lx, Ly)
            if r < cutoff
                push!(tmp, exp(im * thetas[m] - r / R0))
                # push!(tmp_density,exp(-r/R0))
            end
        end
        fine_grid[i, j] = angle(mean(tmp))
        # fine_grid_density[i,j] = sum(tmp_density)
    end
    return fine_grid
end


## Time Evolution
construct_cell_list(system) = construct_cell_list(get_pos(system), system.N, system.Lx, system.Ly, system.R0)
function construct_cell_list(pos::Vector{Tuple{T,T}}, N::Int, Lx::Int, Ly::Int, R0::Number)::Tuple{Vector{Int},Matrix{Int}} where {T<:AbstractFloat}
    nb_cells_x = ceil(Int,Lx/R0)
    nb_cells_y = ceil(Int,Ly/R0)
    head = -ones(Int, nb_cells_x, nb_cells_y) # head[i,j] contains the index of the first particle in cell (i,j). -1 if empty
    list = -ones(Int, N) # list[n] contains the index of the particle to which particle n points. -1 if it points to no one
    for n in 1:N
        cellx, celly = floor(Int,pos[n][1]/R0) + 1, floor(Int,pos[n][2]/R0) + 1 # cell to which particle n belongs
        list[n] = head[cellx, celly]
        head[cellx, celly] = n
    end
    return list, head
end

function get_neighbouring_cells(cellx::Int, celly::Int, nb_cells_x::Int, nb_cells_y::Int)::Vector{Vector{Int}}
    # In the end, this function is quite fast, it contributes +3ms for N=1E4 particles
    should_take_mod = (cellx == 1) || (cellx == nb_cells_x) || (celly == 1) || (celly == nb_cells_y)
    if should_take_mod
        neighbouring_cells = Vector{Int}[[cellx, celly], [cellx, mod1(celly + 1, nb_cells_y)], [mod1(cellx + 1, nb_cells_x), celly], [cellx, mod1(celly - 1, nb_cells_y)], [mod1(cellx - 1, nb_cells_x), celly], [mod1(cellx + 1, nb_cells_x), mod1(celly + 1, nb_cells_y)], [mod1(cellx - 1, nb_cells_x), mod1(celly - 1, nb_cells_y)], [mod1(cellx - 1, nb_cells_x), mod1(celly + 1, nb_cells_y)], [mod1(cellx + 1, nb_cells_x), mod1(celly - 1, nb_cells_y)]]
    else
        neighbouring_cells = Vector{Int}[[cellx, celly], [cellx, celly + 1], [cellx + 1, celly], [cellx, celly - 1], [cellx - 1, celly], [cellx + 1, celly + 1], [cellx - 1, celly - 1], [cellx - 1, celly + 1], [cellx + 1, celly - 1]]
    end

    #= For non-integer L/R0, part of the last cell is empty. 
    This implies that some of the agents of the first row
    interact not only with the last row but also with the before-last row.=# 
    if cellx == 1 push!(neighbouring_cells, [nb_cells_x-1, celly], [nb_cells_x-1, mod1(celly+1, nb_cells_y)], [nb_cells_x-1, mod1(celly-1, nb_cells_y)]) end
    if celly == 1 push!(neighbouring_cells, [cellx, nb_cells_y-1], [mod1(cellx+1, nb_cells_x), nb_cells_y-1], [mod1(cellx-1, nb_cells_x), nb_cells_y-1]) end    
    # Note that this problem does not exist for the last row, which only interacts with the first row and not the second one 
    
    return neighbouring_cells

    # Code below as fast but less clear
    # if should_take_mod
    #     neighbouring_cells = Tuple{Int,Int}[ (cellx,celly) , (cellx,mod1(celly+1,nb_cells_y)) , (mod1(cellx+1,nb_cells_x),celly) , (cellx,mod1(celly-1,nb_cells_y)) , (mod1(cellx-1,nb_cells_x),celly) , (mod1(cellx+1,nb_cells_x),mod1(celly+1,nb_cells_y)) ,  (mod1(cellx-1,nb_cells_x),mod1(celly-1,nb_cells_y)) , (mod1(cellx-1,nb_cells_x),mod1(celly+1,nb_cells_y)) , (mod1(cellx+1,nb_cells_x),mod1(celly-1,nb_cells_y))]
    # else
    #     neighbouring_cells = Tuple{Int,Int}[ (cellx,celly) , (cellx,celly+1) , (cellx+1,celly) , (cellx,celly-1) , (cellx-1,celly) , (cellx+1,celly+1) ,  (cellx-1,celly-1) , (cellx-1,celly+1) , (cellx+1,celly-1)]
    # end
    # return CartesianIndex.(neighbouring_cells)

    # Code below is 2x slower
    # neighbouring_cells = [[cellx,celly]] .+ [[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]]
    # neighbouring_cells = [mod1.(neighbouring_cells[i],nb_cells_x) for i in 1:9]
end


get_list_neighbours(system::System) = get_list_neighbours(get_pos(system), system.N, system.Lx, system.Ly, system.R0)
function get_list_neighbours(pos::Vector{Tuple{T,T}}, N::Int, Lx::Int, Ly::Int, R0::Number) where {T<:AbstractFloat}
    ind_neighbours = Vector{Vector{Int}}(undef, N)
    list, head = construct_cell_list(pos, N, Lx, Ly, R0)
    nb_cells_x, nb_cells_y = size(head)
    R02 = R0^2

    for n in 1:N
        cellx, celly = Int(div(pos[n][1], R0)) + 1, Int(div(pos[n][2], R0)) + 1 # cell to which particle n belongs
        poscur = pos[n]
        ind_neighbours_n = Int[]
        sizehint!(ind_neighbours_n, 10) # average number of neighbours is ~ 4ρ
        for (i, j) in get_neighbouring_cells(cellx, celly, nb_cells_x, nb_cells_y)
            # After some investigation, the line above is not slow, what takes time is the cell algo itself
            next = head[i, j]
            if next ≠ -1
                if 0 < dist2(poscur, pos[next], Lx, Ly) ≤ R02
                    push!(ind_neighbours_n, next)
                end
                while list[next] ≠ -1
                    if 0 < dist2(poscur, pos[list[next]], Lx, Ly) ≤ R02
                        push!(ind_neighbours_n, list[next])
                    end
                    next = list[next]
                end
            end
        end
        ind_neighbours[n] = ind_neighbours_n
    end
    return ind_neighbours
end

function get_number_neighbours(system::System)
    ind_neighbours = get_list_neighbours(get_pos(system), system.N, system.Lx, system.Ly, system.R0)
    # number_of_neighbours = [length(ind_neighbours[n]) for n in eachindex(ind_neighbours)]
    return length.(ind_neighbours)
end

function evolve!(system::System, tmax::Number)
    v0,dt = system.v0, system.dt
    if v0 == 0 && !system.phonons
        ind_neighbours = get_list_neighbours(system)
        if system.periodic
            while system.t < tmax
                system.t += dt
                update_thetas!(system, ind_neighbours)
            end
        else
            while system.t < tmax
                system.t += dt
                update_thetas_nonperiodic!(system, ind_neighbours)
            end
        end
    elseif v0 == 0 && system.phonons
        # If phonons are present, the system is necessarily periodic 
        # (has been checked in the constructor)
        while system.t < tmax
            system.t += dt
            update_positions_phonons!(system)
            ind_neighbours = get_list_neighbours(system)
            update_thetas!(system, ind_neighbours)
        end
    else
        while system.t < tmax
            system.t += dt
            update_positions!(system)
            ind_neighbours = get_list_neighbours(system)
            update_thetas!(system, ind_neighbours)
        end
    end
    return system
end

function check_agents_in_box(system)
    Lx, Ly = system.Lx, system.Ly
    for n in 1:system.N
        pos = system.agents[n].pos
        if pos[1] < 0 || pos[1] ≥ Lx || pos[2] < 0 || pos[2] ≥ Ly
            error("Agent $n is out of the box at time $(system.t), position $(pos), Lx = $(system.Lx), Ly = $(system.Ly)")
        end
    end
end

function update!(system::System)
    update_positions!(system)
    ind_neighbours = get_list_neighbours(system)
    update_thetas!(system, ind_neighbours)
end

mod0(x::T, y::Int) where {T<:Real} = (m = mod(x, y); ifelse(m == y, zero(m), m)) 
function update_positions!(system::System{T}) where {T<:AbstractFloat}
    v0, dt = system.v0, system.dt
    Lx, Ly = system.Lx, system.Ly
    for n in 1:system.N
        agent = system.agents[n]
        x = agent.pos[1] + v0 * dt * cos(agent.psi)
        y = agent.pos[2] + v0 * dt * sin(agent.psi)
        x = round(x, digits=5)
        y = round(y, digits=5)
        agent.pos = (mod0(x,Lx), mod0(y,Ly))
        #= Use of custom modulo function above is required because
            if a particle passes from +ε to -ε in the x or y direction, 
            with ε < EPSILON(AbstractFloat), then the result of the standard  
            mod functuion is contra-intuitive. This issue has already been debated, 
            cf. https://github.com/JuliaLang/julia/issues/36310. 
            Since it appears that a simple fix does not exist, 
            I use the suggested workaround. 
            Adds 200µs for an update of 10^4 particles, to be compared to 3.750 ms in total (+2%)
            =#
    end
end

function update_positions_phonons!(system::System{T}) where {T<:AbstractFloat}
    Lx, Ly = system.Lx, system.Ly
    A = system.phonon_amplitude
    k = system.phonon_k
    omega = system.phonon_omega
    omega_t = omega*system.t
    for n in 1:system.N
        agent = system.agents[n]
        x = agent.pos0[1] + A*cos(k*agent.pos0[1] - omega_t)
        x = round(x, digits=5)
        y = agent.pos0[2] + A*cos(k*agent.pos0[2] - omega_t)
        y = round(y, digits=5)
        agent.pos = (mod0(x,Lx),mod0(y,Ly))
    end
end
function update_thetas!(system::System{FT}, ind_neighbours::Vector{Vector{Int}}) where {FT<:AbstractFloat}
    dt, T = system.dt, system.T
    for n in 1:system.N
        agent = system.agents[n]
        if length(ind_neighbours[n]) > 0
            thetas_neighbours = [system.agents[i].theta for i in ind_neighbours[n]]
            sumsin = 0.0
            for theta in thetas_neighbours
                sumsin += sin(theta - agent.theta)
            end
            agent.theta += dt * (agent.omega + sumsin) + sqrt(2T * dt) * randn()
        else
            agent.theta += dt * agent.omega + sqrt(2T * dt) * randn()
        end
    end
end
function update_thetas_nonperiodic!(system::System{FT}, ind_neighbours::Vector{Vector{Int}}) where {FT<:AbstractFloat}
    dt, T = system.dt, system.T
    R0 = system.R0
    Lx,Ly = system.Lx, system.Ly
    for n in 1:system.N
        agent = system.agents[n]
        x,y = agent.pos
        cellx = ceil(Int, x/R0)
        celly = ceil(Int, y/R0)
        nb_cells_x = ceil(Int, Lx/R0)
        nb_cells_y = ceil(Int, Ly/R0)
        on_the_edge = cellx == 1 || cellx == nb_cells_x || celly == 1 || celly == nb_cells_y
        if !on_the_edge
            if length(ind_neighbours[n]) > 0
                thetas_neighbours = [system.agents[i].theta for i in ind_neighbours[n]]
                sumsin = 0.0
                for theta in thetas_neighbours
                    sumsin += sin(theta - agent.theta)
                end
                agent.theta += dt * (agent.omega + sumsin) + sqrt(2T * dt) * randn()
            else
                agent.theta += dt * agent.omega + sqrt(2T * dt) * randn()
            end
        else 
            # do nothing, leave the edges unchanged thoughout the simulation
        end
    end
end

## Measurements
polarOP(system::System) = polarOP(get_thetas(system))
function polarOP(thetas::Vector{T}) where {T<:AbstractFloat}
    tmp = Complex(0)
    for theta in thetas
        tmp += exp(1im * theta)
    end
    return norm(tmp) / length(thetas), angle(tmp)
end

nematicOP(system::System) = nematicOP(get_thetas(system))
function nematicOP(thetas::Vector{T}) where {T<:AbstractFloat}
    tmp = Complex(0)
    for theta in thetas
        tmp += exp(2im * theta)
    end
    return norm(tmp) / length(thetas), angle(tmp)
end


corr(system::System; dr=system.R0, algo="fast") = corr(get_pos(system), get_thetas(system), system.N, system.Lx, system.Ly, dr)
function corr(pos::Vector{Tuple{T,T}}, thetas::Vector{T}, N, Lx, Ly, dr; algo="fast")::Vector{T} where {T<:AbstractFloat}
    if algo == "slow"
        return corr_slow(pos, thetas, N, Lx, Ly, dr)
    elseif algo == "fast"
        return corr_fast(pos, thetas, N, Lx, Ly, dr)
    end
end


function corr_slow(pos::Vector{Tuple{T,T}}, thetas::Vector{T}, N, Lx, Ly, dr)::Vector{T} where {T<:AbstractFloat}
    #= Iterates over all (unique) pairs of particles (hence i in j+1:N), and computes the correlation function
    Regroup the results in bins of size dr, and compute the average correlation in each bin. =#
    Lmin = min(Lx, Ly)
    C = [T[] for i in 1:round(Int, Lmin / 2 / dr)]
    for j in 1:N, i in j+1:N
        d = dist(pos[i], pos[j], Lx, Ly)
        if d ≤ round(Int, Lmin / 2)
            push!(C[min(ceil(Int, d / dr), length(C))], cos(thetas[i] - thetas[j]))
        end
    end

    Cavg = [mean(C[i]) for i in eachindex(C)]
    return vcat(1, Cavg) # 1 to represent C(0,t)
end

function corr_fast(pos::Vector{Tuple{T,T}}, thetas::Vector{T}, N, Lx, Ly, dr)::Vector{T} where {T<:AbstractFloat}
    #= Same as before, but stops the computation when the number of pairs in all bins exceeds M =#
    M = 50 # benchmarked so that corr_fast is quantitatively equivalent to corr_slow
    Lmin = min(Lx, Ly)
    C = [T[] for i in 1:round(Int, Lmin / 2 / dr)]
    ms = zeros(Int, round(Int, Lmin / 2 / dr))

    for j in 1:N, i in j+1:N
        d = dist(pos[i], pos[j], Lx, Ly)
        if 0 < d ≤ round(Int, Lmin / 2)
            ind = min(ceil(Int, d / dr), length(C))
            push!(C[ind], cos(thetas[i] - thetas[j]))
            ms[ind] += 1
        end
        if sum(ms .> M) == length(ms)
            break
        end
    end

    Cavg = [mean(C[i]) for i in eachindex(C)]
    return vcat(1, Cavg) # 1 to represent C(0,t)
end

## Methods for defects


function arclength(theta1::T, theta2::T)::T where {T<:AbstractFloat}
    #= This function returns the signed arclength on the unit trigonometric circle .
    Clockwise        >> sign -
    Counterclockwise >> sign +
    WARNING
    Note that the inputs thetas need to lie within [0,π] or [0,2π], depending on the symmetry of the model =#
    dtheta = theta2 - theta1
    dtheta_abs = abs(theta2 - theta1)

    shortest_unsigned_arclength = min(2π - dtheta_abs, dtheta_abs)
    if dtheta_abs ≤ π
        signe = sign(dtheta)
    else
        signe = -sign(dtheta)
    end
    return signe * shortest_unsigned_arclength
end

function get_neighbours(thetas::Matrix{<:T}, i::Int, j::Int)::Vector{T} where {T<:AbstractFloat}
    Lx, Ly = size(thetas)
    # convention depuis la droite et sens trigo
    bulk = is_in_bulk(i, j, Lx, Ly)
    if bulk
        jm = j - 1
        jp = j + 1
        imm = i - 1
        ip = i + 1
    else
        jm = mod1(j - 1, Ly)
        jp = mod1(j + 1, Ly)
        imm = mod1(i - 1, Lx)
        ip = mod1(i + 1, Lx)
    end

    @inbounds angles =
        [thetas[i, jp],
            thetas[imm, j],
            thetas[i, jm],
            thetas[ip, j]]

    return angles
end

function get_angles_plaquette(thetas::Matrix{<:T}, i::Int, j::Int)::Vector{T} where {T<:AbstractFloat}
    Lx, Ly = size(thetas)
    # By convention, the 4 angles are taken from the topleft corner and counterclockwise
    bulk = is_in_bulk(i, j, Lx, Ly)
    if bulk
        jp = j + 1
        ip = i + 1
    else
        jp = mod1(j + 1, Ly)
        ip = mod1(i + 1, Lx)
    end

    @inbounds angles_plaquette =
        [thetas[i, j],
            thetas[ip, j],
            thetas[ip, jp],
            thetas[i, jp]]

    return angles_plaquette
end

is_on_border(i::Int, j::Int, Lx::Int, Ly::Int) = (i == 1) || (j == 1) || (i == Lx) || (j == Ly)
is_in_bulk(i::Int, j::Int, Lx::Int, Ly::Int) = !is_on_border(i, j, Lx, Ly)

function get_vorticity(thetas_mat::Matrix{T}, i::Int, j::Int, Lx::Int, Ly::Int)::T where {T<:AbstractFloat}
    # Note : thetas_mat = mod.(cg(thetas),2π) # coarsegrained thetas
    angles_corners = get_angles_plaquette(thetas_mat, i, j)
    perimeter_covered = 0.0
    nnn = length(angles_corners)
    for i in 1:nnn
        perimeter_covered += arclength(angles_corners[i], angles_corners[mod1(i+1,nnn)])
    end
    charge = round(perimeter_covered / 2π, digits=1)
    return charge
end


function number_defects(system::System)
    defects = spot_defects(system)
    return length(defects[1]) + length(defects[2])
end
number_defectsP(system::System) = length(spot_defects(system)[1])
number_defectsN(system::System) = length(spot_defects(system)[2])
number_positive_defects = number_defects_pos = number_defectsP
number_negative_defects = number_defects_neg = number_defectsN

function spot_defects(system::System{T}) where {T<:AbstractFloat}  
    vortices_plus = Tuple{T,T,T}[]
    vortices_minus = Tuple{T,T,T}[]
    R0 = system.R0

    thetasmod = mod.(cg(system), 2π) # Coarse-graining 
    Lx,Ly = size(thetasmod)

    if system.periodic 
        range_i = 1:Lx
        range_j = 1:Ly
    else
        range_i = 3:Lx-2 # to prevent artifical detection of defects on the borders (2:L-1 was not enough in some L,R0 situations)
        range_j = 3:Ly-2
    end
    
    for i in range_i
        for j in range_j
            q = get_vorticity(thetasmod, i, j, Lx, Ly)
            x,y = (R0).*(i+0.5 ,j+0.5) # 
            # NB1 : since defects live on the dual lattice, the coordinates are shifted by 0.5
            # NB2 : defects are localized in the coarse grained lattice, so we need to multiply by R0 to recover the real coordinates
            if q > +0.1
                push!(vortices_plus, (x,y,q)) # we want to keep ±½ defects, and not rounding errors
            elseif q < -0.1
                push!(vortices_minus, (x,y,q))
            end
        end
    end
    
    return vortices_plus, vortices_minus
end


function relax!(thetas::Matrix{FT}, trelax=0.3) where {FT<:AbstractFloat}
    t = 0.0
    dt = 1E-2
    # T = 0.05
    Lx, Ly = size(thetas)

    thetas_old = copy(thetas)

    while t < trelax
        t += dt
        for j in 1:Ly
            for i in 1:Lx
                θ = thetas_old[i, j]
                angle_neighbours = get_neighbours(thetas_old, i, j)
                thetas[i, j] += dt * sum(sin, angle_neighbours .- θ)
                # thetas[i,j] += dt*sum(sin,angle_neighbours .- θ) + sqrt(2T*dt)*randn(FT)
            end
        end
    end

    return thetas
end


function merge_duplicates(list, Lx, Ly; radius=3)
    #= In this list, there might be doubles/triples (2/3 locations for the
    same physical vortex). We thus seek for numerically identified vortices
    which are neighbours and with the same charge to delete them. =#
    pos = [list[i][1:2] for i in eachindex(list)]
    charge = [list[i][3] for i in eachindex(list)]
    dealt_with = falses(length(pos))

    merged_duplicates = []
    for i in 1:length(pos)
        if !dealt_with[i]
            tmp = []
            for j in i:length(pos) # includes defect "i"
                if dist(pos[i], pos[j], Lx, Ly) ≤ radius && charge[i] == charge[j]
                    dealt_with[j] = true
                    push!(tmp, pos[j])
                end
            end
            mean_loc_defect = mean_N_positions(tmp, Lx, Ly, true)
            push!(merged_duplicates, (mean_loc_defect[1], mean_loc_defect[2], charge[i]))
        end
    end
    return merged_duplicates
end

## In vivo Defect Tracking
mutable struct Defect
    id::Int
    charge::Number
    pos::Vector{Tuple{Number,Number}}
    annihilation_time::Union{Float64,Nothing}
    creation_time::Float64
    id_annihilator::Union{Int,Nothing}
end
Defect(; id, charge, loc, t) = Defect(id, charge, [loc], nothing, t, nothing)
first_loc(d::Defect) = d.pos[1]
last_loc(d::Defect) = d.pos[end]
function update_position!(d::Defect, new_loc)
    push!(d.pos, new_loc)
end


mutable struct DefectTracker
    defectsP::Vector{Defect} # the id of a defect is its index in this vector
    defectsN::Vector{Defect} # so there is a (+)defect with id=1 AND and a (-)defect with id=1
    current_time::Float64 # latest update time (by convention, the creation time of the whole data structure = 0)

    function DefectTracker(pos, thetas, N, Lx, Ly, t) # constructor
        vortices, antivortices = spot_defects(pos, thetas, N, Lx, Ly)
        defectsP = [Defect(id=i, charge=vortices[i][3], loc=vortices[i][1:2], t=t) for i in each(vortices)]
        defectsN = [Defect(id=i, charge=antivortices[i][3], loc=antivortices[i][1:2], t=t) for i in each(antivortices)]
        new(defectsP, defectsN, t)
    end

    function DefectTracker(system, t) # constructor
        vortices, antivortices = spot_defects(system)
        defectsP = [Defect(id=i, charge=vortices[i][3], loc=vortices[i][1:2], t=t) for i in each(vortices)]
        defectsN = [Defect(id=i, charge=antivortices[i][3], loc=antivortices[i][1:2], t=t) for i in each(antivortices)]
        new(defectsP, defectsN, t)
    end
end
number_defectsP(dft::DefectTracker) = length(dft.defectsP)
number_defectsN(dft::DefectTracker) = length(dft.defectsN)
number_defects(dft::DefectTracker) = length(dft.defectsN) + length(dft.defectsN)
number_active_defectsP(dft::DefectTracker) = count(isnothing, [d.annihilation_time for d in dft.defectsP])
number_active_defectsN(dft::DefectTracker) = count(isnothing, [d.annihilation_time for d in dft.defectsN])
number_active_defects(dft::DefectTracker) = number_active_defectsN(dft) + number_active_defectsP(dft)


function ID_active_defects(dft::DefectTracker)
    activeP = Int[]
    for i in 1:number_defectsP(dft)
        if dft.defectsP[i].annihilation_time == nothing
            push!(activeP, i)
        end
    end
    activeN = Int[]
    for i in 1:number_defectsN(dft)
        if dft.defectsN[i].annihilation_time == nothing
            push!(activeN, i)
        end
    end
    return activeP, activeN
end

function add_defect!(dft::DefectTracker; charge, loc)
    if charge > 0
        push!(dft.defectsP, Defect(id=1 + number_defectsP(dft), charge=charge, loc=loc, t=dft.current_time))
    else
        push!(dft.defectsN, Defect(id=1 + number_defectsN(dft), charge=charge, loc=loc, t=dft.current_time))
    end
end

function pair_up_hungarian(dt::DefectTracker, new, old, Lx, Ly, charge::String)
    # charge can be either "+" or "-"
    distance_matrixx = distance_matrix(new, old, Lx, Ly) # m_new lignes, m_old colonnes
    proposal = hungarian(distance_matrixx)[1] # length = length(new)
    assignment = copy(proposal) # because it will be modified in the next for loop

    #= A few comments :
    1. The proposal is the match between indices of the vectors
    new,old while the assignment matches actual IDs of the DefectTracker.
    2. If cost_matrix is a NxM matrix (workers x jobs), the output of hungarian(cost_matrix)
    is a Nx1 vector containing the assignments of each of the N workers to the indices of the jobs (1 < indices < M).
    A "0" means the worker has not been assigned to any job, hence the assignment[i] ≠ 0 condition below.
    =#
    if charge == "+"
        for i in eachindex(assignment)
            for j in 1:number_defectsP(dt)
                if assignment[i] ≠ 0 && dt.defectsP[j].annihilation_time == nothing && last_loc(dt.defectsP[j]) == old[proposal[i]]
                    # the first condition is a protection against the creation case, where the Hungarian algo matches the newly created vortex to 0
                    # the second condition ensures that one only considers currently living vortices and not vortices now annihilated
                    assignment[i] = j
                    break # breaks innerloop only
                end
            end
        end
    else # if charge  == "-"
        for i in eachindex(assignment)
            for j in 1:number_defectsN(dt)
                if assignment[i] ≠ 0 && dt.defectsN[j].annihilation_time == nothing && last_loc(dt.defectsN[j]) == old[proposal[i]]
                    assignment[i] = j
                    break
                end
            end
        end
    end
    return assignment
end

function find_closest_before_annihilation(dt, Lx, Ly, old_loc_defect)
    distance = Inf
    ID_antidefect = -1 # dummy
    for i in each(dt.defectsN)
        isnothing(dt.defectsN[i].annihilation_time) ? annihilation_time_defect = nothing : annihilation_time_defect = round(dt.defectsN[i].annihilation_time, digits=2)
        if annihilation_time_defect == round(dt.current_time, digits=2) # it has just annihilated
            tmp = dist(old_loc_defect, last_loc(dt.defectsN[i]), Lx, Ly)
            if tmp < distance
                distance = tmp
                ID_antidefect = i
            end
        end
    end
    if ID_antidefect == -1
        #= Comments :
        This error seems to occur when a defect disappears without any proper
        annihilation event (leaving a different number of defects + and - ) .
        This is why the algo cannot find the annihilator.
        Adding a small relaxation at T=0 before spotting defects seems to
        resolve the issue.
        =#
        println("Something has gone wrong... probably a defect has disappeared without any proper
        annihilation event.")
        println("Number of defects = ", number_active_defectsP(dft), " + ", number_active_defectsN(dft))
        return -1, (-1, -1)
    else
        return ID_antidefect, last_loc(dt.defectsN[ID_antidefect])
    end
end

# function delete_defect(dft::DefectTracker, id::Int, charge::String)
#     if charge == "+"
#         popat!(dft.defectsP, id)
#     else
#         popat!(dft.defectsN, id)
#     end
#     decrease_annihilation_ids!(dft, id, charge)
#     return dft
# end

# function decrease_annihilation_ids!(dft::DefectTracker, id::Int, charge::String)
#     # the id and the charge are those of the defect that was deleted
#     if charge == "+"
#         for def in dft.defectsN
#             if def.id_annihilator > id
#                 def.id_annihilator -= 1
#             end
#         end
#     else
#         for def in dft.defectsP
#             if def.id_annihilator > id
#                 def.id_annihilator -= 1
#             end
#         end
#     end
#     return dft
# end


function annihilate_defects(dt::DefectTracker, ids_annihilated_defects, Lx, Ly)
    for i in ids_annihilated_defects
        old_loc_vortex = last_loc(dt.defectsP[i])
        ID_antivortex, old_loc_antivortex = find_closest_before_annihilation(dt, Lx, Ly, old_loc_vortex)

        dt.defectsP[i].id_annihilator = ID_antivortex
        try
            dt.defectsN[ID_antivortex].id_annihilator = i
        catch
            dt.defectsN[ID_antivortex].id_annihilator = -1
        end

        estimate = mean_2_positions(old_loc_vortex, old_loc_antivortex, Lx, Ly)
        update_position!(dt.defectsP[i], estimate)
        try
            update_position!(dt.defectsN[ID_antivortex], estimate)
        catch
        end
    end
    return dt
end

function track!(dft::DefectTracker, system::System, times::AbstractVector;verbose=false, stop_no_defects=true)
    # Since the DefectTracker has already been initialized, discard t=0 if it exists
    if times[1] == 0
        times = times[2:end]
    end
    for token in each(times)
        evolve!(system, times[token])
        
        try
            update_DefectTracker!(dft, system)
        catch e
            println(e)
            println("Previous DefectTracker saved instead and immediate return.")
            return dft, system
        end
        if stop_no_defects && number_active_defects(dft) == 0
            println("Simulation stopped, there is no defects left.")
            break
        end
        if verbose 
            println("t = ", round(system.t, digits=1), " & n(t) = ", number_active_defectsP(dft), " + ", number_active_defectsN(dft))
        end
        if number_active_defectsN(dft) ≠ number_active_defectsP(dft)
            # p = plot_thetas(system,size=(512,512),
            # particles = true, 
            # nb_neighbours = false, 
            # defects = true)
            # title!()
            # display(p)
            println("There is a problem, the number of (+) and (-) defects is not the same.")
        end
    end
    return dft, system # times[end] is the last time of the simulation
end

function update_DefectTracker!(dft::DefectTracker, system::System) 
    dft.current_time = system.t
    
    Lx, Ly = system.Lx, system.Ly
    vortices_new, antivortices_new = spot_defects(system)

    # if BC == "periodic" @assert length(vortices_new) == length(antivortices_new) && length(vortices_old) == length(antivortices_old) end
    locP_old = [last_loc(dft.defectsP[i]) for i in eachindex(dft.defectsP)]
    locN_old = [last_loc(dft.defectsN[i]) for i in eachindex(dft.defectsN)]
    chargeP_old = [dft.defectsP[i].charge for i in eachindex(dft.defectsP)]
    chargeN_old = [dft.defectsN[i].charge for i in eachindex(dft.defectsN)]

    locP_new = [vortices_new[i][1:2] for i in eachindex(vortices_new)]
    locN_new = [antivortices_new[i][1:2] for i in eachindex(antivortices_new)]
    chargeP_new = [vortices_new[i][3] for i in eachindex(vortices_new)]
    chargeN_new = [antivortices_new[i][3] for i in eachindex(antivortices_new)]

    Np_new, Np_old = length(locP_new), length(locP_old)
    Nn_new, Nn_old = length(locN_new), length(locN_old)
    N_old = Np_old + Nn_old
    N_new = Np_new + Nn_new

    # Special simple cases to deal with upstream
    if N_new == N_old == 0 # do nothing, otherwise, "reducing over empty collection blablabla"

    elseif Nn_new == Nn_old == 0 && Np_new == Np_old > 0 # there are only (+) defects and no creation/annihilation
        assignment_vortices = pair_up_hungarian(dft, locP_new, locP_old, Lx, Ly, "+")
        for i in 1:Np_new
            update_position!(dft.defectsP[assignment_vortices[i]], locP_new[i])
        end

    elseif Np_new == Np_old == 0 && Nn_new == Nn_old > 0 # there are only (-) defects and no creation/annihilation
        assignment_antivortices = pair_up_hungarian(dft, locN_new, locN_old, Lx, Ly, "-")
        for i in 1:Nn_new
            update_position!(dft.defectsN[assignment_antivortices[i]], locN_new[i])
        end

    elseif N_new > 0 && N_old == 0
        for i in 1:Np_new
            add_defect!(dft, charge=chargeP_new[i], loc=locP_new[i])
        end
        for i in 1:Nn_new
            add_defect!(dft, charge=chargeN_new[i], loc=locN_new[i])
        end

    elseif N_new == 0 && N_old > 0 # (+)(-) >> plus rien
        id_just_annihilated_defectP, id_just_annihilated_defectM = ID_active_defects(dft) # seek for not yet annihilated defects

        for i in id_just_annihilated_defectP
            dft.defectsP[i].annihilation_time = dft.current_time
        end
        for i in id_just_annihilated_defectM
            dft.defectsN[i].annihilation_time = dft.current_time
        end
        dft = annihilate_defects(dft::DefectTracker, id_just_annihilated_defectP, Lx, Ly)

    elseif Np_new > 0 && Np_old > 0 && Nn_old > 0 && Nn_new == 0  # (+)(+)(-) >> (+) par exemple
        assignment_vortices = pair_up_hungarian(dft, locP_new, locP_old, Lx, Ly, "+")
        # Update living vortices. NB : the annihilated vortex is absent from the assignment vector : proceed without the condition "≠ 0"
        for i in eachindex(assignment_vortices)
            update_position!(dft.defectsP[assignment_vortices[i]], locP_new[i])
        end
        # Identify annihilated defects
        ID_annihilated_vortices = []
        ID_annihilated_antivortices = []
        for i in 1:number_defectsP(dft)
            if i ∉ assignment_vortices && dft.defectsP[i].annihilation_time == nothing # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
                push!(ID_annihilated_vortices, i)
            end
        end
        ID_annihilated_antivortices = ID_active_defects(dft)[2] # in this special case, there is no antivortices left in the "new" timestep

        for i in ID_annihilated_vortices
            dft.defectsP[i].annihilation_time = dft.current_time
        end
        for i in ID_annihilated_antivortices
            dft.defectsN[i].annihilation_time = dft.current_time
        end
        dft = annihilate_defects(dft, ID_annihilated_vortices, Lx, Ly)

    elseif Nn_new > 0 && Nn_old > 0 && Np_old > 0 && Np_new == 0  # (+)(-)(-) >> (-) par exemple
        assignment_antivortices = pair_up_hungarian(dft, locN_new, locN_old, Lx, Ly, "-")
        # Update living antivortices. NB : the annihilated antivortex is absent from the assignment vector : proceed without the condition "≠ 0"
        for i in eachindex(assignment_antivortices)
            update_position!(dft.defectsN[assignment_antivortices[i]])
        end
        # Identify annihilated defects
        ID_annihilated_vortices = []
        ID_annihilated_antivortices = []
        for i in 1:number_defectsN(dft)
            if i ∉ assignment_antivortices && dft.defectsN[i].annihilation_time == nothing # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
                push!(ID_annihilated_antivortices, i)
            end
        end
        ID_annihilated_vortices = ID_active_defects(dft)[1] # in this special case, there is no vortices left in the "new" timestep

        for i in ID_annihilated_vortices
            dft.defectsP[i].annihilation_time = dft.current_time
        end
        for i in ID_annihilated_antivortices
            dft.defectsN[i].annihilation_time = dft.current_time
        end

        dft = annihilate_defects(dft, ID_annihilated_vortices, Lx, Ly)
    else # end of special cases

        # GENERAL TREATMENT
        assignment_vortices = pair_up_hungarian(dft, locP_new, locP_old, Lx, Ly, "+")
        assignment_antivortices = pair_up_hungarian(dft, locN_new, locN_old, Lx, Ly, "-")

        # CASE 1 : no creation, no annihilation : simply update the data structure
        if N_new == N_old
            for i in 1:Np_new
                update_position!(dft.defectsP[assignment_vortices[i]], locP_new[i])
            end
            for i in 1:Nn_new
                update_position!(dft.defectsN[assignment_antivortices[i]], locN_new[i])
            end

            # CASE 2 : creation !
        elseif N_new > N_old
            # Take care of the newly created defects
            ind_created_vortex = findall(iszero, assignment_vortices) # newly created vortex -> the assignment vector contains a 0
            loc_created_vortex = vortices_new[ind_created_vortex]
            for j in each(loc_created_vortex)
                add_defect!(dft, charge=chargeP_new[j], loc=loc_created_vortex[j][1:2])
            end

            ind_created_antivortex = findall(iszero, assignment_antivortices)
            loc_created_antivortex = antivortices_new[ind_created_antivortex]
            for j in each(loc_created_antivortex)
                add_defect!(dft, charge=chargeN_new[j], loc=loc_created_antivortex[j][1:2])
            end

            # Update the ancient defects' positions
            for i in eachindex(assignment_vortices)
                if assignment_vortices[i] ≠ 0 # avoid newly created defects
                    update_position!(dft.defectsP[assignment_vortices[i]], locP_new[i])
                end
            end
            for i in eachindex(assignment_antivortices)
                if assignment_antivortices[i] ≠ 0 # avoid newly created defects
                    update_position!(dft.defectsN[assignment_antivortices[i]], locN_new[i])
                end
            end

            # CASE 3 : annihilation !
        elseif N_new < N_old
            # Update living vortices. NB : the annihilated vortex is absent from the assignment vector : proceed without the condition "≠ 0"
            for i in eachindex(assignment_vortices)
                update_position!(dft.defectsP[assignment_vortices[i]], locP_new[i])
            end
            for i in eachindex(assignment_antivortices)
                update_position!(dft.defectsN[assignment_antivortices[i]], locN_new[i])
            end

            # Identify annihilated defects
            ID_annihilated_vortices = []
            ID_annihilated_antivortices = []
            for i in 1:number_defectsP(dft)
                if i ∉ assignment_vortices && dft.defectsP[i].annihilation_time == nothing # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
                    dft.defectsP[i].annihilation_time = dft.current_time # from now on, defects that have just annihilated have annihilation_time == t
                    push!(ID_annihilated_vortices, i)
                end
            end
            for i in 1:number_defectsN(dft)
                if i ∉ assignment_antivortices && dft.defectsN[i].annihilation_time == nothing
                    dft.defectsN[i].annihilation_time = dft.current_time
                    push!(ID_annihilated_antivortices, i)
                end
            end
            if length(ID_annihilated_antivortices) >= length(ID_annihilated_vortices)
                dft = annihilate_defects(dft, ID_annihilated_vortices, Lx, Ly)
            else
                dft = annihilate_defects(dft, ID_annihilated_antivortices, Lx, Ly)
            end

        end # end of general treatment
    end # end of special cases & general treatment
    return dft
end

## Defects Analysis : MSD
function MSD(dfts::Union{Vector{DefectTracker},Vector{Union{Missing,DefectTracker}}}, Lx, Ly)
    indices = [] # indices of dft defined (is simulation not finished, dfts[i] == missing)
    for i in 1:length(dfts)
        if !ismissing(dfts[i])
            push!(indices, i)
        end
    end
    maxlength = maximum([maximum([length(d.pos) for d in vcat(dft.defectsP, dft.defectsN)]) for dft in dfts[indices]])
    MSD_P = NaN * zeros(length(indices), maxlength)
    MSD_N = NaN * zeros(length(indices), maxlength)
    MSD_all = NaN * zeros(length(indices), maxlength)
    for i in 1:length(indices)
        msd_all, msd_p, msd_n = MSD(dfts[indices[i]], Lx, Ly)
        MSD_P[i, 1:length(msd_p)] = msd_p
        MSD_N[i, 1:length(msd_n)] = msd_n
        MSD_all[i, 1:length(msd_all)] = msd_all
    end

    MSD_P_avg = nanmean(MSD_P, 1)[1, :]
    MSD_N_avg = nanmean(MSD_N, 1)[1, :]
    MSD_all_avg = nanmean(MSD_all, 1)[1, :]

    return MSD_all_avg, MSD_P_avg, MSD_N_avg
end

function MSD(dft::DefectTracker, Lx, Ly, maxlength=nothing)
    nP = number_defectsP(dft)
    nN = number_defectsN(dft)
    # tmin,tmax = t_bounds(dft) # (tmin,tmax) = timestamps of (first defect creation , last defect annihilation)

    # hasfield(typeof(model),:dt) ? dummy_dt = model.dt : dummy_dt = 1
    if isnothing(maxlength)
        maxlength = maximum([length(d.pos) for d in vcat(dft.defectsP, dft.defectsN)])
    end
    # Compute the SD
    SD_P = NaN * zeros(nP, maxlength)
    SD_N = NaN * zeros(nN, maxlength)
    for n in 1:nP
        defect = dft.defectsP[n]
        tmp = square_displacement(defect, Lx, Ly)
        SD_P[n, 1:length(tmp)] = tmp
    end
    for n in 1:nN
        defect = dft.defectsN[n]
        tmp = square_displacement(defect, Lx, Ly)
        SD_N[n, 1:length(tmp)] = tmp
    end

    # Now average to compute the MSD
    MSD_P = nanmean(SD_P, 1)[1, :]
    MSD_N = nanmean(SD_N, 1)[1, :]
    MSD_all = nanmean(hcat(MSD_P, MSD_N), 2)[:]

    return MSD_all, MSD_P, MSD_N
end

function square_displacement(d::Defect, Lx, Ly)
    loc_t0 = first_loc(d)
    return [dist(loc, loc_t0, Lx, Ly) for loc in d.pos] .^ 2
end

## Defects Analysis : Distance between defects
function interdefect_distance(defect1::Defect, defect2::Defect, Lx, Ly)
    # TODO take care of case with creation and/or annihilation time different.
    # So far, this care is left to the user...
    # @assert defect1.creation_time == defect2.creation_time
    # @assert defect1.annihilation_time == defect2.annihilation_time
    tmax = min(length(defect1.pos), length(defect2.pos))
    R = [dist(defect1.pos[t], defect2.pos[t], Lx, Ly) for t in 1:tmax]
    return R
end

function mean_distance_to_annihilator(dfts::Union{Vector{DefectTracker},Vector{Union{Missing,DefectTracker}}}, Lx, Ly)
    indices = [] # indices of dft defined (is simulation not finished, dfts[i] == missing)
    for i in 1:length(dfts)
        if !ismissing(dfts[i])
            push!(indices, i)
        end
    end
    Rs = [mean_distance_to_annihilator(dfts[indices[n]], Lx, Ly) for n in 1:length(indices)]
    Rs_matrix = vector_of_vector2matrix(Rs)
    return nanmean(Rs_matrix, 2)[:, 1]
end

function mean_distance_to_annihilator(dft::DefectTracker, Lx, Ly)
    nP = number_defectsP(dft)
    Rs = [distance_to_annihilator(dft, n, Lx, Ly) for n in 1:nP]
    Rs_matrix = vector_of_vector2matrix(Rs)
    return nanmean(Rs_matrix, 2)[:, 1], nanstd(Rs_matrix, 2)[:, 1]
end

function distance_to_annihilator(dft::DefectTracker, id1::Int, Lx, Ly; reversed=true)
    if isnothing(dft.defectsP[id1].id_annihilator) # not yet annihilated
        return [NaN]
    else
        R = interdefect_distance(dft, dft.defectsP[id1], dft.defectsN[dft.defectsP[id1].id_annihilator], Lx, Ly)
        if reversed
            reverse!(R)
        end
        return R
    end
end

## Auxiliary functions
function mean_2_positions(pos1, pos2, Lx, Ly, should_take_mod::Bool=true)
    a, b = pos1
    x, y = pos2

    dx = (x - a) #; dx = min(dx,L-dx)
    dy = (y - b) #; dy = min(dy,L-dy)

    if should_take_mod
        if abs(Lx - dx) < abs(dx)
            dx = -(Lx - dx)
        end
        if abs(Ly - dy) < abs(dy)
            dy = -(Ly - dy)
        end
        # return mod1.((a,b) .+ 0.5.*(dx,dy),L)
        return (mod1(a .+ 0.5 * dx, Lx), mod1(b .+ 0.5 * dy, Ly))
    else
        return (a, b) .+ 0.5 .* (dx, dy)
    end
end
# l = 100
# mean_2_positions((50,50),(60,60),l) == (55,55)
# mean_2_positions((10,10),(90,90),l) == (100,100)
# mean_2_positions((49,66),(51,61),l) == (50.0, 63.5)

function mean_N_positions(vec_pos, Lx, Ly, should_take_mod::Bool=true)
    averaged_pos = vec_pos[1]
    for i in 2:length(vec_pos)
        averaged_pos = mean_2_positions(averaged_pos, vec_pos[i], Lx, Ly, should_take_mod)
    end
    return averaged_pos
end

function smooth(x, c=1)
    len = length(x)
    result = zeros(len)
    result[1:c] .= NaN
    result[end-c+1:end] .= NaN
    for i in 1+c:len-c
        result[i] = mean(x[i-c:i+c])
    end
    return result
end

function corr_length(C::Vector{T}, rs=1:length(C); seuil=exp(-1))::T where {T<:AbstractFloat} # from a time series, returns the correlation length ()
    i_after = findfirst(x -> x < seuil, C)
    if i_after ≠ nothing && i_after > 1
        # Linear interpolation
        i_befor = i_after - 1
        r_after = rs[i_after]
        r_befor = rs[i_befor]
        c_after = C[i_after]
        c_befor = C[i_befor]
        ξ = (seuil * (r_after - r_befor) - (c_befor * r_after - r_befor * c_after)) / (c_after - c_befor)
    else
        ξ = NaN
    end
    return ξ
end

function energy(system)
    energy = 0.0
    list_neighbours = get_list_neighbours(system)
    thetas = get_thetas(system)
    for n in 1:system.N
        neighbours_n = list_neighbours[n]
        thetas_neighbours = thetas[neighbours_n]
        for theta in thetas_neighbours
            energy -= cos(thetas[n] - theta)
        end
    end
    return energy
end




########################## Miscellaneous ##########################

function remove_negative(input)
    array = Float64.(copy(input))
    for i in 1:length(array)
        if array[i] ≤ 0
            array[i] = NaN
        end
    end
    return array
end

logspace(x1, x2, n; digits=1) = unique!(round.([10.0^y for y in range(log10(x1), log10(x2), length=n)], digits=digits))


function prinz(z)
    ss = "$(round(Int,z)) seconds"
    mm = "$(round(z/60,digits=2)) minutes"
    hh = "$(round(z/3600,digits=2)) hours"
    dd = "$(round(z/86400,digits=2)) days"
    if z < 60
        println("Runtime : $ss = $mm")
    elseif z < 3600
        println("Runtime : $mm = $hh")
    else
        println("Runtime : $hh = $dd")
    end
    return z
end

function hrun(runtimes,nbins=20)
    mean_runtime = nanmean(runtimes)
    if mean_runtime < 60
        unit = "seconds"
    elseif mean_runtime < 3600
        unit = "minutes"
        runtimes = runtimes/60
    elseif mean_runtime < 86400
        unit = "hours"
        runtimes = runtimes/3600
    else 
        unit = "days"
        runtimes = runtimes/86400
    end
    h = histogram(runtimes, nbins=nbins, title=unit)
    return h
end

each = eachindex

nanmean(x) = mean(filter(!isnan, x))
nanmean(x, y) = mapslices(nanmean, x, dims=y)
nanstd(x) = std(filter(!isnan, x))
nanstd(x, y) = mapslices(nanstd, x, dims=y)

function vector_of_vector2matrix(x::Vector{Vector{T}}) where {T<:AbstractFloat}
    maxlength = maximum([length(x[i]) for i in eachindex(x)])
    matrice = NaN * zeros(T, maxlength, length(x))
    for i in eachindex(x)
        matrice[1:length(x[i]), i] = x[i]
    end
    return matrice
end
vecTuple2matrix(vec::Vector{Tuple{F,F}}) where {F<:AbstractFloat} = F.(vector_of_vector2matrix([[vec[n][1], vec[n][2]] for n in 1:length(vec)]))
# vecTuple2matrix(get_pos(system))


dist(a, b, Lx, Ly) = sqrt(dist2(a, b, Lx, Ly))
function dist2(a, b, Lx, Ly)  # euclidian distance with Periodic BCs
    dx = abs(a[1] - b[1])
    dx = min(dx, Lx - dx)
    dy = abs(a[2] - b[2])
    dy = min(dy, Ly - dy)
    return dx^2 + dy^2
end

function distance_matrix(new, old, Lx, Ly)
    m_new, m_old = length(new), length(old)
    distance_matrix = zeros(m_new, m_old)
    for j in 1:m_old
        for i in 1:m_new
            distance_matrix[i, j] = dist(new[i], old[j], Lx, Ly)
        end
    end
    return distance_matrix
end