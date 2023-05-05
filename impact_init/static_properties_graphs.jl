cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian, LambertW, Graphs, Compose, GraphPlot
include("../methods.jl")
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);

#=
For the different graphs : 
	- Square and Triangular regular lattices
	- Random
	- RSA
	- Poissin Disc Sampling
Investigate their static properties :
	- Degree distribution
	- Clustering coefficient
	- Percolation/existence GCC/ number of components / number of neighbours versus R0
	- Structure factor / exponent of the variance of the number of neighbours
Visualization :
	- Components (colorized)
	- Number of neighbours
=#

# get_list_neighbours(system) to obtain the list of the neighbours of each particle
# length.(get_list_neighbours(system)) to obtain the number of neighbours of each particle

function system2graph(system)
	graph = SimpleGraph(system.N)
	neighbours = get_list_neighbours(system)
	for i in 1:system.N
		for j in neighbours[i]
			add_edge!(graph, i, j)
		end
	end
	return graph
end

function ids_connected_component(g)
	components = connected_components(g)
	d_connected_component = 
	IDs_connected_component = zeros(Int, nv(g))
	for i in each(components)
		IDs_connected_component[components[i]] .= i
	end
	return IDs_connected_component
end

##
include("../parameters.jl")
lattice_type = "square"
lattice_type = "random"
lattice_type = "pds"
params_init[:init_pos] = lattice_type
param[:params_init] = params_init




## --------------- Systems Generation  --------------- ##
include("../parameters.jl")
lattice_types = ["square", "random", "rsa", "pds"]
labels = ["Square", "Random", "RSA", "PDS"]

systems = []
graphs = []
for i in each(lattice_types)
	params_init[:init_pos] = lattice_types[i]
	param[:params_init] = params_init
	system = System(param)
	push!(systems, system)
	push!(graphs, system2graph(system))
end

## --------------- Distributions --------------- ##
# p1 = plot(yaxis=:log, legend=:topright,xlabel="Number of neighbours", ylabel="Frequency")
# vline!([8],c=1)
# for i in 2:length(systems)
# 	data = length.(get_list_neighbours(systems[i]))
# 	h = fit(Histogram, data, 1:20)
# 	ww = h.weights/length(data)
# 	for i in each(ww)
# 		if ww[i] == 0
# 			ww[i] = NaN
# 		end
# 	end
# 	println(ww)
# 	plot!(h.edges[1][2:end], ww,marker=2,c=i,label=labels[i])
# end
# p1
# savefig("impact_init/figures/static_properties_graphs/degree_distribution.png")
## --------------- Visualisation nnn and components --------------- ##
msss = 3
p2_subs = Array{Any}(undef, 3, length(lattice_types))
# gplot
for i in each(graphs)
	p2_subs[1,i] = gplot(graphs[i])
end
# nnn
for i in each(systems)
	pos = get_pos(systems[i])
	tmp = plot(aspect_ratio=1)
	scatter!([el[1] for el in pos], [el[2] for el in pos], 
	marker_z=length.(get_list_neighbours(systems[i])), 
	markerstrokewidth=0, markersize=msss,size=(500,400),
	c=cgrad([:black, :blue, :green, :orange, :red]),title=labels[i])
	p2_subs[2,i] = tmp
end
# components
for i in each(systems)
	pos = get_pos(systems[i])
	tmp = plot(aspect_ratio=1)
	scatter!([el[1] for el in pos], [el[2] for el in pos], 
	marker_z=ids_connected_component(graphs[i]), 
	markerstrokewidth=0, markersize=msss,size=(500,400),
	c=cgrad([:black, :blue, :green, :orange, :red]))
	p2_subs[3,i] = tmp
end
p2 = plot(p2_subs[2,1], p2_subs[2,2], p2_subs[2,3], p2_subs[2,4])
# p2_subs[3,1], p2_subs[3,2], p2_subs[3,3], p2_subs[3,4], 
# layout=(2,length(lattice_types)), size=(2000,800))


# p2 = plot(p2_subs[2,:]..., layout=(1,length(lattice_types)), size=(2000,400))
# title!("Number of neighbours")
# savefig("impact_init/figures/static_properties_graphs/vizu_nnn_cc.png")
p2 = plot(p2_subs[2,1], p2_subs[2,2], p2_subs[2,3], p2_subs[2,4],layout=(2,2), size=(1000,800))
savefig(p2,"impact_init/figures/static_properties_graphs/vizu_nnn.png")

## --------------- Influence R0 on Connected Components --------------- ##
include("../parameters.jl")
R0ss = collect(1.0:0.01:1.5)
R = 100
nccs = zeros(length(R0ss),length(systems),R)
fraction_gcc = zeros(length(R0ss),length(systems),R)
for i in each(systems)
	println("i = $i/$(length(systems))")
	for j in each(R0ss)
		Threads.@threads for r in 1:R
			param[:R0] = R0ss[j]
			params_init[:init_pos] = lattice_types[i]
			param[:params_init] = params_init		
			system = System(param)
			g = system2graph(system)
			seuil_cluster = 3
			nccs[j,i,r] = sum(length.(connected_components(g)) .≥ seuil_cluster)
			fraction_gcc[j,i,r] = length.(connected_components(g))[1]/nv(g)		
		end
	end
end
nccs_avg = mean(nccs,dims=3)[:,:,1]
fraction_gcc_avg = mean(fraction_gcc,dims=3)[:,:,1]
##
p3a=plot(xlabel=L"R_0",ylabel="Number of Components",yaxis=:log)
for i in each(systems)
	plot!(R0ss,nccs_avg[:,i],label=labels[i],m=true,ms=2)
end
p3a

p3b=plot(xlabel=L"R_0",ylabel="Fraction GCC",yaxis=:log,legend=:bottomright)
for i in each(systems)
	plot!(R0ss,fraction_gcc_avg[:,i],label=labels[i],m=true,ms=2)
end
p3b

p3 = plot(p3a,p3b,layout=(1,2),size=(800,400))
# savefig("impact_init/figures/static_properties_graphs/ncc_fcc_R0s.png")

## --------------- Influence R0 on mean number of neighbours --------------- ##
include("../parameters.jl")
R0ss = collect(1.0:0.1:2)
R = 10
nnns = zeros(length(R0ss),length(systems),R)
for i in each(systems)
	println("i = $i/$(length(systems))")
	for j in each(R0ss)
		println("j = $j/$(length(R0ss))")
		Threads.@threads for r in 1:R
			param[:R0] = R0ss[j]
			params_init[:init_pos] = lattice_types[i]
			param[:params_init] = params_init		
			system = System(param)
			g = system2graph(system)
			nnns[j,i,r] = mean(length.(get_list_neighbours(system)))
		end
	end
end
nnns_avg = mean(nnns,dims=3)[:,:,1]
##
p4=plot(xlabel=L"R_0^2",ylabel="Number of Neighbours",uaxis=:log)
for i in each(systems)
	plot!(R0ss.^2,nnns_avg[:,i],label=labels[i],m=true,ms=2)
end
yticks!(0:2:12)
# plot!(R0ss.^2,pi*R0ss.^2 .-pi .+2.5,c=:black,)

p4
savefig("impact_init/figures/static_properties_graphs/mean_nnn_R0s.png")
## --------------- Structure Factor --------------- ##
function structure_factor(system)
	pos = get_pos(system)
	N = system.N
	Lx,Ly = system.Lx, system.Ly
	qs = collect(0.0:0.1:10)
	S = zeros(ComplexF32,length(qs))
	for i in 1:N
		pos_i = pos[i]
		for j in i:N
			distance  = dist(pos_i,pos[j],Lx,Ly)
			S += 2exp.(-im*distance*qs)
		end
	end
	S = S/N
	return qs,S
end

include("../parameters.jl")
system = System(param)
system.N
qs,S = structure_factor(system)
plot(qs,imag(S))
@btime structure_factor(system)


## ChatGPT
using FFTW

function compute_structure_factor(lattice_points, N, grid_size)
    # Input: lattice_points - an array of N points in the lattice with coordinates (x_i, y_i)
    #        N - total number of points in the lattice
    #        grid_size - size of the 2D grid for FFT, should be a power of 2
    # Output: structure_factor - the computed structure factor of the lattice

    # Scale the lattice points to fit within the grid size
    scaled_points = ceil.(Int, lattice_points .* grid_size)
    
    # Create a 2D grid of size grid_size x grid_size and initialize with zeros
    grid = zeros(grid_size, grid_size)
    
    # Loop through all the scaled lattice points and increment the corresponding grid cell by 1
    for i in 1:N
        x = scaled_points[i, 1]
        y = scaled_points[i, 2]
        grid[x, y] += 1
    end
    
    # Perform 2D FFT on the grid
    fft_grid = fftshift(fft(grid))
    
    # Compute the squared magnitude of the FFT result
    power_spectrum = abs2.(fft_grid)
    
    # Normalize the power spectrum
    structure_factor = power_spectrum / (grid_size^2)
    
    return structure_factor
end

# Example usage:
# Generate random lattice points in [0,1]
N = 1000
lattice_points = rand(N, 2)

# Set the grid size for FFT
grid_size = 256

# Compute structure factor
structure_factorr = compute_structure_factor(lattice_points, N, grid_size)
