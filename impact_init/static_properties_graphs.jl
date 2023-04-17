cd("D:/Documents/Research/projects/kuramoto_ballistic")
using JLD2, StatsBase, Distributions, LinearAlgebra, Parameters, Random, BenchmarkTools, Hungarian, LambertW, Graphs, Compose, GraphPlot
include("../methods.jl")
using Plots, ColorSchemes, LaTeXStrings
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
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
lattice_type = "square_lattice"
lattice_type = "random"
lattice_type = "pds"
params_init[:init_pos] = lattice_type
param[:params_init] = params_init
system = System(param)
g = system2graph(system)
# gplot(g)
connected_components(g)
seuil_cluster = 
sum(length.(connected_components(g)) .≥ seuil_cluster)
##
histogram(log.(length.(connected_components(g))))

plot(connected_components(g))

collect(keys(degree_histogram(g)))
meann = round(mean(length.(get_list_neighbours(system))),digits=2)
varr = round(var(length.(get_list_neighbours(system))),digits=2)
histogram(length.(get_list_neighbours(system)),bins=5, title="Mean : $meann ± $varr")

collect(values(degree_histogram(g)))

components(1:N)

## --------------- Systems Generation  --------------- ##
include("../parameters.jl")
lattice_types = ["square_lattice", "random", "rsa", "pds"]
labels = ["Square lattice", "Random", "RSA", "PDS"]

systems = []
graphs = []
for i in each(lattice_types)
	params_init[:init_pos] = lattice_types[i]
	param[:params_init] = params_init
	system = System(param)
	push!(systems, system)
	push!(graphs, system2graph(system))
end

# --------------- Distributions --------------- ##
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
mss = 1
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
	markerstrokewidth=0, markersize=mss,size=(500,400),
	c=cgrad([:black, :blue, :green, :orange, :red]))
	p2_subs[2,i] = tmp
end
# components
for i in each(systems)
	pos = get_pos(systems[i])
	tmp = plot(aspect_ratio=1)
	scatter!([el[1] for el in pos], [el[2] for el in pos], 
	marker_z=ids_connected_component(graphs[i]), 
	markerstrokewidth=0, markersize=mss,size=(500,400),
	c=cgrad([:black, :blue, :green, :orange, :red]))
	p2_subs[3,i] = tmp
end
# p2 = plot(p2_subs[2,1], p2_subs[2,2], p2_subs[2,3], p2_subs[2,4],
# p2_subs[3,1], p2_subs[3,2], p2_subs[3,3], p2_subs[3,4], 
# layout=(2,length(lattice_types)), size=(2000,800))

# p2_subs[2,4]

p2 = plot(p2_subs[2,:]..., layout=(1,length(lattice_types)), size=(2000,400))
title!("Number of neighbours")
# savefig("impact_init/figures/static_properties_graphs/vizu_nnn_cc.png")


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

## --------------- Structure Factor --------------- ##
