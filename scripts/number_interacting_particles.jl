cd("D:/Documents/Research/projects/kuramoto_ballistic")
    using JLD2,StatsBase,Distributions,LinearAlgebra,Parameters,Random,BenchmarkTools,Hungarian
    include("../methods.jl")
    global const R0 = 1
    using Plots,ColorSchemes,LaTeXStrings
    gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    plot()
&

#= The file aims at computing the number of different particles
a test particles interacts with over time. Underlying idea: understanding
the relation between moving particles and MF approximation.
=#
N = Int(4E2)
    rho = 4
    L = round(Int,sqrt(N/rho))

    pos,~,psis,~ = initialisation(N,L,0)
    v0 = 1

    Niter = Int(1E3)
    ID_interactions = [Int[] for i in 1:N]
    number_interactions = zeros(N,Niter)

    z = @elapsed for t in 1:Niter
        # Update position
        for n in eachindex(psis)
            pos[:,n] += v0*0.5E-1*[cos(psis[n]),sin(psis[n])]
        end
        pos = mod.(pos,L)

        nb_cells_1D = Int(div(L,R0)) + 1
        head = -ones(Int,nb_cells_1D,nb_cells_1D) # head[i,j] contains the index of the first particle in cell (i,j). -1 if empty
        list = -ones(Int,N) # list[n] contains the index of the particle to which particle n points. -1 if it points to no one
        for n in 1:N
            cellx,celly = Int(div(pos[1,n],R0)) + 1 , Int(div(pos[2,n],R0)) + 1 # cell to which particle n belongs
            list[n] = head[cellx,celly]
            head[cellx,celly] = n
        end

        for n in 1:N
            poscur = pos[:,n]

            cellx,celly = Int(div(pos[1,n],R0)) + 1 , Int(div(pos[2,n],R0)) + 1 # cell to which particle n belongs
            should_take_mod = (cellx == 1) || (cellx == nb_cells_1D) || (celly == 1) || (celly == nb_cells_1D)
            if should_take_mod
                neighbouring_cells = Vector{Int}[ [cellx,celly] , [cellx,mod1(celly+1,nb_cells_1D)] , [mod1(cellx+1,nb_cells_1D),celly] , [cellx,mod1(celly-1,nb_cells_1D)] , [mod1(cellx-1,nb_cells_1D),celly] , [mod1(cellx+1,nb_cells_1D),mod1(celly+1,nb_cells_1D)] ,  [mod1(cellx-1,nb_cells_1D),mod1(celly-1,nb_cells_1D)] , [mod1(cellx-1,nb_cells_1D),mod1(celly+1,nb_cells_1D)] , [mod1(cellx+1,nb_cells_1D),mod1(celly-1,nb_cells_1D)]]
            else
                neighbouring_cells = Vector{Int}[ [cellx,celly] , [cellx,celly+1] , [cellx+1,celly] , [cellx,celly-1] , [cellx-1,celly] , [cellx+1,celly+1] ,  [cellx-1,celly-1] , [cellx-1,celly+1] , [cellx+1,celly-1]]
            end

            ind_neighbours = Int[]
            for (i,j) in neighbouring_cells
                next = head[i,j]
                if next ≠ -1
                    if 0 < dist(poscur,pos[:,next],L) < R0 push!(ind_neighbours,next) end
                    while list[next] ≠ -1
                        if 0 < dist(poscur,pos[:,list[next]],L) < R0 push!(ind_neighbours,list[next]) end
                        next = list[next]
                    end
                end
            end
            push!(ID_interactions[n],ind_neighbours...)
            unique!(ID_interactions[n])
            number_interactions[n,t] = length(ID_interactions[n])
        end
    end
    prinz(z)
    nintavg = mean(number_interactions,dims=1)[1,:]
    plot!(nintavg[1:end]/N,m=:circle)
    # plot!(x->2/pi*atan(x/28))

plot(nintavg[1:end]/N,m=:circle)
plot!(x->2/pi*atan(x/2100*4),c=:black)



## Solving the dynamics of the fully connected Kuramoto model
N = Int(1E4)
function complexOP(thetas::Vector{Float64})
    z = mean(exp.(im*thetas))
    return abs(z),angle(z)
end

dt = 1E-1 ; tmax = 500
times = logspace(dt,tmax,40,digits=2)
rs = zeros(length(times))
thetas = 2pi*rand(N)
t = 0 ; token = 1
z = @elapsed while t < tmax
    t += dt
    r,psi = complexOP(thetas)
    thetas += dt*0.05r*sin.(psi .- thetas)
    if t ≥ times[token]
        rs[token] = r
        token = min(token+1,length(times))
    end
end
prinz(z)
plot(times,rs,uaxis=:log,m=true)
    plot!(x->1-exp(-x/10000),uaxis=:log)
    plot!(x->2/pi*atan(x/1000))
histogram(mod.(thetas,2pi),bins=50,normalize=true)

&
