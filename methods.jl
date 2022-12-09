"Copyright (c) 2022 Y.Rouzaire All Rights Reserved."
logspace(x1, x2, n;digits=1) = unique!(round.([10.0 ^y for y in range(log10(x1), log10(x2), length=n)],digits=digits))

function prinz(z)
    println("Runtime : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes = $(round(z/3600,digits=2)) hours.")
    return z
end

each = eachindex

nanmean(x) = mean(filter(!isnan,x))
nanmean(x,y) = mapslices(nanmean,x,dims=y)
nanstd(x) = std(filter(!isnan,x))
nanstd(x,y) = mapslices(nanstd,x,dims=y)

function vector_of_vector2matrix(x::Vector{Vector{T}}) where T<:AbstractFloat
    maxlength = maximum([length(x[i]) for i in eachindex(x)])
    matrice = NaN*zeros(T,maxlength,length(x))
    for i in eachindex(x)
        matrice[1:length(x[i]),i] = x[i]
    end
    return matrice
end


function dist(a::Vector{T},b::Vector{T},Lx,Ly) where T<:AbstractFloat  # euclidian distance with Periodic BCs
    dx = abs(a[1] - b[1]) ; dx = min(dx,Lx-dx)
    dy = abs(a[2] - b[2]) ; dy = min(dy,Ly-dy)
    return sqrt(dx^2 + dy^2)
end

function dist(a::Tuple{T,T},b::Tuple{T,T},Lx,Ly) where T<:Number  # euclidian distance with Periodic BCs
    dx = abs(a[1] - b[1]) ; dx = min(dx,Lx-dx)
    dy = abs(a[2] - b[2]) ; dy = min(dy,Ly-dy)
    return sqrt(dx^2 + dy^2)
end

function dist(a,b,Lx,Ly)  # euclidian distance with Periodic BCs
    dx = abs(a[1] - b[1]) ; dx = min(dx,Lx-dx)
    dy = abs(a[2] - b[2]) ; dy = min(dy,Ly-dy)
    return sqrt(dx^2 + dy^2)
end


function distance_matrix(new,old,Lx,Ly)
    m_new,m_old = length(new),length(old)
    distance_matrix = zeros(m_new,m_old)
    for j in 1:m_old
        for i in 1:m_new
            distance_matrix[i,j] = dist(new[i],old[j],Lx,Ly)
        end
    end
    return distance_matrix
end

function determine_dt(T,σ,v0,N,rho)
    estimation_max = sqrt(2log(N)) * σ # from math.stackexchange.com/questions/89030, simple answer of Sivaraman
    arbitrary_coeff =  π/20
    return minimum(
            [1/5/v0 ,
            arbitrary_coeff/estimation_max ,
            arbitrary_coeff/π/rho ,
            arbitrary_coeff^2*π/4/T])
end

## Initialisation and evolution
function initialisation(N,Lx,Ly,σ,params=["disordered"];float_type=Float32)
    pos = rand(2,N)
    pos[1,:] *= Lx
    pos[2,:] *= Ly

    psis = 2π*rand(N)
    omegas = σ*randn(N)

    if params[1] in ["disordered","hightemp"]
        thetas = 2π*rand(N)
    elseif params[1] in ["ordered","lowtemp"]
        thetas = zeros(N)
    elseif params[1] in ["isolated","single"] # params[2] = q the charge
        thetas = zeros(N)
        for n in 1:N
            x,y = pos[:,n]
            thetas[n] = params[2] * atan(Ly/2 - y,Lx/2 - x)
        end
    elseif params[1] == "pair" # params[2] is R0 the separation distance
        # should I implement the different pairs ? source and sink will not act the same way for instance
        thetas = zeros(N)
        for n in 1:N
            x,y = pos[:,n]
            thetas[n] = atan(Ly/2 - y,Lx/2 - x + params[2]/2) - atan(Ly/2 - y,Lx/2 - x - params[2]/2) # +1 defect by convention
        end
    elseif params[1] == "1Dwave"
        # then params[2] is "vertical"/"horizontal"
        # then params[3] is the number of waves
        thetas = zeros(N)
        if params[2] == "horizontal"
            for n in 1:N thetas[n] = mod(params[3]*2pi*pos[1,n]/Lx,2pi) end
        elseif params[2] == "vertical"
            for n in 1:N thetas[n] = mod(params[3]*2pi*pos[2,n]/Ly,2pi) end
        end
    elseif params[1] == "2Dwave"
        # then params[2] is the number of horizontal waves
        # then params[3] is the number of vertical waves
        thetas = zeros(N)
        for n in 1:N thetas[n] = mod(π .+ π/2*(sin(params[2]*pos[1,n]*2pi/Lx)+sin(params[3]*pos[2,n]*2pi/Lx)),2pi) end
    end
    return float_type.(pos),float_type.(thetas),float_type.(psis),float_type.(omegas)
end

function update(pos::Matrix{FT},thetas::Vector{FT},psis::Vector{FT},omegas::Vector{FT},T::Number,v0::Number,N::Int,Lx::Int,Ly::Int,dt::Number) where FT<:AbstractFloat
    pos_new    = zeros(FT,2,N)
    thetas_new = zeros(FT,N)

    # Update position
    for n in eachindex(psis)
        pos_new[:,n] = pos[:,n] + v0*dt*[cos(psis[n]),sin(psis[n])]
    end
    pos_new[1,:] = mod.(pos_new[1,:],Lx)
    pos_new[2,:] = mod.(pos_new[2,:],Ly)

    ## List construction
    nb_cells_x = Int(div(Lx,R0)) + 1
    nb_cells_y = Int(div(Ly,R0)) + 1
    head = -ones(Int,nb_cells_x,nb_cells_y) # head[i,j] contains the index of the first particle in cell (i,j). -1 if empty
    list = -ones(Int,N) # list[n] contains the index of the particle to which particle n points. -1 if it points to no one
    for n in 1:N
        cellx,celly = Int(div(pos[1,n],R0)) + 1 , Int(div(pos[2,n],R0)) + 1 # cell to which particle n belongs
        list[n] = head[cellx,celly]
        head[cellx,celly] = n
    end

    for n in 1:N
        poscur = pos[:,n]

        cellx,celly = Int(div(pos[1,n],R0)) + 1 , Int(div(pos[2,n],R0)) + 1 # cell to which particle n belongs
        should_take_mod = (cellx == 1) || (cellx == nb_cells_x) || (celly == 1) || (celly == nb_cells_y)
        if should_take_mod
            neighbouring_cells = Vector{Int}[ [cellx,celly] , [cellx,mod1(celly+1,nb_cells_y)] , [mod1(cellx+1,nb_cells_x),celly] , [cellx,mod1(celly-1,nb_cells_y)] , [mod1(cellx-1,nb_cells_x),celly] , [mod1(cellx+1,nb_cells_x),mod1(celly+1,nb_cells_y)] ,  [mod1(cellx-1,nb_cells_x),mod1(celly-1,nb_cells_y)] , [mod1(cellx-1,nb_cells_x),mod1(celly+1,nb_cells_y)] , [mod1(cellx+1,nb_cells_x),mod1(celly-1,nb_cells_y)]]
        else
            neighbouring_cells = Vector{Int}[ [cellx,celly] , [cellx,celly+1] , [cellx+1,celly] , [cellx,celly-1] , [cellx-1,celly] , [cellx+1,celly+1] ,  [cellx-1,celly-1] , [cellx-1,celly+1] , [cellx+1,celly-1]]
        end

        ind_neighbours = Int[]
        for (i,j) in neighbouring_cells
            next = head[i,j]
            if next ≠ -1
                if 0 < dist(poscur,pos[:,next],Lx,Ly) < R0 push!(ind_neighbours,next) end
                while list[next] ≠ -1
                    if 0 < dist(poscur,pos[:,list[next]],Lx,Ly) < R0 push!(ind_neighbours,list[next]) end
                    next = list[next]
                end
            end
        end

        thetas_neighbours = thetas[ind_neighbours]

        # Sync Theta
        theta = thetas[n]
        if length(thetas_neighbours) > 0
            thetas_new[n] = theta + dt * omegas[n] + dt * sum(sin,thetas_neighbours .- theta) + sqrt(2T*dt)*randn()
        else
            thetas_new[n] = theta + dt * omegas[n] + sqrt(2T*dt)*randn()
        end
    end
    return pos_new,thetas_new
end

## Measurements
function polarOP(thetas::Vector{T}) where T<:AbstractFloat
    tmp = Complex(0)
    for theta in thetas
        tmp += exp(1im*theta)
    end
    return norm(tmp)/length(thetas),angle(tmp)
end

function nematicOP(thetas::Vector{T}) where T<:AbstractFloat
    tmp = Complex(0)
    for theta in thetas
        tmp += exp(2im*theta)
    end
    return norm(tmp)/length(thetas),angle(tmp)
end

function corr(pos::Matrix{T},thetas::Vector{T},N,Lx,Ly,dr;algo="fast")::Vector{T} where T<:AbstractFloat
    if     algo == "slow" return corr_slow(pos,thetas,N,Lx,Ly,dr)
    elseif algo == "fast" return corr_fast(pos,thetas,N,Lx,Ly,dr)
    end
end

function corr_slow(pos::Matrix{T},thetas::Vector{T},N,Lx,Ly,dr)::Vector{T} where T<:AbstractFloat
    # Construct matrix of distances
    Lmin = min(Lx,Ly)
    C = [T[] for i in 1:round(Int,Lmin/2/dr)]
    # distances = zeros(N,N)
    # for j in 1:N , i in j+1:N
    #     distances[i,j] = dist(pos[:,i],pos[:,j],L)
    # end
    for j in 1:N , i in j+1:N
        d = dist(pos[:,i],pos[:,j],Lx,Ly)
        if d ≤ round(Int,Lmin/2) push!(C[min(ceil(Int,d/dr),length(C))],cos(thetas[i] - thetas[j])) end
    end

    Cavg = [mean(C[i]) for i in eachindex(C)]
    return vcat(1,Cavg) # 1 to represent C(0,t)
end

function corr_fast(pos::Matrix{T},thetas::Vector{T},N,Lx,Ly,dr)::Vector{T} where T<:AbstractFloat
    M = 50
    Lmin = min(Lx,Ly)
    C = [T[] for i in 1:round(Int,Lmin/2/dr)]
    ms = zeros(Int,round(Int,Lmin/2/dr))

    for j in 1:N , i in j+1:N
        d = dist(pos[:,i],pos[:,j],Lx,Ly)
        if d ≤ round(Int,Lmin/2)
            ind = min(ceil(Int,d/dr),length(C))
            push!(C[ind],cos(thetas[i] - thetas[j]))
            ms[ind] += 1
        end
        if sum(ms .> M) == length(ms) break end
    end

    Cavg = [mean(C[i]) for i in eachindex(C)]
    return vcat(1,Cavg) # 1 to represent C(0,t)
end

## Methods for defects
function mod1_2D(xx::Tuple{T,T},Lx::Int,Ly::Int) where T<:Number
    return (mod1(xx[1],Lx) , mod1(xx[2],Ly))
end

function cg(pos::Matrix{T},thetas::Vector{T},N,Lx,Ly) where T<:AbstractFloat
    mesh_size = R0
    cutoff = 5R0 # for contributions

    ## Cell List construction
    nb_cells_x = Int(div(Lx,mesh_size)) + 1
    nb_cells_y = Int(div(Ly,mesh_size)) + 1
    head = -ones(Int,nb_cells_x,nb_cells_y) # head[i,j] contains the index of the first particle in cell (i,j). -1 if empty
    list = -ones(Int,N) # list[n] contains the index of the particle to which particle n points. -1 if it points to no one
    for n in 1:N
        cellx,celly = Int(div(pos[1,n],mesh_size)) + 1 , Int(div(pos[2,n],mesh_size)) + 1 # cell to which particle n belongs
        list[n] = head[cellx,celly]
        head[cellx,celly] = n
    end

    LLx = round(Int,Lx/mesh_size)
    LLy = round(Int,Ly/mesh_size)
    fine_grid = NaN*zeros(LLx,LLy)
    fine_grid_density = NaN*zeros(LLx,LLy)
    for i in 1:LLx, j in 1:LLy # scan fine_grid cells
        center_finegrid_cell = T.([(i-0.5)*mesh_size,(j-0.5)*mesh_size])

        # find cell from coarse mesh correponding to cell i,j belonging to fine_grid
        cellx,celly = Int(div(center_finegrid_cell[1],mesh_size)) + 1 , Int(div(center_finegrid_cell[2],mesh_size)) + 1
        a = round(Int,cutoff/mesh_size)
        neighbouring_cells = vec([(ii,jj) for ii in -a:a, jj in -a:a])
        neighbouring_cells = [mod1_2D(el .+ (cellx,celly),nb_cells_x,nb_cells_y) for el in neighbouring_cells]

        # get indices of particles within those cells (cells belonging to the coarse mesh)
        ind_neighbours = Int[]
        for (i,j) in neighbouring_cells
            next = head[i,j]
            if next ≠ -1
                if dist(center_finegrid_cell,pos[:,next],Lx,Ly) < cutoff push!(ind_neighbours,next) end
                while list[next] ≠ -1
                    if dist(center_finegrid_cell,pos[:,list[next]],Lx,Ly) < cutoff push!(ind_neighbours,list[next]) end
                    next = list[next]
                end
            end
        end

        # compute contributions from those "neighbouring" particles
        tmp = Complex[]
        tmp_density = Float64[]
        for m in ind_neighbours
            r = dist(center_finegrid_cell,pos[:,m],Lx,Ly)
            if r < cutoff
                push!(tmp,exp(im*thetas[m]-r/R0))
                # push!(tmp_density,exp(-r/R0))
            end
        end
        fine_grid[i,j] = angle(mean(tmp))
        # fine_grid_density[i,j] = sum(tmp_density)
    end
    return fine_grid
    # return fine_grid,fine_grid_density
end

function arclength(theta1::T,theta2::T)::T where T<:AbstractFloat
    #= This function returns the signed arclength on the unit trigonometric circle .
    Clockwise        >> sign -
    Counterclockwise >> sign +
    WARNING
    Note that the inputs thetas need to lie within [0,π] or [0,2π], depending on the symmetry of the model =#
    dtheta = theta2 - theta1
    dtheta_abs = abs(theta2 - theta1)

    shortest_unsigned_arclength = min(2π-dtheta_abs,dtheta_abs)
    if dtheta_abs ≤ π
        signe = sign(dtheta)
    else
        signe = -sign(dtheta)
    end
    return signe*shortest_unsigned_arclength
end

function get_neighbours(thetas::Matrix{<:T},i::Int,j::Int,bulk::Bool=false)::Vector{T} where T<:AbstractFloat
    Lx,Ly = size(thetas)
    # convention depuis la droite et sens trigo
    if bulk
        jm  = j-1
        jp  = j+1
        imm = i-1
        ip  = i+1
    else
        jm  = mod1(j-1,Ly)
        jp  = mod1(j+1,Ly)
        imm = mod1(i-1,Lx)
        ip  = mod1(i+1,Lx)
    end

    @inbounds angles =
       [thetas[i,jp],
        thetas[imm,j],
        thetas[i,jm],
        thetas[ip,j]]

    return angles
end

is_on_border(i::Int,j::Int,Lx::Int,Ly::Int) = (i == 1) || (j == 1) || (i == Lx) || (j == Ly)
is_in_bulk(i::Int,j::Int,Lx::Int,Ly::Int) = !is_on_border(i,j,Lx,Ly)

function get_vorticity(thetas_mat::Matrix{T},i::Int,j::Int,Lx::Int,Ly::Int)::T where T<:AbstractFloat
    # Note : thetas_mat = mod.(thetas,2π)
    angles_corners = get_neighbours(thetas_mat,i,j,is_in_bulk(i,j,Lx,Ly))
    perimeter_covered = 0.0
    for i in 1:length(angles_corners)-1
        perimeter_covered += arclength(angles_corners[i],angles_corners[i+1])
    end
    if !isempty(angles_corners)
        perimeter_covered += arclength(angles_corners[end],angles_corners[1])
        charge = round(perimeter_covered/2π,digits=1)
    else charge = NaN
    end
    return charge
end

number_defects(pos,thetas,N,Lx,Ly) = sum(length,spot_defects(pos,thetas,N,Lx,Ly))
number_defectsP(pos,thetas,N,Lx,Ly) = length(spot_defects(pos,thetas,N,Lx,Ly)[1])
number_defectsN(pos,thetas,N,Lx,Ly) = length(spot_defects(pos,thetas,N,Lx,Ly)[2])
function spot_defects(pos::Matrix{T},thetas::Vector{T},N,Lx,Ly) where T<:AbstractFloat
    vortices_plus  = Tuple{Int,Int,T}[]
    vortices_minus = Tuple{Int,Int,T}[]

    thetasmod = mod.(cg(pos,thetas,N,Lx,Ly),2π)
    relax!(thetasmod)

    for i in 1:Lx
        for j in 1:Ly
            q = get_vorticity(thetasmod,i,j,Lx,Ly)
            if     q > + 0.1 push!(vortices_plus,(i,j,q)) # we want to keep ±½ defects, and not rounding errors
            elseif q < - 0.1 push!(vortices_minus,(i,j,q))
            end
        end
    end
    vortices_plus_no_duplicates  = merge_duplicates(vortices_plus,Lx,Ly)
    vortices_minus_no_duplicates = merge_duplicates(vortices_minus,Lx,Ly)

    return vortices_plus_no_duplicates,vortices_minus_no_duplicates
end

function relax!(thetas::Matrix{FT},trelax=0.3) where FT<:AbstractFloat
    t = 0.0
    dt = 1E-2
    # T = 0.05
    Lx,Ly = size(thetas)

    thetas_old = copy(thetas)

    while t<trelax
        t += dt
        for j in 1:Ly
            for i in 1:Lx
                θ = thetas_old[i,j]
                angle_neighbours = get_neighbours(thetas_old,i,j,is_in_bulk(i,j,Lx,Ly))
                thetas[i,j] += dt*sum(sin,angle_neighbours .- θ)
                # thetas[i,j] += dt*sum(sin,angle_neighbours .- θ) + sqrt(2T*dt)*randn(FT)
            end
        end
    end

    return thetas
end


function merge_duplicates(list,Lx,Ly;radius=3)
    #= In this list, there might be doubles/triples (2/3 locations for the
    same physical vortex). We thus seek for numerically identified vortices
    which are neighbours and with the same charge to delete them. =#
    pos    = [list[i][1:2] for i in each(list)]
    charge = [list[i][3]   for i in each(list)]
    dealt_with = falses(length(pos))

    merged_duplicates = []
    for i in 1:length(pos)
        if !dealt_with[i]
            tmp = []
            for j in i:length(pos) # includes defect "i"
                if dist(pos[i],pos[j],Lx,Ly) ≤ radius && charge[i] == charge[j]
                    dealt_with[j] = true
                    push!(tmp,pos[j])
                end
            end
            mean_loc_defect = mean_N_positions(tmp,Lx,Ly,true)
            push!(merged_duplicates,(mean_loc_defect[1],mean_loc_defect[2],charge[i]))
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
Defect(;id,charge,loc,t) = Defect(id,charge,[loc],nothing,t,nothing)
first_loc(d::Defect) = d.pos[1]
last_loc(d::Defect)  = d.pos[end]
function update_position!(d::Defect,new_loc)
    push!(d.pos,new_loc)
end


mutable struct DefectTracker
    defectsP::Vector{Defect} # the id of a defect is its index in this vector
    defectsN::Vector{Defect} # so there is a (+)defect with id=1 AND and a (-)defect with id=1
    current_time::Float64 # latest update time (by convention, the creation time of the whole data structure = 0)

    function DefectTracker(pos,thetas,N,Lx,Ly,t) # constructor
        vortices,antivortices = spot_defects(pos,thetas,N,Lx,Ly)
        defectsP = [Defect(id=i,charge=vortices[i][3],loc=vortices[i][1:2],t=t) for i in each(vortices)]
        defectsN = [Defect(id=i,charge=antivortices[i][3],loc=antivortices[i][1:2],t=t) for i in each(antivortices)]
        new(defectsP,defectsN,t)
    end
end
number_defectsP(dt::DefectTracker) = length(dt.defectsP)
number_defectsN(dt::DefectTracker) = length(dt.defectsN)
number_defects(dt::DefectTracker)  = length(dt.defectsN)  + length(dt.defectsN)
number_active_defectsP(dt::DefectTracker) = count(isnothing,[d.annihilation_time for d in dt.defectsP])
number_active_defectsN(dt::DefectTracker) = count(isnothing,[d.annihilation_time for d in dt.defectsN])
number_active_defects(dt::DefectTracker)  = number_active_defectsN(dt) + number_active_defectsP(dt)

function ID_active_defects(dt::DefectTracker)
    activeP = Int[]
    for i in 1:number_defectsP(dt)
        if dt.defectsP[i].annihilation_time == nothing push!(activeP,i) end
    end
    activeN = Int[]
    for i in 1:number_defectsN(dt)
        if dt.defectsN[i].annihilation_time == nothing push!(activeN,i) end
    end
    return activeP,activeN
end

function add_defect!(dt::DefectTracker;charge,loc)
    if charge > 0 push!(dt.defectsP,Defect(id=1+number_defectsP(dt),charge=charge,loc=loc,t=dt.current_time))
    else          push!(dt.defectsN,Defect(id=1+number_defectsN(dt),charge=charge,loc=loc,t=dt.current_time))
    end
end

function pair_up_hungarian(dt::DefectTracker,new,old,Lx,Ly,charge::String)
    # charge can be either "+" or "-"
    distance_matrixx = distance_matrix(new,old,Lx,Ly) # m_new lignes, m_old colonnes
    proposal         = hungarian(distance_matrixx)[1] # length = length(new)
    assignment       = copy(proposal) # because it will be modified in the next for loop

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

function find_closest_before_annihilation(dt,Lx,Ly,old_loc_defect)
    distance = Inf ; ID_antidefect = -1 # dummy
    for i in each(dt.defectsN)
        isnothing(dt.defectsN[i].annihilation_time) ? annihilation_time_defect = nothing : annihilation_time_defect = round(dt.defectsN[i].annihilation_time,digits=2)
        if annihilation_time_defect == round(dt.current_time,digits=2) # it has just annihilated
            tmp = dist(old_loc_defect,last_loc(dt.defectsN[i]),Lx,Ly)
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
        println("Number of defects = ",number_active_defectsP(dft)," + ",number_active_defectsN(dft))
        return -1,(-1,-1)
    else
        return ID_antidefect,last_loc(dt.defectsN[ID_antidefect])
    end
end

function delete_defect(dft::DefectTracker,id::Int,charge::String)
    if charge == "+"
        popat!(dft.defectsP,id)
    else
        popat!(dft.defectsN,id)
    end
    decrease_annihilation_ids!(dft,id,charge)
    return dft
end

function decrease_annihilation_ids!(dft::DefectTracker,id::Int,charge::String)
    # the id and the charge are those of the defect that was deleted
    if charge == "+"
        for def in dft.defectsN
            if def.id_annihilator > id
                def.id_annihilator -= 1
            end
        end
    else
        for def in dft.defectsP
            if def.id_annihilator > id
                def.id_annihilator -= 1
            end
        end
    end
    return dft
end


function annihilate_defects(dt::DefectTracker,ids_annihilated_defects,Lx,Ly)
    for i in ids_annihilated_defects
        old_loc_vortex = last_loc(dt.defectsP[i])
        ID_antivortex,old_loc_antivortex = find_closest_before_annihilation(dt,Lx,Ly,old_loc_vortex)

        dt.defectsP[i].id_annihilator = ID_antivortex
        try dt.defectsN[ID_antivortex].id_annihilator = i
        catch
            dt.defectsN[ID_antivortex].id_annihilator = -1
        end

        estimate = mean_2_positions(old_loc_vortex,old_loc_antivortex,Lx,Ly)
        update_position!(dt.defectsP[i],estimate)
        try update_position!(dt.defectsN[ID_antivortex],estimate)
        catch ; end
    end
    return dt
end

function update_and_track!(dft::DefectTracker,pos::Matrix{FT},thetas::Vector{FT},psis::Vector{FT},omegas::Vector{FT},T::Number,v0::Number,N::Int,Lx::Int,Ly::Int,dt::Number,t::Number,times::AbstractVector) where FT<:AbstractFloat
    token = 1
    while t < tmax
        t += dt
        pos,thetas = update(pos,thetas,psis,omegas,T,v0,N,Lx,Ly,dt)
        if t ≥ times[token]
            if number_active_defects(dft) == 0
                println("Simulation stopped, there is no defects left.")
                break
            end
            println("t = ",round(t,digits=1)," & n(t) = ",number_active_defectsP(dft)," + ",number_active_defectsN(dft))
            update_DefectTracker!(dft,pos,thetas,N,Lx,Ly,t)
            token = min(token+1,length(times))
        end
    end
    return dft,pos,thetas,t
end

function update_DefectTracker!(dt::DefectTracker,pos::Matrix{T},thetas::Vector{T},N,Lx,Ly,t) where T<:AbstractFloat
    dt.current_time = t
    vortices_new,antivortices_new = spot_defects(pos,thetas,N,Lx,Ly)

    # if BC == "periodic" @assert length(vortices_new) == length(antivortices_new) && length(vortices_old) == length(antivortices_old) end
    locP_old    = [last_loc(dt.defectsP[i]) for i in each(dt.defectsP)]
    locN_old    = [last_loc(dt.defectsN[i]) for i in each(dt.defectsN)]
    chargeP_old = [dt.defectsP[i].charge    for i in each(dt.defectsP)]
    chargeN_old = [dt.defectsN[i].charge    for i in each(dt.defectsN)]

    locP_new    = [vortices_new[i][1:2]     for i in each(vortices_new)]
    locN_new    = [antivortices_new[i][1:2] for i in each(antivortices_new)]
    chargeP_new = [vortices_new[i][3]       for i in each(vortices_new)]
    chargeN_new = [antivortices_new[i][3]   for i in each(antivortices_new)]

    Np_new,Np_old = length(locP_new),length(locP_old)
    Nn_new,Nn_old = length(locN_new),length(locN_old)
    N_old = Np_old + Nn_old
    N_new = Np_new + Nn_new

    # Special simple cases to deal with upstream
    if N_new == N_old == 0 # do nothing, otherwise, "reducing over empty collection blablabla"

    elseif Nn_new == Nn_old == 0 && Np_new == Np_old > 0 # there are only (+) defects and no creation/annihilation
        assignment_vortices = pair_up_hungarian(dt,locP_new,locP_old,Lx,Ly,"+")
        for i in 1:Np_new update_position!(dt.defectsP[assignment_vortices[i]],locP_new[i]) end

    elseif Np_new == Np_old == 0 && Nn_new == Nn_old > 0 # there are only (-) defects and no creation/annihilation
        assignment_antivortices = pair_up_hungarian(dt,locN_new,locN_old,Lx,Ly,"-")
        for i in 1:Nn_new update_position!(dt.defectsN[assignment_antivortices[i]],locN_new[i]) end

    elseif N_new > 0 && N_old == 0
        for i in 1:Np_new add_defect!(dt,charge=chargeP_new[i],loc=locP_new[i]) end
        for i in 1:Nn_new add_defect!(dt,charge=chargeN_new[i],loc=locN_new[i]) end

    elseif N_new == 0 && N_old > 0 # (+)(-) >> plus rien
        id_just_annihilated_defectP,id_just_annihilated_defectM = ID_active_defects(dt) # seek for not yet annihilated defects

        for i in id_just_annihilated_defectP dt.defectsP[i].annihilation_time = dt.current_time end
        for i in id_just_annihilated_defectM dt.defectsN[i].annihilation_time = dt.current_time end
        dt = annihilate_defects(dt::DefectTracker,id_just_annihilated_defectP,L)

    elseif Np_new > 0 && Np_old > 0 && Nn_old > 0 && Nn_new == 0  # (+)(+)(-) >> (+) par exemple
        assignment_vortices = pair_up_hungarian(dt,locP_new,locP_old,Lx,Ly,"+")
        # Update living vortices. NB : the annihilated vortex is absent from the assignment vector : proceed without the condition "≠ 0"
        for i in eachindex(assignment_vortices) update_position!(dt.defectsP[assignment_vortices[i]],locP_new[i]) end
        # Identify annihilated defects
        ID_annihilated_vortices = [] ; ID_annihilated_antivortices = []
        for i in 1:number_defectsP(dt)
            if i ∉ assignment_vortices && dt.defectsP[i].annihilation_time == nothing # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
                push!(ID_annihilated_vortices,i)
            end
        end
        ID_annihilated_antivortices = ID_active_defects(dt)[2] # in this special case, there is no antivortices left in the "new" timestep

        for i in ID_annihilated_vortices     dt.defectsP[i].annihilation_time = dt.current_time end
        for i in ID_annihilated_antivortices dt.defectsN[i].annihilation_time = dt.current_time end
        dt = annihilate_defects(dt,ID_annihilated_vortices,Lx,Ly)

    elseif Nn_new > 0 && Nn_old > 0 && Np_old > 0 && Np_new == 0  # (+)(-)(-) >> (-) par exemple
        assignment_antivortices = pair_up_hungarian(dt,locN_new,locN_old,Lx,Ly,"-")
        # Update living antivortices. NB : the annihilated antivortex is absent from the assignment vector : proceed without the condition "≠ 0"
        for i in eachindex(assignment_antivortices) update_position!(dt.defectsN[assignment_antivortices[i]]) end
        # Identify annihilated defects
        ID_annihilated_vortices = [] ; ID_annihilated_antivortices = []
        for i in 1:number_defectsN(dt)
            if i ∉ assignment_antivortices && dt.defectsN[i].annihilation_time == nothing # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
                push!(ID_annihilated_antivortices,i)
            end
        end
        ID_annihilated_vortices = ID_active_defects(dt)[1] # in this special case, there is no vortices left in the "new" timestep

        for i in ID_annihilated_vortices     dt.defectsP[i].annihilation_time = dt.current_time end
        for i in ID_annihilated_antivortices dt.defectsN[i].annihilation_time = dt.current_time end

        dt = annihilate_defects(dt,ID_annihilated_vortices,Lx,Ly)
    else # end of special cases

    # GENERAL TREATMENT
        assignment_vortices     = pair_up_hungarian(dt,locP_new,locP_old,Lx,Ly,"+")
        assignment_antivortices = pair_up_hungarian(dt,locN_new,locN_old,Lx,Ly,"-")

        # CASE 1 : no creation, no annihilation : simply update the data structure
        if N_new == N_old
            for i in 1:Np_new update_position!(dt.defectsP[assignment_vortices[i]],locP_new[i]) end
            for i in 1:Nn_new update_position!(dt.defectsN[assignment_antivortices[i]],locN_new[i]) end

        # CASE 2 : creation !
    elseif N_new > N_old
            # Take care of the newly created defects
            ind_created_vortex = findall(iszero,assignment_vortices) # newly created vortex -> the assignment vector contains a 0
            loc_created_vortex = vortices_new[ind_created_vortex]
            for j in each(loc_created_vortex) add_defect!(dt,charge=chargeP_new[j],loc=loc_created_vortex[j][1:2]) end

            ind_created_antivortex = findall(iszero,assignment_antivortices)
            loc_created_antivortex = antivortices_new[ind_created_antivortex]
            for j in each(loc_created_antivortex) add_defect!(dt,charge=chargeN_new[j],loc=loc_created_antivortex[j][1:2]) end

            # Update the ancient defects' positions
            for i in eachindex(assignment_vortices)
                if assignment_vortices[i] ≠ 0 # avoid newly created defects
                    update_position!(dt.defectsP[assignment_vortices[i]],locP_new[i])
                end
            end
            for i in eachindex(assignment_antivortices)
                if assignment_antivortices[i] ≠ 0 # avoid newly created defects
                    update_position!(dt.defectsN[assignment_antivortices[i]],locN_new[i])
                end
            end

        # CASE 3 : annihilation !
    elseif N_new < N_old
             # Update living vortices. NB : the annihilated vortex is absent from the assignment vector : proceed without the condition "≠ 0"
             for i in eachindex(assignment_vortices)     update_position!(dt.defectsP[assignment_vortices[i]],locP_new[i]) end
             for i in eachindex(assignment_antivortices) update_position!(dt.defectsN[assignment_antivortices[i]],locN_new[i]) end

            # Identify annihilated defects
            ID_annihilated_vortices = [] ; ID_annihilated_antivortices = []
            for i in 1:number_defectsP(dt)
                if i ∉ assignment_vortices && dt.defectsP[i].annihilation_time == nothing # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
                    dt.defectsP[i].annihilation_time = dt.current_time # from now on, defects that have just annihilated have annihilation_time == t
                    push!(ID_annihilated_vortices,i)
                end
            end
            for i in 1:number_defectsN(dt)
                if i ∉ assignment_antivortices && dt.defectsN[i].annihilation_time == nothing
                    dt.defectsN[i].annihilation_time = dt.current_time
                    push!(ID_annihilated_antivortices,i)
                end
            end
            if length(ID_annihilated_antivortices) >= length(ID_annihilated_vortices)
                dt = annihilate_defects(dt,ID_annihilated_vortices,Lx,Ly)
            else
                dt = annihilate_defects(dt,ID_annihilated_antivortices,Lx,Ly)
            end

        end # end of general treatment
    end # end of special cases & general treatment
    return dt
end

## Defects Analysis : MSD
function MSD(dfts::Union{Vector{DefectTracker},Vector{Union{Missing,DefectTracker}}},Lx,Ly)
    indices = [] # indices of dft defined (is simulation not finished, dfts[i] == missing)
    for i in 1:length(dfts)
        if !ismissing(dfts[i]) push!(indices,i) end
    end
    maxlength = maximum([maximum([length(d.pos) for d in vcat(dft.defectsP,dft.defectsN)]) for dft in dfts[indices]])
    MSD_P   = NaN*zeros(length(indices),maxlength)
    MSD_N   = NaN*zeros(length(indices),maxlength)
    MSD_all = NaN*zeros(length(indices),maxlength)
    for i in 1:length(indices)
        msd_all, msd_p, msd_n = MSD(dfts[indices[i]],Lx,Ly)
        MSD_P[i,1:length(msd_p)] = msd_p
        MSD_N[i,1:length(msd_n)] = msd_n
        MSD_all[i,1:length(msd_all)] = msd_all
    end

    MSD_P_avg = nanmean(MSD_P,1)[1,:]
    MSD_N_avg = nanmean(MSD_N,1)[1,:]
    MSD_all_avg = nanmean(MSD_all,1)[1,:]

    return MSD_all_avg,MSD_P_avg,MSD_N_avg
end

function MSD(dft::DefectTracker,Lx,Ly,maxlength=nothing)
    nP = number_defectsP(dft)
    nN = number_defectsN(dft)
    # tmin,tmax = t_bounds(dft) # (tmin,tmax) = timestamps of (first defect creation , last defect annihilation)

    # hasfield(typeof(model),:dt) ? dummy_dt = model.dt : dummy_dt = 1
    if isnothing(maxlength)
        maxlength = maximum([length(d.pos) for d in vcat(dft.defectsP,dft.defectsN)])
    end
    # Compute the SD
    SD_P = NaN*zeros(nP,maxlength)
    SD_N = NaN*zeros(nN,maxlength)
    for n in 1:nP
        defect = dft.defectsP[n]
        tmp = square_displacement(defect,Lx,Ly)
        SD_P[n,1:length(tmp)] = tmp
    end
    for n in 1:nN
        defect = dft.defectsN[n]
        tmp = square_displacement(defect,Lx,Ly)
        SD_N[n,1:length(tmp)] = tmp
    end

    # Now average to compute the MSD
    MSD_P = nanmean(SD_P,1)[1,:]
    MSD_N = nanmean(SD_N,1)[1,:]
    MSD_all = nanmean(hcat(MSD_P,MSD_N),2)[:]

    return MSD_all, MSD_P, MSD_N
end

function square_displacement(d::Defect,Lx,Ly)
    loc_t0 = first_loc(d)
    return [dist(loc,loc_t0,Lx,Ly) for loc in d.pos] .^ 2
end

## Defects Analysis : Distance between defects
function interdefect_distance(dft::DefectTracker,defect1::Defect,defect2::Defect,Lx,Ly)
    # TODO take care of case with creation and/or annihilation time different.
    # So far, this care is left to the user...
    # @assert defect1.creation_time == defect2.creation_time
    # @assert defect1.annihilation_time == defect2.annihilation_time
    tmax = min(length(defect1.pos),length(defect2.pos))
    R = [dist(defect1.pos[t],defect2.pos[t],Lx,Ly) for t in 1:tmax]
    return R
end

function mean_distance_to_annihilator(dfts::Union{Vector{DefectTracker},Vector{Union{Missing,DefectTracker}}},Lx,Ly)
    indices = [] # indices of dft defined (is simulation not finished, dfts[i] == missing)
    for i in 1:length(dfts)
        if !ismissing(dfts[i]) push!(indices,i) end
    end
    Rs = [mean_distance_to_annihilator(dfts[indices[n]],Lx,Ly) for n in 1:length(indices)]
    Rs_matrix = vector_of_vector2matrix(Rs)
    return nanmean(Rs_matrix,2)[:,1]
end

function mean_distance_to_annihilator(dft::DefectTracker,Lx,Ly)
    nP = number_defectsP(dft)
    Rs = [distance_to_annihilator(dft,n,Lx,Ly) for n in 1:nP]
    Rs_matrix = vector_of_vector2matrix(Rs)
    return nanmean(Rs_matrix,2)[:,1],nanstd(Rs_matrix,2)[:,1]
end

function distance_to_annihilator(dft::DefectTracker,id1::Int,Lx,Ly;reversed=true)
    if isnothing(dft.defectsP[id1].id_annihilator) # not yet annihilated
        return [NaN]
    else
        R = interdefect_distance(dft,dft.defectsP[id1],dft.defectsN[dft.defectsP[id1].id_annihilator],Lx,Ly)
        if reversed reverse!(R) end
        return R
    end
end

## Auxiliary functions
function mean_2_positions(pos1,pos2,Lx,Ly,should_take_mod::Bool=true)
    a,b = pos1 ; x,y = pos2

    dx = (x - a) #; dx = min(dx,L-dx)
    dy = (y - b) #; dy = min(dy,L-dy)

    if should_take_mod
        if abs(Lx-dx) < abs(dx) dx = -(Lx-dx) end
        if abs(Ly-dy) < abs(dy) dy = -(Ly-dy) end
        # return mod1.((a,b) .+ 0.5.*(dx,dy),L)
        return (mod1(a .+ 0.5*dx,Lx) , mod1(b .+ 0.5*dy,Ly))
    else return (a,b) .+ 0.5.*(dx,dy)
    end
end
# l = 100
# mean_2_positions((50,50),(60,60),l) == (55,55)
# mean_2_positions((10,10),(90,90),l) == (100,100)
# mean_2_positions((49,66),(51,61),l) == (50.0, 63.5)

function mean_N_positions(vec_pos,Lx,Ly,should_take_mod::Bool=true)
    averaged_pos = vec_pos[1]
    for i in 2:length(vec_pos)
        averaged_pos = mean_2_positions(averaged_pos,vec_pos[i],Lx,Ly,should_take_mod)
    end
    return averaged_pos
end

function smooth(x,c=1)
    len = length(x)
    result = zeros(len)
    result[1:c] .= NaN
    result[end-c+1:end] .= NaN
    for i in 1+c:len-c
        result[i] = mean(x[i-c:i+c])
    end
    return result
end

function corr_length(C::Vector{T},rs=1:length(C);seuil=exp(-1))::T where T<:AbstractFloat # from a time series, returns the correlation length ()
    i_after = findfirst(x->x<seuil,C)
    if i_after ≠ nothing && i_after > 1
        # Linear interpolation
        i_befor = i_after - 1
        r_after = rs[i_after]
        r_befor = rs[i_befor]
        c_after = C[i_after]
        c_befor = C[i_befor]
        ξ = (seuil*(r_after-r_befor) -  (c_befor*r_after - r_befor*c_after))/(c_after-c_befor)
    else
    ξ = NaN
    end
    return ξ
end

function remove_negative(input)
    array = Float64.(copy(input))
    for i in 1:length(array)
        if array[i] ≤ 0 array[i] = NaN end
    end
    return array
end

## Plotting methods
# import Plots.plot
# function plot(pos,thetas,N,Lx,Ly;particles=false,vertical=false,size=(512,512),defects=false,title="")
#     cols = cgrad([:black,:blue,:green,:orange,:red,:black])
#     if particles
#         if vertical
#             p1 = scatter((pos[1,:],pos[2,:]),marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/Lx,size=size,aspect_ratio=Ly/Lx,xlims=(0,Lx),ylims=(0,Ly))
#             thetas_cg = cg(pos,thetas,N,Lx,Ly)
#             p2 = heatmap(mod.(thetas_cg,2pi)',clims=(0,2pi),c=cols,size=size,aspect_ratio=Ly/Lx)
#             if defects
#                 defects_p,defects_m = spot_defects(pos,thetas,N,Lx,Ly)
#                 locP = [defects_p[i][1:2] for i in each(defects_p)]
#                 locN = [defects_m[i][1:2] for i in each(defects_m)]
#                 highlight_defects!(p2,Lx,Ly,locP,locN)
#             end
#             return plot(p1,p2,layout=(2,1),size=(size[1],2*size[2]),title=title)
#         else
#             p1 = scatter((pos[1,:],pos[2,:]),marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/Lx,size=size,aspect_ratio=Ly/Lx,xlims=(0,Lx),ylims=(0,Ly))
#             thetas_cg = cg(pos,thetas,N,Lx,Ly)
#             p2 = heatmap(mod.(thetas_cg,2pi)',clims=(0,2pi),c=cols,size=size,aspect_ratio=Ly/Lx)
#             if defects
#                 defects_p,defects_m = spot_defects(pos,thetas,N,Lx,Ly)
#                 locP = [defects_p[i][1:2] for i in each(defects_p)]
#                 locN = [defects_m[i][1:2] for i in each(defects_m)]
#                 highlight_defects!(p2,Lx,Ly,locP,locN)
#             end
#             return plot(p1,p2,size=(2*size[1],size[2]),title=title)
#         end
#     else
#         thetas_cg = cg(pos,thetas,N,Lx,Ly)
#         p2 = heatmap(mod.(thetas_cg,2pi)',clims=(0,2pi),c=cols,size=size,aspect_ratio=Ly/Lx,title=title)
#         if defects
#             defects_p,defects_m = spot_defects(pos,thetas,N,Lx,Ly)
#             locP = [defects_p[i][1:2] for i in each(defects_p)]
#             locN = [defects_m[i][1:2] for i in each(defects_m)]
#             highlight_defects!(p2,Lx,Ly,locP,locN)
#         end
#         return p2
#     end
# end
#
# function highlight_defects!(p,Lx,Ly,defects_p,defects_m,symbP=:circle,symbM=:utriangle)
#     for defect in defects_p
#         scatter!((defect), m = (1.5, 1., symbP,:transparent, stroke(1.2, :grey85)))
#     end
#     for defect in defects_m
#         scatter!((defect), m = (1.5, 1., symbM,:transparent, stroke(1.2, :grey85)))
#     end
#     xlims!((1,Lx))
#     ylims!((1,Ly))
#     return p
# end
