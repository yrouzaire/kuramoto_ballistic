## Investigation C(r,t) when \sigma = 0
@unpack rhos,Ts,Ns,v0s,tmax,sigmas,times,P,C,runtimes = load("data/OP_C_sigma_0.jld")
histogram(runtimes/3600,bin=20)
P
Pavg = nanmean(P,7)
Cavg = nanmean(C,8)
p=plot()
    for i in 1:length(Ts)
        # plot!(0:0.5:round(Int,sqrt(Ns)/2),remove_negative(Cavg[:,1,1,1,1,1,t,1]),yaxis=:log)
        plot!(times[2:end]*(Ts[i]),remove_negative(Pavg[1,1,i,1,1,2:end,1]),axis=:log)
    end
    plot!(x->2E-2(x/log(40x))^0.5,c=:black)
    p

p=plot()
    for i in 1:1:length(Ts)
        plot!(0:0.5:round(Int,sqrt(Ns)/2),remove_negative(Cavg[:,1,1,i,1,1,end,1]),yaxis=:log)
    end
    # xlims!(0,10)
    # ylims!(0.2,1)
    p

p=plot(ylims=(5E-3,1.2))
    for i in 1:10:length(times)
        plot!(0:0.5:round(Int,sqrt(Ns)/2),remove_negative(Cavg[:,1,1,1,1,1,i,1]),yaxis=:log)
    end
    p

## Study of pairs of defects when \sigma = 0
gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
N = Int(1E4)
    rho = 1
    L = round(Int,sqrt(N/rho)) ; R0 = 1
    T = 0.1
    v0 = 1
    σ = 0

    tmax = 10 ; dt = determine_dt(T,σ,v0,N,rho)

t = 0.0
params_init = ["pair",L/2]
pos,vel_angles,thetas,omegas = initialisation(N,L,σ,params_init)
    scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/L,size=(512,512),xlims=(0,L),ylims=(0,L))

z = @elapsed while t<10
    t += dt
    global pos,thetas = update(pos,vel_angles,thetas,omegas,T,v0,N,L,rho,R0,σ,dt)
end
    prinz(z)
    scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/L,size=(512,512),xlims=(0,L),ylims=(0,L))

# Coarse graining
cgmat,cgmat_interpolated = cg(pos,thetas,N,L,2R0)
heatmap(mod.(cgmat,2π)',color=cols)
heatmap(mod.(cgmat_interpolated,2π)',color=cols)
exp(1im*NaN)

(1-exp(-0.5))*1.35

##  Analysis dependency Order on vo for sigma>0
@unpack Ns,rhos,Ts,v0s,sigmas,tmax,times,P,runtimes = load("data/onset_order_v0_rhos.jld")
@unpack Ns,rhos,Ts,v0s,sigmas,tmax,times,P,runtimes = load("data/onset_order_v0_rhos_normalized_by_rho.jld")
histogram(runtimes/3600,bins=20,title="Médiane = $(round(median(filter(!isnan,runtimes)/3600),digits=2))",xlabel="Time [hours]")
R = size(P)[end]
Pavg = nanmean(P,7) ; Pstd = nanstd(P,7) # Ns,rhos,Ts,v0s,sigmas,times,1
# L = sqrt(Ns/rhos)
# Find critical sigma
p=plot(xlabel="σ",ylabel="P",legend=:outerright,size=(500,400))
    for i in 5:4:length(v0s) #each(v0s)
        plot!(sigmas/v0s[i],Pavg[1,4,1,i,:,50,1],rib=0Pstd[1,1,1,i,:,50,1],label=L"v_0 = "*string(round(v0s[i],digits=2)))
        # plot!(sigmas,Pavg[1,4,1,i,:,50,1],rib=0Pstd[1,1,1,i,:,50,1],label=L"v_0 = "*string(round(v0s[i],digits=2)))
    end
    plot!(sigmas,0.3*ones(length(sigmas)),c=:grey,line=:dash)
    p
# savefig("figures\\dependency_P_v0.pdf")
&
# p=plot(xlabel=L"v_0",ylabel=L"σ_c",legend=:outerright,size=(600,400))
    tt = 15
    # for rhoo in length(rhos)
    for rhoo in 2:length(rhos)
    seuil = 0.3
    sigmacs = NaN*ones(length(v0s),R)
    for i in each(v0s) , r in 1:R
        ind = findfirst(x->x < seuil,P[1,rhoo,1,i,:,tt,r])
        try
            # println("ind = $ind")
            if ind > 1
                aa = (P[1,rhoo,1,i,ind,tt,r] - P[1,rhoo,1,i,ind-1,tt,r])/(sigmas[ind] - sigmas[ind-1])
                b = P[1,rhoo,1,i,ind-1,tt,r] - sigmas[ind-1]*aa
                sigmacs[i,r] = (seuil - b)/aa
            elseif ind == 1
                sigmacs[i,r] = 0
            end
        catch e ;
        end
    end
    plot!(v0s,nanmean(sigmacs,2)[:,1],rib=0nanstd(sigmacs,2)[:,1],m=false,label="ρ = $(rhos[rhoo])",line=:dash)
    end
    p
# savefig("figures\\sigmacrit_vs_v0.png")
# start = 19
    # plot!(v0s[start:end],0.53sqrt.(v0s[start:end] .- v0s[start]),c=:black)
    # annotate!(1.5,0.15,text(L"v_{0,c} = "*string(round(v0s[start],digits=2)),10))
    # v_crit = 0.37
    # plot!(x->0.5*sqrt(x-v_crit),c=:black,lw=1.8)
    # annotate!(1.3,0.15,text(L"\sigma_c \sim \sqrt{v_0 - v_{0,c}}"),10)
    # annotate!(1.3,0.07,text(L"v_{0,c} = "*string(v_crit),10))
    # lens!([0.3,0.45], [-0.02, 0.1], inset = (1, bbox(0.02, 0.0, 0.4, 0.4)))
# savefig("figures/critical_sigma_v0.png")

##
#= To do : some screenshots for \sigma = 0
 in the ordered phase to show that there is no spin waves
 at the onset of disorder to show that there are defects
 =#

# No spin waves
N = Int(1E4)
    rho = 0.444
    L = round(Int,sqrt(N/rho))
    R0 = 1
    T = 0.001
    v0 = 30
    σ = 0.5
    tmax = 100 ; dt = determine_dt(T,σ,v0,N,rho)

R = 10
z = @elapsed for r in 1:R
    pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
    t = 0.0
    while t<tmax
        t += dt
        pos,thetas = update(pos,vel_angles,thetas,omegas,T,v0,N,L,rho,R0,σ,dt)
    end
    scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=325/L,size=(512,512),xlims=(0,L),ylims=(0,L),title="P=$(round(polarOP(thetas)[1],digits=2))")
    savefig("figures/looking_for_spin_waves/$r.png")
    println("Simulation $r/$R done.")
end
prinz(z)

gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=400/L,size=(512,512),xlims=(0,L),ylims=(0,L),title="P=$(round(polarOP(thetas)[1],digits=2))")

# Presence of defects at the onset of order
N = Int(1E5)
    rho = 1
    L = round(Int,sqrt(N/rho))
    R0 = 1
    T = 0.4
    v0 = 2
    σ = 0.
    tmax = 50 ; dt = determine_dt(T,σ,v0,N,rho)

R = 1
z = @elapsed for r in 1:R
    pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
    t = 0.0
    while t<tmax
        t += dt
        pos,thetas = update(pos,vel_angles,thetas,omegas,T,v0,N,L,rho,R0,σ,dt)
    end
    titre = "P=$(round(polarOP(thetas)[1],digits=2))"
    scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=325/L,size=(512,512),xlims=(0,L),ylims=(0,L),title=titre)
    savefig("figures/presence_defects_onset_order/$r.png")
    println("Simulation $r/$R done.")
end
prinz(z)

## Movies
gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

function movies(params,every,tmax,dt)
    rho,T,v0,N,L,R0,σ = params

    pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
    println("N = $N")
    times = 1:every*dt:round(Int,tmax)
    Q = zeros(length(times)) ; P = zeros(length(times))
    anim = @animate for i in 1:length(times)
        println("$(round(i/length(times)*100,digits=2)) %")
        for j in 1:every pos,thetas = update(pos,vel_angles,thetas,omegas,T,v0,N,L,rho,R0,σ,dt) end
        P[i] = polarOP(thetas)[1]
        Q[i] = nematicOP(thetas)[1]
        titre = "P=$(round(P[i],digits=2)) , Q=$(round(Q[i],digits=2))"
        p1 = scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/L,size=(512,512),xlims=(0,L),ylims=(0,L),title=titre,aspect_ratio=1)
        p2 = plot(ylims=(0,1),size=(512,512))
            plot!(times,P,c=1,label="Polar",legend=:topleft,rib=0)
            plot!(times,Q,c=2,label="Nematic",rib=0)
        p = plot(p1,p2,size=(1024,512))
        end
    return anim
end

# Controlled parameters
rho = 1
    T = 0.4 # temperature for angle diffusion
    v0 = 1 # norm of individual velocities
    σ = 0

    # Other parameters
    N = Int(5E4)
    L = round(Int,sqrt(N/rho))
    R0 = 1

    params = Any[rho,T,v0,N,L,R0,σ] # any to avoid N being interpreted as a Float
    dt = determine_dt(T,σ,v0,N,rho)

pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/L,size=(512,512),xlims=(0,L),ylims=(0,L),aspect_ratio=1,title="Test")

every = 10 ; tmax = 2500
z = @elapsed anim = movies(params,every,tmax,dt)
prinz(z)
# mp4(anim,"films/film_sigma_$(σ)_T_$(T)_vo_$(v0)_alpha_$(α).mp4",fps=20)
# mp4(anim,"films/close_to_transition_sigma_$(σ)_T_$(T)_vo_$(v0)_rho_$(rho).mp4",fps=35)
mp4(anim,"films/TEST.mp4",fps=35)
## Sanity check: final screenshot for σ = 0
gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
N = Int(5E4)
    rho = 1
    L = round(Int,sqrt(N/rho)) ; R0 = 1
    T = 0.05
    v0 = 25
    σ = 0

    tmax = 250 ; dt = determine_dt(T,σ,v0,N,rho)

t = 0.0
pos,vel_angles,thetas,omegas = initialisation(N,L,σ)

z = @elapsed while t<tmax
    t += dt
    global pos,thetas = update(pos,vel_angles,thetas,omegas,T,v0,N,L,rho,R0,σ,dt)
end
prinz(z)

scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/L,size=(512,512),xlims=(0,L),ylims=(0,L))

dr = 0.5
z = @elapsed c = corr_fast(pos,thetas,N,L,dr)
# c_long = corr(pos,thetas,N,L,dr)
plot(0:dr:round(Int,L/2),c)

## Identification of influence of \sigma
#= We use the same thermal realisation, the same initial
distribution of angles and the same distribution of intrinsic frequencies.
The only variable is the value of \sigma. First, let's try to visualize the
difference with movies.
Eventually, I reckon the supposed reentrant behaviour was in fact a silly issue
of dt too large, leading to OP = 0 for low sigmas. =#
using Random
gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()


function update_seed(pos::Matrix{Float64},vel_angles::Vector{Float64},thetas::Vector{Float64},omegas::Vector{Float64},T,v0,N,L,rho,R0,σ,dt,seed)
    pos_new    = zeros(size(pos))
    thetas_new = zeros(size(thetas))

    Random.seed!(seed)
    vec_noise = randn(N)

    # Update position
    for n in eachindex(vel_angles)
        pos_new[:,n] = pos[:,n] + v0*dt*[cos(vel_angles[n]),sin(vel_angles[n])]
    end
    pos_new = mod.(pos_new,L)

    ## List construction
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
        neighbouring_cells = [ [cellx,celly] , [cellx,mod1(celly+1,nb_cells_1D)] , [mod1(cellx+1,nb_cells_1D),celly] , [cellx,mod1(celly-1,nb_cells_1D)] , [mod1(cellx-1,nb_cells_1D),celly] , [mod1(cellx+1,nb_cells_1D),mod1(celly+1,nb_cells_1D)] ,  [mod1(cellx-1,nb_cells_1D),mod1(celly-1,nb_cells_1D)] , [mod1(cellx-1,nb_cells_1D),mod1(celly+1,nb_cells_1D)] , [mod1(cellx+1,nb_cells_1D),mod1(celly-1,nb_cells_1D)]]
        ind_neighbours = []
        for (i,j) in neighbouring_cells
            next = head[i,j]
            if next ≠ -1
                if next ≠ n  # do not include yourself in your neighbourhood
                    if dist(poscur,pos[:,next],L) < R0
                        push!(ind_neighbours,next)
                    end
                end
                while list[next] ≠ -1
                    if dist(poscur,pos[:,list[next]],L) < R0
                        push!(ind_neighbours,list[next])
                    end
                    next = list[next]
                end
            end
        end

        thetas_neighbours = thetas[ind_neighbours]

        # Sync Theta
        if length(thetas_neighbours) > 0
            thetas_new[n] = thetas[n] + dt * omegas[n] + dt * sum(sin.((thetas_neighbours .- thetas[n]))) + sqrt(2T*dt)*vec_noise[n]
        else
            thetas_new[n] = thetas[n] + dt * omegas[n] + sqrt(2T*dt)*vec_noise[n]
        end

    end
    return pos_new,thetas_new
end

N = Int(1E4)
    rho = 2
    L = round(Int,sqrt(N/rho)) ; R0 = 1
    T = 0.25
    v0 = 0.4
    σ = 0.1

    tmax = 30 ; every = 10 ; dt = 0.1
    times = 0:every*dt:round(Int,tmax)


# Initialisation
Random.seed!(1) ; pos = L*rand(2,N)
    Random.seed!(2) ; vel_angles = 2π*rand(N)
    Random.seed!(3) ; thetas = 2π*rand(N)
    Random.seed!(4) ; omegas = σ*randn(N)

z = @elapsed anim = @animate for i in 1:length(times)
    println("$(round(i/length(times)*100,digits=2)) %")
    for j in 1:every
        global pos,thetas = update_seed(pos,vel_angles,thetas,omegas,T,v0,N,L,rho,R0,σ,dt,i+j)
    end
    scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=275/L,size=(512,512),xlims=(0,L),ylims=(0,L),title="P=$(round(polarOP(thetas)[1],digits=2))")
end
prinz(z)
# mp4(anim,"films/reentrant_behaviour/N4_BIS_sigma=$(σ)_rho=$(rho)_T=$(T).mp4",fps=25)

## Order Disorder Transition
@unpack Ns,rhos,Ts,v0s,sigmas,tmax,times,P,runtimes = load("data/order_disorder_transition.jld")
histogram(runtimes/3600,bin=20)
Pavg = nanmean(P,7)
p=plot(xlabel="T/ρ",ylabel="P",ylims=(0,1))
    for i in each(rhos), j in each(sigmas)
        plot!(Ts/rhos[i],Pavg[1,i,:,1,j,10,1],rib=0,label="ρ = $(rhos[i]), σ = $(sigmas[j])")
        # plot!(Ts,Pavg[1,i,:,1,j,end-1,1],label="ρ = $(rhos[i]), σ = $(sigmas[j])")
    end
    p
# savefig("figures\\order_disorder_transition.pdf")

## Analysis dependency Tc on vo (Tc ~ \alpha for all v0 ?)
@unpack Ns,alphas,Ts,v0s,sigmas,tmax,times,P,runtime = load("data/dependency_Tc_on_v0_2.jld")
Pavg = nanmean(P,7) ; Pstd = nanstd(P,7) # Ns,alphas,Ts,v0s,sigmas,times,1
heatmap(Ts,alphas,Pavg[1,:,:,30,1,end,1],clims=(0,1),c=cgrad([:red,:orange, :green]),size=(500,400),xlabel=L"v_0",ylabel="α",colorbar_title="P",title="t=$(times[end])")

p0 = plot(ylabel="Order Parameter",xlabel="T")
    plot!(Ts,Pavg[1,5,:,30,1,end,1],m=true)
    plot!(Ts,0.25*ones(length(Ts)),c=:grey,line=:dash)
    plot!(0.415ones(51),collect(0:0.01:0.5),c=:grey,line=:dash)

tt = 101
    seuil = 0.25
    Tcs = zeros(length(alphas),length(v0s))
    for i in each(v0s) , a in each(alphas)
        ind = findfirst(x->x < seuil,Pavg[1,a,:,i,1,end,1])
        # Tcs[a,i] = Ts[ind]
        aa = (Pavg[1,a,ind,i,1,tt,1] - Pavg[1,a,ind-1,i,1,tt,1])/(Ts[ind] - Ts[ind-1])
        b = Pavg[1,a,ind-1,i,1,tt,1] - Ts[ind-1]*aa
        Tcs[a,i] = (seuil - b)/aa
        # Tcs[a,i] = Ts[ind-1] + (Ts[ind]-Ts[ind-1])*(Pavg[1,a,ind-1,i,1,end,1] - seuil)/(Pavg[1,a,ind,i,1,end,1] - Pavg[1,a,ind-1,i,1,end,1])
    end
    # p1=plot(legend=:topleft,xlabel=L"v_0/L",ylabel=L"T_c",xlims=(-0,0.16),ylims=(0,0.7),size=(400,400))
p1=plot(legend=:topright,xlabel=L"v_0",ylabel=L"T_c",size=(400,400),ylims=(0,.5))
    for a in each(alphas)
        plot!(v0s*sqrt(pi*Ns[1]/alphas[a]),Tcs[a,:],label="ρ = $(round(alphas[a]/π,digits=2))",rib=0.0)
    end
    p1


p2=plot(legend=:bottomright,xlabel=L"v_0",ylabel=L"T_c/ρ",size=(400,400))
    for a in each(alphas)
        plot!(v0s*sqrt(pi*Ns[1]/alphas[a]),Tcs[a,:]/alphas[a]*pi,label="ρ = $(round(alphas[a]/π,digits=2)) ↔ L = $(round(Int,1/sqrt(alphas[a]/π/Ns)))",rib=0.0)
    end
    aa = 1.35 ; bb = 0.5 ;
    plot!(x->aa*(1-exp(-bb*x)),c=:black,line=:dash,lw=2)
    annotate!(12,1.025,text(L"f(x) = a\,(1-e^{-bx})",10))
    annotate!(12,0.9,text("a = $aa , b = $bb ",8))
    p2

plot(p0,p1,p2,layout=(1,3),size=(1200,400))
# savefig("figures\\dependency_Tc_v0.pdf")

## Analysis
@unpack Ns,alphas,Ts,v0s,sigmas,tmax,times,P,runtime = load("data/dependency_alpha_on_v0.jld")
Pavg = nanmean(P,7) ; Pstd = nanstd(P,7) # Ns,alphas,Ts,v0s,sigmas,times,1
heatmap(v0s,alphas,Pavg[3,:,2,:,1,end,1],clims=(0,1),c=cgrad([:red,:orange, :green]),size=(500,400),xlabel=L"v_0",ylabel="α",colorbar_title="P",title="t=$(times[end])")

v0_alpha_film = @animate for i in eachindex(times)
    # heatmap(v0s,alphas,Pavg[1,:,1,:,1,i,1],c=cgrad([:red,:orange, :green]),size=(500,400),xlabel=L"v_0",ylabel="α",colorbar_title="P",title="t=$(times[i])")
    heatmap(v0s,alphas,Pavg[3,:,1,:,1,i,1],clims=(0,1),c=cgrad([:red,:orange, :green]),size=(500,400),xlabel=L"v_0",ylabel="α",colorbar_title="P",title="t=$(times[i])")
end
mp4(v0_alpha_film,"films/phase_diagram_v0_alpha_T0.mp4",fps=20)

heatmap(times,alphas,Pavg[1,:,1,1,1,:,1])
heatmap(times,alphas,P[1,:,1,1,1,:,5])

## Influence of N on Phase Diag
#= We shall take 2-3 slices of the phase diag and see how they
evolve as a function of N (for T = 0, 0.2, 0.4) =#
Ns = Int.(2 .^ [7,9,10,11,12])
    alphas = 1
    Ts     = [0,0.2,0.4]
    v0s    = 0.2
    sigmas = 0:0.1:1.5

    tmax = 10 ; times = 0:0.2:tmax

    P   = zeros(length(Ns),length(alphas),length(Ts),length(v0s),length(sigmas),length(times))
    m = 0 ; M = length(Ns)*length(alphas)*length(Ts)*length(v0s)*length(sigmas)

z = @elapsed for n in eachindex(Ns) , i in eachindex(alphas) , j in eachindex(Ts) , k in eachindex(v0s) , q in eachindex(sigmas)
    N = Ns[n] ; α = alphas[i] ; T = Ts[j] ; v0 = v0s[k] ; σ = sigmas[q]
    global m += 1 ; println("Simulation $m/$M. ")

    L = 1
    R0 = L*sqrt(N*π/α)
    rho = N/L^2
    dt = 1E-2
    # dt = determine_dt(T,σ,v0,R0,L)
    params = Any[α,T,v0,N,L,rho,R0,σ]

    pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
    P[n,i,j,k,q,1] = polarOP(thetas)[1]
    t = 0.0 ; token = 2
    while t < tmax
        t += dt
        pos,thetas = update(pos,vel_angles,thetas,omegas,params,dt)
        if t > times[token]
            op =  polarOP(thetas)[1]
            P[n,i,j,k,q,token] = op
            # if op > 0.9
            #     P[n,i,j,k,q,token:end] .= 1
            #     println("Early stopping at t = $(times[token])")
            #     break
            # end
            token = min(token+1,length(times))
        end
    end
end # scanned parameters
prinz(z)

## Complexity of algo should be O(N)
Ns = round.(Int,2.0 .^ collect(4:12))
α = 1
    T     = 0.1
    v0    = 0.2
    σ = 0
    R = 10
    tmax = 2

    runtimes = zeros(length(Ns),R)

for n in eachindex(Ns)
    N = Ns[n]
    L = 1
    R0 = L*sqrt(α/π/N)
    rho = N/L^2
    println("N = $N started, R0 = $R0.")

    params = Any[α,T,v0,N,L,rho,R0,σ]
    dt = 0.01

    Threads.@threads for r in 1:R
        pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
        t = 0.0 ; token = 2
        z = @elapsed while t < tmax
            t += dt
            pos,thetas = update(pos,vel_angles,thetas,omegas,α,T,v0,N,L,rho,R0,σ,dt)
        end
        runtimes[n,r] = z
    end # realisations
end # scanned parameters

runtimesavg = mean(runtimes,dims=2)
plot!(Ns,runtimesavg[:,1],axis=:log)
# plot!(Ns,Ns .^ 1 * 8E-3,axis=:log,c=:black)



#=
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
.
=#

## Benchmark corr fast
N = 1000
    α = 1
    T  = 0.
    v0 = 1
    σ  = 0.
    R0 = 1
    R  = 10
    tmax = 10 ; times = 0:tmax-0.1:tmax

    dr = R0/2
    L = round(Int,sqrt(pi*N/α)); rho = N/L^2
    param = Any[α,T,v0,N,L,rho,R0,σ]

    pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
    t = 0.0 ; dt = 0.01
while t < 50
        t += dt
        pos,thetas = update(pos,vel_angles,thetas,omegas,param,dt)
    end

C = corr(pos,thetas,param,dr)
    plot(C)
Cfast = corr_fast(pos,thetas,param,dr)
    plot!(Cfast)

plott(pos,thetas)

function plott(pos,thetas,alpha=1)
    N = length(thetas)
    L = sqrt(pi*N/alpha)
    p = scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=500/L,size=(512,512),xlims=(0,L),ylims=(0,L))
    return p
end


## Efficiency of function corr
Ns = Int.(2.0 .^ collect(5:15))
    α = 1
    T  = 0.
    v0 = 1
    σ  = 0.
    R0 = 1
    R  = 10
    tmax = 10 ; times = 0:tmax-0.1:tmax

    dr = R0/2

    runtimes = zeros(length(Ns))
    runtimesfast = zeros(length(Ns))

for n in eachindex(Ns)
    N = Ns[n] ; L = round(Int,sqrt(pi*Ns[n]/α)) ; rho = N/L^2
    println("N = $N")
    params = Any[α,T,v0,Ns[n],L,rho,R0,σ]
    pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
    # runtimes[n] = @elapsed corr(pos,thetas,params,dr)
    runtimesfast[n] = @elapsed corr_fast(pos,thetas,params,dr)
end


# plot(Ns,runtimes,m=true,axis=:log)
plot(Ns,runtimesfast,m=true,axis=:log)
    plot!(Ns,Ns .^2 * 1E-6,m=true,axis=:log)

## Strong Finite Size Effects ? Close to the transitions yes
Ns = Int.(2.0 .^ collect(5:15))
    α = 1
    T  = 0.
    v0 = 1
    σ  = 0.
    R0 = 1
    R  = 1
    tmax = 10 ; times = 0:tmax-0.1:tmax

    P = zeros(length(Ns),R)
    dr = R0/2 ; Crt = Array{Vector{Float64},2}(undef,length(Ns),R)

z = @elapsed for n in eachindex(Ns)
    N = Ns[n]
    R0 = 1
    L = round(Int,sqrt.(π*N/α)*R0)
    rho = N/L^2
    # v0 = v0*L

    params = Any[α,T,v0,N,L,rho,R0,σ]
    dt = determine_dt(T,σ,v0,R0,L)

    println("N = $N started, with L = $L , dt = $dt")
    # for r in 1:R
    Threads.@threads for r in 1:R
        pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
        t = 0.0
        while t < tmax
            t += dt
            pos,thetas = update(pos,vel_angles,thetas,omegas,params,dt)
        end

        P[n,r] = polarOP(thetas)[1]
        Crt[n,r] = corr(pos,thetas,params,dr)
    end # realisations
end # scan parameters
prinz(z)

Pavg = mean(P,dims=2)
Cavg = Array{Vector{Float64},1}(undef,length(Ns))
    for n in eachindex(Ns)
        Cavg[n] = mean([Crt[n,r] for r in 1:R])
    end

p=plot(legend=:bottomright,ylims=(0,1))
    for n in eachindex(Ns)
        rr = collect(0:dr:round(Int,sqrt(R0^2*pi*Ns[n]/α)/2))
        plot!(rr,Cavg[n],label="N = $(Ns[n])")
    end
    p


## Analysis results cluster
@unpack alphas,Ts,v0s,tmax,times,P,runtime,Ns,R0,sigmas = load("data/phasediag_sig_v0_alpha=1_T0.jld")
P
collect(sigmas)
collect(v0s)
Pavg = mean(P,dims=7)
# Movies
phdiag_sig_vo = @animate for i in eachindex(times)
    heatmap(sigmas,v0s,Pavg[1,1,1,:,:,i,1],c=cgrad([:red,:orange, :green]),size=(500,400),clims=(0,1),xlabel="σ",ylabel="v_0/L",colorbar_title="P",title="t=$(times[i])")
end
mp4(phdiag_sig_vo,"films/phase_diagram_sig_v0_alpha=1_T0.mp4",fps=15)


## Analysis of simulations from cluster
@unpack alphas,Ts,v0s,tmax,times,P,runtime,N,L,R0,sigmas = load("data/phasediag_large_sigmas_Ts_alpha=1_v0=0.2.jld")
v0s
Ts
sigmas
P             #  (alphas),(Ts),(v0s),(sigmas),(times),R
Pavg = mean(P,dims=6)
Pavg = nanmean(P,6)

# Movies
anim_varying_time = @animate for i in eachindex(times)
    heatmap(sigmas,Ts,Pavg[1,:,1,:,i,1],c=cgrad([:red,:orange, :green]),size=(500,400),clims=(0,1),xlabel="σ",ylabel="T",colorbar_title="P",title="t=$(round(times[i],digits=1))")
end
mp4(anim_varying_time,"films/phase_diagram_large_sigT_varying_time_vo_0.2.mp4",fps=20)

# heatmap(sigmas,Ts,Pavg[1,:,1,:,end,1],c=cgrad([:red,:orange, :green]),size=(500,400),clims=(0,1),xlabel="σ",ylabel="T",colorbar_title="P",title="t=$(round(times[end],digits=1))")
heatmap(sigmas,Ts,Pavg[1,:,1,:,end,1],c=cgrad([:red,:orange, :green]),size=(500,400),clims=(0,1),xlabel="σ",ylabel="T",colorbar_title="P")#,title="t=$(round(times[end],digits=1))")
savefig("figures\\phdiag_sigmaT.png")


## Strong Finite Size Effects ? Close to the transitions yes
Ns = Int.(2.0 .^ collect(5:9))
    α = 1
    T  = 0.
    v0 = 0.2
    σ  = 0.
    R0 = 1
    R  = 10
    tmax = 10 ; times = 0:tmax-0.1:tmax

    P = zeros(length(Ns),length(times),R)
    dr = R0/2 ; Crt = Array{Vector{Float64},3}(undef,length(Ns),length(times),R)

z = @elapsed for n in eachindex(Ns)
    N = Ns[n]
    R0 = 1
    L = round(Int,sqrt.(π*N/α)*R0)
    rho = N/L^2
    v0 = v0*L
    println("N = $N started, with L = $L .")

    params = Any[α,T,v0,N,L,rho,R0,σ]
    dt = determine_dt(T,σ,v0,R0,L)

    # for r in 1:R
    Threads.@threads for r in 1:R
        pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
        P[n,1,r] = polarOP(thetas)[1]
        Crt[n,1,r] = [NaN]
        t = 0.0 ; token = 2
        while t < tmax
            t += dt
            pos,thetas = update(pos,vel_angles,thetas,omegas,params,dt)
            if t > times[token]
                if r == 1 println("Hello") end
                op =  polarOP(thetas)[1]
                c = corr(pos,thetas,params,dr)
                Crt[n,token,r] = c
                P[n,token,r] = op
                token = min(token+1,length(times))
            end
        end
    end # realisations
end # scan parameters
prinz(z)

Pavg = mean(P,dims=3)
Cavg = Array{Vector{Float64},2}(undef,length(Ns),length(times))
    for n in eachindex(Ns) , t in eachindex(times)
        # Cavg[n,t] = nanmean.([Crt[n,t,r] for r in 1:R])
        Cavg[n,t] = mean([Crt[n,t,r] for r in 1:R])
    end

# p=plot(legend=:bottomright,ylims=(0,1))
#     for n in eachindex(Ns)
#         plot!(times,Pavg[n,:,1],label="N = $(Ns[n])")
#     end
#     p
p=plot(legend=:bottomright,ylims=(0,1))
    for n in eachindex(Ns)
        rr = collect(0:dr:round(Int,sqrt(R0^2*pi*Ns[n]/α)/2))
        plot!(rr,Cavg[n,end],label="N = $(Ns[n])")
    end
    p



## First test for C(r,t)
α = 1 # average number of influenced neighbours
    T = 0.1 # temperature for angle diffusion
    v0 = .1 # norm of individual velocities
    σ = 0

    # Other parameters
    N = Int(1E4)
    R0 = 1
    rho = N/L^2
    L = round(Int,sqrt(π*N/α)*R0) # interaction radius
    println("L = $L")

    params = Any[α,T,v0*L,N,L,rho,R0,σ] # any to avoid N being interpreted as a Float
    dt = 5E-3 ; tmax = 2 ; times = 0:1:tmax
    dr = R0/2
    # determine_dt(T,σ,v0*L,R0,L)

    Crt = zeros(round(Int,L/2/dr)+1,length(times))

pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
t = 0.0 ; token = 2
z = @elapsed while t < tmax
    t += dt
    pos,thetas = update(pos,vel_angles,thetas,omegas,params,dt)
    if t > times[token]
        crt =  corr(pos,thetas,params,dr)
        Crt[:,token] = crt
        token = min(token+1,length(times))
    end
end

scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=500/L,size=(512,512),xlims=(0,L),ylims=(0,L))
cg_theta = cg(pos,thetas,params,dr=3R0)
heatmap(cg_theta',color=cols,clims=(0,2pi))


rr = collect(0:dr:round(Int,L/2))
p=plot()
    for t in 1:length(times)
        plot!(rr,Crt[:,t])
    end
    p

## Time Complexity of the code
#= Note, if R0 is very small, the number of cells in
the cell list method explodes (rounding issues ?), thus critically slowing down the execution.
Resolved by fixing R0 = 1 instead of fixing L = 1 =#
Ns = round.(Int,2.0 .^ collect(4:15))
α = 1
    T     = 0.4
    v0    = 0.2
    σ = 0.8
    R = 10
    tmax = 2

    runtimes = zeros(length(Ns),R)

for n in eachindex(Ns)
    N = Ns[n]
    R0 = 1
    L = sqrt.(π*N/α)*R0 # interaction radius
    rho = N/L^2
    # R0 = 1/10 # interaction radius
    println("N = $N started, L = $L.")

    params = Any[α,T,v0,N,L,rho,R0,σ]
    dt = determine_dt(T,σ,v0,R0,L)

    Threads.@threads for r in 1:R
        pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
        t = 0.0 ; token = 2
        z = @elapsed while t < tmax
            t += dt
            pos,thetas = update(pos,vel_angles,thetas,omegas,params,dt)
        end
        runtimes[n,r] = z
    end # realisations
end # scanned parameters

runtimesavg = mean(runtimes,dims=2)
plot(Ns,runtimesavg[:,1],axis=:log)
plot!(Ns,Ns .^ 2 / 1E5,axis=:log,c=:black)
plot!(Ns,Ns .^ 1 / 1E2,axis=:log,c=:black)


# Phase Diagram
alphas = 1
    Ts     = 0:0.05:0.5
    v0s    = 0.2
    sigmas = 0:0.05:0.3
    R = 30

    tmax = 20 ; times = 0:0.1:tmax

    P = zeros(length(alphas),length(Ts),length(v0s),length(sigmas),length(times),R)
    m = 0 ; M = length(alphas)*length(Ts)*length(v0s)*length(sigmas)


z = @elapsed for i in eachindex(alphas) , j in eachindex(Ts) , k in eachindex(v0s) , q in eachindex(sigmas)
    α = alphas[i] ; T = Ts[j] ; v0 = v0s[k] ; σ = sigmas[q]
    m += 1 ; println("Simulation $m/$M : α = $(α) , T = $(T) , v0 = $(v0) . ")

    N = Int(1E3)
    L = 1
    rho = N/L^2
    R0 = sqrt.(α/π/N)*L # interaction radius

    params = Any[α,T,v0,N,L,rho,R0,σ]
    dt = determine_dt(T,σ,v0,R0,L)

    Threads.@threads for r in 1:R
        pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
        P[i,j,k,q,1,r] = polarOP(thetas)[1]
        t = 0.0 ; token = 2
        while t < tmax
            t += dt
            pos,thetas = update(pos,vel_angles,thetas,omegas,params,dt)
            if t > times[token]
                op =  polarOP(thetas)[1]
                P[i,j,k,q,token,r] = op
                token = min(token+1,length(times))
                # if op > .95
                #     P[i,j,k,token:end,r] .= 1
                #     println("Early stopping.")
                #     break
                # end
            end
        end
    end # realisations
end # scanned parameters
prinz(z)

# JLD2.jldsave("data\\firstdata_polarOP.jld2";alphas,Ts,v0s,P,R,tmax,times,dt)
Pavg = mean(P,dims=6)
heatmap(Ts,sigmas,Pavg[1,:,1,:,end,1],c=cgrad([:red,:green]),clims=(0,1))



# anim = @animate for j in eachindex(Ts)
#     plot(xlabel="α",ylabel=L"v_0",title="T=$(Ts[j])")
#     heatmap!(alphas,v0s,Pavg[:,j,:,;,end,1]',c=cgrad([:red,:green]),clims=(0,1))
# end
# mp4(anim,"films/phase_diagram_sliding_T.mp4",fps=1)

# @unpack alphas,Ts,v0s,P,R,tmax,times,dt = JLD2.load("data\\firstdata_polarOP.jld2")
# Pavg = mean(P,dims=5)
# p=plot(xlabel="t",ylabel="OP(t)",legend=:outerright,size=(500,400),ylims=(0,1))
#     for i in eachindex(alphas)
#         plot!(times,Pavg[i,11,1,:,1],label=string(alphas[i]))
#     end
#     p
&


## Analysis of simulations from cluster
@unpack alphas,Ts,v0s,tmax,times,P,runtime,N,L,R0,init_thetas,init_alphas,init_omegas,Var = load("data/phasediag_Var0.jld")
v0s
Ts
P             #  (alphas),(Ts),(v0s),(times)
Pavg = mean(P,dims=5)

# Movies
anim_varyingT = @animate for i in eachindex(Ts)
    heatmap(v0s,π*alphas,Pavg[:,i,:,end-1,1],c=cgrad([:red,:orange, :green]),size=(500,400),clims=(0,1),xlabel="v",ylabel="α",colorbar_title="P",title="T=$(Ts[i])")
end
anim_varyingv0 = @animate for i in eachindex(v0s)
    heatmap(Ts,π*alphas,Pavg[:,:,i,end-1,1],c=cgrad([:red,:orange, :green]),size=(500,400),clims=(0,1),xlabel="T",ylabel="α",colorbar_title="P",title="v_0=$(v0s[i])")
end
anim_varyingalpha = @animate for i in eachindex(alphas)
    heatmap(v0s,Ts,Pavg[i,:,:,end-1,1],c=cgrad([:red,:orange, :green]),size=(500,400),clims=(0,1),xlabel="v",ylabel="T",colorbar_title="P",title="α=$(round(π*alphas[i],digits=1))")
end

mp4(anim_varyingT,"films/phase_diagram_varyingT_coarse.mp4",fps=2)
mp4(anim_varyingv0,"films/phase_diagram_varyingv0_coarse.mp4",fps=1)
mp4(anim_varyingalpha,"films/phase_diagram_varyingalpha_coarse.mp4",fps=3)

## First Movies of Ballistic Motion
gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

# Controlled parameters
    α = 2 # average number of influenced neighbours
    T = 0.5 # temperature for angle diffusion
    v0 = 0.1 # norm of individual velocities
    σ = 1.5

    # Other parameters
    N = Int(1E4)
    L = 1
    rho = N/L^2
    R0 = sqrt(α/pi/N) # interaction radius
    println("R0 = $R0")

    params = Any[α,T,v0,N,L,rho,R0,σ] # any to avoid N being interpreted as a Float
    dt = 1E-2

# pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
#     scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=300/L,size=(512,512),xlims=(0,L),ylims=(0,L))

every = 10 ; tmax = 50
z = @elapsed anim = movies(params,every,tmax,dt) ; prinz(z)
mp4(anim,"films/N1E4_sigma_$(σ)_T_$(T)_alpha_$(α).mp4")
