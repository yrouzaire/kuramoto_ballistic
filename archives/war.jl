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
    determine_dt(T,σ,v0*L,R0,L)

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
    T = 0.1 # temperature for angle diffusion
    v0 = 5 # norm of individual velocities
    σ = 0

    # Other parameters
    N = Int(1E3)
    L = 1
    rho = N/L^2
    R0 = 1
    L = round(Int,sqrt(π*N/α)*R0) # interaction radius
    println("L = $L")

    params = Any[α,T,v0,N,L,rho,R0,σ] # any to avoid N being interpreted as a Float
    dt = 1E-2

# pos,vel_angles,thetas,omegas = initialisation(N,L,σ)
#     scatter(pos[1,:],pos[2,:],marker_z = mod.(thetas,2pi),color=cols,clims=(0,2pi),ms=300/L,size=(512,512),xlims=(0,L),ylims=(0,L))

every = 10 ; tmax = 30
z = @elapsed anim = movies(params,every,tmax,dt) ; prinz(z)
mp4(anim,"films/demos_R0fixed/film_sigma_$(σ)_T_$(T)_vo_$(v0)_alpha_$(α).mp4",fps=30)
