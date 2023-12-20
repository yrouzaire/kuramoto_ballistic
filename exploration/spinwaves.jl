# cd("D:/Documents/Research/projects/kuramoto_ballistic")
    using JLD2,StatsBase,Distributions,LinearAlgebra,Parameters,Random,BenchmarkTools,Hungarian, JLD2
    include("../methods.jl")
    global const R0 = 1
    using Plots,ColorSchemes,LaTeXStrings
    pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    plot()
&
    
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
plot()

## ------------------------------ Cluster data analysis ------------------------------ ##
filename = "data/proba_spinwaves_scan_phase_space2.jld2"
filename = "data/proba_spinwaves_scan_phase_space_N4000_complement.jld2"
@load filename R_per_core Rtot R all_nb_detected_spinwave all_times_detected_spinwave all_Ps_detected_spinwave all_thetas_detected_spinwave all_pos_detected_spinwave sigmas v0s tmax times p_threshold init_pos init_theta Ntarget rho T aspect_ratio runtimes
all_nb_detected_spinwave_complement = all_nb_detected_spinwave
Rtot_complement = Rtot

filename = "data/proba_spinwaves_scan_phase_space_N4000.jld2"
# filename = "data/proba_spinwaves.jld2"
@load filename R_per_core Rtot R all_nb_detected_spinwave all_times_detected_spinwave all_Ps_detected_spinwave all_thetas_detected_spinwave all_pos_detected_spinwave sigmas v0s tmax times p_threshold init_pos init_theta Ntarget rho T aspect_ratio runtimes
proba_spinwave = all_nb_detected_spinwave / Rtot
hrun(runtimes)

Ntarget 
Rtot
tmax
sigmas
v0s

#
proba_spinwave = all_nb_detected_spinwave / Rtot 
proba_spinwave = 0all_nb_detected_spinwave / Rtot + all_nb_detected_spinwave_complement / Rtot_complement


all_nb_detected_spinwave
all_times_detected_spinwave
all_Ps_detected_spinwave
all_thetas_detected_spinwave
all_pos_detected_spinwave

## ---------------- Probas ---------------- ##
# p=plot(xlabel=L"v_0",ylabel=L"\mathbb{P}"*"(spinwave)", size=(500,400))
# for j in each(sigmas)
#     plot!(p,v0s,proba_spinwave[:,j],label="σ = $(sigmas[j])",c=j,rib=0,m=:circle)
# end
# p
# # savefig("figures/proba_spinwave.png")


## ---------------- Heatmap Probas ---------------- ##
## ---------------- Heatmap Probas ---------------- ##
## ---------------- Heatmap Probas ---------------- ##
## ---------------- Heatmap Probas ---------------- ##
sigmas
v0s
colss = cgrad([:black, :red, :orange, :gold])
plot(xaxis=:log, size=(470,400))
heatmap!(v0s, sigmas, 100*proba_spinwave',clims=(0,10), c=colss,colorbartitle=L"\mathbb{P}" * "(spinwave) [%]")
plot!(v0s, x -> 1 / 2 * max(0, x - (0.23)^2), c=:white, lw=0.8)
ylims!(-0.001,0.4)
xlims!(minimum(v0s), 1.18maximum(v0s))
xticks!([1E-2, 1E-1, 1E-0], [L"10^{-2}", L"10^{-1}", L"10^{0}"])
# xticks!([1E-2, 2E-2, 3E-2, 4E-2, 5E-2, 6E-2, 7E-2, 8E-2, 9E-2, 1E-1, 2E-1, 3E-1, 4E-1, 5E-1, 6E-1, 7E-1, 8E-1, 9E-1, 1E-0, 2, 3], 
#     [L"10^{-2}", "", "", "", "", "", "", "", "", L"10^{-1}", "", "", "", "", "", "", "", "", L"10^{0}", "", ""])
# critere : is_green_region = sigma < 1 / 2 * max(0, sqrt(v0) - 0.25) , at rho=1
annotate!((0.05, 0.98), text(L"\sigma", 17, :left, :top, :white))
annotate!((0.96, 0.03), text(L"v_0", 17, :right, :bottom, :white))
title!(L"ρ = 1, N = 10^4")
# savefig(pwd()*"/figures_paper/proba_spinwave.svg")


## ---------------- Detection times ---------------- ##
## ---------------- Detection times ---------------- ##
## ---------------- Detection times ---------------- ##
## ---------------- Detection times ---------------- ##

mean_time_detected_spinwave = NaN*zeros(length(v0s),length(sigmas))
for i in each(v0s), j in each(sigmas)
    is_empty = isempty(all_times_detected_spinwave[i,j])
    if !is_empty mean_time_detected_spinwave[i,j] = mean(all_times_detected_spinwave[i,j]) end
end
p = plot(xlabel=L"v_0", ylabel=L"\langle t \rangle" * "(spinwave)", size=(500, 400))
for j in each(sigmas)
    plot!(p, v0s, mean_time_detected_spinwave[:, j], label="σ = $(sigmas[j])", c=j, rib=0, m=:circle)
end
p


p = plot(xlabel=L"v_0", ylabel=L"\langle t_{creation} \rangle", size=(450, 400), yaxis=:log)
plot!(p, v0s, nanmean(mean_time_detected_spinwave, 2)[:, 1], rib=0, m=:circle)
plot!(v0s, x->8.5E2exp(-x/1.25), c=:black, lw=0.8, label=L"\exp(-v_0/1.25)")
p


## ---------------- Profiles : dependence on sigma and v0 ---------------- ##
## ---------------- Profiles : dependence on sigma and v0 ---------------- ##
## ---------------- Profiles : dependence on sigma and v0 ---------------- ##
## ---------------- Profiles : dependence on sigma and v0 ---------------- ##

# probas_renormalized = proba_spinwave 


# savefig("figures/proba_spinwave.png")


## ---------------- Plot all spinwaves ---------------- ##
## ---------------- Plot all spinwaves ---------------- ##
## ---------------- Plot all spinwaves ---------------- ##
## ---------------- Plot all spinwaves ---------------- ##

gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

for ind_v0 in each(v0s), ind_sig in each(sigmas)
    if length(all_times_detected_spinwave[ind_v0,ind_sig]) > 0
        reall = rand(1:length(all_times_detected_spinwave[ind_v0,ind_sig]))
        pp = Vector(undef,reall)
        for r in 1:reall
            pp[r] = scatter(all_pos_detected_spinwave[ind_v0, ind_sig][r],
                marker_z=mod.(all_thetas_detected_spinwave[ind_v0, ind_sig][r], 2π),
                c=cols, markersize=2, aspect_ratio=1, legend=false, axis=false)
        end
        # scatter(all_pos_detected_spinwave[ind_v0,ind_sig][reall],
        #     marker_z=mod.(all_thetas_detected_spinwave[ind_v0,ind_sig][reall],2π),
        #     c=cols,markersize=2,aspect_ratio=1,legend=false,
        #     title = "v,σ=$(v0s[ind_v0]),$(sigmas[ind_sig])  P=$(all_Ps_detected_spinwave[ind_v0,ind_sig][reall])")

        nb_cols = 3
        p=plot(pp..., layout=(round(Int, length(pp) / nb_cols, RoundUp), nb_cols), size=(400 * nb_cols, 400 * reall / nb_cols))
        display(p)
    end
end

## ---------------- Plot all spinwaves for given parameters ---------------- ##
## ---------------- Plot all spinwaves for given parameters ---------------- ##
## ---------------- Plot all spinwaves for given parameters ---------------- ##
## ---------------- Plot all spinwaves for given parameters ---------------- ##
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);

v0s
sigmas
ind_v0 = 20
ind_sig = 6 

reall = rand(1:length(all_times_detected_spinwave[ind_v0,ind_sig]))
pp = Vector(undef,reall)
for r in 1:reall
    pp[r] = scatter(all_pos_detected_spinwave[ind_v0, ind_sig][r],
        marker_z=mod.(all_thetas_detected_spinwave[ind_v0, ind_sig][r], 2π),
        c=cols, markersize=2, aspect_ratio=1, legend=false, axis=false)
end
# scatter(all_pos_detected_spinwave[ind_v0,ind_sig][reall],
#     marker_z=mod.(all_thetas_detected_spinwave[ind_v0,ind_sig][reall],2π),
#     c=cols,markersize=2,aspect_ratio=1,legend=false,
#     title = "v,σ=$(v0s[ind_v0]),$(sigmas[ind_sig])  P=$(all_Ps_detected_spinwave[ind_v0,ind_sig][reall])")

nb_cols = 3
p=plot(pp..., layout=(round(Int, length(pp) / nb_cols, RoundUp), nb_cols), size=(400 * nb_cols, 400 * reall / nb_cols))
display(p)






## --------------- Development P-Criteria && n == 0 to detect spin waves --------------- ##
## --------------- Development P-Criteria && n == 0 to detect spin waves --------------- ##
## --------------- Development P-Criteria && n == 0 to detect spin waves --------------- ##
## --------------- Development P-Criteria && n == 0 to detect spin waves --------------- ##
include("../parameters.jl")
# Fixed important params 
Ntarget = Int(1E4)
aspect_ratio = 1
rho = 1
T = 0.1
N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)
init_pos = "random"
init_theta = "hightemp"

# Scanned params
p_thresholds = [0.3,0.5]
sigmas = [0]
v0s = collect(0.5:0.5:5) 
v0s = [5]

R_per_core = 1 
tmax = 1E2 # max time
times = [tmax]

m = 0 
M = length(v0s)*length(sigmas)*R_per_core
nb_detected_spinwave = zeros(length(v0s),length(sigmas), length(p_thresholds))
plot_detected_spinwave = [Any[] for i in each(v0s), j in each(sigmas), k in each(p_thresholds)]

z = @elapsed for i in each(v0s)
    for j in each(sigmas)
        v0 = v0s[i]
        sigma = sigmas[j]

        params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
        params_phonons = Dict(:phonons => phonons, :phonon_amplitude => phonon_amplitude, :phonon_k => phonon_k, :phonon_omega => phonon_omega)
        param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
        :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
        :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init, :params_phonons => params_phonons)

        for r in 1:R_per_core
            global m += 1
            println("v0 = $v0, σ = $sigma, r = $r : Simu #$m/$M")
            system = System(param)

            for tt in eachindex(times)
                evolve!(system, times[tt])
                P = round(polarOP(system)[1],digits=2)
                n = number_defects(system)
                for k in each(p_thresholds)
                    if P < p_thresholds[k] && n == 0
                        println("v0 = $v0, σ = $sigma, P = $P < $(p_thresholds[k])")
                        p=plot_thetas(system,particles=true,nb_neighbours=true,defects=true,vertical=true,size=(512,512))
                        title!("v0 = $v0, σ = $sigma, P = $P, t = $(times[tt])")
                        push!(plot_detected_spinwave[i,j,k],p)
                        nb_detected_spinwave[i,j,k] += 1
                    end
                end
            end
        end
    end
end
prinz(z) 
# prinz(1000/40*10*3*z) 

# filename = "data/looking_for_spinwaves/proba_N$(Ntarget)_rho$(rho)_v$(v0s[1])_$(init_pos)_σ$(sigmas[1]).jld2"
# @save filename proba_spinwave nb_detected_spinwave plot_detected_spinwave sigmas v0s R tmax times p_threshold 

##
styless = [:solid,:dash,:dot]
p=plot(xlabel=L"v_0",ylabel=L"\mathbb{P}"*"(spinwave)")
for j in each(sigmas), k in each(p_thresholds)
    plot!(p,v0s,proba_spinwave[:,j,k],label="σ = $(sigmas[j]), p* = $(p_thresholds[k])",
    line=styless[k],c=j)
end
p

## Visualisation of the spinwave 
plot_detected_spinwave[1,1,1][1]


## --------------------- Old work ---------------------
## Generation of many realisations for a given parameter set
Ns = Int.(1E4)
    rhos = [1]
    Ts     = 0.1
    # v0s    = [0.0,0.01,0.016,0.025,0.04,0.063,0.1,0.16,0.25,0.4,0.63,1,1.6,2.5,4,6.3]
    # v0s    = [0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.5,1,2,3,5]
    v0s    = [5]
    sigmas = [0.2]
    inits  = ["disordered"]
    R = 50

    tmax = 50; #times_log = logspace(0.1,tmax,1)

    P  = zeros(length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),R)#,length(times_log))
    n  = zeros(length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),R)#,length(times_log))
    # @assert length(Ns) == 1
    # @assert length(rhos) == 1
    # dr = 0.5 ; C = zeros(length(0:dr:round(Int,sqrt(Ns)/2)),length(Ns),length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),length(times_log))

    pos_saved  = zeros(Float16,2,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),R)
    thetas_saved = zeros(Float16,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),R)
    psis_saved   = zeros(Float16,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),R)
    omegas_saved = zeros(Float16,Ns[1],length(rhos),length(Ts),length(v0s),length(sigmas),length(inits),R)

    m = 0 ; M = length(P)

z = @elapsed for a in eachindex(Ns)
    for i in eachindex(rhos) , j in eachindex(Ts) , k in eachindex(v0s) , q in eachindex(sigmas) , ini in eachindex(inits) , r in 1:R
    N = Ns[a] ; rho = rhos[i] ; T = Ts[j] ; v0 = v0s[k] ; σ = sigmas[q] ; init = inits[ini]

    L = round(Int,sqrt(N/rho))
    dt = determine_dt(T,σ,v0,N,rho)
    # params = Any[T,v0,N,L,rho,σ]
    # m += 1 ; println("Simu $m/$M")

    Random.seed!(r+100)
    pos,thetas,omegas,psis = initialisation(N,L,σ,[init])
    t = 0.0 ; token = 1
    while t < tmax
        t += dt
        ind_neighbours = get_list_neighbours(pos,N,L,L)
        pos,thetas = update(pos,thetas,omegas,psis,ind_neighbours,T,v0,N,L,dt)
        # if t ≥ times_log[token]
        #     P[a,i,j,k,q,ini,token]   = polarOP(thetas)[1]
        #     n[a,i,j,k,q,ini,token]   = number_defects(pos,thetas,N,L)
        #     # C[:,a,i,j,k,q,ini,token] = corr_fast(pos,thetas,N,L,dr)
        #     token = min(token+1,length(times_log))
        # end
    end
        P[a,i,j,k,q,ini,r]   = polarOP(thetas)[1]
        n[a,i,j,k,q,ini,r]   = number_defects(pos,thetas,N,L)


    pos_saved[:,:,i,j,k,q,ini,r]  = pos
    thetas_saved[:,i,j,k,q,ini,r] = thetas
    psis_saved[:,i,j,k,q,ini,r]   = psis
    omegas_saved[:,i,j,k,q,ini,r] = omegas

    p=plot(pos,thetas,N,L)
    title!("r=$r, σ=$σ, v=$v0")
    display(p)

   end
end
prinz(z)

# jldsave("data/looking_for_spinwaves/proba_N$(Ns[1])_rho$(rhos[1])_v$(v0s)_$(inits[1])_σ$(sigmas).jld2";R,P,n,pos_saved,thetas_saved,psis_saved,omegas_saved,Ns,rhos,sigmas,Ts,v0s,inits,tmax)
jldsave("data/looking_for_spinwaves/N$(Ns[1])_rho$(rhos[1])_v$(v0s[1])_$(inits[1])_σ$(sigmas[1]).jld2";R,P,n,pos_saved,thetas_saved,psis_saved,omegas_saved,Ns,rhos,sigmas,Ts,v0s,inits,tmax)

ro = 1
    sig= 1
    vo = 1
    plot(vec(P[1,ro,1,vo,sig,1,:]),ylims=(0,1))

## Make movies from promising seeds
promising_seeds = [112,133,122,127,109,101]
N = Int.(1E4) ; rho = 1 ; T = 0.1 ; v0 = 5 ; sigma = 0.2
L = round(Int,sqrt(N/rho));
Tmax = 5000
dt = 0.04
z = @elapsed for seeed in promising_seeds
    Random.seed!(seeed)
    pos,thetas,omegas,psis = initialisation(N,L,sigma,["disordered"])
    animation = @animate for i in 1:25:Tmax
        for j in 1:25 pos,thetas = update(pos,thetas,omegas,psis,T,v0,N,L,dt) end
        p=plot(pos,thetas,N,L)
    end
    mp4(animation,"films/looking_for_spinwaves/promising_seeds/$(seeed).mp4")
end
prinz(z)

## Investigation
filename = "data/looking_for_spinwaves/SPINWAVE_N10000_rho1_v5_disordered_σ0.2.jld2" ;
@unpack R,P,n,pos_saved,thetas_saved,psis_saved,omegas_saved,Ns,rhos,sigmas,Ts,v0s,inits,tmax = load(filename)
N = Ns[1] ; L = round(Int,sqrt(N/rhos[1])) ; σ = sigmas[1]
rea = 4
    pos    = Float32.(pos_saved[:,:,1,1,1,1,1,rea])
    thetas = Float32.(mod.(thetas_saved[:,1,1,1,1,1,rea] + 0randn(N),2π))
    psis   = Float32.(psis_saved[:,1,1,1,1,1,rea])
    omegas = Float32.(omegas_saved[:,1,1,1,1,1,rea])
    plot(pos,thetas,Ns,sqrt(Ns/rhos[1]),particles=false)

# savefig("figures/spinwaveR$(rea)_N10000_rho1_v5_σ0.2.png")
psis   = Float32.(2π*rand(N))
omegas = Float32.(2*randn(N))
for i in 1:500 pos,thetas = update(pos,thetas,omegas,psis,0.1,0.01,N,L,0.05) end
    plot(pos,thetas,N,L,particles=false)
savefig("figures/spinwaveR$(rea)_perturb1.png")
