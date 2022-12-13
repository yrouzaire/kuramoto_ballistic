cd("D:/Documents/Research/projects/kuramoto_ballistic")
    using JLD2,StatsBase,Distributions,LinearAlgebra,Parameters,Random,BenchmarkTools,Hungarian
    include("../methods.jl")
    global const R0 = 1
    using Plots,ColorSchemes,LaTeXStrings
    pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    plot()
&

filename = "data/nature_phase_impact_rho_N1E4_tmax1E4.jld2"
@load filename Ps Cs ns runtimes Ts Ns v0s rhos sigmas times_log tmax comments R

length(v0s)
length(sigmas)
histogram(runtimes/3600/24*2/3*10*6/9*5,bins=20)
histogram(runtimes/3600/24,bins=20)
Ps_avg = nanmean(Ps,8) # N rho T v0 sigma t init R
Ps_std = nanstd(Ps,8) # N rho T v0 sigma t init R
ns_avg = nanmean(ns,8) # N rho T v0 sigma t init R
real_completed = []
for r in 1:R
    try
        Cs[r]
        push!(real_completed,r)
    catch ;
    end
end
Cs_avg = mean(Cs[real_completed])
xis_avg = NaN*zeros(length(rhos),length(v0s),length(sigmas),length(times_log))
for p in each(rhos), i in each(v0s) , j in each(sigmas) , t in each(times_log)
    rr = 0:0.5:round(sqrt(Ns)/rhos[p])
    xis_avg[p,i,j,t] = corr_length(Cs_avg[1,p,1,i,j,1,t],rr)
end

## Recall the phase space
v0s
heatmap(v0s[2:end],sigmas,Ps_avg[1,end,1,2:end,:,1,end,1]',xaxis=:log,c=cgrad([:red,:orange,:green]),clims=(0,1))

## Order or disorder at XY model (v = 0, \sigma = 0)?
p=plot(xlabel="ρ",ylabel="P(v = σ = 0)")
    plot!(rhos,Ps_avg[1,:,1,1,1,1,end,1],m=true,line=true)

p=plot(xlabel="r",ylabel=L"C(r,t_{max})",axis=:log,legend=:bottomleft)#,ylims=(1E-2,1.1))
    ind_sig = 1
    ind_vel = 1
    for i in 1:2:length(rhos)
        rr = 0.5:0.5:round(sqrt(Ns)/rhos[i]/2,RoundUp)
        plot!(rr,remove_negative(Cs_avg[1,i,1,ind_vel,ind_sig,1,end][2:length(rr)+1]),label="ρ = $(rhos[i])",rib=0)
    end
    p

savefig("figures\\nature_phase_diagram/D_loglin_C_several_rhos_v$(v0s[ind_vel])_σ$(sigmas[ind_sig]).png")
&

p=plot(xlabel="t",ylabel=L"ξ(t)",axis=:log)
    ind_sig = 2
    ind_vel = 3
    for i in 1:1:length(rhos)
        plot!(times_log,xis_avg[i,ind_vel,ind_sig,:])
    end
    plot!(x->3sqrt(x/log(100x)),c=:black)
    p

times_log[end]
rhos
p=plot(xlabel="r",ylabel=L"C(r,t)",axis=:log)
    ind_sig = 1
    ind_vel = 1
    ind_rho = 4
    rr = 0.5:0.5:round(sqrt(Ns)/rhos[ind_rho]/2,RoundUp)
    for t in 25:length(times_log)
        plot!(rr,remove_negative(Cs_avg[1,ind_rho,1,ind_vel,ind_sig,1,t][2:length(rr)+1]))
    end
    p
savefig("figures\\nature_phase_diagram/D_loglin_C_several_t_v$(v0s[ind_vel])_σ$(sigmas[ind_sig])_rho2.png")


## Correlation function at final time for different v0
ind_sig = 1
    p=plot(xlabel="r",ylabel="C(r,∞)",legend=:outerright,size=(550,400),yaxis=:log,title="σ = $(sigmas[ind_sig])",ylims=(1E-2,1.2))
    for i in 1:length(v0s)
        plot!(0.5:0.5:50,remove_negative(Cs_avg[2:end,1,1,1,i,ind_sig,1,end,1]),label="v = $(v0s[i])",rib=0)
    end
    savefig("figures\\transition_SR-LR_sigma$(sigmas[ind_sig]).png")
    p

## Correlation function at final time for different sigmas, at given v0
ind_vel = 1
    ind_rho =
    p=plot(xlabel="r",ylabel="C(r,∞)",legend=:outerright,size=(550,400),yaxis=:log,title="v = $(v0s[ind_vel])",ylims=(1E-3,1.2))
    for i in 1:length(sigmas)
        plot!(0.5:0.5:round(Int,sqrt(Ns/rhos[ind_rho])/2,RoundDown),remove_negative(Cs_avg[1,ind_rho,1,ind_vel,i,1,end,1][3:end]),label="σ = $(sigmas[i])",rib=0)
    end
    # savefig("figures\\transition_SR-LR_sigma$(sigmas[ind_sig]).png")
    p

## Correlation length over time for different v0
ind_sig = 1
    p=plot(xlabel="t",ylabel="ξ(t)",legend=:outerright,size=(550,400),axis=:log,title="σ = $(sigmas[ind_sig])")
    for i in 1:length(v0s)
        plot!(times_log,xis_avg[i,ind_sig,:],label="v = $(v0s[i])",rib=0)
    end
    plot!(x->sqrt(x),c=:black)
    savefig("figures\\xit_SR-LR_sigma$(sigmas[ind_sig]).png")
    p

## Correlation functions over time
ind_v0 = 8
    ind_sigma = 2
    p=plot(xlabel="r",ylabel="C(r,t)",legend=false,axis=:log,ylims=(1E-3,1.2))
    for tt in each(times_log)
        plot!(0.5:0.5:50,remove_negative(Cs_avg[2:end,1,1,1,ind_v0,ind_sigma,1,tt,1]))
    end
    p

## Correlation functions over phase space
styles = [:solid,:dash]
    p=plot(xlabel="r",ylabel="C(r,t)",legend=false,axis=:log,ylims=(1E-3,1.2))
    for ind_v0 in each(v0s)
        for ind_sigma in each(sigmas)
            plot!(0.5:0.5:50,remove_negative(Cs_avg[2:end,1,1,1,ind_v0,ind_sigma,1,end,1]),c=ind_v0,line=styles[ind_sigma])
        end
    end
    p

## Correlation length over time
p=plot(xlabel="t",ylabel="ξ(t)",legend=:outerright,axis=:log,size=(600,400))#,ylims=(1E-3,1.2))
    for ind_v0 in 1:length(v0s)
        plot!([NaN,NaN],rib=0,label="v = $(v0s[ind_v0])")
        for ind_sigma in each(sigmas)
            plot!(times_log,xis_avg[ind_v0,ind_sigma,:],c=ind_v0,line=styles[ind_sigma])
        end
    end
    plot!(x->sqrt(x/log(x)),c=:black)
    plot!(x->sqrt(x),c=:grey)
    hline!([exp(-1)*100],c=:grey,line=:dot)

## XY Collapse Scaling ?
rr = collect(0.5:0.5:50)
ind_v0 = 1
    ind_sigma = 1
    p=plot(xlabel="r/ξ(t)",ylabel=L"r^{-η}\,C(r/ξ,t)",legend=false,axis=:log,ylims=(1E-3,1.2))
    for tt in 15:2:length(times_log)
        # plot!(rr ./xis_avg[ind_v0,ind_sigma,tt],rr .^(-0.1/2pi) .* remove_negative(Cs_avg[2:end,1,1,1,ind_v0,ind_sigma,1,tt,1]))
        plot!(rr ,remove_negative(Cs_avg[2:end,1,1,1,ind_v0,ind_sigma,1,tt,1]))
    end
    p

## XY Relation Length and number of vortices ?
rr = collect(0.5:0.5:50)
    p=plot(xlabel="t",ylabel=L"ξ.\sqrt{n}/L",legend=false,xaxis=:log)
    for ind_v0 in 1:length(v0s)
        # plot!([NaN,NaN],rib=0,label="v = $(v0s[ind_v0])")
        for ind_sigma in 1#each(sigmas)
            # plot!(times_log,.(v0s[ind_v0]*times_log).^(-0.01) .* xis_avg[ind_v0,ind_sigma,:] .* sqrt.(ns_avg[1,1,1,ind_v0,ind_sigma,1,:,1])/sqrt(Ns[1]))
            plot!(times_log,xis_avg[ind_v0,ind_sigma,:] .* sqrt.(ns_avg[1,1,1,ind_v0,ind_sigma,1,:,1])/sqrt(Ns[1]))
        end
    end
    p

## Visual aspect of realisations
@load filename thetas_saveds pos_saveds
thetas_saved = thetas_saveds
pos_saved = pos_saveds
ind_sig = 2
    ind_vel = 1
    ind_rho = 11
    r = 5
    L = round(Int,sqrt(Ns/rhos[ind_rho]))
    plot(pos_saved[:,:,ind_rho,1,ind_vel,ind_sig,1,r],thetas_saved[:,ind_rho,1,ind_vel,ind_sig,1,r],Ns,L,L,particles=true)
