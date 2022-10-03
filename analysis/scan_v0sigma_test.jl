cd("D:/Documents/Research/projects/kuramoto_ballistic")
    using JLD2,StatsBase,Distributions,LinearAlgebra,Parameters,Random,BenchmarkTools,Hungarian
    include("../methods.jl")
    global const R0 = 1
    using Plots,ColorSchemes,LaTeXStrings
    pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    plot()
&

filename = "data/scan_v0_sigma_N1E4_rho1_tmax100000.0.jld2"
@load filename Ps Cs ns runtimes Ts Ns v0s rhos sigmas times_log tmax comments R

histogram(runtimes/3600/24,bins=20)
Ps_avg = nanmean(Ps,7) # N rho T v0 sigma t R
Ps_std = nanstd(Ps,7) # N rho T v0 sigma t R
ns_avg = nanmean(ns,7) # N rho T v0 sigma t R
Cs_avg = nanmean(Cs,8)

# Over time
lss = [:solid,:dash]
p=plot(legend=false)#:bottomright)
    for i in each(v0s) , j in each(sigmas)
        plot!(times_log,Ps_avg[1,1,1,i,j,:,1],axis=:log,rib=0Ps_std[1,1,1,i,j,:,1],
        label="v = $(v0s[i]), σ = $(sigmas[j])",c=i,line=lss[j])
    end
    plot!(x->1E-2sqrt(x))
    plot!(x->0.8)
    p


p=plot(legend=false)#:bottomright)
    for i in 2:length(v0s) , j in 2# each(sigmas)
        plot!(remove_negative(Cs_avg[:,1,1,1,i,j,end,1]),axis=:log,rib=0Ps_std[1,1,1,i,j,:,1],
        label="v = $(v0s[i]), σ = $(sigmas[j])",c=i,line=lss[j])
    end
    # plot!(x->1E-2sqrt(x))
    # plot!(x->0.8)
    p


p=plot(legend=false)#:bottomright)
    for i in each(v0s) , j in each(sigmas)
        plot!(times_log,ns_avg[1,1,1,i,j,:,1].+1,axis=:log,label="v = $(v0s[i]), σ = $(sigmas[j])",c=i,line=lss[j])
    end
    p
