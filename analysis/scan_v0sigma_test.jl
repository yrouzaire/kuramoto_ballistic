cd("D:/Documents/Research/projects/kuramoto_ballistic")
    using JLD2,StatsBase,Distributions,LinearAlgebra,Parameters,Random,BenchmarkTools,Hungarian
    include("../methods.jl")
    global const R0 = 1
    using Plots,ColorSchemes,LaTeXStrings
    pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    plot()
&

filename = "data/scan_for_discussion_Parisa_v0_sigma_rho_N1E3_tmax2E3.jld2"
filename = "data/XIN_during_coarsening_v0_sigma_N1E4_tmax2E3.jld2"
filename = "data/scan_v0_sigma_rho1_N1E4_tmax1E4.jld2"
@load filename Ps Cs ns runtimes Ts Ns v0s rhos sigmas times_log tmax comments R

histogram(runtimes/3600/24,bins=20)
Ps_avg = nanmean(Ps,7) # N rho T v0 sigma t init R
Ps_std = nanstd(Ps,7) # N rho T v0 sigma t init R
ns_avg = nanmean(ns,7) # N rho T v0 sigma t init R
Cs_avg = nanmean(Cs,9)
xis_avg = zeros(length(v0s),length(sigmas),length(times_log))
for i in each(v0s) , j in each(sigmas) , t in each(times_log)
    xis_avg[i,j,t] = corr_length(Cs_avg[:,1,1,1,i,j,1,t,1],0:0.5:50)
end


plot(times_log,Ps_avg[1,1,1,end,1,:,1],axis=:log,rib=0)


# Heatmap of P at final time
p=plot(xlabel=L"v_0",ylabel=L"\sigma",colorbartitle="P",size=(500,400))
    heatmap!(v0s[2:end],sigmas,Ps_avg[1,2,1,2:end,:,1,end,1]',xaxis=:log,c=cgrad([:red,:orange,:green]),clims=(0,1))
    # plot!(x->0.25(x-0.1)^0.5,xlims=(v0s[2],v0s[end]),ylims=(0,sigmas[end]),c=:black)
savefig("figures/heatmap_sigmav0_rho2.png")

# Animation of Heatmap over time
anim = @animate for t in each(times_log)
    p_low=plot(xlabel=L"v_0",ylabel=L"\sigma",title="t = $(times_log[t])")
        heatmap!(v0s[2:end],sigmas,Ps_avg[1,2,1,2:end,:,1,t,1]',xaxis=:log,c=cgrad([:red,:orange,:green]),size=(512,512),clims=(0,1))
    p_high=plot(xlabel=L"v_0",ylabel=L"\sigma",title="t = $(times_log[t])")
        heatmap!(v0s[2:end],sigmas,Ps_avg[1,2,1,2:end,:,2,t,1]',xaxis=:log,c=cgrad([:red,:orange,:green]),size=(512,512),clims=(0,1))
    plot(p_low,p_high,size=(1024,512))
end
mp4(anim,"figures/animation_phasespace_v0σ_N1E3_rho2.mp4",fps=10)

# Heatmap of N at final time
p=plot(xlabel=L"v_0",ylabel=L"\sigma",colorbartitle=L"\log_{10}(n+1)",size=(500,400))
    heatmap!(v0s[2:end],sigmas,log10.(ns_avg[1,2,1,2:end,:,1,end,1]' .+1),xaxis=:log,c=cgrad([:green,:orange,:red]))
    # plot!(x->0.25(x-0.1)^0.5,xlims=(v0s[2],v0s[end]),ylims=(0,sigmas[end]),c=:black)
savefig("figures/heatmapN_sigmav0_rho2.png")

# Animation of Heatmap over time
anim = @animate for t in each(times_log)
    p_low=plot(xlabel=L"v_0",ylabel=L"\sigma",title="t = $(times_log[t])")
    heatmap!(v0s[2:end],sigmas,log10.(ns_avg[1,1,1,2:end,:,1,t,1]' .+1),xaxis=:log,c=cgrad([:green,:orange,:red]),clims=(0,log10(maximum(ns_avg))))
    p_high=plot(xlabel=L"v_0",ylabel=L"\sigma",title="t = $(times_log[t])")
    heatmap!(v0s[2:end],sigmas,log10.(ns_avg[1,1,1,2:end,:,2,t,1]' .+1),xaxis=:log,c=cgrad([:green,:orange,:red]),clims=(0,log10(maximum(ns_avg))))
    plot(p_low,p_high,size=(1024,512))
end
mp4(anim,"figures/animation_phasespaceN_v0σ_N1E3.mp4",fps=10)

# Relation between P(t) and n(t) ?
ind = 11
    plot(v0s[2:end],Ps_avg[1,1,1,2:end,ind,1,end,1].*(ns_avg[1,1,1,2:end,ind,1,end,1].+1),uaxis=:log)

p=plot(xlabel="t",ylabel="P(t) / n(t)",legend=false)
    for indv in 7# each(v0s)
        for inds in each(sigmas)
            if ns_avg[1,1,1,indv,inds,1,end,1] > 2
                plot!(times_log,Ps_avg[1,1,1,indv,inds,1,:,1].*(ns_avg[1,1,1,indv,inds,1,:,1]),uaxis=:log,label="v=$(v0s[indv])")
            end
        end
    end
    p
# savefig("figures/PsurN_diffv0_onesigma.png")
# savefig("figures/PsurN_diffsigma_onev0.png")

# Determination of critical sigma as a function of v0
sigmac = zeros(length(v0s))
seuil = 0.2
tmp = findfirst(x->x>seuil,Ps_avg[1,1,1,end,:,1,end,1])

for i in each(v0s)
    tmp = findfirst(x->x<seuil,Ps_avg[1,1,1,i,:,1,end,1])
    if isnothing(tmp)
        sigmac[i] = sigmas[end]
    else
        sigmac[i] = sigmas[tmp]
    end
end
plot(v0s[2:end],sigmac[2:end],xaxis=:log)
    plot!(x->0.5sqrt(x))

# Over time
lss = [:solid,:dash,:dot,:solid]
p=plot(legend=false)#:bottomright)
    for i in each(v0s)
        for j in each(sigmas)
            plot!(times_log,Ps_avg[1,1,1,i,j,1,:,1],axis=:log,rib=0,
            label="v = $(v0s[i]), σ = $(sigmas[j])",c=i,line=lss[j])
        end
    end
    plot!(x->1E-2x^0.25,c=:black)
    # plot!(x->0.8)
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
    for i in 3 #each(v0s)
        for j in each(sigmas)
            plot!(times_log,ns_avg[1,1,1,i,j,:,1].+1,axis=:log,label="v = $(v0s[i]), σ = $(sigmas[j])",c=i,line=lss[j])
        end
    end
    plot!(times_log[7:end],5E2*log.(times_log[7:end]) ./ times_log[7:end],c=:black)
    # plot!(x->1E2log(x)/x)

    p

## No hysteresis : lowtemp and hightemp init lead to similar results.
Ps_avg = nanmean(Ps,8) # N rho T v0 sigma t init R

indss = [(1,1),(4,10),(10,4),(10,7)]
    p=plot(xlabel="t",ylabel="P(t)",legend=:bottomright)
    for i in each(indss)
        indv,inds = indss[i]
        plot!(times_log,Ps_avg[1,1,1,indv,inds,1,:,1],line=:solid,c=i,axis=:log)
        plot!(times_log,Ps_avg[1,1,1,indv,inds,2,:,1],line=:dash,c=i,axis=:log)
    end
    plot!([NaN,NaN],label="Ordered",line=:solid,c=:grey)
    plot!([NaN,NaN],label="Disordered",line=:dash,c=:grey)
    plot!(times_log[5:end-5],4E-2sqrt.(times_log[5:end-5]),c=:black,label="√t")
    p
savefig("figures/no_difference_lowtemp_hightemp_loglog.png")

## Does xi \sqrt(n) ~ cst hold ?
# D'un côté, pas hyper clair mais oui, ca a l'air constant (on arrive au terme du coarsening)
# De l'autre, non pas constant ca augmente lentement mais surement...
p=plot()
    for i in 2# each(v0s)
        for j in 1#each(sigmas)
            plot!(times_log,xis_avg[i,j,:],m=true,axis=:log)
            plot!(times_log,sqrt.(ns_avg[1,1,1,i,j,1,:,1]),m=true,axis=:log)
            plot!(times_log,xis_avg[i,j,:] .* sqrt.(ns_avg[1,1,1,i,j,1,:,1]),m=true,axis=:log)
        end
    end
    p
    plot!(x->x^0.5)
