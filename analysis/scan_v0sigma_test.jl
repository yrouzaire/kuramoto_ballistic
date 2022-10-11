cd("D:/Documents/Research/projects/kuramoto_ballistic")
    using JLD2,StatsBase,Distributions,LinearAlgebra,Parameters,Random,BenchmarkTools,Hungarian
    include("../methods.jl")
    global const R0 = 1
    using Plots,ColorSchemes,LaTeXStrings
    pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    plot()
&

filename = "data/bigscan_v0_sigma_N1E3_tmax1E4.jld2"
filename = "data/scan_v0_sigma_rho1_N1E4_tmax1E4.jld2"
filename = "data/scan_for_discussion_Parisa_v0_sigma_rho_N1E3_tmax2E3.jld2"
@load filename Ps ns runtimes Ts Ns v0s rhos sigmas times_log tmax comments R

histogram(runtimes/3600/24,bins=20)
Ps_avg = nanmean(Ps,8) # N rho T v0 sigma t init R
Ps_std = nanstd(Ps,8) # N rho T v0 sigma t init R
ns_avg = nanmean(ns,8) # N rho T v0 sigma t init R
Cs_avg = nanmean(Cs,9)


# Heatmap of P at final time
p=plot(xlabel=L"v_0",ylabel=L"\sigma")
    heatmap!(v0s[2:end],sigmas,Ps_avg[1,1,1,2:end,:,1,end,1]',xaxis=:log,c=cgrad([:red,:orange,:green]),clims=(0,1))
    # plot!(x->0.25(x-0.1)^0.5,xlims=(v0s[2],v0s[end]),ylims=(0,sigmas[end]),c=:black)
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
p=plot(xlabel=L"v_0",ylabel=L"\sigma")
    heatmap!(v0s[2:end],sigmas,log10.(ns_avg[1,2,1,2:end,:,1,end,1]' .+1),xaxis=:log,c=cgrad([:green,:orange,:red]))
    # plot!(x->0.25(x-0.1)^0.5,xlims=(v0s[2],v0s[end]),ylims=(0,sigmas[end]),c=:black)

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

p=plot(legend=false)
    for indv in each(v0s)
        for inds in 10# each(sigmas)
            if ns_avg[1,1,1,indv,inds,1,end,1] > 2
                plot!(times_log,Ps_avg[1,1,1,indv,inds,1,:,1].*(ns_avg[1,1,1,indv,inds,1,:,1]),uaxis=:log,label="v=$(v0s[indv])")
            end
        end
    end
    p

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
lss = [:solid,:dash,:dot]
p=plot(legend=false)#:bottomright)
    for i in each(v0s)
        for j in each(sigmas)
            plot!(times_log,Ps_avg[1,1,1,i,j,1,:,1],axis=:log,rib=0Ps_std[1,1,1,i,j,1,:,1],
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
