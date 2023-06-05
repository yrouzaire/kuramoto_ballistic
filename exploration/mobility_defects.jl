cd("D:/Documents/Research/projects/kuramoto_ballistic")
    using JLD2,StatsBase,Distributions,LinearAlgebra,Parameters,Random,BenchmarkTools,Hungarian
    include("../methods.jl")
    global const R0 = 1
    using Plots,ColorSchemes,LaTeXStrings
    pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    plot()
&
20*20*10000*128/8/1E9

## ------------------------------ Cluster data analysis ------------------------------ ##
filename = "data/mobility_defects_sigma_v0.jld2"
@load filename sigmas v0s Ts all_xy_pos all_xy_neg all_rr all_times_collision R_per_core Rtot params_init Ntarget R0 q init_theta init_pos aspect_ratio times tmax comments rhoc runtimes
hrun(runtimes) 
hrun(runtimes*2*10/4/3*5) 
Rtot
tmax
times
# times = collect(0:1:tmax)

all_xy_pos[1,1,1,1] # v, sigma, T, r 
L = round(Int,sqrt(Ntarget))

# Compute fluctuations in the y-direction
xs_pos = NaN*zeros(length(v0s), length(sigmas), length(Ts),length(times),Rtot)
xs_neg = NaN*zeros(length(v0s), length(sigmas), length(Ts),length(times),Rtot)
ys_pos = NaN*zeros(length(v0s), length(sigmas), length(Ts),length(times),Rtot)
ys_neg = NaN*zeros(length(v0s), length(sigmas), length(Ts),length(times),Rtot)
for i in each(v0s)
    for j in each(sigmas)
        for k in each(Ts)
            for r in 1:Rtot 
                for tt in 1:length(all_xy_pos[i,j,k,r])
                    xs_pos[i,j,k,tt,r] = all_xy_pos[i,j,k,r][tt][1]
                    xs_neg[i,j,k,tt,r] = all_xy_neg[i,j,k,r][tt][1]
                    ys_pos[i,j,k,tt,r] = all_xy_pos[i,j,k,r][tt][2]
                    ys_neg[i,j,k,tt,r] = all_xy_neg[i,j,k,r][tt][2]
                end
            end
        end
    end
end

## Analysis rr(t*)
all_rr_reverse = NaN*zeros(length(v0s),length(sigmas),length(Ts),length(times),Rtot)
all_rr_ = NaN*zeros(length(v0s),length(sigmas),length(Ts),length(times),Rtot)
for i in each(v0s)
    for j in each(sigmas)
        for k in each(Ts)
            for r in 1:Rtot 
                data = all_rr[i,j,k,r]
                ll = length(data)
                all_rr_reverse[i,j,k,1:ll,r] = reverse(data)
                all_rr_[i,j,k,1:ll,r] = (data)
                all_rr_[i,j,k,1+ll:end,r] .= 0
            end
        end
    end
end
all_rr_reverse_avg = nanmean(all_rr_reverse,5)[:,:,:,:,1]
all_rr_avg = nanmean(all_rr_,5)[:,:,:,:,1]

## Ys +
p=plot(xlabel=L"\sqrt{t}",ylabel="std(y+)")
for i in each(v0s)
    for j in 1#each(sigmas)
        for k in each(Ts)
            plot!(sqrt.(times),[std(ys_pos[i,j,k,t,:]) for t in 1:length(times)])
        end
    end
end
# xlims!(-1,maximum(filter(isnan,vec(ys_pos))))
p

## Ys -
p=plot(xlabel=L"\sqrt{t}",ylabel="std(y-)")
for i in each(v0s)
    for j in 1#each(sigmas)
        for k in each(Ts)
            plot!(sqrt.(times),[std(ys_neg[i,j,k,t,:]) for t in 1:length(times)])
        end
    end
end
p

## Ys +/-
p=plot(xlabel=L"\sqrt{t}",ylabel="std(y-)")
for i in each(v0s)
    for j in 1#each(sigmas)
        for k in each(Ts)
            plot!(sqrt.(times),[std(vcat(ys_pos[i,j,k,t,:],ys_neg[i,j,k,t,:])) for t in 1:length(times)])
        end
    end
end
p

## Fit the slope of variance in 
first_index = 7 # robust to any change
slopes = NaN*zeros(length(v0s), length(sigmas), length(Ts))

for i in each(v0s)
    for j in each(sigmas)
        for k in each(Ts)
            data = smooth([var(vcat(ys_pos[i,j,k,t,:],ys_neg[i,j,k,t,:])) for t in 1:length(times)],0)
            last_index = findfirst(isnan,data)
            if isnothing(last_index) last_index = length(data) end
            slopes[i,j,k] = (data[last_index-1] - data[first_index])/((times[last_index-1]) - (times[first_index]))
        end
    end
end

plot(v0s, slopes[:,1,1],m=true)
plot!(v0s, 0.08v0s .+ 0.03,c=:black)

## times_collision
for i in 1:length(v0s)
    p=histogram(all_times_collision[i,1,1,:],bins=50)
    title!("i = $i, v0 = $(v0s[i])")
    display(p)
end

## R(t) 
p=plot(xlabel=L"t", ylabel=L"R(t)", uaxis=:log,ylims=(0,30))
for i in each(v0s)
    for j in each(sigmas)
        for k in each(Ts)
            plot!(times[3:end],all_rr_avg[i,j,k,3:end])
        end
    end
end
p

# Add the theoretical prediction for the fastest curve
mean_annihilation_time = mean(all_times_collision[end,1,1,:])
# plot!(0:1:mean_annihilation_time-10,x->3.8exp(0.5*(1+lambertw((mean_annihilation_time-x)*2/exp(1)))),c=:black)
plot!(0:1:mean_annihilation_time,x->28*sqrt((-x+mean_annihilation_time)/mean_annihilation_time),c=:black)

## R(t*) 
using LambertW
p=plot(xlabel=L"t^*", ylabel=L"R(t^*)", axis=:log)
for i in each(v0s)
    for j in 1#each(sigmas)
        for k in each(Ts)
            # plot!(times[2:end],all_rr_reverse_avg[i,j,k,2:end])
            plot!((v0s[i])*times[2:end],all_rr_reverse_avg[i,j,k,2:end])
            # plot!(v0s[i]*times[2:end],(0.05+v0s[i])^-0.5*all_rr_reverse_avg[i,j,k,2:end])
        end
    end
end
mu = 1
D = 2
plot!(times[2:end],x->exp(0.5*lambertw(2π*x/mu)),c=:black)
plot!(times[2:end],x->D*sqrt(x),c=:black, line=:dash)
p

