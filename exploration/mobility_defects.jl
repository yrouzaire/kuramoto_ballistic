cd("D:/Documents/Research/projects/kuramoto_ballistic")
    using JLD2,StatsBase,Distributions,LinearAlgebra,Parameters,Random,BenchmarkTools,Hungarian
    include("../methods.jl")
    global const R0 = 1
    using Plots,ColorSchemes,LaTeXStrings
    pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    plot()
&

## ------------------------------ Cluster data analysis ------------------------------ ##
filename = "data/mobility_defects_sigma_v0.jld2"
@load filename sigmas v0s Ts all_xy_pos all_xy_neg all_rr all_times_collision R_per_core Rtot params_init Ntarget R0 q init_theta init_pos aspect_ratio times tmax comments rhoc runtimes
hrun(runtimes) 

sigmas
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
p=plot(xlabel=L"t",ylabel="Var(y+)")
for i in each(v0s)
    for j in 1#each(sigmas)
        for k in each(Ts)
            plot!(times,[var(ys_pos[i,j,k,t,:]) for t in 1:length(times)])
        end
    end
end
xlims!(-10,400)
p

## Ys -
p=plot(xlabel=L"t",ylabel="Var(y-)")
for i in each(v0s)
    for j in 1#each(sigmas)
        for k in each(Ts)
            plot!(times,[var(ys_neg[i,j,k,t,:]) for t in 1:length(times)])
        end
    end
end
xlims!(-10,400)
p

## Ys +/-
p=plot(xlabel=L"t",ylabel="Var(y±)")
for i in each(v0s)
    for j in 1#each(sigmas)
        for k in each(Ts)
            plot!((times),[var(vcat(ys_pos[i,j,k,t,:],ys_neg[i,j,k,t,:])) for t in 1:length(times)])
        end
    end
end
xlims!(-10,400)
p

# ## Fit the slope of variance in 
# first_index = 4 # robust to any change
# slopes = NaN*zeros(length(v0s), length(sigmas), length(Ts))

# for i in each(v0s)
#     for j in each(sigmas)
#         for k in each(Ts)
#             data = smooth([var(vcat(ys_pos[i,j,k,t,:],ys_neg[i,j,k,t,:])) for t in 1:length(times)],0)
#             last_index = findfirst(isnan,data)
#             if isnothing(last_index) last_index = length(data) end
#             slopes[i,j,k] = (data[last_index-1] - data[first_index])/((times[last_index-1]) - (times[first_index]))
#         end
#     end
# end

# plot(v0s, slopes[:,1,1],m=true)
# plot!(v0s, 0.02v0s .+ 0.13,c=:black)

## Fit the slope of variance in 
xx = rand(10,10)
var(xx[1,:])
ind0 = 50
mobilities = NaN*zeros(length(v0s), length(sigmas), length(Ts), length(times)-ind0) # Var(y±) = mobility*v0*t
for i in each(v0s)
    for j in each(sigmas)
        for k in each(Ts)
            for tt in 1+ind0:length(times)
                var_plus = var(filter(!isnan,ys_pos[i,j,k,tt,:]))
                var_minus = var(filter(!isnan,ys_neg[i,j,k,tt,:]))
                mobilities[i,j,k,tt-ind0] = 0.5*(var_minus + var_plus)/(times[tt])
            end
        end
    end
end
mobilities
mobilities_avg = nanmean(mobilities,(4))[:,:,:,1]
##
plot(v0s,mobilities_avg[1:end,1,2],m=true)
plot!(v0s, 0.05v0s .+ 0.,c=:black)

## times_collision
for i in 1:length(v0s)
    p=histogram(all_times_collision[i,1,1,:],bins=50)
    title!("i = $i, v0 = $(v0s[i])")
    display(p)
end

## R(t) 
p=plot(xlabel=L"t", ylabel=L"R(t)", uaxis=:log)#,ylims=(0,33))
for i in each(v0s)
    for j in each(sigmas)
        for k in 1#each(Ts)
            plot!(times[3:end],all_rr_avg[i,j,k,3:end])
        end
    end
end
p

# Add the theoretical prediction for the fastest curve
# mean_annihilation_time = mean(all_times_collision[end,1,1,:])
# # plot!(0:1:mean_annihilation_time-10,x->3.8exp(0.5*(1+lambertw((mean_annihilation_time-x)*2/exp(1)))),c=:black)
# plot!(0:1:mean_annihilation_time,x->28*sqrt((-x+mean_annihilation_time)/mean_annihilation_time),c=:black)

## R(t*) 
ind_T = 2
using LambertW
p=plot(xlabel=L"t^*", ylabel=L"R(t^*)", xaxis=:log)
for i in 2:length(v0s)
    for j in 1#each(sigmas)
        for k in ind_T
            plot!(times[2:end],all_rr_reverse_avg[i,j,k,2:end])
        end
    end
end
xlims!(0.9,1E3)
xticks!([1,10,100,1000],[L"10^0",L"10^1",L"10^2",L"10^3"])

pcollapse=plot(xlabel=L"\sqrt{v_0}t^*", ylabel=L"R(t^*)/\sqrt{v_0}", xaxis=:log, legend=:topleft)
for i in 2:length(v0s)
    for j in 1#each(sigmas)
        for k in ind_T
            plot!(sqrt(v0s[i])*times[2:end],1/sqrt(v0s[i])*all_rr_reverse_avg[i,j,k,2:end], 
            # plot!(sqrt(0.45*Ts[ind_T]+0.5v0s[i])*times[2:end],1/sqrt(v0s[i])*all_rr_reverse_avg[i,j,k,2:end], 
            label=L"v_0 = "*"$(v0s[i])", rib=0)
        end
    end
end
mu = 1/2
plot!(times[2:700],x->exp(0.5*lambertw(2π*x/mu)),c=:black)
xlims!(0.9,1E3)
ylims!(0,40)
xticks!([1,10,100,1000],[L"10^0",L"10^1",L"10^2",L"10^3"])
pcollapse

plot(p,pcollapse,layout=(1,2),size=(800,400))
# savefig("figures/mobility_defects_v0_T$(Ts[ind_T]).pdf")



## ------------------------------ Annihilation times ------------------------------ 
## ------------------------------ Annihilation times ------------------------------ 
## ------------------------------ Annihilation times ------------------------------ 
## ------------------------------ Annihilation times ------------------------------ 

filename = "data/mobility_defects_sigma_v0.jld2"
@load filename sigmas v0s Ts all_xy_pos all_xy_neg all_rr all_times_collision R_per_core Rtot params_init Ntarget R0 q init_theta init_pos aspect_ratio times tmax comments rhoc runtimes

Ts
avg_times_collision = mean(all_times_collision, dims=4)[:, :, :, 1]
std_times_collision = std(all_times_collision, dims=4)[:, :, :, 1]

##
plot(yaxis=:log)
plot!(v0s, std_times_collision[:, 1, ind_T] ./ avg_times_collision[:, 1, ind_T], m=true)
# plot!(v0s, avg_times_collision[:, 1, ind_T], rib=std_times_collision[:, 1, ind_T], m=true)

##

inds = [1, 5, length(v0s)]
phistogram2 = plot(size=(400, 400), uaxis=:log10, legend=false, legend_title=L"v_0")
for ind in inds
    data = all_times_collision[ind, 1, ind_T, :]
    avg = mean(data)
    histogram!(data / avg, bins=30, lw=0.2, label=L"0.5", normalize=true)
end
# ylims!(0, 93)
# xticks!(1:3, [L"10^{1}", L"10^{2}", L"10^{3}"])
# xlims!(1, 4)
annotate!((0.5, 0.9), text("Distribution of " * L"\tau", 11, :center, :bottom, :black))
annotate!((0.93, 0.02), text(L"\tau", 13, :center, :bottom, :black))

