cd("D:/Documents/Research/projects/kuramoto_ballistic")
    using JLD2,StatsBase,Distributions,LinearAlgebra,Parameters,Random,BenchmarkTools,Hungarian
    include("../methods.jl")
    global const R0 = 1
    using Plots,ColorSchemes,LaTeXStrings
    gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    plot()
&
#= This file investigates the possible finite size effect
of the critical velocity with N.
=#

Ns = round.(Int,logspace(1E3,1E5,10,digits=0))
rho = 1
T = 0.1
v0s    = logspace(0.01,2,20,digits=2)
sigmas = [0,0.2]
R = 10

tmax = 1E3

vc = zeros(length(Ns),length(sigmas),R)
seuil = 0.5
z = @elapsed for r in 1:R
    println("r=$r/$R")
    for j in each(sigmas)
        σ = sigmas[j]
        for n in each(Ns)
            N = Ns[n] ; L = round(Int,sqrt(N/rho))
            P = zeros(length(v0s))
            for i in each(v0s)
                v0 = v0s[i]
                dt = determine_dt(T,σ,v0,N,rho)
                t = 0. ; pos,thetas,psis,omegas = initialisation(N,L,σ,["hightemp"])
                while t<tmax
                    t += dt
                    pos,thetas = update(pos,thetas,psis,omegas,T,v0,N,L,dt)
                end
                P[i] = polarOP(thetas)[1]
            end
            tmp = findfirst(x->x>seuil,P)
            if tmp == nothing tmp = length(v0s) end
            vc[n,j,r] = v0s[tmp]
        end
    end
end
prinz(z)
vc_avg = nanmean(vc,3)
p=plot(xlabel="N",ylabel=L"v_c")
    for j in each(sigmas)
        plot!(Ns,vc_avg[:,j,1],xaxis=:log,m=true)
    end
    p
