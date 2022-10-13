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

Ns = round.(Int,logspace(5E2,1E4,10,digits=0))
rho = 1
T = 0.1
v0s    = reverse(logspace(0.01,2,20,digits=2))
sigmas = [0]
R = 5

tmax = 1E3

vc = zeros(length(Ns),length(sigmas),R)
seuil_break = 0.8
seuil_crit = 0.3
times_break = range(1, stop=tmax, length=20)
z = @elapsed for r in 1:R
    for j in each(sigmas)
        σ = sigmas[j]
        for n in each(Ns)
            N = Ns[n] ; L = round(Int,sqrt(N/rho))
            for i in each(v0s)
                v0 = v0s[i]
                println("r=$r/$R, σ = $σ, N = $N, v0 = $v0")
                dt = determine_dt(T,σ,v0,N,rho)
                t = 0. ; pos,thetas,psis,omegas = initialisation(N,L,σ,["hightemp"])
                token = 1
                while t<tmax
                    t += dt
                    pos,thetas = update(pos,thetas,psis,omegas,T,v0,N,L,dt)
                    if t ≥ times_break[token]
                        if polarOP(thetas)[1] > seuil_break
                            println("Simu stopped because P>$(seuil_break) : v0 > vc")
                            break # stop this simulation, pass to the next v0
                        end
                        token = min(token+1,length(times_break))
                    end
                end
                if polarOP(thetas)[1] < seuil_crit
                    vc[n,j,r] = v0
                    println("Critical velocity because P<$(seuil_crit) : v0 < vc")
                    break # you found the critical velocity
                end
            end
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
