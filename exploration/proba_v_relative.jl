cd("D:/Documents/Research/projects/kuramoto_ballistic")
    using JLD2,StatsBase,Distributions,LinearAlgebra,Parameters,Random,BenchmarkTools,Hungarian
    include("../methods.jl")
    global const R0 = 1
    using Plots,ColorSchemes,LaTeXStrings
    gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
    plot()
&
N = Int(1E4)
psis1 = 2pi * rand(N)
psis2 = 2pi * rand(N)
phis = (psis2 - psis1)


histogram(phis, bins=100)
histogram(mod.(phis, 2pi), bins=100, normalize=true)
histogram(vrel.(1, phis), bins=100, normalize=true)
histogram(1 ./vrel.(1, phis), bins=100, normalize=false)
xlims!(0, 10)

vrel(v, phi)=v*(1-cos(phi))