# Physical parameters 
Ntarget = Int(1E4)
aspect_ratio = 1
rho = 1
T = 0.1
sigma = 0.0
v0 = 0
R0 = 1.95

N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)

# Initialisation parameters
init_pos = "random"
init_theta = "hightemp"
r0 = 20#round(Int,Lx/2)
q = 1.0

# Time parameters
tmax = 1E2
times = collect(0:1:tmax) # linear time
# times = logspace(1,tmax,30,digits=1) # log time


params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
:rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
:N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)
