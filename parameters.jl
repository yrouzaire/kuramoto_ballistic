Ntarget = Int(1E3)
aspect_ratio = 1
rho = 1
T = 0.0
sigma = 1.0
v0 = 1.0
R0 = 1

N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)

# Initialisation parameters
init_pos = "random"
init_theta = "pair"
r0 = 20.0
q = 1.0
params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio, :rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0, :N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init)
