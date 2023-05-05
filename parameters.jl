# Physical parameters 
Ntarget = Int(1E3)
aspect_ratio = 1
rho = 1
T = 0.1
sigma = 0.0
v0 = 0
R0 = 2

N, Lx, Ly = effective_number_particle(Ntarget, rho, aspect_ratio)

# Initialisation parameters
init_pos = "square"
init_theta = "hightemp"
r0 = round(Int,Lx/2)
q = 1.0

# Phonons parameters, for immobile particles (v = 0) only
phonons = false
if phonons @assert v0 == 0.0 "Phonons only make sense for immobile particles! " end
if phonons @assert init_theta ≠ "single" "Phonons only make sense for PBC! " end
if phonons @assert aspect_ratio == 1 "Phonons only implemented for square box ! (for now) " end
phonon_amplitude = 1
phonon_k = 1*(2π/Lx) # wavenumber
phonon_omega = 0 # "frequency" (up to a factor 2π)


# Time parameters
tmax = 1E2
times = collect(0:1:tmax) # linear time
# times = logspace(1,tmax,30,digits=1) # log time


params_init = Dict(:init_pos => init_pos, :init_theta => init_theta, :r0 => r0, :q => q)
params_phonons = Dict(:phonons => phonons, :phonon_amplitude => phonon_amplitude, :phonon_k => phonon_k, :phonon_omega => phonon_omega)
param = Dict(:Ntarget => Ntarget, :aspect_ratio => aspect_ratio,
:rho => rho, :T => T, :R0 => R0, :sigma => sigma, :v0 => v0,
:N => N, :Lx => Lx, :Ly => Ly, :params_init => params_init, :params_phonons => params_phonons)
