#Load important packages

using LinearAlgebra
using NLsolve

debug = false;
const D = 4
include("flux-and-jacobian-functions.jl")
include("inversion.jl")
include("flux-functions.jl")
include("initial-data.jl")
include("choques_utils.jl")
include("auxfunctions.jl")
include("metric.jl")


initial_data_type = :chichon_2D

#Grid parameters
N = 50
M = 50
N_Fields = 14
Lx = 100.0
Ly = 100.0
lengthpars = (Lx, Ly, N, M, N_Fields, mink_metric!)




#initial data vector
initial_data = ones(N,M,2*N_Fields);


C0 = 10
C1 = 10
C2 = 10

χ = (1.0, 1.0)
C = (C0, C1, C2)

tol = 1e-15
iter_max = 10

par_inidat, par_source, par_flux, par_inv, par_ev = parameters_creator(χ, C, lengthpars, tol, iter_max)


println("Starting Dissipative Run...")

dx = Lx/N
dt = dx*0.1

tf = 50.0

create_initial_data(initial_data_type, initial_data, par_inidat)
println("initialized.")
    
prob = ODEProblem(evolution!, initial_data, (0.0,tf), (par_inv, par_ev));
println("Starting to solve...")

sol = solve(prob,SSPRK33(),dt=dt,saveat=0.01*tf);