using StaticArrays
using Symbolics


using OrdinaryDiffEq
using Plots
using FileIO
using IJulia
using Plots
using JLD2
using Printf


include("inversion_ext.jl")
include("choques_utils.jl") # all functions needed for evolution
include("Flux_function_ext.jl")
include("initial_data.jl")

N_FIELDS = 10
name = :chichonpaper_wide
out_name = "chichon_paper_long" 
#name = :square
#out_name = "square" 
L=200.0 # Length of the integration region
M=2000 # number of space-points in the discretization
tf=100.0 # final time
dx = L/M
dt = 0.1 * dx # this depends on the maximal propagation speed
tol = 10^(-10) # error tolerance in inversion function (Newton-Raphson)
iter_max = 400  # maximum number of NR iterations
χ = [-1.0,1.0,-1.0]# [1.,1.,0.06]# the equation parameters   #par_f
C = [-10; 10; -10] # new source parámeters     1/κ, 1/λ, 1/η   #par_s
x = [-L/2+(i-1)*dx for i in 1:M];

#we initialize empty arrays
par_mem = auxvectors(N_FIELDS);
wenoz_auxvectors = createWENOZvectors(N_FIELDS);

par_source = (χ, C) # parameters to use on Is! 
par_inv = (χ, tol, iter_max, N_FIELDS, M, f2c!, Jf2c) # parameters for the inversion
par = (par_source, par_inv, 1.0/dx, N_FIELDS, M, Flux_imp!, Speed_max, Is!, c_to_f!, par_mem);

scheme = :wenoz
boundary = :free

evolution! = construct_evolution_function(conservative_scheme = scheme, boundary = boundary)

# χ[3] = χ2

χ = [-1.0, 1.0, -1.0]

C0 = -10*1e3
C1 = 10*1e3
C2 = -10*1e3
C = [C0, C1, C2]*1.0

bigfilename = "Resultados/chichon_runs_C0variable.jld"

jldopen(bigfilename, "w") do file
    if "wenoz/" in keys(file)
        println("test")
    end
                
end

function checktree(dictionary, keylist)
    #recursive dictionary
    indict = dictionary
    for key in keylist
        if key in keys(indict)
            indict = indict[key]
        else
            println("key not in dict, exiting")
            return false
            break
        end
    end
    println("key in dict!")
    return true
end

χ[1] = -1.0  #χ0
χ[2] = 1.0   #χ1
χ[3] = -1.0  #χ2

χ0 = χ[1]
χ1 = χ[2]
χ2 = χ[3]
C = [-10, 10, -1]*300
C_0 = C[1]
C_1 = C[2]
C_2 = C[3]
scale = 
#for scale in [31,32,33,34,35,36,37,38,39]
dt = dx
#out_name = "square"

sol = 0

C0 = -10.0
C1 = 10.0
C2 = -10.0
C = [C0,C1,C2]

χ[1] = -1.0  #χ0
χ[2] = 0.0   #χ1
χ[3] = -1.0  #χ2

println("Starting Euler Run...")

jldopen(bigfilename, "a+") do file
        tf = 100
        #χ2 = 1.0    
        χ[3] = χ2
        C = [C0, C1, C2]
                
        par_source = (χ, C) # parameters to use on the equations 
        par_inv = (χ, tol, iter_max, N_FIELDS, M, f2c!, Jf2c) # parameters for the inversion

        u_i = zeros(M,N_FIELDS) # flux variables for the initial data (these are used for the flux function)
        du = zeros(M,N_FIELDS) # for the rhs
        par_inidat = χ, C, N_FIELDS, M, dx, Euler, f2c!, Jf2c , JA, Initialize_data!
        nu = 0.0
        r1 = 0.0
        tau = 0.0
        mu0 = -6.0
        δmu = 3.0
        filenamelist = [string(scheme),out_name, "euler"]
        
        #filenamelist = ["euler"]
        if checktree(file, filenamelist)
            println("filename already pressent, avoiding")
        else
            paraux = (mu0, δmu, nu, r1, tau)
            println("initializing data...")
            create_initial_data_pablo!(name, x, u_i, par_inidat, paraux)
            println("initialized.")
            con_0 = view(u_i, :, 1:N_FIELDS÷2)
            flu_0 = view(u_i, :, N_FIELDS÷2+1:N_FIELDS)
            θ = 1.0
            par = (par_source, par_inv, 1.0/dx, N_FIELDS, M, Flux_imp!, Speed_max, Is!, c_to_f!, par_mem);

            prob = ODEProblem(evolution!,u_i,(0.0,tf),par);
            println("Starting to solve...")
            #solmp5 = solve(prob,SSPRK33(),dt=dt,saveat=0.01*tf);
            sol = solve(prob,SSPRK33(),dt=dt,saveat=0.01*tf);

            println("Solved")
            file_name = ""
            for element in filenamelist
                file_name = file_name*"/"*element
            end

            file[file_name] = sol
            file["x_array"] = x
            file["dx"] = dx
            file["dt"] = dt
        end
end


χ[1] = -1.0  #χ0
χ[2] = 1.0   #χ1
χ[3] = -1.0  #χ2



#=
println("Starting to modify C2")

N_FIELDS = 10
name = :chichonpaper_wide
out_name = "chichon_paper_long" 
#name = :square
#out_name = "square" 
L=200.0 # Length of the integration region
M=2000 # number of space-points in the discretization
tf=100.0 # final time
dx = L/M
dt = 0.1 * dx # this depends on the maximal propagation speed
tol = 10^(-10) # error tolerance in inversion function (Newton-Raphson)
iter_max = 40  # maximum number of NR iterations
χ = [-1.0,1.0,-1.0]# [1.,1.,0.06]# the equation parameters   #par_f
C = [-10; 10; -10] # new source parámeters     1/κ, 1/λ, 1/η   #par_s
x = [-L/2+(i-1)*dx for i in 1:M];

#we initialize empty arrays
par_mem = auxvectors(N_FIELDS);
wenoz_auxvectors = createWENOZvectors(N_FIELDS);

par_source = (χ, C) # parameters to use on Is! 
par_inv = (χ, tol, iter_max, N_FIELDS, M, f2c!, Jf2c) # parameters for the inversion
par = (par_source, par_inv, 1.0/dx, N_FIELDS, M, Flux_imp!, Speed_max, Is!, c_to_f!, par_mem);

scheme = :wenoz
boundary = :free

χ2 = -1.0
χ1 = 1.0
χ[2] = χ1

C0 = -1000
C1 = 1000
C2 = -1000


jldopen(bigfilename, "a+") do file
    for C2 in [-1000,-100,-10,-1, 0]
        tf = 100
        #χ2 = 1.0    
        χ[3] = χ2
        C = [C0, C1, C2]
        
        println("C = ", C)
        println("χ = ", χ)
        par_source = (χ, C) # parameters to use on the equations 
        par_inv = (χ, tol, iter_max, N_FIELDS, M, f2c!, Jf2c) # parameters for the inversion

        u_i = zeros(M,N_FIELDS) # flux variables for the initial data (these are used for the flux function)
        du = zeros(M,N_FIELDS) # for the rhs
        par_inidat = χ, C, N_FIELDS, M, dx, Euler, f2c!, Jf2c , JA, Initialize_data!
        nu = 0.0
        r1 = 0.0
        tau = 0.0
        mu0 = -6.0
        δmu = 3.0
        filenamelist = [string(scheme),out_name,"χ2="*@sprintf("%.2lf", χ2),"C0="*@sprintf("%.4lf", C[1]),"C1="*@sprintf("%.4lf", C[2]),"C2 ="*@sprintf("%.4lf", C[3]),"nu="*@sprintf("%.2lf", nu),"r1="*@sprintf("%.2lf", r1),"tau="*@sprintf("%.2lf", tau)]
        
        #filenamelist = ["euler"]
        if checktree(file, filenamelist)
            println("filename already pressent, avoiding")
            continue
        end
        paraux = (mu0, δmu, nu, r1, tau)
        println("initializing data...")
        create_initial_data_pablo!(name, x, u_i, par_inidat, paraux)
        println("initialized.")
        con_0 = view(u_i, :, 1:N_FIELDS÷2)
        flu_0 = view(u_i, :, N_FIELDS÷2+1:N_FIELDS)
        θ = 1.0
        par = (par_source, par_inv, 1.0/dx, N_FIELDS, M, Flux_imp!, Speed_max, Is!, c_to_f!, par_mem);
        
        prob = ODEProblem(evolution!,u_i,(0.0,tf),par);
        println("Starting to solve...")
        sol = solve(prob,SSPRK33(),dt=dt,saveat=0.01*tf);

        println("Solved")
        file_name = ""
        for element in filenamelist
            file_name = file_name*"/"*element
        end

        file[file_name] = sol
        
    end
end



println("Starting to modify C1")

N_FIELDS = 10
name = :chichonpaper_wide
out_name = "chichon_paper_long" 
#name = :square
#out_name = "square" 
L=200.0 # Length of the integration region
M=2000 # number of space-points in the discretization
tf=100.0 # final time
dx = L/M
dt = 0.1 * dx # this depends on the maximal propagation speed
tol = 10^(-10) # error tolerance in inversion function (Newton-Raphson)
iter_max = 40  # maximum number of NR iterations
χ = [-1.0,1.0,-1.0]# [1.,1.,0.06]# the equation parameters   #par_f
C = [-10; 10; -10] # new source parámeters     1/κ, 1/λ, 1/η   #par_s
x = [-L/2+(i-1)*dx for i in 1:M];

#we initialize empty arrays
par_mem = auxvectors(N_FIELDS);
wenoz_auxvectors = createWENOZvectors(N_FIELDS);

par_source = (χ, C) # parameters to use on Is! 
par_inv = (χ, tol, iter_max, N_FIELDS, M, f2c!, Jf2c) # parameters for the inversion
par = (par_source, par_inv, 1.0/dx, N_FIELDS, M, Flux_imp!, Speed_max, Is!, c_to_f!, par_mem);

scheme = :wenoz
boundary = :free

χ2 = -1.0
χ1 = 1.0
χ[2] = χ1

C0 = -1000
C1 = 1000
C2 = -1000


jldopen(bigfilename, "a+") do file
    for C1 in [1000,100,10,1, 0]
        tf = 100
        #χ2 = 1.0    
        χ[3] = χ2
        C = [C0, C1, C2]
        
        println("C = ", C)
        println("χ = ", χ)
        par_source = (χ, C) # parameters to use on the equations 
        par_inv = (χ, tol, iter_max, N_FIELDS, M, f2c!, Jf2c) # parameters for the inversion

        u_i = zeros(M,N_FIELDS) # flux variables for the initial data (these are used for the flux function)
        du = zeros(M,N_FIELDS) # for the rhs
        par_inidat = χ, C, N_FIELDS, M, dx, Euler, f2c!, Jf2c , JA, Initialize_data!
        nu = 0.0
        r1 = 0.0
        tau = 0.0
        mu0 = -6.0
        δmu = 3.0
        filenamelist = [string(scheme),out_name,"χ2="*@sprintf("%.2lf", χ2),"C0="*@sprintf("%.4lf", C[1]),"C1="*@sprintf("%.4lf", C[2]),"C2 ="*@sprintf("%.4lf", C[3]),"nu="*@sprintf("%.2lf", nu),"r1="*@sprintf("%.2lf", r1),"tau="*@sprintf("%.2lf", tau)]
        
        #filenamelist = ["euler"]
        if checktree(file, filenamelist)
            println("filename already pressent, avoiding")
            continue
        end
        paraux = (mu0, δmu, nu, r1, tau)
        println("initializing data...")
        create_initial_data_pablo!(name, x, u_i, par_inidat, paraux)
        println("initialized.")
        con_0 = view(u_i, :, 1:N_FIELDS÷2)
        flu_0 = view(u_i, :, N_FIELDS÷2+1:N_FIELDS)
        θ = 1.0
        par = (par_source, par_inv, 1.0/dx, N_FIELDS, M, Flux_imp!, Speed_max, Is!, c_to_f!, par_mem);
        
        prob = ODEProblem(evolution!,u_i,(0.0,tf),par);
        println("Starting to solve...")
        sol = solve(prob,SSPRK33(),dt=dt,saveat=0.01*tf);

        println("Solved")
        file_name = ""
        for element in filenamelist
            file_name = file_name*"/"*element
        end

        file[file_name] = sol
        
    end
end

=#

println("Starting to modify C0")

N_FIELDS = 10
name = :chichonpaper_wide
out_name = "chichon_paper_long" 
#name = :square
#out_name = "square" 
L=200.0 # Length of the integration region
M=2000 # number of space-points in the discretization
tf=100.0 # final time
dx = L/M
dt = 0.1 * dx # this depends on the maximal propagation speed
tol = 10^(-10) # error tolerance in inversion function (Newton-Raphson)
iter_max = 40  # maximum number of NR iterations
χ = [-1.0,1.0,-1.0]# [1.,1.,0.06]# the equation parameters   #par_f
C = [-10; 10; -10] # new source parámeters     1/κ, 1/λ, 1/η   #par_s
x = [-L/2+(i-1)*dx for i in 1:M];

#we initialize empty arrays
par_mem = auxvectors(N_FIELDS);
wenoz_auxvectors = createWENOZvectors(N_FIELDS);

par_source = (χ, C) # parameters to use on Is! 
par_inv = (χ, tol, iter_max, N_FIELDS, M, f2c!, Jf2c) # parameters for the inversion
par = (par_source, par_inv, 1.0/dx, N_FIELDS, M, Flux_imp!, Speed_max, Is!, c_to_f!, par_mem);

scheme = :wenoz
boundary = :free

χ2 = -1.0
χ1 = 1.0
χ[2] = χ1

C0 = -1000
C1 = 1000
C2 = -1000


jldopen(bigfilename, "a+") do file
    for C0 in [-1000,-100,-10,-1, 0]
        tf = 100
        #χ2 = 1.0    
        χ[3] = χ2
        C = [C0, C1, C2]
        
        println("C = ", C)
        println("χ = ", χ)
        par_source = (χ, C) # parameters to use on the equations 
        par_inv = (χ, tol, iter_max, N_FIELDS, M, f2c!, Jf2c) # parameters for the inversion

        u_i = zeros(M,N_FIELDS) # flux variables for the initial data (these are used for the flux function)
        du = zeros(M,N_FIELDS) # for the rhs
        par_inidat = χ, C, N_FIELDS, M, dx, Euler, f2c!, Jf2c , JA, Initialize_data!
        nu = 0.0
        r1 = 0.0
        tau = 0.0
        mu0 = -6.0
        δmu = 3.0
        filenamelist = [string(scheme),out_name,"χ2="*@sprintf("%.2lf", χ2),"C0="*@sprintf("%.4lf", C[1]),"C1="*@sprintf("%.4lf", C[2]),"C2 ="*@sprintf("%.4lf", C[3]),"nu="*@sprintf("%.2lf", nu),"r1="*@sprintf("%.2lf", r1),"tau="*@sprintf("%.2lf", tau)]
        
        #filenamelist = ["euler"]
        if checktree(file, filenamelist)
            println("filename already pressent, avoiding")
            continue
        end
        paraux = (mu0, δmu, nu, r1, tau)
        println("initializing data...")
        create_initial_data_pablo!(name, x, u_i, par_inidat, paraux)
        println("initialized.")
        con_0 = view(u_i, :, 1:N_FIELDS÷2)
        flu_0 = view(u_i, :, N_FIELDS÷2+1:N_FIELDS)
        θ = 1.0
        par = (par_source, par_inv, 1.0/dx, N_FIELDS, M, Flux_imp!, Speed_max, Is!, c_to_f!, par_mem);
        
        prob = ODEProblem(evolution!,u_i,(0.0,tf),par);
        println("Starting to solve...")
        sol = solve(prob,SSPRK33(),dt=dt,saveat=0.01*tf);

        println("Solved")
        file_name = ""
        for element in filenamelist
            file_name = file_name*"/"*element
        end

        file[file_name] = sol
        
    end
end

