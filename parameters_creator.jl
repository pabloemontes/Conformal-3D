include("metric-creator.jl")
include("inversion.jl")
include("flux-functions.jl")
#=NOT AXISYMMETRIC.=#
function parameters_creator(chi, C, lengthpars, tol, iter_max, Mass, a)

    χ1_num, χ2_num = chi
    C0, C1, C2 = C


    Lx, Ly, Lz, N, M, O = lengthpars;

    println("creating metric")
    g_ab_table, gᵃᵇ_table = metric_creator(lengthpars, Mass, a)
    println("metric created")
    println("creating christoffels")
    Γ_table = christoffel_creator(lengthpars, Mass, a, (zeros(4),zeros(4),zeros(4,4),zeros(4,4,4)))
    println("christoffels created")
    auxvectors = [zeros(4),zeros(4),zeros(4),zeros(4)]
    auxmatrices = [zeros(4,4), zeros(4,4), zeros(4,4)]

    aux_residual = zeros(14)
    aux_jacobian = zeros(14,14)

    println("I am in parameters creator")
    par_source = (C0, C1, C2, gᵃᵇ_table, (zeros(4),zeros(4)), (zeros(4,4),zeros(4,4),zeros(4,4)), Γ_table)
    #println(par_source[end])
    par_flux = (χ1_num, χ2_num, auxvectors, auxmatrices, g_ab_table, gᵃᵇ_table)
    println("hello!")
    par_inv = (par_flux, tol, iter_max, abs2c_res!, Jabs2c!, aux_residual, aux_jacobian, g_ab_table, gᵃᵇ_table);
    par_inidat = (lengthpars, par_flux)

println("I get here!")
    auxvecs_evolution = []
    for i in 1:21
        push!(auxvecs_evolution, zeros(28))
    end
    auxvecs_evolution = Tuple(auxvecs_evolution)

    T = zeros(4,4)
    A = zeros(4,4,4)

    par_ev = (lengthpars, par_flux, par_source, Fx!, Fy!, Fz!, Is!, Speed_max, auxvecs_evolution, T, A);
    return par_inidat, par_source, par_flux, par_inv, par_ev
    
end

#=
function parameters_creator_axisymmetric(chi, C, lengthpars, tol, iter_max, Mass)

    χ1_num, χ2_num = chi
    C0, C1, C2 = C


    Lx, Ly, N, M = lengthpars;

    println("creating metric")
    g_ab_table, gᵃᵇ_table= metric_creator(lengthpars, Mass)
    println("metric created")
    println("creating christoffels")
    Γ_table = christoffel_creator(lengthpars, Mass, (zeros(4),zeros(4),zeros(4,4),zeros(4,4,4)))
    println("christoffels created")
    auxvectors = [zeros(4),zeros(4),zeros(4),zeros(4)]
    auxmatrices = [zeros(4,4), zeros(4,4), zeros(4,4)]

    aux_residual = zeros(14)
    aux_jacobian = zeros(14,14)

    par_source = (C0, C1, C2, gᵃᵇ_table, (zeros(4),zeros(4)), (zeros(4,4),zeros(4,4),zeros(4,4)), Γ_table)
    par_flux = (χ1_num, χ2_num, auxvectors, auxmatrices, g_ab_table, gᵃᵇ_table)
    par_inv = (par_flux, tol, iter_max, abs2c_res!, Jabs2c!, aux_residual, aux_jacobian, g_ab_table, gᵃᵇ_table);
    par_inidat = (lengthpars, par_flux)


    auxvecs_evolution = []
    for i in 1:41
        push!(auxvecs_evolution, zeros(28))
    end
    auxvecs_evolution = Tuple(auxvecs_evolution)

    T = zeros(4,4)
    A = zeros(4,4,4)

    par_ev = (lengthpars, par_flux, par_source, Fx!, Fz!, Is!, Speed_max, auxvecs_evolution, T, A);

    return par_inidat, par_source, par_flux, par_inv, par_ev
    
end
=#