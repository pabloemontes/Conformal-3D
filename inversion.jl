
#import Pkg; 
#Pkg.activate("Conf_Fluids")
#Pkg.add("Symbolics")
#using Symbolics
#using StaticArrays
using LinearAlgebra

include("auxfunctions.jl")
include("auxfunctions-varcreation.jl")
include("auxmacros-varcreation.jl")
include("trace-free-utils.jl")
#If you are lost, some of the definitions of things are in TandA calculations.

function abs2c!(convars, absvars, idx_x, idx_y, idx_z, p)
    
    χ_1, χ_2, auxvectors, auxmatrices, g_ab_table, gᵃᵇ_table = p
  
    g_ab = @view g_ab_table[idx_x, idx_y, idx_z, :, :]
    gᵃᵇ = @view gᵃᵇ_table[idx_x, idx_y, idx_z, :, :]

    lᵃ, l_a, ηᵃ, ξᵃ = auxvectors
    ξ_ab, ξᵃᵇ, Mᵃᵇ = auxmatrices

    ξ_a = absvars[1:4]
    ξ_ab_vec = absvars[5:end]

    #We fill ξ_ab
    ξ_ab_fun!(ξ_ab_vec, ξ_ab)
    #Vectors
    rise_index!(ξᵃ, gᵃᵇ, ξ_a)
    ξᵃᵇ_fun!(ξᵃᵇ, gᵃᵇ, ξ_ab)
    mul!(lᵃ, ξᵃᵇ, ξ_a)
    lower_index!(l_a, g_ab, lᵃ)
    mul!(ηᵃ, ξᵃᵇ, l_a)
    #Tensors
    Mᵃᵇ_fun!(Mᵃᵇ, g_ab, ξᵃᵇ)

    #Scalars
    μ = ξᵃ'ξ_a
    ν = ξᵃ'l_a
    ψ_1 = ψ_1_function(ξᵃᵇ, ξ_ab)
    ψ_2 = lᵃ'l_a
   
    Conservative_fun!(convars, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab)
    
    return convars;
    #We set the dummy fluxes of the abstract variables to zero.
end


function abs2c_res!(residual, convars, absvars, idx_x, idx_y, idx_z, p)
    abs2c!(residual, absvars, idx_x, idx_y, idx_z, p)
    residual .-= convars;
    return residual;
end




function Jabs2c!(Jacobian, convars, absvars, idx_x, idx_y, idx_z, p)
    
    χ_1, χ_2, auxvectors, auxmatrices, g_ab_table, gᵃᵇ_table = p
  
    g_ab = @view g_ab_table[idx_x, idx_y, idx_z, :, :]
    gᵃᵇ = @view gᵃᵇ_table[idx_x, idx_y, idx_z, :, :]

    lᵃ, l_a, ηᵃ, ξᵃ = auxvectors
    ξ_ab, ξᵃᵇ, Mᵃᵇ = auxmatrices


    ξ_a = absvars[1:4]
    ξ_ab_vec = absvars[5:end]

    #We fill ξ_ab
    ξ_ab_fun!(ξ_ab_vec, ξ_ab)
    #Vectors
    rise_index!(ξᵃ, gᵃᵇ, ξ_a)
    ξᵃᵇ_fun!(ξᵃᵇ, gᵃᵇ, ξ_ab)
    mul!(lᵃ, ξᵃᵇ, ξ_a)
    lower_index!(l_a, g_ab, lᵃ)
    mul!(ηᵃ, ξᵃᵇ, l_a)
    #Tensors
    Mᵃᵇ_fun!(Mᵃᵇ, g_ab, ξᵃᵇ)

    #Scalars
    μ = ξᵃ'ξ_a
    ν = ξᵃ'l_a
    ψ_1 = ψ_1_function(ξᵃᵇ, ξ_ab)
    ψ_2 = lᵃ'l_a
   
    #Indices of the conservative/flux variables
    #Energy-momentum
    Jacobian_fun!(Jacobian, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab)
      
    return Jacobian;
end


function NR_step!(F!, Jac!, aux_jacobian, residual, convars, absvars, idx_x, idx_y, idx_z, p)
    χ_1, χ_2, auxvectors, auxmatrices, g_ab_table, gᵃᵇ_table = p
    #g_ab = @view g_ab_table[idx_x, idx_y, idx_z, :, :]
    #gᵃᵇ = @view gᵃᵇ_table[idx_x, idx_y, idx_z, :, :]
    F!(residual,convars,absvars, idx_x, idx_y, idx_z, p)
    Jac!(aux_jacobian, convars,absvars, idx_x, idx_y, idx_z, p)
    residual .= aux_jacobian\residual
    #remove_trace(residual[5:end], g_ab, gᵃᵇ)
    absvars .-= residual
end


function c2abs!(convars, absvars, idx_x, idx_y, idx_z, p)
    pars, tol, iter_max, F!, Jac!, residual, aux_jacobian = p 
    residual .= 1.0 #residual
    iter = 1
    #while F(flu[j,:],con[j,:], χ)'*F(flu[j,:],con[j,:], χ) > tol && iter < iter_max # para F que no cambia valores
    while maximum(abs.(residual)) > tol && iter < iter_max
        NR_step!(F!, Jac!, aux_jacobian, residual, convars, absvars, idx_x, idx_y, idx_z, pars)
        iter = iter + 1
        #abs2c!(convars, absvars, pars)
    end
end


#=
@variables a_jvar[1:14], c_jvar[1:14], r_jvar[1:14]
@variables χ_jvar
@variables auxmatrix1_jvar[1:4, 1:4]

aa_jvar = Symbolics.scalarize(a_jvar)
cc_jvar = Symbolics.scalarize(c_jvar)
rr_jvar = Symbolics.scalarize(r_jvar)
aauxmatrix1_jvar = Symbolics.scalarize(auxmatrix1_jvar)

p_jvar = [χ_jvar, auxmatrix1_jvar]
pp_jvar = [χ_jvar, aauxmatrix1_jvar]

JS_alt = Symbolics.jacobian(abs2c!(rr_jvar,cc_jvar,aa_jvar,pp_jvar),aa_jvar);

J_exp_alt = Symbolics.build_function(JS_alt, r_jvar, c_jvar, a_jvar, p_jvar);
Jabs2c = eval(J_exp_alt[1]);
=#