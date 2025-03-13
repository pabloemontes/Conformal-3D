using Symbolics
using Latexify
using LinearAlgebra
using NLsolve

#Definiendo dimensión del problema
const D = 4
const SymTensorSize = D*(D+1)÷2;
const Dlarge = D + SymTensorSize

include("variablescreator.jl")
include("auxfunctions.jl")
include("auxmacros.jl")
include("auxmacros-varcreation.jl")
include("auxfunctions-varcreation.jl")



debug = false
@trace_free_derivatives(true)
## Generating function
@variables χ_1, χ_2
@variables μ, ν, ψ_1, ψ_2

chi_cero(μ) = -1/μ
chi_uno(μ, χ_1) = χ_1/(μ^3)
chi_dos_1(μ, χ) = -χ_2/μ^(3)
chi_dos_2(μ, χ) = 12*χ_2/μ^(4)
chi_dos_3(μ, χ) = -24*χ_2/μ^(5);
chi(μ, ν, ψ_1, ψ_2, ψ_3, χ_1, χ_2) = chi_cero(μ) + ν*chi_uno(μ, χ_1) + ψ_1*chi_dos_1(μ, χ_2) + ψ_2*chi_dos_2(μ, χ_2)+ ψ_3*chi_dos_3(μ, χ_2)
μ_fun_eval = μ_fun(ξ_a...)
ν_fun_eval = ν_fun(ξ_full...)
ψ_1_fun_eval = ψ_1_fun(ξ_ab_vec...)
ψ_2_fun_eval = ψ_2_fun(ξ_full...)
ψ_3_fun_eval = ν_fun_eval^2


debug && println("Derivadas de chi...")
chi_expr = chi(μ_fun_eval, ν_fun_eval, ψ_1_fun_eval, ψ_2_fun_eval, ψ_3_fun_eval, χ_1, χ_2)
#derivada de chi respecto de ξa
dchi_dcsia = Symbolics.gradient(chi_expr, ξ_a, simplify = true)

debug && println("Calculando Tab y Aabc...")
#Tab and Aabc
Tab = [Symbolics.derivative(dchi_dcsia[i], ξ_a[j]) for i in 1:D, j in 1:D]
Aabc = [Symbolics.derivative(dchi_dcsia[i], ξ_ab[j,k]) for i in 1:D, j in 1:D, k in 1:D];

#Unique indices for a symmetric matrix
non_trivial_indexes = []
for i in 1:4
    for j in i:4
        push!(non_trivial_indexes, (i,j))
    end
end

#Expression for Conservative Variables
ConservativeVariables = [Tab[1,:]..., [Aabc[1,i,j] for (i,j) in non_trivial_indexes]...];


#Jacobian of the Conservative variables

@trace_free_derivatives(false)  #We need the full derivatives, not the trace-free ones
JacobianConservativeVariables = zeros(14,14)*Num(0)
for i in 1:14
    for j in 1:14
        JacobianConservativeVariables[i,j] = Symbolics.derivative(ConservativeVariables[i], ξ_full[j])
    end
end

#Right now, Tab, Aabc, ConservativeVariables and JacobianConservativeVariables have terms like
#ξᵃ(ξ_a...), lᵃ(ξ_a..., ξ_ab_vec...)..., and we want to change it to ξᵃ, lᵃ...

replace_dictionary = Dict([μ_fun(ξ_a...) => μ])
replace_dictionary[ν_fun(ξ_full...)] = ν
replace_dictionary[ψ_1_fun(ξ_ab_vec...)] = ψ_1
replace_dictionary[ψ_2_fun(ξ_full...)] = ψ_2

for i in 1:D
    replace_dictionary[lᵃ_funs[i](ξ_full...)] = lᵃ[i]
    replace_dictionary[ξᵃ_funs[i](ξ_a...)] = ξᵃ[i]
    replace_dictionary[ηᵃ_funs[i](ξ_full...)] = ηᵃ[i]
end
for i in 1:D
    for j in i:D
        replace_dictionary[ξᵃᵇ_funs[i,j](ξ_ab_vec...)] = ξᵃᵇ[i,j]
        replace_dictionary[Mᵃᵇ_funs[i,j](ξ_ab_vec...)] = Mᵃᵇ[i,j]
        replace_dictionary[g.uu[i,j]] = gᵃᵇ[i,j]
        replace_dictionary[g.dd[i,j]] = g_ab[i,j]
    end
end

Tab_simple = Symbolics.simplify(substitute(Tab, replace_dictionary));
Aabc_simple = Symbolics.simplify(substitute(Aabc, replace_dictionary));
Conservative_simple = Symbolics.simplify(substitute(ConservativeVariables, replace_dictionary));
Jacobian_simple = Symbolics.simplify(substitute(JacobianConservativeVariables, replace_dictionary));


debug && println("Construyendo funciones...")

Tab_fun = eval(Symbolics.build_function(Tab_simple, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab)[1])
Aabc_fun = eval(Symbolics.build_function(Aabc_simple, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab)[1])
Conservative_fun! = eval(Symbolics.build_function(Conservative_simple, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab)[2])
Jacobian_fun! = eval(Symbolics.build_function(Jacobian_simple, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab)[2])


#We also want the functions as an array of functions instead of a function of arrays
Tab_functions = Array{Function, 2}(undef, 4,4)
Aabc_functions = Array{Function, 3}(undef, 4,4,4)
for i in 1:4
    for j in 1:4
        Tab_functions[i,j] = eval(Symbolics.build_function(Tab_simple[i,j], ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ,g_ab))
    end
end

for i in 1:4
    for j in 1:4
        for k in 1:4
            Aabc_functions[i,j,k] = eval(Symbolics.build_function(Aabc_simple[i,j,k], ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ,g_ab))
        end
    end
end

debug && println("Cargando inversión y flujos...")
#We define invertion and flux functions
include("inversion.jl")
include("flux-functions.jl")
debug && println("Done");