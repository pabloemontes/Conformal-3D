using Symbolics
using Latexify
using LinearAlgebra
include("auxfunctions-varcreation.jl")
include("auxmacros-varcreation.jl")
include("myregister.jl")


if (@isdefined debug) == false
    debug = false
end

#Definiciones de ξa
debug && println("Definiendo ξᵃ y ξ_a")
@my_covector ξ
@my_vector ξ

#Definiciones de ξab
debug && println("Definiendo ξᵃᵇ y ξ_ab")
@my_symmetric_tensor_upup ξ
@my_symmetric_tensor_dndn ξ
@my_symmetric_tensor_upup M

#Construcción de vectores auxiliares para ξab
#=Definimos también ξ_ab_vec, que tiene los elementos únicos de ξ_ab,
y ξ_full, que tiene a ξ_a y a ξ_ab_vec=#
ξ_ab_vec_symbols = vec_triu(ξ_ab_symbols) #only the upper triangular form of ξ_ab_symbols as a vector
ξ_ab_vec = vec_triu(ξ_ab) #only the upper triangular form of ξ_ab as a vector
ξ_full_symbols = [ξ_a_symbols..., ξ_ab_vec_symbols...] #concatenation of ξ_a_symbols and ξ_ab_vec_symbols
ξ_full = [ξ_a..., ξ_ab_vec...] #concatenation of ξ_a and ξ_ab_vec


#Definición de la métrica
debug && println("Definiendo la métrica")
@my_symmetric_tensor_dndn g
@my_symmetric_tensor_upup g

mutable struct tensor
    dd
    uu
    ud
end

g = tensor(g_ab, gᵃᵇ, Diagonal(ones(D)))
#=
g_mink_dd = Diagonal(ones(D))
g_mink_dd[1,1] = -1
g_mink_uu = copy(g_mink_dd)
g_mink_ud = Diagonal(ones(D))
g = tensor(g_mink_dd, g_mink_uu, g_mink_ud)
=#

#Variables auxiliares
debug && println("Definiendo variables auxiliares")
@my_vector η
@my_vector l
@variables μ, ν, ψ_1, ψ_2

#Registrando funciones necesarias
for f in (:μ_fun, ξᵃ_fun_symbols...) #Functions that depend on ξ_a
    @eval @my_register_symbolic $f($(ξ_a_symbols...))
end
for f in (:ψ_1_fun, ξ_ab_fun_symbols..., ξᵃᵇ_fun_symbols..., Mᵃᵇ_fun_symbols...) #Functions that depend on ξ_ab
    @eval @my_register_symbolic $f($(ξ_ab_vec_symbols...))
end
for f in (:ν_fun, :ψ_2_fun, ηᵃ_fun_symbols..., lᵃ_fun_symbols...) #Functions that depend on bot ξ_a and ξ_ab
    @eval @my_register_symbolic $f($(ξ_full_symbols...))
end
@eval lᵃ_funs = [$(lᵃ_fun_symbols...)]
@eval ηᵃ_funs = [$(ηᵃ_fun_symbols...)]
@eval ξᵃ_funs = [$(ξᵃ_fun_symbols...)]
@eval ξ_ab_funs = reshape([$(ξ_ab_fun_symbols...)], D, D);
@eval ξᵃᵇ_funs = reshape([$(ξᵃᵇ_fun_symbols...)], D, D);
@eval Mᵃᵇ_funs = reshape([$(Mᵃᵇ_fun_symbols...)], D, D);