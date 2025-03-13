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
X1j = [Tab[2,:]..., [Aabc[2,i,j] for (i,j) in non_trivial_indexes]...];
X2j = [Tab[3,:]..., [Aabc[3,i,j] for (i,j) in non_trivial_indexes]...];
X3j = [Tab[4,:]..., [Aabc[4,i,j] for (i,j) in non_trivial_indexes]...];


#Jacobian of the Conservative variables

@trace_free_derivatives(false)  #We need the full derivatives, not the trace-free ones
JacobianConservativeVariables = Symbolics.jacobian(ConservativeVariables, ξ_full)
J1 = Symbolics.simplify(Symbolics.jacobian(X1j, ξ_full))
J2 = Symbolics.simplify(Symbolics.jacobian(X2j, ξ_full))
J3 = Symbolics.simplify(Symbolics.jacobian(X3j, ξ_full))


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
J1_simple = Symbolics.simplify(substitute(J1, replace_dictionary));
J2_simple = Symbolics.simplify(substitute(J2, replace_dictionary));
J3_simple = Symbolics.simplify(substitute(J3, replace_dictionary));


debug && println("Construyendo funciones...")
outfile = "flux-and-jacobian-functions.jl"
write(outfile, "#Funciones T, A, Conservative_simple y Jacobian_simple=#\n")


open(outfile,"a") do outfile

    Tab_fun = Symbolics.build_function(Tab_simple, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab)[1]
    Aabc_fun = Symbolics.build_function(Aabc_simple, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab)[1]
    Conservative_fun! = Symbolics.build_function(Conservative_simple, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab)[2]
    Jacobian_fun! = Symbolics.build_function(Jacobian_simple, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab)[2]
    J1_fun = eval(Symbolics.build_function(J1_simple, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab)[2])
    J2_fun = eval(Symbolics.build_function(J2_simple, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab)[2])
    J3_fun = eval(Symbolics.build_function(J3_simple, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab)[2])

    write(outfile, "\n#=Tab as a function that receives ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab and returns a matrix.=#\n")
    write(outfile, "Tab_fun = " * string(Tab_fun))

    write(outfile, "\n#=Aabc as a function that receives ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab and returns a matrix.=#\n")
    write(outfile, "Aabc_fun = " * string(Aabc_fun))
    
    write(outfile, "\n#=Function that receives cons, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab and stores the values of
    the conservative variables in cons.=#\n")
    write(outfile, "Conservative_fun! = " * string(Conservative_fun!))

    write(outfile, "\n#=Function that receives J, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab and stores the values of
    the jacobian of the conservative variables in cons.=#\n")
    write(outfile, "Jacobian_fun! = " * string(Jacobian_fun!))
    write(outfile, "\n#=Function that receives J, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab and stores the values of
    the jacobian of X¹ᴬ variables in cons.=#\n")
    write(outfile, "J1_fun! = " * string(Jacobian_fun!))
    write(outfile, "\n#=Function that receives J, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab and stores the values of
    the jacobian of X²ᴬ variables in cons.=#\n")
    write(outfile, "J2_fun! = " * string(Jacobian_fun!))
    write(outfile, "\n#=Function that receives J, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ, g_ab and stores the values of
    the jacobian of X³ᴬ variables in cons.=#\n")
    write(outfile, "J3_fun! = " * string(Jacobian_fun!))


    #We also want the functions as an array of functions instead of a function of arrays
    
    
    write(outfile, "\n#=Now we create two arrays of functions with all the indexes of Tab and Aabc=#\n")
    
    write(outfile, "\nTab_functions = Array{Function, 2}(undef, 4,4)\n")
    write(outfile, "\nAabc_functions = Array{Function, 3}(undef, 4,4,4)\n")
    Tab_functions = Array{Function, 2}(undef, 4,4)
    Aabc_functions = Array{Function, 3}(undef, 4,4,4)
    for i in 1:4
        for j in 1:4
            write(outfile, "\nTab_functions[$i,$j] = " * string(Symbolics.build_function(Tab_simple[i,j], ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ,g_ab)))
            write(outfile, "\nTab_functions$i$j = " * string(Symbolics.build_function(Tab_simple[i,j], ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ,g_ab)))
        end
    end

    for i in 1:4
        for j in 1:4
            for k in 1:4
                write(outfile, "\nAabc_functions[$i,$j,$k] = " * string(Symbolics.build_function(Aabc_simple[i,j,k], ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ,g_ab)))
                write(outfile, "\nAabc_functions$i$j$k = " * string(Symbolics.build_function(Aabc_simple[i,j,k], ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ, Mᵃᵇ, gᵃᵇ,g_ab)))
            end
        end
    end

    write(outfile, "\n function Tab_fun!(T, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab)\n")
    for i in 1:4
        for j in i:4
            write(outfile, "T[$i,$j] = Tab_functions$i$j(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab)\n")
            (i!=j) && write(outfile, "T[$j,$i] = T[$i,$j]\n")
        end
    end
    write(outfile, "end\n")


    write(outfile, "\n function Aabc_fun!(A, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab)\n")
    for i in 1:4
        for j in i:4
            for k in j:4
            write(outfile, "A[$i,$j,$k] = Aabc_functions$i$j$k(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab)\n")
            (i!=j) && write(outfile, "A[$j,$i,$k] = A[$i,$j,$k]\n")
            (i!=k) && write(outfile, "A[$k,$j,$i] = A[$i,$j,$k]\n")
            (j!=k) && write(outfile, "A[$i,$k,$j] = A[$i,$j,$k]\n")
            ((i!=j) && (j!=k)) && write(outfile, "A[$k,$i,$j] = A[$i,$j,$k]\n")
            ((i!=j) && (j!=k)) && write(outfile, "A[$j,$k,$i] = A[$i,$j,$k]\n")
            end
        end
    end
    write(outfile, "end\n")
            

    
end