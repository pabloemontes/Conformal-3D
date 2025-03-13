if (@isdefined debug) == false
    println("Debug was not defined, defining global variable debug = false")
    debug = false
end
#===========================================================================
Derivadas de μ
∂μ/∂ξ_a = 2ξᵃ
============================================================================#
using Symbolics
@register μ_fun(x)
for i in 1:D
    return_expression = :(2*$(ξᵃ_fun_symbols[i])(args...))
    debug == true && println(return_expression)
    @eval Symbolics.derivative(::typeof(μ_fun), args::NTuple{$D,Any}, ::Val{$i}) = $return_expression
end

#===========================================================================#

#===========================================================================
Derivadas de ν
∂ν/∂ξ_a = 2lᵃ
∂ν/∂ξ_ab = ξᵃξᵇ
============================================================================#
#Derivadas de ν respecto de ξ_a
for i in 1:D
    return_expression = :(2*$(lᵃ_fun_symbols[i])(args...))
    (debug == true) && println(return_expression)
    @eval Symbolics.derivative(::typeof(ν_fun), args::NTuple{$(D+D*(D+1)÷2),Any}, ::Val{$i}) = $return_expression
end

#Derivadas de ν respecto de ξ_ab
begin
    local k = D
    for i in 1:D
        for j in i:D
            k = k+1
            return_expression = :(($(ξᵃ_fun_symbols[i]))(args[1:$D]...)*($(ξᵃ_fun_symbols[j]))(args[1:$D]...))
            if i != j
                return_expression = Expr(:call, *, 2, return_expression)
            end
            
            (debug == true) && println(return_expression)
            @eval Symbolics.derivative(::typeof(ν_fun), args::NTuple{$(D+D*(D+1)÷2),Any}, ::Val{$k}) = $return_expression
        end
    end
end

#============================================================================#

### Derivadas de lᵃ

#lᵃ = ξᵃᵇξ_b

#∂lᵃ/∂ξ_b = ξᵃᵇ

for i in 1:D
    for j in 1:D
        return_expression = :($(ξᵃᵇ_fun_symbols[i,j])(args[$(D+1):end]...))
        (debug == true) && println(return_expression)
        @eval Symbolics.derivative(::typeof($(lᵃ_fun_symbols[i])), args::NTuple{$(D+D*(D+1)÷2),Any}, ::Val{$j}) = $return_expression
    end
end

#∂lᵃ/∂ξ_bc = (1/2)(ξᶜgᵃᵇ + ξᵇgᵃᶜ)

#Derivadas sin traza de lᵃ respecto de ξ_bc
for i in 1:D
    n = D
    for j in 1:D
        for k in j:D
            n = n+1
            ξᵃ_i = ξᵃ_fun_symbols[i]
            ξᵃ_j = ξᵃ_fun_symbols[j]
            ξᵃ_k = ξᵃ_fun_symbols[k]
            return_expression = :(0.5*($(ξᵃ_k)(args[1:$D]...)*$(g.uu[i,j]) + $(ξᵃ_j)(args[1:$D]...)*$(g.uu[i,k])))
            if j != k
                return_expression = Expr(:call, *, 2, return_expression)
            end
            (debug == true) && println(return_expression)
            @eval Symbolics.derivative(::typeof($(lᵃ_fun_symbols[i])), args::NTuple{$(D+D*(D+1)÷2),Any}, ::Val{$n}) = $return_expression
        end
    end
end

#================================================================================#
#Derivadas de ηᵃ

#ηᵃ = ξᵃᵇl_b

#Derivadas de ηᵃ respecto de ξ_b

#∂ηᵃ/∂ξ_b = ξᵃᶜξᵇᵈg_cd = Mᵃᵇ

for i in 1:D
    for j in 1:D
        #return_expression = Expr(:call, :+, [:($(ξᵃᵇ_fun_symbols[i,k])(args[$(D+1):end]...)*$(ξᵃᵇ_fun_symbols[j,l])(args[$(D+1):end]...)*$(g.dd[k,l]) ) for k in 1:D, l in 1:D]...)
        return_expression = :($(Mᵃᵇ_fun_symbols[i,j])(args[$(D+1):end]...))
        (debug == true) && println(return_expression)
        @eval Symbolics.derivative(::typeof($(ηᵃ_fun_symbols[i])), args::NTuple{$(D+D*(D+1)÷2),Any}, ::Val{$j}) = $return_expression
    end
end

#Derivadas de ηᵃ respecto de ξ_bc

#∂ηᵃ/∂ξ_bc = (1/2)*(gᵃᵇlᶜ+gᵃᶜlᵇ+ξᵃᵇξᶜ+ξᵃᶜξᵇ)


for i in 1:D
    n = D
    for j in 1:D
        for k in j:D
            n = n+1
            ξᵃ_i = :($(ξᵃ_fun_symbols[i])(args[1:$D]...))
            ξᵃ_j = :($(ξᵃ_fun_symbols[j])(args[1:$D]...))
            ξᵃ_k = :($(ξᵃ_fun_symbols[k])(args[1:$D]...))
            lᵃ_i = :($(lᵃ_fun_symbols[i])(args...))
            lᵃ_j = :($(lᵃ_fun_symbols[j])(args...))
            lᵃ_k = :($(lᵃ_fun_symbols[k])(args...))
            ξᵃᵇ_ij = :($(ξᵃᵇ_fun_symbols[i,j])(args[$(D+1):end]...))
            ξᵃᵇ_ik = :($(ξᵃᵇ_fun_symbols[i,k])(args[$(D+1):end]...))
            
            return_expression = :(0.5*($(g.uu[i,j])*$(lᵃ_k) + $(g.uu[i,k])*$(lᵃ_j) + $(ξᵃᵇ_ij)*$(ξᵃ_k) + $(ξᵃᵇ_ik)*$(ξᵃ_j)))
            if j != k
                return_expression = Expr(:call, *, 2, return_expression)
            end
            (debug == true) && println(return_expression)
            @eval Symbolics.derivative(::typeof($(ηᵃ_fun_symbols[i])), args::NTuple{$(D+D*(D+1)÷2),Any}, ::Val{$n}) = $return_expression
        end
    end
end

#===================================================================================#
#Derivadas de ξᵃ

#∂ξᵃ/∂ξ_b = gᵃᵇ

for i in 1:D
    for j in 1:D
        return_expression = :($(g.uu[i,j]))
        (debug == true) && println(return_expression)
        @eval Symbolics.derivative(::typeof($(ξᵃ_fun_symbols[i])), args::NTuple{$D, Any}, ::Val{$j}) = $return_expression
    end
end

#===================================================================================#
#Derivadas de ξᵃᵇ

#∂ξᵃᵇ/∂ξ_cd = (1/2)*(gᵃᶜgᵇᵈ + gᵃᵈgᵇᶜ)

for i in 1:D
    for j in i:D
        m_00000001 = 0
        for k in 1:D
            for l in k:D
                m_00000001 = m_00000001+1
                return_expression = :(0.5*($(g.uu[i,k])*$(g.uu[j,l]) + $(g.uu[i,l])*$(g.uu[j,k])))
                if k != l
                    return_expression = Expr(:call, *, 2, return_expression)
                end
                (debug == true) && println(return_expression)
                @eval Symbolics.derivative(::typeof($(ξᵃᵇ_fun_symbols[i,j])), args::NTuple{$(D*(D+1)÷2), Any}, ::Val{$m_00000001}) = $return_expression
            end
        end
    end
end


#===================================================================================#
#Derivadas de Mᵃᵇ

#Mᵃᵇ = ξᵃᵉξᵇ_e

#∂Mᵃᵇ/∂ξ_cd = (1/2)(gᵃᶜξᵇᵈ + gᵃᵈξᵇᶜ + gᵇᶜξᵃᵈ + gᵇᵈξᵃᶜ)

for i in 1:D
    for j in i:D
        m_00000001 = 0
        for k in 1:D
            for l in k:D
                m_00000001 = m_00000001+1

                ξᵃᵇ_ij = :($(ξᵃᵇ_fun_symbols[i,j])(args...))
                ξᵃᵇ_ik = :($(ξᵃᵇ_fun_symbols[i,k])(args...))
                ξᵃᵇ_il = :($(ξᵃᵇ_fun_symbols[i,l])(args...))
                ξᵃᵇ_jk = :($(ξᵃᵇ_fun_symbols[j,k])(args...))
                ξᵃᵇ_jl = :($(ξᵃᵇ_fun_symbols[j,l])(args...))
                return_expression = :(0.5*($(g.uu[i,k])*$ξᵃᵇ_jl + $(g.uu[i,l])*$ξᵃᵇ_jk + $(g.uu[j,k])*$ξᵃᵇ_il + $(g.uu[j,l])*$ξᵃᵇ_ik)) 
                if k != l
                    return_expression = Expr(:call, *, 2, return_expression)
                end
                (debug == true) && println(return_expression)
                @eval Symbolics.derivative(::typeof($(Mᵃᵇ_fun_symbols[i,j])), args::NTuple{$(D*(D+1)÷2), Any}, ::Val{$m_00000001}) = $return_expression
            end
        end
    end
end



#====================================================================================#
#Derivadas de ψ_1

#ψ_1 = ξᵃᵇξ_ab

#∂ψ_2/∂ξ_ab = 2ξᵃᵇ  (which should already be trace free)

#Derivadas de ψ_1 respecto de ξ_ab (already trace free)
begin
    local k = 0
    for i in 1:D
        for j in i:D
            k = k+1
            return_expression = :(2*$(ξᵃᵇ_fun_symbols[i,j])(args...))
            if i != j
                return_expression = Expr(:call, *, 2, return_expression)
            end
            (debug == true) && println(return_expression)
            @eval Symbolics.derivative(::typeof(ψ_1_fun), args::NTuple{$(D*(D+1)÷2), Any}, ::Val{$k}) = $return_expression
        end
    end
end


#====================================================================================#
#Derivadas de ψ_2

#ψ_2 = lᵃl_a    

#∂ψ_2/∂ξ_a = 2ηᵃ  

#Derivadas de ψ_2 respecto de ξ_a
for i in 1:D
    return_expression = :(2*$(ηᵃ_fun_symbols[i])(args...))
    (debug == true) && println(return_expression)
    @eval Symbolics.derivative(::typeof(ψ_2_fun), args::NTuple{$(D+D*(D+1)÷2), Any}, ::Val{$(i)}) = $return_expression
end

#∂ψ_2/∂ξ_ab = ξᵃlᵇ + ξᵇlᵃ


#Derivadas de ψ_2 respecto de ξ_ab sin traza
begin
    local k = D
    for i in 1:D
        for j in i:D
            k = k+1
            return_expression = :((($(lᵃ_fun_symbols[i])(args...)) * ($(ξᵃ_fun_symbols[j])(args[1:$D]...)) + ($(lᵃ_fun_symbols[j])(args...) * $(ξᵃ_fun_symbols[i])(args[1:$D]...))))
            if i != j
                return_expression = Expr(:call, *, 2, return_expression)
            end
            (debug == true) && println(return_expression)
            @eval Symbolics.derivative(::typeof(ψ_2_fun), args::NTuple{$(D+D*(D+1)÷2), Any}, ::Val{$(k)}) = $return_expression
        end
    end
end
