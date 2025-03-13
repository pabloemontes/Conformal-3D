macro my_covector(var)
    varname = string(var)
    var_a_strings = Symbol(string(var)*"_a_strings")
    var_a_symbols = Symbol(string(var)*"_a_symbols")
    var_a_fun_symbols = Symbol(string(var)*"_a_fun_symbols")
    var_a = Symbol(string(var)*"_a") 
    var_i_strings = [string(var)*"_$i" for i in 0:(D-1)]
    var_i_symbols = Symbol.(var_i_strings)
    var_i_fun_symbols = Symbol.(var_i_strings.*"_fun")
    endmessage = "$var_a, $var_a_symbols, $var_a_strings, $var_a_fun_symbols and $(join(var_i_symbols, ", ")) created"
    quote
        $(esc(var_a_strings)) = $(var_i_strings) #We define D strings of the form ξ_i, with i ∈ 0:(D-1)
        $(esc(var_a_symbols)) = $(var_i_symbols)      #We create a vector of symbols created from those strings
        @eval @variables $(var_i_symbols...)
        $(esc(var_a)) = [$(var_i_symbols...)]
        $(esc(var_a_fun_symbols)) = $(var_i_fun_symbols)
        debug && println($endmessage)
    end
end

macro my_vector(var)
    varname = string(var)
    var_a_strings = Symbol(string(var)*"ᵃ_strings")
    var_a_symbols = Symbol(string(var)*"ᵃ_symbols")
    var_a_fun_symbols = Symbol(string(var)*"ᵃ_fun_symbols")
    var_a = Symbol(string(var)*"ᵃ") 
    var_i_strings = [string(var)*superscriptnumber(i) for i in 0:(D-1)]
    var_i_symbols = Symbol.(var_i_strings)
    var_i_fun_symbols = Symbol.(var_i_strings.*"_fun")
    endmessage = "$var_a, $var_a_symbols, $var_a_fun_symbols $var_a_strings and $(join(var_i_symbols, ", ")) created"
    quote
        $(esc(var_a_strings)) = $(var_i_strings) #We define D strings of the form ξ_i, with i ∈ 0:(D-1)
        $(esc(var_a_symbols)) = $(var_i_symbols)      #We create a vector of symbols created from those strings
        $(esc(var_a_fun_symbols)) = $(var_i_fun_symbols)
        @eval @variables $(var_i_symbols...)
        $(esc(var_a)) = [$(var_i_symbols...)]
        debug && println($endmessage)
    end
end

macro my_tensor_upup(var)
    varname = string(var)
    var_ab_strings = Symbol(string(var)*"ᵃᵇ_strings")
    var_ab_symbols = Symbol(string(var)*"ᵃᵇ_symbols")
    var_ab_fun_symbols = Symbol(string(var)*"ᵃᵇ_fun_symbols")
    var_ab = Symbol(string(var)*"ᵃᵇ") 
    var_ij_strings = [string(var)*superscriptnumber(i)*superscriptnumber(j) for i in 0:(D-1), j in 0:(D-1)]
    var_ij_symbols = Symbol.(var_ij_strings)
    var_ij_fun_symbols = Symbol.(var_ij_strings.*"_fun")
    endmessage = "$var_ab, $var_ab_symbols, $var_ab_strings, $var_ab_fun_symbols, and $(join(var_ij_symbols, ", ")) created"
    quote
        $(esc(var_ab_strings)) = $(var_ij_strings) #We define D strings of the form ξ_i, with i ∈ 0:(D-1)
        $(esc(var_ab_symbols)) = $(var_ij_symbols)      #We create a vector of symbols created from those strings
        $(esc(var_ab_fun_symbols)) = $(var_ij_fun_symbols)
        @eval @variables $(var_ij_symbols...)
        $(esc(var_ab)) = reshape([$(var_ij_symbols...)], D, D)
        debug && println($endmessage)
    end
end


macro my_tensor_updn(var)
    varname = string(var)
    var_ab_strings = Symbol(string(var)*"ᵃ_b_strings")
    var_ab_symbols = Symbol(string(var)*"ᵃ_b_symbols")
    var_ab_fun_symbols = Symbol(string(var)*"ᵃ_b_fun_symbols")
    var_ab = Symbol(string(var)*"ᵃ_b") 
    var_ij_strings = [string(var)*superscriptnumber(i)*"_$(j)" for i in 0:(D-1), j in 0:(D-1)]
    var_ij_symbols = Symbol.(var_ij_strings)
    var_ij_fun_symbols = Symbol.(var_ij_strings.*"_fun")
    endmessage = "$var_ab, $var_ab_symbols, $var_ab_strings, $var_ab_fun_symbols, and $(join(var_ij_symbols, ", ")) created"
    quote
        $(esc(var_ab_strings)) = $(var_ij_strings) #We define D strings of the form ξ_i, with i ∈ 0:(D-1)
        $(esc(var_ab_symbols)) = $(var_ij_symbols)      #We create a vector of symbols created from those strings
        $(esc(var_ab_fun_symbols)) = $(var_ij_fun_symbols)
        @eval @variables $(var_ij_symbols...)
        $(esc(var_ab)) = reshape([$(var_ij_symbols...)], D, D)
        debug && println($endmessage)
    end
end


macro my_tensor_dndn(var)
    varname = string(var)
    var_ab_strings = Symbol(string(var)*"_ab_strings")
    var_ab_symbols = Symbol(string(var)*"_ab_symbols")
    var_ab_fun_symbols = Symbol(string(var)*"_ab_fun_symbols")
    var_ab = Symbol(string(var)*"_ab") 
    var_ij_strings = [string(var)*"_$(i)$(j)" for i in 0:(D-1), j in 0:(D-1)]
    var_ij_symbols = Symbol.(var_ij_strings)
    var_ij_fun_symbols = Symbol.(var_ij_strings.*"_fun")
    endmessage = "$var_ab, $var_ab_symbols, $var_ab_strings, $var_ab_fun_symbols and $(join(var_ij_symbols, ", ")) created"
    quote
        $(esc(var_ab_strings)) = $(var_ij_strings) #We define D strings of the form ξ_i, with i ∈ 0:(D-1)
        $(esc(var_ab_symbols)) = $(var_ij_symbols)      #We create a vector of symbols created from those strings
        $(esc(var_ab_fun_symbols)) = $(var_ij_fun_symbols)
        @eval @variables $(var_ij_symbols...)
        $(esc(var_ab)) = reshape([$(var_ij_symbols...)], D, D)
        debug && println($endmessage)
    end
end

macro my_symmetric_tensor_upup(var)
    varname = string(var)
    var_ab_strings = Symbol(string(var)*"ᵃᵇ_strings")
    var_ab_symbols = Symbol(string(var)*"ᵃᵇ_symbols")
    var_ab_fun_symbols = Symbol(string(var)*"ᵃᵇ_fun_symbols")
    var_ab = Symbol(string(var)*"ᵃᵇ") 
    var_ij_strings = [string(var)*superscriptnumber(min(i,j))*superscriptnumber(max(i,j)) for i in 0:(D-1), j in 0:(D-1)]
    var_ij_symbols = Symbol.(var_ij_strings)
    var_ij_fun_symbols = Symbol.(var_ij_strings.*"_fun")
    endmessage = "$var_ab, $var_ab_symbols, $var_ab_strings, $var_ab_fun_symbols, and $(join(var_ij_symbols, ", ")) created"
    quote
        $(esc(var_ab_strings)) = $(var_ij_strings) #We define D strings of the form ξ_i, with i ∈ 0:(D-1)
        $(esc(var_ab_symbols)) = $(var_ij_symbols)      #We create a vector of symbols created from those strings
        $(esc(var_ab_fun_symbols)) = $(var_ij_fun_symbols)      #We create a vector of symbols created from those strings
        @eval @variables $(var_ij_symbols...)
        $(esc(var_ab)) = Symmetric(reshape([$(var_ij_symbols...)], D, D))
        debug && println($endmessage)
    end
end
macro my_symmetric_tensor_dndn(var)
    varname = string(var)
    var_ab_strings = Symbol(string(var)*"_ab_strings")
    var_ab_symbols = Symbol(string(var)*"_ab_symbols")
    var_ab_fun_symbols = Symbol(string(var)*"_ab_fun_symbols")
    var_ab = Symbol(string(var)*"_ab") 
    var_ij_strings = [string(var)*"_$(min(i,j))$(max(i,j))" for i in 0:(D-1), j in 0:(D-1)]
    var_ij_symbols = Symbol.(var_ij_strings)
    var_ij_fun_symbols = Symbol.(var_ij_strings.*"_fun")
    endmessage = "$var_ab, $var_ab_symbols, $var_ab_strings, $var_ab_fun_symbols and $(join(var_ij_symbols, ", ")) created"
    quote
        $(esc(var_ab_strings)) = $(var_ij_strings) #We define D strings of the form ξ_i, with i ∈ 0:(D-1)
        $(esc(var_ab_symbols)) = $(var_ij_symbols)      #We create a vector of symbols created from those strings
        $(esc(var_ab_fun_symbols)) = $(var_ij_fun_symbols)
        @eval @variables $(var_ij_symbols...)
        $(esc(var_ab)) = Symmetric(reshape([$(var_ij_symbols...)], D, D))
        debug && println($endmessage)
    end
end


macro my_symmetric_tensor_updn(var)
    varname = string(var)
    var_ab_strings = Symbol(string(var)*"ᵃ_b_strings")
    var_ab_symbols = Symbol(string(var)*"ᵃ_b_symbols")
    var_ab_fun_symbols = Symbol(string(var)*"ᵃ_b_fun_symbols")
    var_ab = Symbol(string(var)*"ᵃ_b") 
    var_ij_strings = [string(var)*superscriptnumber(i)*"_$(j)" for i in 0:(D-1), j in 0:(D-1)]
    var_ij_symbols = Symbol.(var_ij_strings)
    var_ij_fun_symbols = Symbol.(var_ij_strings.*"_fun")
    endmessage = "$var_ab, $var_ab_symbols, $var_ab_strings, $var_ab_fun_symbols, and $(join(var_ij_symbols, ", ")) created"
    quote
        $(esc(var_ab_strings)) = $(var_ij_strings) #We define D strings of the form ξ_i, with i ∈ 0:(D-1)
        $(esc(var_ab_symbols)) = $(var_ij_symbols)      #We create a vector of symbols created from those strings
        $(esc(var_ab_fun_symbols)) = $(var_ij_fun_symbols)
        @eval @variables $(var_ij_symbols...)
        $(esc(var_ab)) = Symmetric(reshape([$(var_ij_symbols...)], D, D))
        debug && println($endmessage)
    end
end
