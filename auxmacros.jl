macro trace_free_derivatives(boolean)
    str = ""
    if boolean
        str = "derivatives-tracefree.jl"
    else
        str = "derivatives.jl"
    end
    quote
        include($str)
    end
end
