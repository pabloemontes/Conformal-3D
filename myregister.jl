macro my_register_symbolic(expr, define_promotion = true, Ts = [])
    if expr.head === :(::)
        ret_type = expr.args[2]
        expr = expr.args[1]
    else
        ret_type = Real
    end

    @assert expr.head === :call

    f = expr.args[1]
    args = expr.args[2:end]

    if f isa Expr && f.head == :(::)
        @assert length(f.args) == 2
    end

    types = map(args) do x
        if x isa Symbol
            :(($Real, $Symbolics.wrapper_type($Real), $SymbolicUtils.Symbolic{<:$Real}))
        elseif Meta.isexpr(x, :(::))
            T = x.args[2]
            :($SymbolicUtils.has_symwrapper($T) ?
              ($T, $SymbolicUtils.Symbolic{<:$T}, $Symbolics.wrapper_type($T),) :
              ($T, $SymbolicUtils.Symbolic{<:$T}))
        else
            error("Invalid argument format $x")
        end
    end

    eval_method = :(@eval function $f($(Expr(:$, :(var"##_register_macro_s"...))),)
                        args = [$(Expr(:$, :(var"##_register_macro_s_syms"...)))]
                        unwrapped_args = map($Symbolics.unwrap, args)
                        res = if !any(x->$SymbolicUtils.issym(x) || $SymbolicUtils.istree(x), unwrapped_args)
                            $f(unwrapped_args...)
                        else
                            $SymbolicUtils.Term{$ret_type}($f, unwrapped_args)
                        end
                        if typeof.(args) == typeof.(unwrapped_args)
                            return res
                        else
                            return $Symbolics.wrap(res)
                        end
                    end)
    verbose = false
    mod, fname = f isa Expr && f.head == :(.) ? f.args : (:(@__MODULE__), QuoteNode(f))
    Ts = Symbol("##__Ts")
    ftype = if f isa Expr && f.head == :(::)
        f.args[end]
    else
        :($typeof($f))
    end
    quote
        $Ts = [Tuple{x...} for x in [typelist for typelist in zip($(types...))]
                if any(x->x <: $SymbolicUtils.Symbolic || Symbolics.is_wrapper_type(x), x)]
        if $verbose
            println("Candidates")
            map(println, $Ts)
        end

        for sig in $Ts
            var"##_register_macro_s" = map(((i,T,),)->Expr(:(::), Symbol("arg", i), T), enumerate(sig.parameters))
            var"##_register_macro_s_syms" = map(x->x.args[1], var"##_register_macro_s")
            $eval_method
        end
        if $define_promotion
            (::$typeof($SymbolicUtils.promote_symtype))(::$ftype, args...) = $ret_type
        end
    end |> esc
end