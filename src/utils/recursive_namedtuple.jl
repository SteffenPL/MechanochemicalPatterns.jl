preprocess(x) = x

function preprocess(x::String) 
    if startswith(x, "julia:") 
        try
            expr = replace(x[7:end], "\n" => ";")
            eval(Meta.parse(expr))
        catch e
            println("Error in parsing: ", x[7:end])
            showerror(stdout, e)
            @error "Unable to compile Julia expressions in parameter file."
        end
    else 
        x
    end
end

function recursive_namedtuple(d::Dict, preprocess = preprocess)
    for (k, v) in d
        if isa(v, Dict)
            d[k] = recursive_namedtuple(v, preprocess)
        end
    end
    return NamedTuple(Symbol(k) => preprocess(v) for (k,v) in d)
end