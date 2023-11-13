preprocess(x) = x
preprocess(x::String) = startswith(x, "julia:") ? eval(Meta.parse(x[7:end])) : x

function recursive_namedtuple(d::Dict, preprocess = preprocess)
    for (k, v) in d
        if isa(v, Dict)
            d[k] = recursive_namedtuple(v, preprocess)
        end
    end
    return NamedTuple(Symbol(k) => preprocess(v) for (k,v) in d)
end