function recursive_namedtuple(d::Dict, process = x -> x)
    for (k, v) in d
        if isa(v, Dict)
            d[k] = recursive_namedtuple(v, process)
        end
    end
    return NamedTuple( Symbol(k) => process(v) for (k,v) in d)
end