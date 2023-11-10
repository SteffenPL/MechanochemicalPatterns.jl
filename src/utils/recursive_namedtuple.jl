function recursive_namedtuple(d::Dict)
    for (k, v) in d
        if isa(v, Dict)
            d[k] = recursive_namedtuple(v)
        end
    end
    return NamedTuple( Symbol(k) => v for (k,v) in d)
end