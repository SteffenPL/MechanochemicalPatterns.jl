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

preprocess(x::Vector{Union{Float64, Int64}}) = preprocess(Float64.(x))

preprocess(x::Vector{Float64}) = SVector{length(x), Float64}(x)
preprocess(x::Vector{Int64}) = SVector{length(x), Int64}(x)
preprocess(x::Vector) = tuple(x...)

function preprocess(x::@NamedTuple{center::SVector{N,Float64}, size::SVector{N,Float64}}) where {N}
    c = SVector{N, Float64}(x.center)
    s = SVector{N, Float64}(x.size)
    min = c - 0.5 .* s
    max = c + 0.5 .* s
    inv_size = 1 ./ s
    return (center = c, size = s, min = min, max = max, inv_size = inv_size)
end

function recursive_namedtuple(d::Dict, preprocess = preprocess)
    for (k, v) in d
        if isa(v, Dict)
            d[k] = recursive_namedtuple(v, preprocess)
        end
    end
    return NamedTuple(Symbol(k) => preprocess(v) for (k,v) in d)
end