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

function preprocess(x::NamedTuple)
    if isempty(x) || !hasproperty(first(x), :ID)
        return x
    else
        data = sort(zip(keys(x), values(x)), by = x -> x[2].ID)
        return NamedTuple(key => preprocess(value) for (key, value) in data)
    end
end

preprocess(x::Vector{Union{Float64, Int64}}) = preprocess(Float64.(x))
preprocess(x::Vector{Float64}) = SVector{length(x), Float64}(x)
preprocess(x::Vector{Int64}) = SVector{length(x), Int64}(x)
preprocess(x::Vector) = Tuple(preprocess.(x))

preprocess(x::@NamedTuple{periodic::T}) where {T} = PeriodicBoundary(0.0)
preprocess(x::@NamedTuple{neumann::Float64}) = NeumannBoundary(x.neumann)
preprocess(x::@NamedTuple{dirichlet::Float64}) = DirichletBoundary(x.dirichlet)

function preprocess(x::@NamedTuple{center::SVector{N,Float64}, size::SVector{N,Float64}}) where {N}
    c = SVector{N, Float64}(x.center)
    s = SVector{N, Float64}(x.size)
    min = c - 0.5 .* s
    max = c + 0.5 .* s
    inv_size = 1 ./ s
    return (center = c, size = s, min = min, max = max, inv_size = inv_size)
end

function preprocess(d::Dict)
    return preprocess(NamedTuple(Symbol(k) => preprocess(v) for (k,v) in d))
end

recursive_namedtuple(d) = preprocess(d)
