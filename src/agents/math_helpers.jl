function dist²(a,b)  # backup for non static vectors
    d = zero(eltype(a))
    for i in eachindex(a)
        d += (a[i] - b[i])^2
    end
    return d
end

dist²(a::SVector,b::SVector) = sum(x -> x^2, a-b)
dist(a,b) = sqrt(dist²(a,b))

dist²(a::SVector) = sum(x -> x^2, a)
dist(a::SVector) = sqrt(dist²(a))

function random_direction(Dim)
    v = SVector{Dim,Float64}(randn(Dim))
    return v ./ norm(v)
end


@inline function safe_normalize(x, d = x)
    n = norm(x)
    if n > 0
        return x/n
    else
        return d
    end
end