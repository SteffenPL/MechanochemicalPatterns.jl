function dist²(a,b)  # backup for non static vectors
    d = zero(eltype(a))
    for i in eachindex(a)
        d += (a[i] - b[i])^2
    end
    return d
end

dist²(a::SVector,b::SVector) = sum(x -> x^2, a-b)
dist(a,b) = sqrt(dist²(a,b))

function random_direction(Dim)
    v = SVector{Dim,Float64}(randn(Dim))
    return v ./ norm(v)
end