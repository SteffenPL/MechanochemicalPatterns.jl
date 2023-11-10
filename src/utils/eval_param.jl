eval_param(x) = x 

function eval_param(x::@NamedTuple{uniform::Vector}) 
    # Example: {uniform = [-1, 2]}
    return x.uniform[1] + rand() * (x.uniform[2] - x.uniform[1])
end

function eval_param(x::@NamedTuple{normal::@NamedTuple{mean::T1, std::T2}}) where {T1 <: Number, T2} 
    # Example: {normal = {mean = 2.0, std = 1.0}}
    return randn() * x.normal.std + x.normal.mean
end

function eval_param(x::@NamedTuple{normal::@NamedTuple{mean::T1, std::T2}}) where {T1 <: AbstractArray,T2} 
    # Example: {normal = {mean = [0,0,75], std = [100,100,15]}}
    r = randn(length(x.normal.mean))
    @. r = r * x.normal.std + x.normal.mean
    return r
end


function eval_param(x::@NamedTuple{uniform_domain::@NamedTuple{domain::T1, fnc::T2}}) where {T1, T2} 
    fnc = eval(Meta.parse("(x,y,z) -> $(x.uniform_domain.fnc)"))
    dom = x.uniform_domain.domain
    d = length(dom.center)

    # generate random points inside the bounds and with fnc(x,y,z) > 0
    r = zeros(Float64, d)
    for i in 1:1000  # to avoid infinite loops
        r .= dom.center .+ (-0.5 .+ rand(d)) .* dom.size
        if Base.invokelatest(fnc, r...) > 0
            return r
        end
    end

    @error "No point founds with (x,y,z) -> $(x.uniform_domain.fnc) > 0 and inside the bounds $(dom)."
    return NaN
end
