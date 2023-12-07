eval_param(p, x, args...) = x 

function eval_param(p, x::Function, args...)
    return x(args...)
end

function eval_param(p, x::@NamedTuple{uniform::Vector}) 
    # Example: {uniform = [-1, 2]}
    return x.uniform[1] + rand() * (x.uniform[2] - x.uniform[1])
end

function eval_param(p, x::@NamedTuple{uniform::String}) 
    # Example: {uniform = "domain"}
    if x.uniform == "domain"
        return p.env.domain.min + rand(typeof(p.env.domain.max)) * (p.env.domain.max - p.env.domain.min)
    else
        @error "Unknown uniform distribution: $(x.uniform)"
    end
end

function eval_param(p, x::@NamedTuple{normal::@NamedTuple{mean::T1, std::T2}}) where {T1 <: Number, T2} 
    # Example: {normal = {mean = 2.0, std = 1.0}}
    return randn() * x.normal.std + x.normal.mean
end

function eval_param(p, x::@NamedTuple{normal::@NamedTuple{mean::T1, std::T2}}) where {T1 <: AbstractArray,T2} 
    # Example: {normal = {mean = [0,0,75], std = [100,100,15]}}
    r = randn(length(x.normal.mean))
    @. r = r * x.normal.std + x.normal.mean
    return r
end
X = 10
function eval_param(p, nt::@NamedTuple{custom_distr::@NamedTuple{domain_fnc::T1, domain_bounds::T2, distr::T3}}, init_vector = MVector{3,Float64}(zeros(3))) where {T1 <: Function, T2, T3 <: Function} 
    dom = nt.custom_distr.domain_bounds

    # generate random points inside the bounds and with fnc(x,y,z) > 0 
    x = init_vector
    for i in 1:1000  # to avoid infinite loops
        # new random point
        nt.custom_distr.distr(x, p)

        # indomain 
        if nt.custom_distr.domain_fnc(x, p) > 0
            if all(2 * abs.(x - dom.center) .< dom.size)  
                return x
            end
        end
    end

    @error "No point in custom domain found."
    return NaN
end


function eval_param(p, nt::@NamedTuple{custom_distr::@NamedTuple{domain_fnc::T1, domain_bounds::T2}}, init_vector = MVector{3,Float64}(zeros(3))) where {T1 <: Function, T2} 

    distr = function (x, dom)
        rand!(x)
        @. x = x * dom.size + dom.center - 0.5 * dom.size
    end

    return eval_param(p, (;custom_distr = (nt.custom_distr..., distr = distr)), init_vector)
end




function get_param(p, cell_type, sym, default = 0.0)
    ct = p.cells.types[cell_type]
    if hasproperty(ct, sym) 
        return ct[sym]
    elseif hasproperty(p.cells, sym)
        return p.cells[sym]
    else
        return default
    end
end
