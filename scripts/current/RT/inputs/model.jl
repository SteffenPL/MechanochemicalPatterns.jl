
function pde!(du, u::ArrayPartition{Float64, Tuple{MT,MT,MT}}, p, t) where {MT <: Array}
    tmp = p.tmp 
    (FGF, Wnt, Dom) = u.x
    (dFGF, dWnt, dDom) = du.x

    dFGF .= 0.0
    dWnt .= 0.0
    dDom .= 0.0

    tmp.x[1] .=  

    laplace!(dFGF, FGF, p.FGF.D, p.dV)
    laplace!(dWnt, Wnt, p.Wnt.D, p.dV)
    laplace!(dDom, Dom, p.Dom.D, p.dV)

    # decay
    @. dFGF -= p.FGF.decay * FGF
    @. dWnt -= p.Wnt.decay * Wnt
    @. dDom -= p.Dom.decay * Dom

    return nothing
end


function pde!(du, u::ArrayPartition{Float64, Tuple{MT}}, p, t) where {MT <: Array}
    st = p.signals.types

    du .= 0.0
    
    laplace!(du.x[1], u.x[1], st[1].D, p.dV)
    
    # decay
    @. du.x[1] -= st[1].decay * u.x[1]

    return nothing
end


function pde!(du, u::ArrayPartition{Float64, Tuple{MT, MT}}, p, t) where {MT <: Array}
    st = p.signals.types

    du .= 0.0
    
    laplace!(du.x[1], u.x[1], st[1].D, p.dV)
    laplace!(du.x[2], u.x[2], st[2].D, p.dV)
    
    # decay
    @. du.x[1] -= st[1].decay * u.x[1]
    @. du.x[2] -= st[2].decay * u.x[2]

    return nothing
end

hill(x, k, n) = x^n / (k^n + x^n)

function pde!(du, u::ArrayPartition{Float64, Tuple{MT, MT, MT}}, p, t) where {MT <: Array}
    st = p.signals.types
    tmp = p.tmp 

    du .= 0.0
    
    # diffusion factor 
    f(x) = x > 0.5 ? 0.1 : 1.0

    @. tmp.x[1] = f(u.x[1])
    Dom = tmp.x[1]

    # Domain, Sox9, Wnt, BMP, FGF8
    laplace!(du.x[1], u.x[1], st[1].D, p.dV, boundaries = p.boundaries[1])
    laplace!(du.x[2], u.x[2], st[2].D, p.dV, boundaries = p.boundaries[2], factor = Dom)
    laplace!(du.x[3], u.x[3], st[3].D, p.dV, boundaries = p.boundaries[3])
    
    # decay
    @. du.x[1] -= st[1].decay * u.x[1]
    @. du.x[2] -= st[2].decay * u.x[2]
    @. du.x[3] -= st[3].decay * u.x[3]

    return nothing
end



function pde!(du, u::ArrayPartition{Float64, Tuple{MT, MT, MT, MT, MT}}, p, t) where {MT <: Array}
    st = p.signals.types
    tmp = p.tmp 

    (;k2, k3, k4, k1FGF, k5, k7, k2FGF, k9) = p.signals

    # Domain, FGF8, Wnt, BMP, Sox9
    
    du .= 0.0
    
    # diffusion factor 
    f(x) = x > 0.5 ? 0.2 : 1.0

    @. tmp.x[1] = f(u.x[1])
    Dom = tmp.x[1]

    # Domain, Sox9, Wnt, BMP, FGF8
    laplace!(du.x[1], u.x[1], st[1].D, p.dV, boundaries = p.boundaries[1])
    laplace!(du.x[2], u.x[2], st[2].D, p.dV, boundaries = p.boundaries[2], factor = Dom)
    laplace!(du.x[3], u.x[3], st[3].D, p.dV, boundaries = p.boundaries[3])
    laplace!(du.x[4], u.x[4], st[4].D, p.dV, boundaries = p.boundaries[4])
    laplace!(du.x[5], u.x[5], st[5].D, p.dV, boundaries = p.boundaries[5])
    
    # decay
    @. du.x[1] -= st[1].decay * u.x[1]
    @. du.x[2] -= st[2].decay * u.x[2]
    @. du.x[3] -= st[3].decay * u.x[3]
    @. du.x[4] -= st[4].decay * u.x[4]
    @. du.x[5] -= st[5].decay * u.x[5]

    return nothing
end