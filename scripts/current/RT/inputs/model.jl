
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



function pde!(du, u::ArrayPartition{Float64, Tuple{MT, MT, MT}}, p, t) where {MT <: Array}
    st = p.signals.types

    du .= 0.0
    
    laplace!(du.x[1], u.x[1], st[1].D, p.dV)
    laplace!(du.x[2], u.x[2], st[2].D, p.dV)
    laplace!(du.x[3], u.x[3], st[3].D, p.dV)
    
    # decay
    @. du.x[1] -= st[1].decay * u.x[1]
    @. du.x[2] -= st[2].decay * u.x[2]
    @. du.x[3] -= st[3].decay * u.x[3]

    return nothing
end