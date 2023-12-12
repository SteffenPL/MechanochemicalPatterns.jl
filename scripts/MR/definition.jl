
# Model parameters 

const BDMGraph = BoundedDegreeDiGraph
const AdhesionGraph = BDMGraph{Int64}

# The state contains all data which is needed to produce the next 
# time step (provided the parameters are known).
@proto mutable struct State
    const X::Vector{SVecD} = SVecD[]  # position 
    const P::Vector{SVecD} = SVecD[]  # polarity
    const colony::Vector{Vector{Int}} = Vector{Int}[]  # segments
    const tsr::Vector{Float64} = Float64[]  # time since reversal
    const frustration::Vector{Float64} = Float64[]
    t::Float64 = 0.0
end

function partialcopy(s, lazy = false)
    return State(; 
        X = copy(s.X), 
        P = copy(s.P), 
        colony = deepcopy(s.colony), 
        clock = deepcopy(s.clock),
        frustration = deepcopy(s.frustration),
        t = s.t
        )
end

# Cache contains auxiliary data which is useful for the implementation, 
# but actually reduandant if one knows the parameters and the state.
@proto mutable struct Cache
    outdated::Bool = true  # for internal use

    # parameters 
    N::Int64 = 0

    # forces 
    const F::Vector{SVecD} = SVecD[]
    const Heads::Vector{SVecD} = SVecD[]
    const prevHeads::Vector{SVecD} = SVecD[]

    # flipping 
    const flipping::BitVector = BitVector([])

    # collision detection 
    const st::BoundedHashTable{Dim,Vector{Int64},Float64,Int64}
    const st_heads::BoundedHashTable{Dim,Vector{Int64},Float64,Int64}
end

# for better printing
Base.show(io::IO, s::State) = @printf io "State(%d cells @ t = %.2fh)" length(s.X) s.t
Base.show(io::IO, c::Cache) = @printf io "Cache(%d cells)" c.N


function segments(segments)
    return ( (segments[i], segments[i+1]) for i in 1:length(segments)-1 )
end

function init_state(p)
    
    n_disks = p.cells.n_disks
    n_bact = p.cells.n_bact

    N = n_bact * n_disks

    # initial cell positions
    X = [ rand(SVecD) .* p.env.domain.size .+ p.env.domain.min for i in 1:N]
    P = [ random_direction(Dim) for i in 1:N]
    
    # segments 
    colony = [collect(i:i+n_disks-1) for i in 1:n_disks:N]

    for i in 1:n_disks:N
        for j in i+1:i+n_disks-1
            X[j] = X[j-1] - 0.1 * 2 * p.cells.R_hard * P[i]
            P[j] = P[i]
        end
    end

    if hasproperty(p.cells.reversal_mechanism, :clock)
        tsr = rand(n_bact) .* p.cells.reversal_mechanism.clock.period
    else 
        tsr = Float64[]
    end

    frustration = zeros(n_bact)

    return State(; X, P, colony, tsr, frustration)
end


function resize_cache!(s, p, cache)
    if cache.N != length(s.X)
        cache.N = length(s.X)
        resize!(cache.st, cache.N)
        resize!(cache.F, cache.N)

        cache.outdated = true
    end

    if length(cache.Heads) != length(s.colony)
        resize!(cache.Heads, length(s.colony))
        resize!(cache.prevHeads, length(s.colony))
        resize!(cache.st_heads, length(s.colony))
        resize!(cache.flipping, length(s.colony))
        
        cache.outdated = true
    end
end

function update_cache!(s, p, cache)
    if cache.outdated  
        # after resize, one might need to update other data in the cache
        cache.outdated = false
    end

end

function init_cache(p, s)
    cd = p.sim.collision_detection
    margin = (0.5 + cd.margin)
    dom = p.env.domain
    st = BoundedHashTable(length(s.X), cd.boxsize, 
                            dom.min - margin .* dom.size, 
                            dom.max + margin .* dom.size)

    st_heads = BoundedHashTable(length(s.colony), cd.boxsize, 
                            dom.min - margin .* dom.size, 
                            dom.max + margin .* dom.size)

    c = Cache(; st, st_heads)
    resize_cache!(s, p, c)
    update_cache!(s, p, c)
    return c
end



function presimulate(p, substeps = 1, bending = p.cells.bending_stiffness)

    p = @set p.sim.dt = p.sim.dt / substeps 
    p = @set p.cells.bending_stiffness = bending
    p_ = @set p.cells.n_disks = 1
    
    s = init_state(p_)
    cache = init_cache(p_, s)

    for k in 2:p.cells.n_disks 
        n_steps =  round(Int64, (2 * p.cells.R_hard) / (p.cells.v_self_prop * p.sim.dt))
        for _ in 1:n_steps
            # update cache (in case of cell division)
            resize_cache!(s, p, cache)
            update_cache!(s, p, cache)
            updatetable!(cache.st, s.X)
        
            reset_forces!(s, p, cache)
            pull_polarities!(s, p, cache)
            compute_bending_forces!(s, p, cache)
            add_self_prop!(s, p, cache)
            
            # add forces
            for i in eachindex(s.X)
                s.X[i] += cache.F[i] * p.sim.dt / p.env.damping
            end
        
            project_non_overlap!(s, p, cache)
            project_bonds!(s, p, cache)
            project_onto_domain!(s, p, cache)
        end

        for (i, seg) in enumerate(s.colony)
            head = seg[1]
            tail = seg[end]
            P = s.P[tail]
            add_disk!(s, p, cache, i, s.X[tail] - 2*p.cells.R_hard*P, P)
        end
    end

    return s, cache
end








# helper functions
@inline function neighbours_bc(p, st, pos, r)
    if p.env.periodic
        return periodic_neighbours(st, pos, r)
    else
        return ( (i, zero(SVecD)) for i in neighbours(st, pos, r) )
    end
end

@inline function wrap(p, dx)
    if p.env.periodic
        dom = p.env.domain
        dx = @. dx - dom.size * round(dx * dom.inv_size)
    else 
        return dx
    end
end

import MechanochemicalPatterns: dist², dist
@inline function dist²(p, a, b)
    return sqrt(sum( z -> z^2, wrap(p, a - b)))
end

dist(p, a, b) = sqrt(dist²(p, a, b))
