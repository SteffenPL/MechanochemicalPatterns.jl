
# Model parameters 

const BDMGraph = BoundedDegreeDiGraph
const AdhesionGraph = BDMGraph{Int64}

# The state contains all data which is needed to produce the next 
# time step (provided the parameters are known).
@proto mutable struct State
    const X::Vector{SVecD} = SVecD[]  # position 
    const P::Vector{SVecD} = SVecD[]  # polarity
    const bonds::AdhesionGraph = BDMGraph(0, 2)  # adhesion bonds
    t::Float64 = 0.0
end

function partialcopy(s, lazy = false)
    return State(; 
        X = copy(s.X), 
        P = copy(s.P), 
        bonds = deepcopy(s.bonds), 
        t = s.t
        )
end

# Cache contains auxiliary data which is useful for the implementation, 
# but actually reduandant if one knows the parameters and the state.
@proto mutable struct Cache
    outdated::Bool = true  # for internal use

    # parameters 
    N::Int = 0

    # forces 
    const F::Vector{SVecD} = SVecD[]

    # collision detection 
    const st::BoundedHashTable{Dim,Vector{Int64},Float64,Int64}
end

# for better printing
Base.show(io::IO, s::State) = @printf io "State(%d cells @ t = %.2fh)" length(s.X) s.t
Base.show(io::IO, c::Cache) = @printf io "Cache(%d cells)" c.N



function init_state(p)
    
    n_disks = p.cells.n_disks
    n_bact = p.cells.n_bact

    N = n_bact * n_disks

    # initial cell positions
    pconst = random_direction(Dim)
    X = [ rand(SVecD) .* p.env.domain.size .+ p.env.domain.min for i in 1:N]
    P = [ p for i in 1:N]
    
    # adhesive network 
    bonds = BDMGraph(N, 2)

    for i in 1:n_disks:N
        for j in i+1:i+n_disks-1
            add_edge!(bonds, j, j-1)
            X[j] = X[j-1] - 0.1 * 2 * p.cells.R_hard * P[i]
            P[j] = P[i]
        end
    end

    # # initialize positions in a smart way 
    # st = BoundedHashTable(N, p.sim.collision_detection.boxsize, 
    #                         p.env.domain.min, p.env.domain.max)

    # n_steps = 10
    # dx = 2 * p.cells.R_hard / n_steps

    # for k_disks in 1:n_disks
    #     N_k = p.cells.n_bact * k_disks      # current number of disks we consider 
    #     X = X[1:k_disks, :]                 # consider only the first k_disks disks
    #     F = zeros(SVecD, k_disks, n_bact)   # force on each disk
    #     updatetable!(st, X)

    #     tmp_state = State(; X, P, bonds)
    #     tmp_cache = (;st, F)

    #     for k_steps in 1:n_steps
    #         # move leader forward 
    #         for i in axes(X,2)
    #             X[1, i] += dx * P[i]
    #         end

    #         # deal with contacts 
    #     end

    # end

    return State(; X, P, bonds)
end

function segment(i)
    n = p.cells.n_disks
    i_head = div(i-1, n) * n + 1
    return i_head:i_head + p.cells.n_disks - 1
end

function resize_cache!(s, p, cache)
    if cache.N != length(s.X)
        cache.N = length(s.X)
        for prop_name in propertynames(cache) 
            prob = getproperty(cache, prop_name)
            if prob isa Vector 
                resize!(prob, cache.N)
            end 
        end 

        resize!(cache.st, cache.N)
        cache.outdated = true
    end
end

function update_cache!(s, p, cache)
    if cache.outdated  
        cache.outdated = false
    end
end

function init_cache(p, s)
    cd = p.sim.collision_detection
    margin = (0.5 + cd.margin)
    dom = p.env.domain
    sht = BoundedHashTable(length(s.X), cd.boxsize, 
                            dom.min - margin .* dom.size, 
                            dom.max + margin .* dom.size)


    c = Cache(; st = sht)
    resize_cache!(s, p, c)
    update_cache!(s, p, c)
    return c
end

@inline function neighbours_bc(p, cache, pos, r)
    if p.env.periodic
        return periodic_neighbours(cache.st, pos, r)
    else
        return ( (i, SVecD(0.0, 0.0)) for i in neighbours(cache.st, pos, r) )
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
