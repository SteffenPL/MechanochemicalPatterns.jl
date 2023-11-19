
# Model parameters 

const BDMGraph = BoundedDegreeMetaGraph
const AdhesionGraph = BDMGraph{Int,Float64,Nothing}

# The state contains all data which is needed to produce the next 
# time step (provided the parameters are known).
Base.@kwdef mutable struct State
    const X::Vector{SVecD} = SVecD[]  # position 
    const P::Vector{SVecD} = SVecD[]  # polarity
    const adh_bonds::AdhesionGraph = BoundedDegreeMetaGraph(0, 10)  # adhesion bonds
    const cell_type::Vector{Int} = Int[]
    const u::Matrix{Float64} = zeros(0,0)
    t::Float64 = 0.0
end


# Cache contains auxiliary data which is useful for the implementation, 
# but actually reduandant if one knows the parameters and the state.
Base.@kwdef mutable struct Cache 
    outdated::Bool = true  # for internal use

    # parameters 
    N::Int = 0
    const R_soft::Vector{Float64} = Float64[]
    const R_hard::Vector{Float64} = Float64[]
    const R_adh::Vector{Float64} = Float64[]
    const R_attract::Vector{Float64} = Float64[]
    const t_divide::Vector{Float64} = Float64[]
    const repulsion_stiffness::Vector{Float64} = Float64[]
    const adhesion_stiffness::Vector{Float64} = Float64[]
    const attraction_stiffness::Vector{Float64} = Float64[]

    # forces 
    const F::Vector{SVecD} = SVecD[]

    # collision detection 
    const st::BoundedHashTable{Dim,Vector{Int64},Float64,Int64}
end

# for better printing
Base.show(io::IO, s::State) = @printf io "State(%d cells @ t = %.2fh)" length(s.X) s.t
Base.show(io::IO, c::Cache) = @printf io "Cache(%d cells)" c.N



function init_state(p)
    N = sum(ct.N for ct in p.cells.types)
    cell_type = [ i for ct in p.cells.types for i in fill(ct.ID, ct.N) ]
    shuffle!(cell_type)

    # initial cell positions
    X = [ SVecD(eval_param(p, get_param(p, cell_type[i], :init_pos))) for i in 1:N ]

    # adhesive network 
    adh_bonds = BDMGraph(N, 10)

    return State(; X, cell_type, adh_bonds)

end

function resize_cache!(s, p, cache)
    cache.N = length(s.X)
    for prop_name in propertynames(cache) 
        prob = getproperty(cache, prop_name)
        if prob isa Vector 
            resize!(prob, cache.N)
        end 
    end 
end

function update_cache!(s, p, cache)
    if cache.outdated  
        for i in eachindex(s.X)
            ct = s.cell_type[i]
            
            update_p(sym) = getproperty(cache, sym)[i] = get_param(p, ct, sym)

            update_p(:R_soft)
            update_p(:R_hard)
            update_p(:R_adh)
            update_p(:R_attract)
            update_p(:repulsion_stiffness)
            update_p(:adhesion_stiffness)
            update_p(:attraction_stiffness)
        end
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

