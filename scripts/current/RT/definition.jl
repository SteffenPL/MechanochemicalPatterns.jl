
# Model parameters 

const BDMGraph = BoundedDegreeMetaGraph
const AdhesionGraph = BDMGraph{Int,Float64,Nothing}

using RecursiveArrayTools
import ComponentArrays 
const CA = ComponentArrays

# The state contains all data which is needed to produce the next 
# time step (provided the parameters are known).
@proto mutable struct State{Dim,ST}
    const X::Vector{SVector{Dim,Float64}} = SVector{Dim,Float64}[]  # position 
    const P::Vector{SVector{Dim,Float64}} = SVector{Dim,Float64}[]  # polarity
    const adh_bonds::AdhesionGraph = BoundedDegreeMetaGraph(0, 10)  # adhesion bonds
    const cell_type::Vector{Int} = Int[]
    const cell_age::Vector{Float64} = Float64[]
    const U::ST = nothing
    t::Float64 = 0.0
end

function partialcopy(s, lazy = false)
    return State(; 
        X = copy(s.X), 
        P = copy(s.P), 
        adh_bonds = deepcopy(s.adh_bonds), 
        cell_type = s.cell_type,  # contant cell types...
        U = lazy ? s.U : copy.(s.U), 
        t = s.t
        )
end

@proto struct CellCache{SVT, SVTD} 
    F::SVT
    dX::SVT
    grads::SVTD
    neighbour_avg::SVT
    neighbour_count::Int

    R_soft::Float64
    R_hard::Float64
    R_adh::Float64
    R_attract::Float64
    t_divide::Float64
    repulsion_stiffness::Float64
    adhesion_stiffness::Float64
    attraction_stiffness::Float64
    biased_adhesion::Float64
    new_adh_rate::Float64
    break_adh_rate::Float64
    run_time::Float64
end

using StructArrays

# Cache contains auxiliary data which is useful for the implementation, 
# but actually reduandant if one knows the parameters and the state.
@proto mutable struct Cache{Dim, ODEP, ODEI, ODEG, SA}

    outdated::Bool = true  # for internal use
    N::Int = 0

    const data::SA

    # ODE problem for the signals
    const ode_prob::ODEP = nothing
    const ode_integrator::ODEI = nothing
    const grid::ODEG = nothing

    # collision detection 
    const st::BoundedHashTable{Dim,Vector{Int64},Float64,Int64}
end

# for better printing
Base.show(io::IO, s::State) = @printf io "State(%d cells @ t = %.2fh)" length(s.X) s.t
Base.show(io::IO, c::Cache) = @printf io "Cache(%d cells)" c.N

function CellDataType(p)
    Dim = dim(p)
    return @NamedTuple begin 
        F::SVector{Dim,Float64}
        dX::SVector{Dim,Float64}
        grad::SVector{Dim,Float64}
        neighbour_avg::SVector{Dim,Float64}
        neighbour_count::Int
        R_soft::Float64
        R_hard::Float64
        R_adh::Float64
        R_attract::Float64
        t_divide::Float64
        repulsion_stiffness::Float64
        adhesion_stiffness::Float64
        attraction_stiffness::Float64
        biased_adhesion::Float64
        new_adh_rate::Float64
        break_adh_rate::Float64
        run_time::Float64
    end
end 

function init_state(p)
    N = sum(ct.N for ct in p.cells.types)

    cell_type = reduce(vcat, fill(ct_index, ct.N) for (ct_index, ct) in enumerate(p.cells.types))
    shuffle!(cell_type)

    # initial cell positions
    X = [ svec(p)(eval_param(p, get_param(p, cell_type[i], :init_pos))) for i in 1:N ]
    P = [ random_direction(dim(p)) for i in 1:N]

    # initial age 
    cell_age = [ p.cells.lifespan * rand() for i in 1:N ]
    
    # adhesive network 
    adh_bonds = BDMGraph(N, 10)

    # signals 
    if hasproperty(p, :signals)
        cts = p.signals.types

        grid = Tuple(LinRange.( p.env.domain.min, p.env.domain.max, p.signals.grid))
        n_signals = length(cts)

        U = map( i -> Array{Float64}(undef, p.signals.grid...), Tuple(1:n_signals))

        for (k, u) in enumerate(U)
            ct = cts[k]
            for I in CartesianIndices(u)
                pos = getindex.(grid, Tuple(I))
                u[I] = eval_param(p, ct.init, pos)
            end
        end
    else
        U = Float64[]
    end

    return State(; X, P, cell_type, cell_age, adh_bonds, U = ArrayPartition(U...), t = 0.0)
end

function resize_cache!(s, p, cache)
    if cache.N != length(s.X)
        cache.N = length(s.X)

        resize!(cache.data, cache.N)
        resize!(cache.st, cache.N)

        cache.outdated = true
    end
end

function update_cache!(s, p, cache)
    if cache.outdated  
        for i in eachindex(s.X)
            ct = s.cell_type[i]
            
            update_p(sym) = getproperty(cache.data, sym)[i] = get_param(p, ct, sym)

            update_p(:R_soft)
            update_p(:R_hard)
            update_p(:R_adh)
            update_p(:R_attract)
            update_p(:repulsion_stiffness)
            update_p(:adhesion_stiffness)
            update_p(:attraction_stiffness)
            update_p(:new_adh_rate)
            update_p(:break_adh_rate)
            update_p(:biased_adhesion)
        end
        cache.outdated = false
    end
end

function init_cache(p, s)
    cd = p.sim.collision_detection
    dom = p.env.domain

    margin = (p.env.periodic ? 0.0 : cd.margin)

    sht = BoundedHashTable(length(s.X), cd.boxsize, 
                            dom.min - margin .* dom.size, 
                            dom.max + margin .* dom.size)

    # function rhs_periodic!(dz, z, p_ode, t)
    #     p_ = p_ode.p
    #     dz .= 0.0

    #     if hasproperty(p_.signals.types, :u)
    #         pu = p_.signals.types.u
    #         laplace_periodic!(dz.u, z.u, pu.D, p_ode.dV)
    #         @. dz.u -= pu.decay * z.u
    #     end

    #     if hasproperty(p_.signals.types, :v)
    #         pv = p_.signals.types.v
    #         laplace_periodic!(dz.v, z.v, pv.D, p_ode.dV)
    #         @. dz.v -= pv.decay * z.v
    #     end
    # end
    
    # function rhs!(dz, z, p_ode, t)
    #     p_ = p_ode.p
    #     dz .= 0.0

    #     if hasproperty(p_.signals.types, :u)
    #         pu = p_.signals.types.u
    #         laplace!(dz.u, z.u, pu.D, p_ode.dV)
    #         @. dz.u -= pu.decay * z.u
    #     end

    #     if hasproperty(p_.signals.types, :v)
    #         pv = p_.signals.types.v
    #         laplace!(dz.v, z.v, pv.D, p_ode.dV)
    #         @. dz.v -= pv.decay * z.v
    #     end
    # end

    if hasproperty(p, :signals)
        
        grid = Tuple(LinRange.( p.env.domain.min, p.env.domain.max, p.signals.grid))
        z0 = s.U
        p_ode = (p..., dV = p.env.domain.size ./ p.signals.grid, tmp = similar(z0))
        ode_prob = ODEProblem(pde!, z0, (0.0, p.sim.t_end), p_ode)
        
        ode_integrator = init(ode_prob, Heun(); 
                    save_everystep=false, 
                    reltol = get(p.signals, :reltol, 1e-3), 
                    abstol = get(p.signals, :abstol, 1e-6))
    else
        ode_prob = nothing
        ode_integrator = nothing
        grid = nothing
    end
    
    data = StructArray{CellDataType(p)}(undef, length(s.X))
    c = Cache(; st = sht, data, ode_prob, ode_integrator, grid)
    
    resize_cache!(s, p, c)
    update_cache!(s, p, c)

    return c
end



# helper functions
@inline function neighbours_bc(p, st, pos, r)
    if p.env.periodic
        return periodic_neighbours(st, pos, r)
    else
        return ( (i, zero(svec(p))) for i in neighbours(st, pos, r) )
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
    return sum( z -> z^2, wrap(p, a - b))
end

dist(p, a, b) = sqrt(dist²(p, a, b))
