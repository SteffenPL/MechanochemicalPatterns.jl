
# Model parameters 

const BDMGraph = BoundedDegreeMetaGraph
const AdhesionGraph = BDMGraph{Int,Float64,Nothing}

import ComponentArrays 
const CA = ComponentArrays

# The state contains all data which is needed to produce the next 
# time step (provided the parameters are known).
@proto mutable struct State{Dim}
    const X::Vector{SVector{Dim,Float64}} = SVector{Dim,Float64}[]  # position 
    const P::Vector{SVector{Dim,Float64}} = SVector{Dim,Float64}[]  # polarity
    const adh_bonds::AdhesionGraph = BoundedDegreeMetaGraph(0, 10)  # adhesion bonds
    const cell_type::Vector{Int} = Int[]
    const cell_age::Vector{Float64} = Float64[]
    const u::Array{Float64,Dim} = zeros(0,0,0)
    const v::Array{Float64,Dim} = zeros(0,0,0)
    t::Float64 = 0.0
end

function partialcopy(s, lazy = false)
    return State(; 
        X = copy(s.X), 
        P = copy(s.P), 
        adh_bonds = deepcopy(s.adh_bonds), 
        cell_type = s.cell_type,  # contant cell types 
        u = lazy ? s.u : copy(s.u), 
        v = lazy ? s.v : copy(s.v), 
        t = s.t
        )
end

# Cache contains auxiliary data which is useful for the implementation, 
# but actually reduandant if one knows the parameters and the state.
@proto mutable struct Cache{Dim, ODEP, ODEI}

    outdated::Bool = true  # for internal use

    # forces 
    const F::Vector{SVector{Dim,Float64}} = SVector{Dim,Float64}[]
    const dX::Vector{SVector{Dim,Float64}} = SVector{Dim,Float64}[]

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
    const new_adh_rate::Vector{Float64} = Float64[]
    const break_adh_rate::Vector{Float64} = Float64[]
    const run_time::Vector{Float64} = Float64[]

    # neighbour avg directions 
    const neighbour_avg::Vector{SVector{Dim,Float64}} = SVector{Dim,Float64}[]
    const neighbour_count::Vector{Int} = Int[]

    # ODE problem for the signals
    const ode_prob::ODEP = nothing
    const ode_integrator::ODEI = nothing


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
    X = [ svec(p)(eval_param(p, get_param(p, cell_type[i], :init_pos))) for i in 1:N ]
    P = [ random_direction(dim(p)) for i in 1:N]

    # initial age 
    cell_age = [ p.cells.lifespan * rand() for i in 1:N ]
    
    # adhesive network 
    adh_bonds = BDMGraph(N, 10)

    # signals 
    if dim(p) == 3
        xs, ys, zs = LinRange.( p.env.domain.min, p.env.domain.max, p.signals.grid )

        u = [Base.invokelatest(p.signals.types.u.init, (x,y,z), p) for x in xs, y in ys, z in zs]
        v = if hasproperty(p.signals.types, :v)
            [Base.invokelatest(p.signals.types.v.init, (x,y,z), p) for x in xs, y in ys, z in zs]
        else
            zeros(0,0,0)
        end
    elseif dim(p) == 2
        xs, ys = LinRange.( p.env.domain.min, p.env.domain.max, p.signals.grid )

        u = [Base.invokelatest(p.signals.types.u.init, (x,y), p) for x in xs, y in ys]
        v = if hasproperty(p.signals.types, :v)
            [Base.invokelatest(p.signals.types.v.init, (x,y), p) for x in xs, y in ys]
        else
            zeros(0,0)
        end
    else
        u = zeros(0,0)
        v = similar(u)
    end

    return State(; X, P, cell_type, cell_age, adh_bonds, u, v)
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
            update_p(:new_adh_rate)
            update_p(:break_adh_rate)
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

    function rhs_periodic!(dz, z, p_ode, t)
        p_ = p_ode.p
        dz .= 0.0

        if hasproperty(p_.signals.types, :u)
            pu = p_.signals.types.u
            laplace_periodic!(dz.u, z.u, pu.D, p_ode.dV)
            @. dz.u -= pu.decay * z.u
        end

        if hasproperty(p_.signals.types, :v)
            pv = p_.signals.types.v
            laplace_periodic!(dz.v, z.v, pv.D, p_ode.dV)
            @. dz.v -= pv.decay * z.v
        end
    end
    
    function rhs!(dz, z, p_ode, t)
        p_ = p_ode.p
        dz .= 0.0

        if hasproperty(p_.signals.types, :u)
            pu = p_.signals.types.u
            laplace!(dz.u, z.u, pu.D, p_ode.dV)
            @. dz.u -= pu.decay * z.u
        end

        if hasproperty(p_.signals.types, :v)
            pv = p_.signals.types.v
            laplace!(dz.v, z.v, pv.D, p_ode.dV)
            @. dz.v -= pv.decay * z.v
        end
    end

    if dim(p) == 3
        xs, ys, zs = LinRange.( p.env.domain.min, p.env.domain.max, p.signals.grid )
        p_ode = (;  p = p,
                    dV = (xs[2]-xs[1], ys[2]-ys[1], zs[2]-zs[1]))
    elseif dim(p) == 2
        xs, ys = LinRange.( p.env.domain.min, p.env.domain.max, p.signals.grid )
        p_ode = (;  p = p,
                    dV = (xs[2]-xs[1], ys[2]-ys[1]))
    end

    z0 = CA.ComponentArray(u = s.u, v = s.v)

    ode_prob = ODEProblem(p.env.periodic ? rhs_periodic! : rhs!, z0, (0.0, p.sim.t_end), p_ode)
    
    ode_integrator = init(ode_prob, Heun(); 
                save_everystep=false, 
                reltol = get(p.signals, :reltol, 1e-3), 
                abstol = get(p.signals, :abstol, 1e-6))

    c = Cache(; st = sht, ode_prob, ode_integrator, F = svec(p)[], dX = svec(p)[], neighbour_avg = svec(p)[])
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
