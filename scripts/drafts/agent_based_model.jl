# Note: all quantities carry units in terms of μm and hours.
using Revise
using TOML, Random, Printf
using GLMakie, StaticArrays, SpatialHashTables, ProtoStructs, ProgressMeter, Accessors
using OrdinaryDiffEq
using MechanochemicalPatterns

config = load_config()
init_makie(config)

const Dim = 3
const SVecD = SVector{Dim, Float64}


# Model parameters 
function load_parameters(fn = "scripts/drafts/parameters.toml")
    p_dict = TOML.parsefile(fn)
    return recursive_namedtuple(p_dict)
end



# The state contains all data which is needed to produce the next 
# time step (provided the parameters are known).
@proto mutable struct State
    const X::Vector{SVecD} = SVecD[]  # position 
    const P::Vector{SVecD} = SVecD[]  # polarity
    const cell_type::Vector{Int} = Int[]
    t::Float64 = 0.0
end

# Cache contains auxiliary data which is useful for the implementation, 
# but actually reduandant if one knows the parameters and the state.
@proto mutable struct Cache 
    outdated::Bool = true  # for internal use

    # parameters 
    N::Int = 0
    const R_soft::Vector{Float64} = Float64[]
    const R_hard::Vector{Float64} = Float64[]
    const t_divide::Vector{Float64} = Float64[]
    const repulsion_stiffness::Vector{Float64} = Float64[]
    const adhesion_stiffness::Vector{Float64} = Float64[]

    # forces 
    const F::Vector{SVecD} = SVecD[]

    # collision detection 
    const st::SpatialHashTable{Dim}
end

# for better printing
Base.show(io::IO, s::State) = @printf io "State(%d cells @ t = %.2fh)" length(s.X) s.t
Base.show(io::IO, c::Cache) = @printf io "Cache(%d cells)" c.N



function init_state(p)
    N = sum(ct.N for ct in p.cells.types)
    cell_type = [ i for ct in p.cells.types for i in fill(ct.ID, ct.N) ]
    shuffle!(cell_type)

    # initial cell positions
    @time X = [ SVecD(eval_param(get_param(p, cell_type[i], :init_pos))) for i in 1:N ]

    return State(; X, cell_type)
end

function init_cache(p, s)
    margin = 1.0
    dom = p.env.domain
    sht = SpatialHashTable((  min = SVecD(dom.center - (0.5 + margin) .* dom.size), 
                        max = SVecD(dom.center + (0.5 + margin) .* dom.size)), 
                        (p.sim.grid...,), 
                        length(s.X))

    return Cache(; st = sht)
end


# begin
#     fig = Figure()
#     ax = LScene(fig[1, 1], scenekw=(ssao=Makie.SSAO(radius = 15.0, blur = 3),))
#     ax.scene.ssao.bias[] = 0.025

#     X_node = Observable(X)

#     meshscatter!(X_node, markersize = 2 * p.cells.R_hard, space = :data, color = cell_type, colormap = [:magenta, :green], ssao=true)
#     display(fig)
# end



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
            update_p(:repulsion_stiffness)
            update_p(:adhesion_stiffness)
        end
        cache.outdated = false
    end
end

function reset_forces!(s, p, cache)
    for i in eachindex(cache.F)
        cache.F[i] = zero(SVecD)
    end
end

function add_noise!(s, p, cache)
    sqrt_dt = sqrt(p.sim.dt)
    for i in eachindex(s.X)
        s.X[i] += randn(3) .* p.cells.sigma * sqrt_dt
    end
end

function repulsion_kernel!(s, p, cache, i, j, Xi, Xj, dij)
    Rij = cache.R_soft[i] + cache.R_soft[j] 
    if dij < Rij 
        k = cache.repulsion_stiffness[i] + cache.repulsion_stiffness[j]
        cache.F[i] += (Rij - dij) * k / dij * (Xi - Xj) 
        cache.F[j] -= (Rij - dij) * k / dij * (Xi - Xj) 
    end
end

function interaction_force_kernel!(s, p, cache, i, j, Xi, Xj, dij)
    repulsion_kernel!(s, p, cache, i, j, Xi, Xj, dij)
end

function compute_interaction_forces!(s, p, cache)
    for i in eachindex(s.X)
        Xi = s.X[i]
        for j in neighbours(cache.st, Xi, 20.0)  # 1:i-1
            if i < j 
                Xj = s.X[j]
                dij = dist(Xi, Xj)
                if dij < 2*p.cells.R_interact 
                    interaction_force_kernel!(s, p, cache, i, j, Xi, Xj, dij)
                end
            end
        end
    end
end


function time_step!(s, p, cache) 
    
    # update cache (in case of cell division)
    update_cache!(s, p, cache)

    # reset forces
    reset_forces!(s, p, cache)

    updateboxes!(cache.st, s.X)

    # update forces
    compute_interaction_forces!(s, p, cache)

    # add forces
    for i in eachindex(s.X)
        s.X[i] += cache.F[i] * p.sim.dt / p.env.damping
    end

    # add noise 
    add_noise!(s, p, cache)
    
end




function simulate(s, p, cache, callbacks = Function[])

    (; dt ) = p.sim

    s = deepcopy(s)
    states = [deepcopy(s)]

    resize_cache!(s, p, cache)

    t_unsaved = 0.0
    n_steps = Int(round(p.sim.t_end / dt))

    prog = Progress(n_steps, 1, "Simulating... ")
    for k_step in 1:n_steps 
        
        time_step!(s, p, cache)
        s.t += dt 

        if t_unsaved > p.sim.saveat
            push!(states, deepcopy(s))
            t_unsaved = 0.0
        end
        t_unsaved += dt

        next!(prog)
    end

    return states
end


p = load_parameters()
s = init_state(p)
cache = init_cache(p, s)

p = @set p.sim.dt = 0.1
@time states = simulate(s, p, cache)


begin
    fig = Figure()
    ax = LScene(fig[1, 1])
    alpha = clamp(20 / length(states[end].X), 0.01, 0.2)

    i_slider = Slider(fig[2, 1], range = 1:length(states), startvalue = 1)
    X_node = @lift states[$(i_slider.value)].X
    ct = @lift states[$(i_slider.value)].cell_type

    meshscatter!(X_node, 
                markersize = p.cells.R_hard, space = :data, 
                color = ct, colormap = [:magenta, :green], transparency = true)


    meshscatter!(X_node, 
                markersize = p.cells.R_soft, space = :data, 
                color = ct, colormap = [(:magenta,alpha), (:green,alpha)], transparency = true)
    display(fig)
end
