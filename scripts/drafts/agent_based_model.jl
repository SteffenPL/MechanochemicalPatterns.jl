# Note: all quantities carry units in terms of μm and hours.
using Revise
using TOML, Random, Printf
using GLMakie, StaticArrays, ProtoStructs, ProgressMeter, Accessors 
using Graphs
using BoundedDegreeGraphs, SpatialHashTables
using OrdinaryDiffEq
using MechanochemicalPatterns

config = load_config()
init_makie(config)

const Dim = 3
const SVecD = SVector{Dim, Float64}


# Model parameters 

const BDMGraph = BoundedDegreeMetaGraph
const AdhesionGraph = BDMGraph{Int,Float64,Nothing}

# The state contains all data which is needed to produce the next 
# time step (provided the parameters are known).
@proto mutable struct State
    const X::Vector{SVecD} = SVecD[]  # position 
    const P::Vector{SVecD} = SVecD[]  # polarity
    const adh_bonds::AdhesionGraph = BoundedDegreeMetaGraph(0, 10)  # adhesion bonds
    const cell_type::Vector{Int} = Int[]
    const u::Matrix{Float64} = zeros(0,0)
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
    const attraction_stiffness::Vector{Float64} = Float64[]

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

    # adhesive network 
    adh_bonds = BDMGraph(N, 10)

    return State(; X, cell_type, adh_bonds)

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
            update_p(:attraction_stiffness)
        end
        cache.outdated = false
    end
end

function init_cache(p, s)
    cd = p.sim.collision_detection
    margin = cd.margin
    dom = p.env.domain
    sht = SpatialHashTable((  
                            min = SVecD(dom.min - margin .* dom.size), 
                            max = SVecD(dom.max + margin .* dom.size)), 
                        (cd.grid...,), 
                        length(s.X))

    c = Cache(; st = sht)
    resize_cache!(s, p, c)
    update_cache!(s, p, c)

    return c
end




function reset_forces!(s, p, cache)
    for i in eachindex(cache.F)
        cache.F[i] = zero(SVecD)
    end
end

function noise_kernel!(s, p, cache)
    sqrt_dt = sqrt(p.sim.dt)
    for i in eachindex(s.X)
        s.X[i] += randn(3) .* p.cells.sigma * sqrt_dt
    end
end

function add_bonds!(s, p, cache, i, j, Xi, Xj, dij)
    if dij < p.cells.R_adh && !has_edge(s.adh_bonds, i, j)
        if rand() < 1.0 - exp(-p.cells.new_adh_rate * p.sim.dt)
            add_edge!(s.adh_bonds, i, j, 0.0)
        end
    end
end

add_bonds!(s,p,cache) = apply_interaction_kernel!(s, p, cache, add_bonds!, p.cells.R_adh)

function apply_interaction_kernel!(s, p, cache, fnc, R)
    for i in eachindex(s.X)
        Xi = s.X[i]
        for j in neighbours(cache.st, Xi, R)
            if i < j 
                Xj = s.X[j]
                dij = dist(Xi, Xj)
                fnc(s, p, cache, i, j, Xi, Xj, dij)
            end
        end
    end
end

function remove_bonds!(s, p, cache)
    bonds = s.adh_bonds
    for e in edges(bonds)
        i, j = src(e), dst(e)
        dij = dist(s.X[i], s.X[j])

        if bonds[i,j] > p.cells.adh_duration || dij > 2*p.cells.R_adh
            rem_edge!(bonds, i, j)
        end
    end
end

function compute_adhesive_forces!(s, p, cache)
    bonds = s.adh_bonds
    for e in edges(bonds)
        i, j = src(e), dst(e)
        Xi, Xj = s.X[i], s.X[j]
        k = cache.adhesion_stiffness[i] + cache.adhesion_stiffness[j]
        cache.F[i] += k * (Xj - Xi)
        cache.F[j] -= k * (Xj - Xi)
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

function attraction_kernel!(s, p, cache, i, j, Xi, Xj, dij)
    Rij = cache.R_soft[i] + cache.R_soft[j] 
    if dij < Rij 
        k = cache.attraction_stiffness[i] + cache.attraction_stiffness[j]
        cache.F[i] -= (Rij - dij) * k / dij * (Xi - Xj) 
        cache.F[j] += (Rij - dij) * k / dij * (Xi - Xj) 
    end
end


function interaction_force_kernel!(s, p, cache, i, j, Xi, Xj, dij)
    repulsion_kernel!(s, p, cache, i, j, Xi, Xj, dij)
    attraction_kernel!(s, p, cache, i, j, Xi, Xj, dij)
end

function compute_interaction_forces!(s, p, cache)
    R_int = 2*p.cells.R_interact
    for i in eachindex(s.X)
        Xi = s.X[i]
        for j in neighbours(cache.st, Xi, R_int)  # 1:i-1
            if i < j 
                Xj = s.X[j]
                dij² = dist²(Xi, Xj)
                if dij² < 4*R_int^2
                    interaction_force_kernel!(s, p, cache, i, j, Xi, Xj, sqrt(dij²))
                end
            end
        end
    end
end


function project_non_overlap!(s, p, cache)
    for i in eachindex(s.X)
        Ri = cache.R_hard[i]
        for j in neighbours(cache.st, s.X[i], 20.0) # 1:i-1
            if i < j 
                dij² = dist²(s.X[i], s.X[j])
                Rij = Ri + cache.R_hard[j]
                if dij² < 4*Rij^2
                    dij = sqrt(dij²)
                    Xij = s.X[j] - s.X[i]
                    s.X[i] += (Rij - dij) / dij * Xij
                    s.X[j] -= (Rij - dij) / dij * Xij
                end
            end
        end
    end
end


function project_onto_domain!(s, p, cache)
    for i in eachindex(s.X)
        R_hard = cache.R_hard[i]
        s.X[i] = clamp.(s.X[i], p.env.domain.min .+ R_hard, p.env.domain.max .- R_hard)
    end     
end


function time_step!(s, p, cache) 


    # update cache (in case of cell division)
    update_cache!(s, p, cache)

    updateboxes!(cache.st, s.X)
    add_bonds!(s, p, cache)
    remove_bonds!(s, p, cache)

    # reset forces
    reset_forces!(s, p, cache)


    # update forces
    compute_adhesive_forces!(s, p, cache)
    compute_interaction_forces!(s, p, cache)

    # add forces
    for i in eachindex(s.X)
        s.X[i] += cache.F[i] * p.sim.dt / p.env.damping
    end

    # add noise 
    add_noise!(s, p, cache)
    
    # deal with constraints
    project_non_overlap!(s, p, cache)
    project_onto_domain!(s, p, cache)
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
    ax = Axis3(fig[1,1], aspect = :data) # LScene(fig[1, 1])
    alpha = clamp(20 / length(states[end].X), 0.01, 0.2)
    colors = [:magenta, :lightgreen]

    xlims!(ax, p.env.domain.min[1], p.env.domain.max[1])
    ylims!(ax, p.env.domain.min[2], p.env.domain.max[2])
    zlims!(ax, p.env.domain.min[3], p.env.domain.max[3])

    i_slider = Slider(fig[2, 1], range = 1:length(states), startvalue = 1)
    X_node = @lift states[$(i_slider.value)].X
    ct = @lift states[$(i_slider.value)].cell_type

    E_node = lift(i_slider.value) do i 
        bonds = states[i].adh_bonds
        @show ne(bonds)
        Tuple{SVecD,SVecD}[ (states[i].X[src(e)], states[i].X[dst(e)]) for e in edges(bonds) ]
    end

    meshscatter!(X_node, 
                markersize = p.cells.R_hard, space = :data, 
                color = ct, colormap = colors, transparency = true)

    meshscatter!(X_node, 
                markersize = p.cells.R_soft, space = :data, 
                color = ct, colormap = (c -> (c,alpha)).(colors), transparency = true)
    
    linesegments!(E_node, color = :darkorange, linewidth = 2, transparency = true)

    display(fig)
end
