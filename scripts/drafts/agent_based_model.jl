# Note: all quantities carry units in terms of μm and hours.
using Revise
using TOML, Random, Printf
using GLMakie, StaticArrays, SpatialHashTables, ProtoStructs, ProgressMeter, Accessors

using MechanochemicalPatterns

config = load_config()
init_makie(config)

const Dim = 3
const SVecD = SVector{Dim, Float64}


# Model parameters 
p_dict = TOML.parsefile("scripts/drafts/parameters.toml")
p = recursive_namedtuple(p_dict)


dom = p.env.domain
SpatialHashTable((  min = SVecD(dom.center - 0.5 .* dom.size), 
                    max = SVecD(dom.center + 0.5 .* dom.size)), 
                    (p.sim.grid...,), 
                    100)

N = sum(ct.N for ct in p.cells.types)
cell_type = [ i for ct in p.cells.types for i in fill(ct.ID, ct.N) ]
shuffle!(cell_type)

# initial cell positions
@time X = [ SVecD(eval_param(get_param(p, cell_type[i], :init_pos))) for i in 1:N ]



begin
    fig = Figure()
    ax = LScene(fig[1, 1], scenekw=(ssao=Makie.SSAO(radius = 15.0, blur = 3),))
    ax.scene.ssao.bias[] = 0.025

    X_node = Observable(X)

    meshscatter!(X_node, markersize = 2 * p.cells.R_hard, space = :data, color = cell_type, colormap = [:magenta, :green], ssao=true)
    display(fig)
end


# The state contains all data which is needed to produce the next 
# time step (provided the parameters are known).
@proto mutable struct State
    const X::Vector{SVecD} = SVecD[]  # position 
    const P::Vector{SVecD} = SVecD[]  # polarity
    const cell_type::Vector{Int} = Int[]
    t::Float64 = 0.0
end

Base.show(io::IO, s::State) = @printf io "State(%d cells @ t = %.2fh)" length(s.X) s.t

# Cache contains auxiliary data which is useful for the implementation, 
# but actually reduandant if one knows the parameters and the state.
@proto mutable struct Cache 
    # parameters 
    N::Int = 0
    const R_soft::Vector{Float64} = Float64[]
    const R_hard::Vector{Float64} = Float64[]
    const t_divide::Vector{Float64} = Float64[]
    const repulsion_stiffness::Vector{Float64} = Float64[]
    const adhesion_stiffness::Vector{Float64} = Float64[]

    # forces 
    const F::Vector{SVecD} = SVecD[]
end

# for better printing
Base.show(io::IO, s::State) = @printf io "State(%d cells @ t = %.2fh)" length(s.X) s.t
Base.show(io::IO, c::Cache) = @printf io "Cache(%d cells)" c.N




function simulate(p, s_init, cache, callbacks = Function[])

    (; dt ) = p.sim

    s = deepcopy(s_init)
    states = [deepcopy(s)]

    t_unsaved = 0.0
    n_steps = Int(round(p.sim.t_end / dt))

    prog = Progress(n_steps, 1, "Simulating... ")
    for k_step in 1:n_steps 
        

        time_step!(p, s, cache)
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


function time_step!(p, s, cache) 
    
    sqrt_dt = sqrt(p.sim.dt)
    for i in eachindex(s.X)
        s.X[i] += randn(3) .* p.cells.sigma * sqrt_dt
    end
end


s_init = State(; X, cell_type)
cache = Cache(;)
p = @set p.sim.dt = 0.001
states = simulate(p, s_init, cache);


begin
    fig = Figure()
    ax = LScene(fig[1, 1], scenekw=(ssao=Makie.SSAO(radius = 15.0, blur = 3),))
    ax.scene.ssao.bias[] = 0.025

    i_slider = Slider(fig[2, 1], range = 1:length(states), startvalue = 1)
    X_node = @lift states[$(i_slider.value)].X

    meshscatter!(X_node, markersize = 2 * p.cells.R_hard, space = :data, color = cell_type, colormap = [:magenta, :green], ssao=true)
    display(fig)
end