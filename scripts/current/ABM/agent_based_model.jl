# Note: all quantities carry units in terms of μm and hours.
using Revise
using TOML, Random, Printf, LinearAlgebra, Random
using GLMakie, StaticArrays, ProtoStructs, ProgressMeter, Accessors, Graphs
using BoundedDegreeGraphs, SpatialHashTables
using OrdinaryDiffEq
using MechanochemicalPatterns

config = load_config()
init_makie(config)
set_theme!(theme_black())

const Dim = 3
const SVecD = SVector{Dim, Float64}
const R_max = 20.0

include("definition.jl")
include("forces.jl")
include("events.jl")
include("signals.jl")
include("simulation.jl")
include("plots.jl")

Random.seed!(4)
begin 
    p = load_parameters()
    s = init_state(p)
    cache = init_cache(p, s)

    fig, s_obs = init_plot(s, p; show_polarities = true)
    display(fig)
    
    states = [deepcopy(s)]
    s = states[end]
    simulate(s, p, cache; callbacks = (update_plot_callback!(fig, s_obs, 0.05),), states = states)
end
add_slider!(fig, s_obs, states, p)
display(fig)
play_animation!(fig, s_obs, states, 10)

include("analysis.jl")

fig, s_obs = init_plot(s, p; show_polarities = true)
display(fig)

record(fig, "big4.mp4", 1:20:length(states); framerate = 30) do i
    update_plot!(fig, s_obs, states[i], p)
end