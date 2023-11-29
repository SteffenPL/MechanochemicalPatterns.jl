# Note: all quantities carry units in terms of μm and hours.
using Revise
using TOML, Random, Printf, LinearAlgebra, Random
using GLMakie, StaticArrays, ProtoStructs, ProgressMeter, Accessors, Graphs
using BoundedDegreeGraphs, SpatialHashTables
using OrdinaryDiffEq
using MechanochemicalPatterns

config = load_config()
init_makie(config)

const Dim = 3
const SVecD = SVector{Dim, Float64}
const R_max = 20.0

include("definition.jl")
include("forces.jl")
include("signals.jl")
include("simulation.jl")
include("plots.jl")

Random.seed!(4)

p = load_parameters()
s = init_state(p)
cache = init_cache(p, s)
states = simulate(s, p, cache)

fig, s_obs = init_plot(states[1], p)
add_slider!(fig, s_obs, states, p)
display(fig)
play_animation!(fig, s_obs, states, 10)
