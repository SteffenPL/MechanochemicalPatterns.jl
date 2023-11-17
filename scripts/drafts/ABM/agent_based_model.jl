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
const R_max = 20.0

include("definition.jl")
include("forces.jl")
include("simulation.jl")



p = load_parameters()
s = init_state(p)
cache = init_cache(p, s)

# p = @set p.sim.dt = 0.1
@profview states = simulate(s, p, cache)

include("plots.jl")