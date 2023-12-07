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

const Dim = 2
const SVecD = SVector{Dim, Float64}

include("definition.jl")
include("forces.jl")
include("simulation.jl")
include("plots.jl")
include("events.jl")