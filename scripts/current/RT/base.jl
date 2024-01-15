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


dim(x::Val{n}) where n = n
dim(p) = dim(p.env.dim)
@inline svec(p) = SVector{dim(p), Float64}


include("definition.jl")
include("forces.jl")
include("events.jl")
include("signals.jl")
include("modulation.jl")
include("simulation.jl")
include("analysis.jl")
include("plots.jl")