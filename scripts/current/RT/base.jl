# Note: all quantities carry units in terms of μm and hours.

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


includet("definition.jl")
includet("forces.jl")
includet("events.jl")
includet("signals.jl")
includet("modulation.jl")
includet("simulation.jl")
includet("analysis.jl")
includet("plots.jl")