# Note: all quantities carry units in terms of μm and hours.
using Revise
using GLMakie, StaticArrays, SpatialHashTables, TOML, ProtoStructs, Printf, Random
using MechanochemicalPatterns

config = load_config()
init_makie(config)

const Dim = 3
const SVecD = SVector{Dim, Float64}


function get_param(p, cell_type, sym, default = 0.0)
    ct = p.cells.types[cell_type]
    if hasproperty(ct, sym) 
        return ct[sym]
    elseif hasproperty(p.cells, sym)
        return p.cells[sym]
    else
        return default
    end
end

# Model parameters 
p_dict = TOML.parsefile("scripts/drafts/parameters.toml")

preprocess(x) = x
preprocess(x::String) = startswith(x, "julia:") ? eval(Meta.parse(x[7:end])) : x

p = recursive_namedtuple(p_dict, preprocess)

# maximal distance between complex interactions

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
    ax = LScene(fig[1, 1], scenekw=(ssao=Makie.SSAO(radius = 15.0, blur = 2),))
    ax.scene.ssao.bias[] = 0.25

    X_node = Observable(X)

    meshscatter!(X_node, markersize = 2 * p.cells.R_hard, space = :data, color = cell_type, colormap = [:magenta, :green], ssao=true)
    display(fig)
end


#=
@proto mutable struct State
    const X::Vector{SVecD} = SVecD[]
    const cell_type::Vector{Int} = Int[]
    t::Float64 = 0.0
end

Base.show(io::IO, s::State) = @printf io "State(%d cells @ t = %.2fh)" length(s.X) s.t

@proto mutable struct Cache 
    x::Int64
    y::Int64
end

s_init = State(; X, cell_type)

function simulate(p, s_init, callbacks = Function[])

    s = deepcopy(s_init)
    states = [deepcopy(s)]

    k = 1 
    t = 0.0 

    while t < t_end 
            
            t += dt 
    
            s = time_step(p, s, cache)
    
            if k % save_every == 0
                push!(states, deepcopy(s))
            end
    
            for callback in callbacks
                callback(p, s, cache)
            end

            for i in 1:N
                for j in 1:N
                    if i != j
                        r = norm(s.X[i] - s.X[j])
                        if r < p.cell.repulsion_radius
                            s.X[i] += (s.X[i] - s.X[j]) / r * p.cell.repulsion_stiffness
                        end
                    end
                end
            end



    
            k += 1
    end

    return states
end


function time_step(p, s, cache) 
    
    for i in 1:eachindex(s.X)
        s.X[i] += randn(3) .* sqrt(2 * p.sim.diffusion * dt)
    end

end


struct CellPlot 

end
=#