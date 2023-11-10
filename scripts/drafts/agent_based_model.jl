# Note: all quantities carry units in terms of μm and hours.
using Revise
using GLMakie, StaticArrays, SpatialHashTables, TOML, ProtoStructs
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
p = recursive_namedtuple(p_dict)

# maximal distance between complex interactions

dom = p.env.domain

SpatialHashTable((  min = SVecD(dom.center - 0.5 .* dom.size), 
                    max = SVecD(dom.center + 0.5 .* dom.size)), 
                    (p.sim.grid...,), 
                    100)

N = sum(ct.N for ct in p.cells.types)
cell_type = rand([1,2], N)

# initial cell positions
X = [ SVecD(eval_param(get_param(p, cell_type[i], :init_pos))...) for i in 1:N ]


begin
    fig = Figure()
    ax = LScene(fig[1, 1], scenekw=(ssao=Makie.SSAO(radius = 15.0, blur = 2),))
    ax.scene.ssao.bias[] = 0.25

    X_node = Observable(X)

    meshscatter!(X_node, markersize = 10.0, space = :data, color = cell_type, colormap = [:magenta, :green], ssao=true)
    fig
end

@proto mutable struct State
    X::Vector{SVecD} = SVecD[]
    cell_type::Vector{Int} = Int[]
    t::Float64 = 0.0
end



s = State()

function init_cache()

end

function simulate(p, s_init, callbacks = nothing)

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
    
            if !isnothing(callbacks)
                for callback in callbacks
                    callback(p, s, cache)
                end
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