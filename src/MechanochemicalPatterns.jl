module MechanochemicalPatterns
    using TOML, Random, StaticArrays, GLMakie

    include("agents/math_helpers.jl")
    
    include("utils/load_config.jl")
    include("utils/recursive_namedtuple.jl")

    include("utils/parameters.jl")

    include("reaction_diffusion/operators.jl")

    export load_config, init_makie
    
    export load_parameters, eval_param, get_param

    export dist, dist²

    export laplace!
end
