module MechanochemicalPatterns
    using TOML, Random, StaticArrays, GLMakie

    include("utils/load_config.jl")
    include("utils/recursive_namedtuple.jl")

    include("utils/parameters.jl")

    include("agents/math_helpers.jl")

    export load_config, init_makie
    
    export load_parameters, eval_param, get_param

    export dist, dist²
end
