module MechanochemicalPatterns
    using TOML, Random, StaticArrays, GLMakie

    include("utils/load_config.jl")
    include("utils/recursive_namedtuple.jl")
    include("utils/eval_param.jl")

    include("agents/math_helpers.jl")

    export load_config, init_makie
    export recursive_namedtuple
    export eval_param, get_param

    export dist, dist²
end
