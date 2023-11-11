module MechanochemicalPatterns
    using TOML, Random, StaticArrays

    include("utils/load_config.jl")
    include("utils/recursive_namedtuple.jl")
    include("utils/eval_param.jl")

    export load_config, init_makie
    export recursive_namedtuple
    export eval_param
end
