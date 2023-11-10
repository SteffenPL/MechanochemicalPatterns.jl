module MechanochemicalPatterns
    using TOML

    include("utils/load_config.jl")
    include("utils/recursive_namedtuple.jl")
    include("utils/eval_param.jl")

    export load_config, init_makie
    export recursive_namedtuple
    export eval_param
end
