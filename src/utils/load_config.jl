function load_config(fn = "scripts/config.toml")
    config = TOML.parsefile("scripts/config.toml")


    return config
end

function init_makie(config)
    if config["makie_inline"] && @isdefined Makie 
        Makie.inline!(true)
    end

    if @isdefined GLMakie
        GLMakie.activate!(ssao=true)
        GLMakie.closeall()
    end
end