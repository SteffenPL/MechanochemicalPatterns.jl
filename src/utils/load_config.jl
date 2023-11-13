function load_config(fn = "scripts/config.toml")
    config = TOML.parsefile("scripts/config.toml")


    return config
end

function init_makie(config)
    if config["makie_inline"]
        Makie.inline!(true)
    end

    GLMakie.activate!(ssao=true)
    GLMakie.closeall()
    println("Use SSAO rendering with GLMakie.")
end