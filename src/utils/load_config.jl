function load_config(fn = "scripts/config.toml")
    if !isfile(fn)
        println("Error: ", fn, " does not exist.")
        println("A possible fix is to copy the template file to ", fn, " and edit it.")
        return nothing
    end

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