include("base.jl")

rng_seeds = MersenneTwister()

function run_sim(fn = "$(@__DIR__)/inputs/parameters_2D.toml"; seed = missing, slider = true)
    
    if ismissing(seed)
        seed = rand(RandomDevice(), 1:1000)
    end

    Random.seed!(seed)
    printstyled("Setting seed to $seed\n", color = :green)
    
    p = load_parameters(fn)
    s = init_state(p)
    cache = init_cache(p, s)
    
    fig, s_obs = init_plot(s, p, cache; 
                                show_polarities = true, 
                                bottom_plots = false, 
                                show_concentration = true,
                                n_peaks = 10)
    display(fig)
        
    states = [deepcopy(s)]
    s = states[end]
    simulate(s, p, cache; callbacks = (update_plot_callback!(fig, s_obs, 0.01),), states = states)
    if slider
        add_slider!(fig, s_obs, states, p)
    end
    display(fig)

    return states, p, cache, s_obs, fig, seed
end


begin 
    last_seed = rand(rng_seeds, 1:1000)
    Random.seed!(last_seed)

    printstyled("Setting seed to $last_seed\n", color = :green)

    p = load_parameters("$(@__DIR__)/inputs/parameters_2D.toml")
    s = init_state(p)
    cache = init_cache(p, s)

    fig, s_obs = init_plot(s, p, cache; show_polarities = true, bottom_plots = true, show_concentration = true, show_v = false)
    display(fig)
    
    states = [deepcopy(s)]
    s = states[end]
    simulate(s, p, cache; callbacks = (update_plot_callback!(fig, s_obs, 0.05),), states = states)
    add_slider!(fig, s_obs, states, p)
end
display(fig)
play_animation!(fig, s_obs, states, 10)

include("analysis.jl")

fig, s_obs = init_plot(s, p, cache; show_polarities = true, bottom_plots = true, show_concentration = true, show_v = false)
display(fig)

record(fig, "example_spirals.mp4", 1:2:length(states); framerate = 30) do i
    update_plot!(fig, s_obs, states[i], p)
end

using Dates
function save_video(states, p, s_obs, fig, seed, folder = "./scripts/current/RT/videos", fn = "./scripts/current/RT/inputs/parameters_2D.toml")

    uuid = "sim_" * get(p.env, :name, "") * " " * Dates.format(now(), dateformat"yyyy-mm-dd HH:MM")
    folder = joinpath(folder, uuid)
    mkpath(folder)

    record(fig, joinpath(folder,"sim.mp4"), 1:2:length(states); framerate = 30) do i
        update_plot!(fig, s_obs, states[i], p)
    end

    cp(fn, joinpath(folder, "parameters.toml"))

    open(joinpath(folder, "seed.txt"), "w") do io
        write(io, string(seed))
    end
end

save_video(states, p, s_obs, last_seed)

function doit(seed = missing)
    states, p, cache, s_obs, fig, seed = run_sim(;seed, slider = false)
    save_video(states, p, s_obs, fig, seed)
    return states, p, cache, s_obs, fig, seed
end