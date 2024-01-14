include("base.jl")


function run_sim(fn = "$(@__DIR__)/inputs/parameters_2D.toml"; seed = missing)
    
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
    add_slider!(fig, s_obs, states, p)
    display(fig)

    return states, cache
end


Random.seed!(4)
begin 
    p = load_parameters("$(@__DIR__)/inputs/parameters_2D.toml")
    s = init_state(p)
    cache = init_cache(p, s)

    fig, s_obs = init_plot(s, p, cache; show_polarities = true, bottom_plots = false, show_concentration = true, show_v = false)
    display(fig)
    
    states = [deepcopy(s)]
    s = states[end]
    simulate(s, p, cache; callbacks = (update_plot_callback!(fig, s_obs, 0.05),), states = states)
end
add_slider!(fig, s_obs, states, p)
display(fig)
play_animation!(fig, s_obs, states, 10)

include("analysis.jl")

fig, s_obs = init_plot(s, p; show_polarities = true)
display(fig)

record(fig, "big4.mp4", 1:20:length(states); framerate = 30) do i
    update_plot!(fig, s_obs, states[i], p)
end