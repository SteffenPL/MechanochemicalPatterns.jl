include("base.jl")


function run_sim(fn = "$(@__DIR__)/inputs/parameters_2D.toml")
    Random.seed!(4)
    p = load_parameters(fn)
    s = init_state(p)
    cache = init_cache(p, s)
    
    fig, s_obs = init_plot(s, p, cache; show_polarities = true, bottom_plots = false, show_concentration = true)
    display(fig)
        
    states = [deepcopy(s)]
    s = states[end]
    simulate(s, p, cache; callbacks = (update_plot_callback!(fig, s_obs, 0.05),), states = states)
    add_slider!(fig, s_obs, states, p)
    display(fig)
end

Random.seed!(4)
begin 
    p = load_parameters("$(@__DIR__)/inputs/parameters_2D.toml")
    s = init_state(p)
    cache = init_cache(p, s)

    fig, s_obs = init_plot(s, p, cache; show_polarities = true, bottom_plots = false)
    display(fig)
    
    states = [deepcopy(s)]
    s = states[end]
    @profview simulate(s, p, cache) #; callbacks = (update_plot_callback!(fig, s_obs, 0.05),), states = states)
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