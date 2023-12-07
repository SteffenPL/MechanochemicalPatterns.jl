include("base.jl")



Random.seed!(4)
begin 
    p = load_parameters("scripts/MR/parameters.toml")
    s, cache = presimulate(p, 1, 1.0)

    states = [deepcopy(s)]
    fig, s_obs = init_plot(s, p)
    display(fig)

    simulate(s, p, cache; callbacks = (update_plot_callback!(fig, s_obs, 0.05),), states = states)

    add_slider!(fig, s_obs, states, p)
end
#play_animation!(fig, s_obs, states, 10)


fig, s_obs = init_plot(s, p)
display(fig)

record(fig, "scripts/MR/videos/example.mp4", 1:20:length(states); framerate = 30) do i
    update_plot!(fig, s_obs, states[i], p)
end