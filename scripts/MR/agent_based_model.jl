include("base.jl")

#to set automatically the reversalmechanism before the simulation
#revm = p.cells.reversal_mechanism 
#p = @set p.cells.reversal_mechanism  = (; random = revm.random)


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

record(fig, "scripts/MR/videos/example2.mp4", 1:2:length(states); framerate = 30) do i
    update_plot!(fig, s_obs, states[i], p)
end